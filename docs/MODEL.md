# Model and methods

Every governing equation, numerical method and data source in the engine, with references. Validation
numbers are in [VALIDATION.md](VALIDATION.md).

Units throughout the constraint-based part are mmol/gDW/h for fluxes and 1/h for growth, as in the
models. Kinetic models carry the units their SBML defines.

## 1. Stoichiometric models

A model has metabolites `i`, reactions `j`, a stoichiometric matrix `S` with `S_ij < 0` for substrates,
bounds `l <= v <= u`, and gene-protein-reaction rules. SBML species with `boundaryCondition="true"` are
exempt from the balance. Reading: SBML Levels 2 and 3 with FBC versions 1 to 3 (`fbc:fluxBound`,
`fbc:lowerFluxBound`, `fbc:geneProductAssociation`), the groups package for subsystems, COBRA-style
`GENE_ASSOCIATION` notes as a fallback. Initial assignments and rules that define flux bounds or
stoichiometries are evaluated by the kinetic compiler of section 4. Writing: SBML Level 3 Version 1,
declaring FBC version 2 and groups only when a model uses them.

### Gene deletions

A rule is a Boolean formula over genes, with `and` joining subunits and `or` joining isozymes. Deleting
genes sets `l_j = u_j = 0` for each reaction whose rule is false (Thiele and Palsson 2010). Reactions
without a rule are never disabled. In BiGG models the placeholder `s0001` stands for spontaneous
reactions; the low-level deletion functions treat it as a gene, matching COBRApy, and the screens in the
desktop application and the notebook leave it out because it cannot be deleted.

## 2. Solvers

Programs have the form `min or max c'x + x'Qx/2`, `a <= A x <= b`, `l <= x <= u`, with integrality on
marked columns. Two backends implement it.

HiGHS keeps matrix entries above 1e-12 (`small_matrix_value`). Its default of 1e-9 removes real enzyme costs from GECKO
models: `eciML1515` has 17 entries at or below it, and solving that file with the default overstates anaerobic growth by
2e-5, because the affected enzymes become free.

**HiGHS** (Huangfu and Hall 2018) solves linear, mixed-integer and convex quadratic programs, through its
C API. Feasibility tolerances are 1e-9, the relative MIP gap 1e-6. Each worker thread keeps one session
and re-solves from the previous basis after a change of bounds or objective.

**Dense simplex** (`lp/simplex.rs`) is a bounded-variable primal simplex written here, for linear programs
only. It exists for WebAssembly, where HiGHS is unavailable, and as an independent check on HiGHS: dense
storage limits it to a few thousand rows. Each row `i` gets a logical variable, `A x - s = 0`, with the
row's bounds moved to `s`. Composite phase 1 minimizes the sum of bound violations of basic variables.
The inverse of the basis is updated by Gauss-Jordan pivots and rebuilt from an LU factorization every 64
iterations. Ratio tests use two passes with a tolerance, preferring the largest pivot (Harris 1973).
After 200 degenerate pivots in a row, pricing and leaving choice switch to Bland's rule (Bland 1977),
which guarantees termination.

## 3. Constraint-based analyses

**FBA** (Orth, Thiele and Palsson 2010): `max c'v` subject to `S v = 0`, `l <= v <= u`.

**Parsimonious FBA** (Lewis et al. 2010): with `z*` the FBA optimum and `f` a retained fraction, minimize
`sum (v+_j + v-_j)` subject to `S (v+ - v-) = 0` and `c'(v+ - v-) >= f z*`, on separate forward and reverse
parts of each flux.

**Flux variability** (Mahadevan and Schilling 2003): minimum and maximum of each `v_j` subject to
`S v = 0`, `l <= v <= u`, `c'v >= f z*`. No slack is added to that last row: some fluxes respond ~1000-fold
to the objective (`ICDHyr` in `e_coli_core` moves 920 units per unit of growth), so 1e-9 of slack is
visible at 1e-6. LPs share a session per worker (Gudmundsson and Thiele 2010).

**Loopless FBA** (Schellenberger, Lewis and Palsson 2011): for each internal reaction a binary `a_j` and a
potential `G_j`, with `-M(1 - a_j)`-style rows `v_j <= u_j a_j`, `v_j >= l_j (1 - a_j)`, `G_j <= -1 + M(1 - a_j)`,
`G_j >= 1 - (M + 1) a_j`, and `N_int' G = 0` over a null-space basis of the internal stoichiometry.

**CycleFreeFlux** (Desouki et al. 2015): given `v`, minimize `sum |w_j|` over internal reactions subject to
`S w = 0`, `w_j` between 0 and `v_j`, and boundary and objective reactions fixed at their values in `v`.

**MOMA** (Segre, Vitkup and Church 2002): minimize `sum (v_j - w_j)^2` over the perturbed network for a
reference `w`, as a convex QP. **Linear MOMA** minimizes `sum |v_j - w_j|`.

**ROOM** (Shlomi, Berkman and Ruppin 2005): minimize `sum y_j`, binary `y_j`, with
`w^u_j = w_j + delta |w_j| + epsilon`, `w^l_j = w_j - delta |w_j| - epsilon`, and
`v_j - y_j (u_j - w^u_j) <= w^u_j`, `v_j - y_j (l_j - w^l_j) >= w^l_j`. Defaults: `delta = 0.03`, `epsilon = 0.001`.

**Production envelope**: on an even grid of a control flux over its feasible range, the minimum and maximum
of a target flux.

**Minimal medium**: minimize total uptake subject to `c'v >= g`. The component-minimizing form first
minimizes the number of open exchanges with binaries, then total uptake among them.

**Dynamic FBA** (Mahadevan, Edwards and Doyle 2002; Varma and Palsson 1994): static optimization. Over a step
`dt` the fluxes are constant, with uptake `v_i >= -vmax_i S_i / (K_i + S_i)` and capped so that a substrate
cannot be exhausted within the step. Then `X(t + dt) = X e^(mu dt)` and
`S_i(t + dt) = S_i + v_i X (e^(mu dt) - 1) / mu`. The scheme is exact while fluxes stay constant and first-order
once they change.

### Flux sampling

Reactions whose FVA range is under 1e-9 are fixed, the remaining fluxes are written `v = x0 + N z` over an
orthonormal null-space basis, and FVA ranges serve as bounds. Keeping an objective fraction adds its
half-space.

**CHRR** (Haraldsdottir et al. 2017): coordinate hit-and-run from the centre of the largest inscribed ball,
found by a linear program, along columns of a rounded basis. *This implementation rounds differently from
the original*: it moves to isotropic position with `T T' = C` for the covariance `C` of samples, over 5 rounds
of `60 d` points (Lovasz and Vempala 2006), rather than using the maximum-volume ellipsoid. Rounding affects
only mixing: the target stays uniform.

**OptGP** (Megchelenbrink, Huynen and Marchiori 2014): chains of artificial-centering hit-and-run from the
centre of FVA-derived warmup points, with projection onto the null space every 100 moves. *A heuristic
sampler*: since artificial centering adapts its directions, it has no exact guarantee of uniformity, so read
its output with the diagnostics below.

**Diagnostics**: split R-hat (Gelman et al. 2013), effective sample size from autocorrelations summed over
Geyer's initial positive sequence, and the Geweke (1992) z-score of the first 10 % against the last 50 %.

### Thermodynamics-based flux analysis

Henry, Broadbelt and Hatzimanikatis (2007). For each reaction with a known standard energy:

    DrG'_j = DrG'0_j + e_j + RT sum_i n_ij ln c_i,      |e_j| <= k sigma_j

with `ln c_i` bounded, binaries `z+_j`, `z-_j` with `z+_j + z-_j <= 1`, and
`v_j <= u_j z+_j`, `v_j >= l_j z-_j`, `DrG'_j <= -eps + M (1 - z+_j)`, `DrG'_j >= eps - M (1 - z-_j)`. Water and
protons do not enter the sum, since their activities are in the transformed standard state. Defaults:
`c` from 1e-6 to 0.02 M, `k = 2`, `M = 1000`, `eps = 1e-3` kJ/mol, `T` from the data.

*Source of the energies*: `DrG'0` is not computed by this crate. `scripts/export_thermo.py` reads component
contribution estimates (Noor et al. 2013) from eQuilibrator (Beber et al. 2022) at pH 7.5, pMg 3, I = 0.25 M,
for reactions within one compartment. Transport across a membrane needs the membrane potential and is left
unconstrained.

### Enzyme constraints

GECKO-style (Sanchez et al. 2017). Enzyme `k` has a usage `e_k` in mmol/gDW, so with turnover numbers in 1/h and
molecular weights in g/mmol:

    v_j <= sum_k kcat+_jk e_k,   v_j >= -sum_k kcat-_jk e_k,   sum_k MW_k e_k <= P,   e_k <= abundance_k

A published GECKO file can also be read as it stands, since it encodes enzymes as metabolites and reactions. Both
routes reach the same optimum on `eciML1515` (`validation/v12_ecmodel.py`). A complex enters as one enzyme whose
weight is what its subunits take per catalytic cycle, `MW = kcat sum_s c_s MW_s`, with `c_s = n_s / kcat`.

Turnover numbers, weights and the pool `P` are inputs. The module holds no biological data, and validation
covers its arithmetic only.

### Metabolic tasks

The format and semantics of RAVEN's `checkTasks` (Agren et al. 2013; Robinson et al. 2020). A task passes when
a steady state exists with its inputs and outputs as the system's whole boundary: an input lets a metabolite's
balance row reach `S v` in `[-UB, -LB]`, an output in `[LB, UB]`. Extra reactions, bound changes and
`ALLMETSIN[c]` tokens are supported.

## 4. Kinetic models

### Compilation

SBML core Levels 2 and 3. Function definitions are inlined, local parameters replaced by their values, and
symbols resolved to slots in a value array. Each species has two slots, the value math sees (a concentration, or
an amount with `hasOnlySubstanceUnits`) and its amount. The state is:

- the amount of each species not otherwise determined;
- the math value of a species with a rate rule or that an algebraic rule determines;
- each compartment, parameter or species reference with a rate rule, an algebraic rule or an event assignment.

All other values are derived, in dependency order: assignment rules, species views `S = n / V`, reaction rates,
`rateOf`. The equations are

    dn_i/dt = conv_i sum_j nu_ij v_j        dx/dt = f(x)  (rate rules)        0 = g(y)  (algebraic rules)

The variable each algebraic rule determines comes from bipartite matching over symbols no other rule or reaction
determines, preferring symbols that no event assigns. Algebraic variables start from Newton iteration and are
restored by Newton iteration after every event.

**Initial values**: compartment sizes, attributes, initial assignments and assignment rules are one system of
definitions, evaluated repeatedly until they stop changing, which resolves them in any order.

**`delay(x, tau)`**: from a stored history of accepted steps, with cubic Hermite interpolation of the state.
Before the start time, the initial definitions are re-evaluated at that earlier time and rates are zero, which
is how the SBML Test Suite treats history before the start. Steps are limited to the smallest constant lag.
This is a method of steps with interpolation, and it is not designed for state-dependent lag discontinuities.

**Fast reactions** (removed from SBML in Level 3 Version 2): *a stand-in*, the singular perturbation of
multiplying the law by `1e8`, with the fast subsystem brought to equilibrium in pseudo-time at the start and after
events. The residual of the equilibrium condition is of the order of `1e-8`.

**Jacobians** come from symbolic differentiation (`kinetic/diff.rs`) of the right-hand side with all derived values
inlined; `piecewise` differentiates branch by branch, step-like operators to zero, and `delay` is treated as
independent of the current state, so its contribution is absent from the Jacobian (which affects only the rate of
convergence of the integrators). Past a size limit of 4 million nodes the integrators fall back on forward
differences.

### Integrators

**Rodas4** (Hairer and Wanner 1996, coefficients from their `rodas.f`): L-stable, stiffly accurate Rosenbrock
method of order 4 with an embedded order-3 estimate, for `M y' = f(t, y)` with a diagonal `M` and index-1 algebraic
equations. Step control follows `rodas.f`: `h_new = h fac`, `fac = clamp(0.9 err^(-1/4), 1/6, 5)`, with `err` the
weighted root mean square of the last stage. Output and trigger search use a cubic Hermite interpolant, linear for
algebraic variables.

**BDF** (Shampine and Reichelt 1997, in the form of SciPy's implementation): numerical differentiation formulas of
order 1 to 5 on a quasi-constant step, with backward differences, `kappa = (0, -0.185, -1/9, -0.0823, -0.0415)`,
Newton iteration on `(M - c J) dy = c f - M (psi + d)`.

Default tolerances: relative `1e-9`, absolute `1e-12`.

### Events

SBML Level 3 semantics. A trigger fires on a false-to-true transition; `initialValue` sets its state before the
start. Delays and, with `useValuesFromTriggerTime`, assignment values are computed at trigger time. A non-persistent
pending event cancels when its trigger goes false. Due events run one at a time in order of priority, with ties
broken by a seeded random generator as the specification requires, and all triggers are re-evaluated after each.
Species values convert with the compartment size before the event, and a species keeps its amount when its
compartment's size jumps.

Each accepted step is scanned at 16 points for a trigger change, which catches a trigger that opens and closes
inside a step. The first change is refined by bisection on the interpolant to a relative width of `1e-11`, and the
step is repeated exactly up to it. A window narrower than a sixteenth of a step can still be missed.

### Steady states

The stoichiometric matrix factors as `N = L Nr`, with `Nr` the independent rows found by full pivoting and
`L = N Nr' (Nr Nr')^-1` (Reder 1988). Dependent species follow from the independent ones through the conserved totals.
The system integrates for a chosen time first, which selects the attractor reached from the initial state, and damped
Newton then solves `Nr v(y(s)) = 0`, running until it stops improving. This finds a fixed point and says nothing about
its stability.

### Metabolic control analysis

With unscaled elasticities `E = dv/dy` from symbolic differentiation (Reder 1988; Heinrich and Schuster 1996):

    Cs = -L (Nr E L)^-1 Nr              Cj = I + E Cs

The solve uses row and column equilibration and two steps of iterative refinement, since enzymes near equilibrium
make `Nr E L` badly scaled. Scaled coefficients are `C^S_ij = Cs_ij v_j / y_i`, `C^J_ij = Cj_ij v_j / v_i`. A flux under
`1e-9` of the largest counts as zero and has no scaled coefficients.

The engine reports residuals of the theorems in the scaled form in which they are stated, `sum_j C^J_ij = 1`,
`sum_j C^S_ij = 0`, `(y_k / v_i) (Cj E L)_ik = 0`, `(y_k / y_i) (Cs E L + L)_ik = 0`, beside the unscaled residuals of
`Cj E L` and `Cs E L + L`, and `Cj K = K`, `Cs K = 0` for the steady-state flux modes `K`.

### Sensitivity and estimation

**Local**: `dJ/dp = Cj dv/dp`, `dy/dp = Cs dv/dp` with `dv/dp` symbolic.

**Global**: first-order Sobol indices from the estimator of Saltelli et al. (2010) and total indices from Jansen's
(1999) formula, with `n (d + 2)` evaluations. A sample that fails to evaluate is dropped from every arm.

**Estimation**: Levenberg-Marquardt on `chi2(p) = sum ((y_model - y_obs) / sigma)^2` in the logarithm of each parameter,
with forward-difference Jacobians over full simulations. **Profile likelihood** (Raue et al. 2009): each grid value of one
parameter re-optimizes the others, and the interval ends where `chi2` exceeds its minimum by 3.84.

### Stochastic simulation

**Direct method** (Gillespie 1977), exact for mass-action networks with constant propensities between firings.
**Tau-leaping** (Cao, Gillespie and Petzold 2006): *an approximation*, with Poisson increments and a leap bound of the
form of their equation 8 with `g = 2` for all species, which is conservative for first-order reactions; it falls back to
direct steps when a leap would cover fewer than 10 expected firings, and halves a leap that would make an amount negative.
With rate rules or events, propensities are treated as constant over steps of at most `T/2000`.

## 5. Outputs

**Escher** (King et al. 2015): reaction data, and map files of schema 1-0-0 with a force-directed layout.

**Atlas overlay** (`rashevsky-flux-overlay/1`): projects absolute fluxes onto the Reactome events of the metabolism map
at `biochemistry.draconitic.com`. A reaction reaches an event through a shared Reactome or Rhea identifier, or, failing
those, through a shared EC number. An event reached by an identifier reports only those reactions; the `match` field of
each event gives the type. *This is a projection for reading fluxes in context*: the model and the map are separate
reconstructions, and an EC number can cover several events and compartments.

**Layouts**: `models/e_coli_core/layout.json` is drawn by hand. Other models use Fruchterman and Reingold (1991) with
gravity, over reactions and non-currency metabolites; this draws a readable network, and its geometry means nothing.

## References

- Agren R et al. (2013) The RAVEN toolbox. PLoS Comput Biol 9:e1002980.
- Beber ME et al. (2022) eQuilibrator 3.0. Nucleic Acids Res 50:D603.
- Bland RG (1977) New finite pivoting rules for the simplex method. Math Oper Res 2:103.
- Cao Y, Gillespie DT, Petzold LR (2006) Efficient step size selection for the tau-leaping simulation method. J Chem Phys 124:044109.
- Desouki AA et al. (2015) CycleFreeFlux. Bioinformatics 31:2159.
- Evans TW, Gillespie CS, Wilkinson DJ (2008) The SBML discrete stochastic models test suite. Bioinformatics 24:285.
- Fruchterman TMJ, Reingold EM (1991) Graph drawing by force-directed placement. Softw Pract Exp 21:1129.
- Gelman A et al. (2013) Bayesian Data Analysis, 3rd ed.
- Geweke J (1992) Evaluating the accuracy of sampling-based approaches. Bayesian Statistics 4.
- Gillespie DT (1977) Exact stochastic simulation of coupled chemical reactions. J Phys Chem 81:2340.
- Gudmundsson S, Thiele I (2010) Computationally efficient flux variability analysis. BMC Bioinformatics 11:489.
- Hairer E, Wanner G (1996) Solving Ordinary Differential Equations II, 2nd ed.
- Haraldsdottir HS et al. (2017) CHRR: coordinate hit-and-run with rounding. Bioinformatics 33:1741.
- Harris PMJ (1973) Pivot selection methods of the Devex LP code. Math Program 5:1.
- Heinrich R, Schuster S (1996) The Regulation of Cellular Systems.
- Henry CS, Broadbelt LJ, Hatzimanikatis V (2007) Thermodynamics-based metabolic flux analysis. Biophys J 92:1792.
- Huangfu Q, Hall JAJ (2018) Parallelizing the dual revised simplex method. Math Program Comput 10:119.
- Jansen MJW (1999) Analysis of variance designs for model output. Comput Phys Commun 117:35.
- Kendall DG (1948) On the generalized birth-and-death process. Ann Math Stat 19:1.
- King ZA et al. (2015) Escher. PLoS Comput Biol 11:e1004321.
- King ZA et al. (2016) BiGG Models. Nucleic Acids Res 44:D515.
- Lewis NE et al. (2010) Omic data from evolved E. coli are consistent with computed optimal growth. Mol Syst Biol 6:390.
- Lovasz L, Vempala S (2006) Simulated annealing in convex bodies and an O*(n^4) volume algorithm. J Comput Syst Sci 72:392.
- Mahadevan R, Edwards JS, Doyle FJ (2002) Dynamic flux balance analysis of diauxic growth. Biophys J 83:1331.
- Mahadevan R, Schilling CH (2003) The effects of alternate optimal solutions. Metab Eng 5:264.
- Megchelenbrink W, Huynen M, Marchiori E (2014) optGpSampler. PLoS ONE 9:e86587.
- Monk JM et al. (2017) iML1515, a knowledgebase that computes Escherichia coli traits. Nat Biotechnol 35:904.
- Noor E et al. (2013) Consistent estimation of Gibbs energy using component contributions. PLoS Comput Biol 9:e1003098.
- Orth JD, Fleming RMT, Palsson BO (2010) Reconstruction and use of microbial metabolic networks. EcoSal Plus 4.
- Orth JD, Thiele I, Palsson BO (2010) What is flux balance analysis? Nat Biotechnol 28:245.
- Raue A et al. (2009) Structural and practical identifiability analysis. Bioinformatics 25:1923.
- Reder C (1988) Metabolic control theory: a structural approach. J Theor Biol 135:175.
- Robinson JL et al. (2020) An atlas of human metabolism. Sci Signal 13:eaaz1482.
- Saltelli A et al. (2010) Variance based sensitivity analysis of model output. Comput Phys Commun 181:259.
- Sanchez BJ et al. (2017) Improving the phenotype predictions of a yeast genome-scale metabolic model. Mol Syst Biol 13:935.
- Schellenberger J, Lewis NE, Palsson BO (2011) Elimination of thermodynamically infeasible loops. Biophys J 100:544.
- Segre D, Vitkup D, Church GM (2002) Analysis of optimality in natural and perturbed metabolic networks. PNAS 99:15112.
- Shampine LF, Reichelt MW (1997) The MATLAB ODE suite. SIAM J Sci Comput 18:1.
- Shlomi T, Berkman O, Ruppin E (2005) Regulatory on/off minimization. PNAS 102:7695.
- Teusink B et al. (2000) Can yeast glycolysis be understood in terms of in vitro kinetics? Eur J Biochem 267:5313.
- Thiele I, Palsson BO (2010) A protocol for generating a high-quality genome-scale metabolic reconstruction. Nat Protoc 5:93.
- Varma A, Palsson BO (1994) Stoichiometric flux balance models quantitatively predict growth and metabolic by-product secretion. Appl Environ Microbiol 60:3724.
