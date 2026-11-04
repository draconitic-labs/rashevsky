# Model and methods

Every governing equation, numerical method and data source in the engine, with references. Validation
numbers are in [VALIDATION.md](VALIDATION.md).

Units throughout the constraint-based part are mmol/gDW/h for fluxes and 1/h for growth, as in the
models. Kinetic models carry the units their SBML defines.

## 1. Stoichiometric models

A model has metabolites `i`, reactions `j`, a stoichiometric matrix `S` with `S_ij < 0` for substrates,
bounds `l <= v <= u`, and gene-protein-reaction rules. SBML species with `boundaryCondition="true"` are
exempt from the balance. A reaction with one balanced metabolite is an exchange, whether written `A <->` or through
a boundary species as `A_b <-> A`, and its uptake direction comes from the coefficient on that metabolite.

Reading: SBML Levels 2 and 3 with FBC versions 1 to 3 (`fbc:fluxBound`, including the operations `less` and `greater`;
`fbc:geneAssociation`; `fbc:lowerFluxBound`; `fbc:geneProductAssociation`), the groups package for subsystems, and the
COBRA conventions before FBC (`LOWER_BOUND`, `UPPER_BOUND` and `OBJECTIVE_COEFFICIENT` among kinetic law parameters,
`GENE_ASSOCIATION` in notes). Stoichiometry expressions, initial assignments and rules that define flux bounds or
stoichiometries are evaluated by the kinetic compiler of section 4. FBC version 3 user constraints and quadratic objective
terms are refused for flux analyses, as is any model using a package it marks as required, other than FBC and groups.
Identifiers are kept as written: a COBRA escape such as `__45__` is not decoded. Namespace bindings last as long as their
element, and a package can arrive under a default namespace.

Writing: SBML Level 3, declaring FBC version 2 and groups only when a model uses them, with annotations and constraints kept.
Version 1 unless the math needs Version 2 (`max`, `min`, `quotient`, `rem`, `implies`, `rateOf`). Level 2 defaults become explicit,
and a Level 2 stoichiometry expression becomes an assignment rule on a named species reference. Content that cannot be written
faithfully, FBC version 3 constraints or a gene rule that does not parse, is an error rather than a different model.

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

Rows that multiply a flux bound by a binary variable (ll-FBA, ROOM, TFA) need finite bounds. An infinite bound is replaced
by the largest flux the network allows, found by linear programming, and a reaction with no such limit is an error.

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

**Retaining an objective**: with `z*` the optimum and `f` the retained fraction, a maximized objective must reach
`min(f z*, z* / f)` and a minimized one must stay under `max(f z*, z* / f)`. For a positive optimum this is the usual
`f z*` and `z* / f`, and it also holds for a negative one. At `f = 0` no row is added.

**Parsimonious FBA** (Lewis et al. 2010): minimize `sum (v+_j + v-_j)` subject to `S (v+ - v-) = 0` and the retained
objective, on separate forward and reverse parts of each flux. The solution's `objective` is this total, as in COBRApy.

**Flux variability** (Mahadevan and Schilling 2003): minimum and maximum of each `v_j` subject to
`S v = 0`, `l <= v <= u` and the retained objective. No slack is added to that last row: some fluxes respond ~1000-fold
to the objective (`ICDHyr` in `e_coli_core` moves 920 units per unit of growth), so 1e-9 of slack is
visible at 1e-6. LPs share a session per worker (Gudmundsson and Thiele 2010).

**Loopless FBA** (Schellenberger, Lewis and Palsson 2011): for each internal reaction a binary `a_j` and a
potential `G_j`, with `-M(1 - a_j)`-style rows `v_j <= u_j a_j`, `v_j >= l_j (1 - a_j)`, `G_j <= -1 + M(1 - a_j)`,
`G_j >= 1 - (M + 1) a_j`, and `N_int' G = 0` over a null-space basis of the internal stoichiometry.

**CycleFreeFlux** (Desouki et al. 2015): given `v`, minimize `sum |w_j|` over internal reactions subject to
`S w = 0`, `w_j` between 0 and `v_j` within its own bounds, and boundary and objective reactions fixed at their values
in `v`.

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
`dt` the fluxes are constant. With `c_i` the exchange's coefficient on its metabolite, uptake is `u_i = c_i v_i`,
limited by `vmax_i S_i / (K_i + S_i)` and never beyond the model's own bounds, except that what is there to take up overrides a
bound that forces uptake. The growth rate `mu` is the flux of the one maximized reaction in the objective. Then
`X(t + dt) = X e^(mu dt)` and `S_i(t + dt) = S_i - u_i I`, with `I = X (e^(mu dt) - 1) / mu`. Every followed species, kinetic or
not, is capped at `u_i I <= S_i`. Since `I` grows with `mu` and `mu` with the cap, the cap that uses up exactly what remains is
found by bisection, which holds biomass to the yield at any step size. The scheme is exact while fluxes stay constant and
first-order once they change. Each step reports its solver status, and an unbounded problem is an error.

### Flux sampling

Reactions whose FVA range is under 1e-9 are fixed, the remaining fluxes are written `v = x0 + N z` over an
orthonormal null-space basis, from the same elimination that selects independent rows, and FVA ranges serve as bounds.
Keeping an objective fraction adds its half-space, for either sense. With no fraction, the objective is free.

**CHRR** (Haraldsdottir et al. 2017): coordinate hit-and-run from the centre of the largest inscribed ball,
found by a linear program, along columns of a rounded basis. *This implementation rounds differently from
the original*: it moves to isotropic position with `T T' = C` for the covariance `C` of samples, over 5 rounds
of `60 d` points (Lovasz and Vempala 2006), rather than using the maximum-volume ellipsoid. Rounding affects
only mixing: the target stays uniform.

**OptGP** (Megchelenbrink, Huynen and Marchiori 2014): chains of artificial-centering hit-and-run from the
centre of FVA-derived warmup points, with projection onto the null space every 100 moves. Moves run along warmup
directions, so if those vertices do not span the polytope, vertices of random objectives are added, and it is an error
if they still do not. *A heuristic sampler*: since artificial centering adapts its directions, it has no exact guarantee of uniformity, so read
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

GECKO-style (Sanchez et al. 2017). Each catalysis `c` (one reaction by one isozyme or complex) has an allocation `a_c` in
mmol/gDW. Enzyme `k` has a usage `e_k`. With turnover numbers in 1/h, molecular weights in g/mmol and `n_ck` copies of
enzyme `k` in the complex of `c`:

    v_j <= sum_c kcat+_c a_c,   v_j >= -sum_c kcat-_c a_c,   sum_c n_ck a_c <= e_k,   sum_k MW_k e_k <= P,   e_k <= abundance_k

Reactions therefore share a promiscuous enzyme, and a complex takes all of its subunits, as in the paper's single
utilization constraint per enzyme. The pool is a limit and not a use, so a second program finds the least protein that reaches the
optimum: `usage` and `pool_used` report that, and it sets the fluxes returned. A published GECKO file can also be read as it stands, since it encodes enzymes as
metabolites and reactions, with coefficients `n_ck / kcat`. Both routes reach the same optimum on `eciML1515`
(`validation/v12_ecmodel.py`).

Turnover numbers, weights and the pool `P` are inputs, and the module holds no biological data.

### Metabolic tasks

The format and semantics of RAVEN's `checkTasks` (Agren et al. 2013; Robinson et al. 2020). A task passes when
a steady state exists with its inputs and outputs as the system's whole boundary: an input lets a metabolite's
balance row reach `S v` in `[-UB, -LB]`, an output in `[LB, UB]`. Extra reactions, bound changes and
`ALLMETSIN[c]` tokens are supported.

## 4. Kinetic models

### Compilation

SBML core Levels 2 and 3. Function definitions are inlined, and symbols resolved to slots in a value array. A local
parameter gets a slot named `reaction.parameter`, so that it can be overridden, fitted and analysed like a global one.
MathML and SBML core are recognized under any namespace prefix. Each species has two slots, the value math sees (a concentration, or
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
definitions, evaluated repeatedly until they stop changing, which resolves them in any order. Algebraic rules alternate
with them, since each can read what the other determines.

**Overrides** (`--set`, `set_value`): rebuild the initial values from the compiled ones, with parameters and sizes first and
species values after, so that the order of overrides does not matter. A symbol that a rule determines cannot be overridden.

**`delay(x, tau)`**: from a stored history of accepted steps, with cubic Hermite interpolation of the state.
Before the start time, the initial definitions are re-evaluated at that earlier time and rates are zero, which
is how the SBML Test Suite treats history before the start. Steps are limited to the smallest constant lag.
This is a method of steps with interpolation, and it is not designed for state-dependent lag discontinuities.

**Fast reactions** (removed from SBML in Level 3 Version 2): *a stand-in*, the singular perturbation of
multiplying the law by `1e8`, with the fast subsystem brought to equilibrium in pseudo-time at the start and after
events. A law whose largest partial derivative at the start is more than 1000 times from one is rescaled so that
this derivative becomes `1e8`, which makes the method independent of the law's size. The residual of the equilibrium
condition is of the order of `1e-8`.

**Jacobians** come from symbolic differentiation (`kinetic/diff.rs`) of the right-hand side with all derived values
inlined; `piecewise` differentiates branch by branch and step-like operators to zero. A model with `delay` uses forward
differences, since the symbolic form leaves `delay` out and would cost a Rosenbrock method its order. So does a model past
a size limit of 4 million nodes.

### Integrators

**Rodas4** (Hairer and Wanner 1996, coefficients from their `rodas.f`): L-stable, stiffly accurate Rosenbrock
method of order 4 with an embedded order-3 estimate, for `M y' = f(t, y)` with a diagonal `M` and index-1 algebraic
equations. Step control follows `rodas.f`: `h_new = h fac`, `fac = clamp(0.9 err^(-1/4), 1/6, 5)`, with `err` the
weighted root mean square of the last stage. Output and trigger search use a cubic Hermite interpolant, linear for
algebraic variables.

**BDF** (Shampine and Reichelt 1997, in the form of SciPy's implementation): numerical differentiation formulas of
order 1 to 5 on a quasi-constant step, with backward differences, `kappa = (0, -0.185, -1/9, -0.0823, -0.0415)`,
Newton iteration on `(M - c J) dy = c f - M (psi + d)`.

Default tolerances: relative `1e-9`, absolute `1e-12`. The dynamics start at time zero whatever the first output time. Far
from time zero, `t + h - t` can differ measurably from `h`, so the BDF difference table is rescaled to the step actually taken.

**Jumps in rates**: a comparison or rounding of an expression of the time alone marks where the right-hand side can jump,
whatever its form: `time > c`, `2 * time > c`, `lt(a, time, b)`, a rule `T = time`, `abs(time - c) < w`, `floor(time / p)`. A
comparison affine in the time by its form has its root solved, wherever that lies. Anything else is examined on a grid of
160,000 points over the run, doubled until the number of jumps stops changing, and refined by bisection; a run that never
settles on a number is an error, and a single window narrower than the finest grid can still be missed. A state that only events
change, such as the time of the next dose, holds still between them, so these times are computed from current values and again
after each event. Integration stops `4e-14 (1 + t)` before each such time, moves across without a step, and restarts just after
it, so that no stage of a step straddles the jump.

Discontinuities that depend on the state are left to step-size control. A step cannot shrink below `4e-14 t`, the limit of what
the arithmetic resolves, and at that size it is taken even over its error limit, since a jump can only be stepped across. Two
hundred such steps in a row, as with a solution that blows up, are an error.

### Events

SBML Level 3 semantics. A trigger fires on a false-to-true transition; `initialValue` sets its state before the
start. Delays and, with `useValuesFromTriggerTime`, assignment values are computed at trigger time. A non-persistent
pending event cancels when its trigger goes false. Due events run one at a time in order of priority, with ties
broken by a seeded random generator as the specification requires, and all triggers are re-evaluated after each.
Species values convert with the compartment size before the event, and a species keeps its amount when its
compartment's size jumps.

A step ends on every time found this way in a trigger, so a window set by the time cannot fall between the points examined. Each accepted step is also scanned at 16 points for a trigger change, which catches a state-driven
trigger that opens and closes inside a step. The first change is refined by bisection on the interpolant to a relative
width of `1e-11`, and the step is repeated exactly up to it. A state-driven window narrower than a sixteenth of a step can
still be missed. A negative delay is an error. `delay` history keeps the value on both sides of an event.

### Steady states

The stoichiometric matrix factors as `N = L Nr`, with `Nr` the independent rows found by full pivoting and
`L = N Nr' (Nr Nr')^-1` (Reder 1988). Dependent species follow from the independent ones through the conserved totals.
The system integrates for a chosen time first, extended past the last time a rate or trigger can jump and past any event still
pending. If events are still firing in the second half of the run, or rates still change with the time, it extends tenfold, up to
twice, and then there is no steady state, which is an error. Rates count as settled when no jump lies ahead on a long grid and they
read the same at times spread over decades, relative to their own size. Damped Newton then solves `Nr v(y(s)) = 0`, running until it
stops improving. A last Newton step that is not small means the iteration is running away, as with `S' = 1 / (1 + S)`, where the
residual falls but no root exists, and that is an error. An event that would fire between the integrated state and the root is an
error too.

The result reports the largest real part of the eigenvalues of the Jacobian over the independent amounts. Newton can reach any
root, so a root that is unstable or has a negative amount extends the integration tenfold, up to twice, before it is reported: an
unstable root is returned with `stable = false`, as on a limit cycle, and one with negative amounts is an error. `delay` reads the
present at a steady state, where those eigenvalues do not decide stability, so a delayed model reports no eigenvalue and counts as
stable only if the integration reached the root, relative to the size of its state. A stoichiometric coefficient or conversion
factor that varies is an error.

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

**Local**: a parameter moves the rates of change at the frozen steady state by `b = N dv/dp + (dN/dp) v`, and it can move the
start. With `d0 = dy0/dp - L dy0_ind/dp`, the part of that move the dynamics cannot undo, such as a conserved total:
`ds/dp = -(Nr E L)^-1 (b + Nr E d0)`, `dy/dp = d0 + L ds/dp`, `dJ/dp = dv/dp + E dy/dp`. The partial derivatives come from a
fourth-order central difference with `h = 1e-3 |p|` through the override mechanism, so that a parameter counts everywhere it
reaches: laws, stoichiometries, compartment sizes, starting values and values assigned from it. In a model with events, or for a
parameter that changes the stoichiometry and with it what is conserved, the conserved totals do not follow from the start alone,
so the sensitivity is a central difference of steady states solved again. A name that does not exist, varies during a simulation or
is set by a rule is an error.

**Global**: first-order Sobol indices from the estimator of Saltelli et al. (2010) and total indices from Jansen's
(1999) formula, with `n (d + 2)` evaluations. The output is a symbol at the steady state or at a time. A sample that fails to
evaluate, or whose steady state is not an attractor, is dropped from every arm. Failures concentrated in one region bias the
indices, so more than 5 % is an error. An output whose variance is of the size of rounding has indices of zero.

**Estimation**: Levenberg-Marquardt on `chi2(p) = sum ((y_model - y_obs) / sigma)^2` in the logarithm of each parameter,
with forward-difference Jacobians over full simulations. **Profile likelihood** (Raue et al. 2009): each grid value of one
parameter re-optimizes the others, and the interval ends where `chi2` exceeds its minimum by 3.84. The grid only brackets each
limit, which bisection on the re-optimized profile then finds, since a sharp likelihood is far from linear across a grid cell.
The walk starts at the lowest point on the profile, which is reported, so that a fit stuck in a local minimum shows. A limit
beyond a failed simulation stays undetermined.

### Stochastic simulation

**Direct method** (Gillespie 1977), exact for networks with constant propensities between firings. A reaction whose firing
would leave a species below zero has propensity zero. Stoichiometry is read at each firing. Event semantics and the handling of
compartment jumps follow the deterministic driver, with priorities and random ties.

**Tau-leaping** (Cao, Gillespie and Petzold 2006): *an approximation*. Reactions within 10 firings of exhausting a reactant are
critical. Over the others, the leap is the `tau'` of their equation 33, over every species a propensity reads, catalysts and
autocatalysts included, with `g_i` from the highest-order reaction reading species `i` and the order from the reactant
stoichiometry, which is the true order for mass-action laws. Below 10 expected firings it
takes 100 direct steps. Otherwise the leap is the smaller of `tau'` and the exponential waiting time of the critical reactions,
clipped to the next stop, with Poisson firings of the others. A leap that makes a species negative is repeated at half of `tau'`,
and the time advances only on acceptance.

The method bounds the relative change of each propensity by `epsilon`, except that one molecule's change is always allowed
(their `max(epsilon x_i / g_i, 1)`). The explicit method's first-order error is therefore of the order of `epsilon` per decay time
(3 % at 0.03 on pure decay), and a low-copy species with fast turnover is biased whatever `epsilon` is: a birth-death species of true
mean 1, beside a large population that sets the leap, reads 1.45 here and 1.81 in GillesPy2's tau-leaping solver on the same system.
Use the direct method for such species.

A propensity may jump at the times of the integrators' analysis above, where a step ends, and is read just past each one, since
`time > 5` is still false at 5. Those times follow the current values, so repeated doses all fire. A propensity that depends on a
rate rule is held constant over steps of at most `T/2000`, and at most a quarter of the rule's own time scale, which is an
approximation. These are errors: negative or undefined propensities, rates that vary continuously with the time, algebraic rules,
`delay`, fast reactions, cascading events that do not settle, amounts or reactant coefficients that are not whole numbers, and runs
whose amounts pass `1e18` or that exceed `1e8` steps.

## 5. Outputs

**Escher** (King et al. 2015): reaction data, and map files of schema 1-0-0 with a force-directed layout.

**Atlas overlay** (`rashevsky-flux-overlay/1`): projects absolute fluxes onto the Reactome events of the metabolism map
at `biochemistry.draconitic.com`. A reaction reaches an event through a shared Reactome (`reactome` or `reactome.reaction`) or
Rhea identifier, or, failing those, through a shared EC number. An event reached by an identifier reports only those reactions,
unless they carry no flux and its EC matches do; the `match` field of each event gives the type. A drawn reaction takes the
value and the line style of its largest event. *This is a projection for reading fluxes in context*: the model and the map are separate
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
