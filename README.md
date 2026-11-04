# Rashevsky

Genome-scale constraint-based modelling and detailed kinetic modelling of metabolism, in one engine.

Which gene knockouts raise product yield? Which enzymes control a pathway's flux? What can a cell do on given nutrients, within thermodynamics? Rashevsky answers these from SBML models through a Rust library, a command line, Python bindings and a desktop workbench.

![The flux map of E. coli central metabolism](docs/screenshots/map.png)

Named after Nicolas Rashevsky (1899-1972), founder of mathematical biophysics.

## What it does

**Constraint-based:** flux balance analysis, parsimonious FBA, flux variability, loopless FBA and CycleFreeFlux, MOMA and ROOM, single and double gene or reaction deletions, synthetic lethal pairs, production envelopes, minimal media, dynamic FBA, flux sampling (CHRR, OptGP) with convergence diagnostics, thermodynamics-based flux analysis, enzyme constraints in the style of GECKO, and metabolic task checks in the format of RAVEN and Human-GEM.

**Kinetic:** SBML core Levels 2 and 3 with rules, events, delays, conversion factors and `rateOf`; stiff integration with Rodas4 and variable-order BDF on analytic Jacobians; steady states with moiety conservation; metabolic control analysis; local and global sensitivity; parameter estimation with profile likelihood; stochastic simulation (Gillespie, tau-leaping).

**Formats:** SBML Level 3 with FBC and groups, read and written. Escher maps and flux data. Flux overlays for the metabolism map at [biochemistry.draconitic.com](https://biochemistry.draconitic.com).

## How far to trust it

Each row is a named published result or another tool on the same input, with the commands to reproduce it in [docs/VALIDATION.md](docs/VALIDATION.md).

| Check | Reference | Result |
| --- | --- | --- |
| Bundled iCH360 | COBRApy on 15 analyses, among them all 63,190 gene pairs | all within tolerance, growth 0.865627 |
| `e_coli_core` growth | 0.8739 1/h, BiGG and COBRApy | 0.873922, equal to COBRApy within 1e-11 |
| Flux variability | COBRApy, tolerance 1e-6 on `e_coli_core`, 1e-5 on iML1515 | largest difference 4e-12 and 1.8e-6 |
| iML1515 gene essentiality | Keio screens, Monk et al. 2017: 93.4 % | 93.48 % of 23,617 phenotypes; 1 row differs from the paper, traced to the model file |
| Human-GEM | 57 essential tasks, Robinson et al. 2020 | 57 of 57, same verdicts as RAVEN |
| Recon3D | COBRApy | same objective within 1e-15, all 2,248 gene deletions within 5e-14 relative |
| Enzyme constraints | Published `eciML1515`, 7 conditions | within 2e-14 of an optimum certified by a dual bound |
| Yeast glycolysis | Teusink et al. 2000, libRoadRunner | same steady state within 1e-14 |
| Control coefficients | COPASI | largest difference 4e-7 for flux, 2e-6 for concentration |
| Control theorems | exact | hold to 6e-12 on five published models |
| SBML semantic suite | SBML Test Suite 3.5.0 | 1,700 of 1,700 applicable cases, with either integrator, and after writing and re-reading each model |
| Stochastic suite | DSMTS | all 39 cases; case 00003 is decided by its exact law, because its kurtosis of 96 makes the suite's variance test depend on the random stream |

## Install

The desktop application, on Linux with WebKitGTK 4.1:

```
sudo dnf install webkit2gtk4.1-devel gtk3-devel cmake clang-libs   # or the equivalents for your system
desktop/scripts/install.sh
```

This adds **Rashevsky** to the applications menu and the desktop. Drop an SBML file on the icon to open it; its right-click menu goes straight to a knockout screen, yeast glycolysis, the human atlas or the validation results.

The command line:

```
cargo install --path crates/rashevsky-cli
```

Python 3.9 or later:

```
pip install maturin && maturin develop --release
```

The build compiles HiGHS from source, and it needs CMake, a C++ compiler and libclang.

## Use

```
rashevsky fba models/ich360/iCH360.xml.gz
rashevsky fba model.xml --ko-gene b4025 --bound EX_o2_e=0:1000
rashevsky fva iML1515.xml.gz --fraction 0.9
rashevsky deletion iML1515.xml.gz --kind gene --order 1
rashevsky envelope model.xml --control EX_o2_e --target EX_ac_e
rashevsky sample model.xml --method chrr --fraction 0.9
rashevsky tasks Human-GEM.xml metabolicTasks_Essential.txt
rashevsky simulate models/teusink/BIOMD0000000064.xml --end 20 --method bdf
rashevsky mca models/teusink/BIOMD0000000064.xml
rashevsky sobol models/teusink/BIOMD0000000064.xml --range KmPFKF6P=0.05:0.2:log --range gR=2:10 --output GLK
rashevsky fit model.xml --start k1=0.5 --data measured.csv --profile k1
rashevsky --help
```

Every command prints JSON, or writes it to the file given by `--out`.

```python
import rashevsky as rv

model = rv.read_sbml("models/ich360/iCH360.xml.gz")
model.fba()["objective"]                                  # 0.8656
model.fba(knockout_genes=["b4025"])["objective"]          # 0.8632
model.fva(fraction=0.9, reactions=["PFK"])                # {"PFK": (low, high)}
model.single_gene_deletion()                              # {gene: growth}

glycolysis = rv.KineticModel.read_sbml("models/teusink/BIOMD0000000064.xml")
glycolysis.simulate(end=20, method="bdf", symbols=["G6P", "ATP"])
glycolysis.mca()["flux_control"]
glycolysis.local_sensitivity(["gR"])["flux"]
glycolysis.sobol([{"name": "gR", "low": 2, "high": 10}], "vGLK", samples=512)["total"]
```

[notebooks/knockout_screen.ipynb](notebooks/knockout_screen.ipynb) ranks every single and double deletion for succinate production.

```rust
use rashevsky::{Model, cb};

let model = Model::read_sbml("models/ich360/iCH360.xml.gz")?;
let growth = cb::fba(&model)?.objective;
```

## Desktop workbench

| | |
| --- | --- |
| ![Deleting a reaction reroutes flux](docs/screenshots/knockout.png) | ![Knockout screen](docs/screenshots/screen.png) |
| ![Control coefficients of yeast glycolysis](docs/screenshots/control.png) | ![Human-GEM on the human metabolism atlas](docs/screenshots/atlas.png) |

Keys: `1` to `8` change workspace, `/` finds a reaction, `k` deletes the selected reaction, `0` fits the map, `Esc` closes the panel. Double-click a reaction on the map to delete it.

## Repository

| Path | Contents |
| --- | --- |
| `crates/rashevsky` | The engine |
| `crates/rashevsky-cli` | The `rashevsky` command |
| `crates/rashevsky-py` | Python bindings |
| `crates/rashevsky-wasm`, `web/` | The WebAssembly demo |
| `desktop/` | The Tauri workbench |
| `models/` | Bundled models |
| `validation/` | Scripts that compare the engine with its references, and their results |
| `docs/MODEL.md` | Every equation and method, with references |
| `docs/VALIDATION.md` | Every validation number, with its command |

## Methods and data

Please cite the method papers for what you use: `docs/MODEL.md` lists them. If you use the bundled or fetched models, cite them as well.

The engine is under the [MIT license](LICENSE), and so is iCH360, the bundled model. Everything in this repository can be used for any purpose, commercial included, under the terms listed in [docs/DATA.md](docs/DATA.md). Models you download from BiGG Models, such as `e_coli_core` and iML1515, come under BiGG's terms: educational, research and non-profit use, with commercial use needing an agreement with UC San Diego.
