# Data and licenses

The engine is under the MIT license. Models and reference data keep their own terms.

## In this repository

| File | Source | Terms |
| --- | --- | --- |
| `models/e_coli_core/e_coli_core.xml.gz` | BiGG Models; Orth et al. 2010 | Educational, research and non-profit use; commercial use needs a license from UC San Diego. See its `NOTICE`. |
| `models/e_coli_core/layout.json` | Original to this project | MIT |
| `models/teusink/BIOMD0000000064.xml` | BioModels; Teusink et al. 2000 | CC0 1.0 |
| `data/thermo/e_coli_core.json` | Computed from eQuilibrator | CC BY 4.0 |
| `data/biochemistry/*` | Atlas of Human Biochemistry; Reactome release 97 | Apache 2.0 with attribution. See its `NOTICE`. |
| `desktop/ui/fonts/*` | Public Sans, JetBrains Mono | SIL Open Font License 1.1 |

The bundled `e_coli_core` makes the desktop application and the WebAssembly demo work on first start.
If you use them commercially, replace it or obtain the license above.

## Fetched by `scripts/fetch_data.py`

`validation/data-manifest.json` lists each file with its source and its SHA-256 digest. A digest that
differs is reported, since a changed upstream file changes the results.

| File | Source | Terms |
| --- | --- | --- |
| `iML1515.xml.gz` | BiGG Models; Monk et al. 2017 | As `e_coli_core` above |
| `Human-GEM.xml`, `metabolicTasks_*.txt` | SysBioChalmers/Human-GEM v2.1.1; Robinson et al. 2020 | CC BY 4.0 |
| `semantic_tests.zip`, `stochastic_tests.zip` | SBML Test Suite 3.5.0 | LGPL |
| `monk_2017_data11.xlsx` | Supplementary Data 11 of Monk et al. 2017 | Journal supplement; read at validation time, not redistributed |
| `BIOMD*.xml` | BioModels | CC0 1.0 |

The desktop application downloads iML1515, iJO1366, Recon3D and Human-GEM only when asked, into
`~/.local/share/rashevsky/models`.
