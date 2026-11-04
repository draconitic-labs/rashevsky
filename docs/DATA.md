# Data and licenses

The engine is under the MIT license. Everything in this repository can be used for any purpose, commercial included, under the terms below. Models fetched from BiGG Models stay out of it, because BiGG's terms do not allow that.

## In this repository

| File | Source | Terms |
| --- | --- | --- |
| `models/ich360/iCH360.xml.gz` | Corrao et al. 2025, repository `marco-corrao/iCH360` | MIT. See its `NOTICE`. |
| `models/ich360/layout.json` | Converted from the authors' Escher map by `scripts/make_ich360_layout.py` | MIT. See the same `NOTICE`. |
| `models/teusink/BIOMD0000000064.xml` | BioModels; Teusink et al. 2000 | CC0 1.0 |
| `data/thermo/ich360.json` | Computed from eQuilibrator | CC BY 4.0 |
| `data/biochemistry/*` | Atlas of Human Biochemistry, Apache 2.0; science from Reactome release 97, whose data is CC0 | See its `NOTICE`. |
| `desktop/ui/fonts/*` | Public Sans, JetBrains Mono | SIL Open Font License 1.1 |

iCH360 makes the desktop application and the WebAssembly demo work on first start. Its license asks for its notice in every copy, so each build carries it: the desktop binary shows it under Models and its packages install it as a file, and the demo serves `NOTICE-iCH360.txt` with a link from the page.

iCH360's authors extracted its reactions from iML1515, which BiGG Models distributes, and curated them against the *E. coli* core model, EcoCyc and KEGG. They release the result under MIT, and that grant is what this repository relies on.

## Fetched by `scripts/fetch_data.py`

`validation/data-manifest.json` lists each file with its source and its SHA-256 digest. A digest that
differs is reported, since a changed upstream file changes the results.

| File | Source | Terms |
| --- | --- | --- |
| `e_coli_core.xml.gz` | BiGG Models; Orth et al. 2010 | BiGG's terms, below |
| `iML1515.xml.gz` | BiGG Models; Monk et al. 2017 | BiGG's terms, below |
| `Recon3D.xml.gz` | BiGG Models; Brunk et al. 2018 | BiGG's terms, below |
| `iCH360_full_map.json` | Corrao et al. 2025, at commit 08dd438 | MIT. The source of `models/ich360/layout.json`. |
| `ecmodels/eciML1515_batch.xml` | SysBioChalmers/ecModels at commit 8c50008; Sanchez et al. 2017, Domenzain et al. 2022 | CC BY 4.0 |
| `Human-GEM.xml`, `metabolicTasks_*.txt` | SysBioChalmers/Human-GEM v2.1.1; Robinson et al. 2020 | CC BY 4.0 |
| `semantic_tests.zip`, `stochastic_tests.zip` | SBML Test Suite 3.5.0 | LGPL |
| `monk_2017_data11.xlsx` | Supplementary Data 11 of Monk et al. 2017 | Journal supplement; read at validation time, not redistributed |
| `BIOMD*.xml` | BioModels | CC0 1.0 |

The desktop application downloads iML1515, iJO1366, Recon3D and Human-GEM only when asked, into
`~/.local/share/rashevsky/models`.

## BiGG's terms

BiGG Models allows use, copying and distribution "for educational, research and non-profit purposes, without fee", and commercial use needs a separate agreement with the University of California: <http://bigg.ucsd.edu/license>. `e_coli_core`, iML1515, iJO1366 and Recon3D come under these terms. No file in this repository derives from them. The desktop application and `scripts/fetch_data.py` download them from BiGG on request, so anyone who does so takes them under BiGG's terms.
