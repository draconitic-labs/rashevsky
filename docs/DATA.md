# Data and licenses

The engine's code is under the MIT license, and it opens any model a user gives it. The data in this repository keeps the terms listed here. One entry carries an open question, which only its rights holder can close: see [iCH360 and BiGG](#ich360-and-bigg).

## In this repository

| File | Source | Terms |
| --- | --- | --- |
| `models/ich360/iCH360.xml.gz` | Corrao et al. 2025, repository `marco-corrao/iCH360` | MIT from its authors, with the question below. See its `NOTICE`. |
| `models/ich360/layout.json` | Converted from the authors' Escher map by `scripts/make_ich360_layout.py` | As the model. See the same `NOTICE`. |
| `models/teusink/BIOMD0000000064.xml` | BioModels; Teusink et al. 2000 | CC0 1.0 |
| `validation/fixtures/biomodels/BIOMD*.xml` | BioModels; pinned validation references | CC0 1.0 |
| `data/thermo/ich360.json` | Computed from eQuilibrator | CC BY 4.0 |
| `data/biochemistry/*` | Atlas of Human Biochemistry, Apache 2.0; science from Reactome release 97, whose data is CC0 | See its `NOTICE`. |
| `desktop/ui/fonts/*` | Public Sans, JetBrains Mono | SIL Open Font License 1.1 |
| `validation/results/ich360_escher_map.json`, `ich360_reaction_data.json` | Written by the engine from iCH360 | As the model |
| `validation/results/human_gem_overlay.svg` | Written by the engine from Human-GEM and the Atlas scene | CC BY 4.0 for Human-GEM; `data/biochemistry/NOTICE` for the scene |
| Other files in `validation/results` | Numbers measured by this project | MIT |

iCH360 makes the desktop application and the WebAssembly demo work on first start. Its notice has to travel with every copy, so each build carries it: the desktop binary shows it under Models and its packages install it as a file, and the demo serves `NOTICE-iCH360.txt` with a link from the page. The licenses of linked software are in `THIRD-PARTY-NOTICES.md`, shipped the same way.

## Fetched by `scripts/fetch_data.py`

`validation/data-manifest.json` lists each file with its source and its SHA-256 digest. A file whose digest differs stops the script, since a changed upstream file changes the results.

The five CC0 BioModels references are copied from the bundled paths in that
manifest and checked against the same pinned digests. This keeps validation
independent of BioModels blocking requests from hosted CI runners. Teusink reuses
`models/teusink/BIOMD0000000064.xml`; the other four files are under
`validation/fixtures/biomodels/`. Other reference data is still fetched under the
terms below.

| File | Source | Terms |
| --- | --- | --- |
| `e_coli_core.xml.gz` | BiGG Models; Orth et al. 2010 | BiGG's terms, below |
| `iML1515.xml.gz` | BiGG Models; Monk et al. 2017 | BiGG's terms, below |
| `Recon3D.xml.gz` | BiGG Models; Brunk et al. 2018 | BiGG's terms, below |
| `iCH360_full_map.json` | Corrao et al. 2025, at commit 08dd438 | As the model. The source of `models/ich360/layout.json`. |
| `ecmodels/eciML1515_batch.xml` | SysBioChalmers/ecModels at commit 8c50008; Sanchez et al. 2017, Domenzain et al. 2022 | CC BY 4.0 |
| `Human-GEM.xml`, `metabolicTasks_*.txt` | SysBioChalmers/Human-GEM v2.1.1; Robinson et al. 2020 | CC BY 4.0 |
| `semantic_tests.zip`, `stochastic_tests.zip` | SBML Test Suite 3.5.0 | LGPL |
| `monk_2017_data11.xlsx` | Supplementary Data 11 of Monk et al. 2017 | Journal supplement; read at validation time, not redistributed |
| `BIOMD*.xml` | BioModels | CC0 1.0 |

The desktop application downloads `e_coli_core`, iML1515, iJO1366, Recon3D and Human-GEM only when asked, over TLS, into `~/.local/share/rashevsky/models`, and refuses a file whose digest differs from the one built into it.

## BiGG's terms

BiGG Models allows use, copying and distribution of "any part of BiGG Models" "for educational, research and non-profit purposes, without fee", provided its notice appears in all copies, and sends commercial use to the University of California, San Diego: <https://bigg.ucsd.edu/license>. `e_coli_core`, iML1515, iJO1366 and Recon3D come under these terms. BiGG now lives at the Technical University of Denmark, where the content inherited from San Diego "remains under the BiGG Models license": <https://bigg.bio/about>. The desktop application and `scripts/fetch_data.py` download these models from BiGG on request, so anyone who does so takes them under BiGG's terms.

## iCH360 and BiGG

No file in this repository is a copy of a BiGG model. iCH360 does descend from one. Its authors extracted its reactions from iML1515 and curated them against the *E. coli* core model, EcoCyc and KEGG: 342 of its 349 reactions match iML1515 in identifier and stoichiometry, 6 share the identifier with a different stoichiometry, and the last is its own biomass reaction. They release the result under MIT.

Their grant cannot give more than they received. Whether BiGG's terms reach a new model that restates reaction stoichiometry has not been settled by anyone here. Reactions are facts of chemistry, and facts are not copyrightable, but a compilation of them can be (*Feist v. Rural*, 499 U.S. 340, 1991). What follows from that:

- Educational, research and non-profit use is covered on either reading. The `NOTICE` of the model reproduces BiGG's terms, which is the condition they set.
- Anyone building iCH360 into a commercial product should first get an answer in writing from UC San Diego (`invent@ucsd.edu`) or from counsel.
- The engine has no such question over it. A commercial user can open a model of their own and leave iCH360 out.
