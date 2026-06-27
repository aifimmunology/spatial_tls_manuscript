# downstream_analysis

Hypothesis-testing analyses and manuscript figures, run on the fully annotated AnnData object
built in `../initial_processing/`. Notebooks are **run in numeric order**; the spatial-distance
notebooks share an object (`adata_distance.h5ad`) written by `02` and read by `03`/`04`.

| Notebook | Purpose | Environment |
|---|---|---|
| `00_sq_interactions` | Squidpy neighborhood interaction scoring, comparing cell–cell interaction enrichment across spatial zones (TLS / adventitia / parenchyma). | `space2` (CPU) |
| `01_DE` | Differential expression in CD4 T cells and between zones (volcano plots, zone marker dotplots); DE functions imported from `utils/de_utils.py`. | `space2` (CPU) |
| `02_distance_tls_bronchi` | Compute signed distance-to-TLS (negative inside) and distance-to-bronchi, derive within-TLS radial distance/categories, and gate spatial regions (IMAP-style). | `space2` (CPU) |
| `03_distance_tls_bronchi_vis` | Visualize the distance analyses — spatial distance maps, CD4 IMAP density plots, expression-over-distance heatmaps, and ligand–receptor spatial co-localization (Jensen–Shannon). | `space2` (CPU) |
| `04_distance_tls_radial_vis` | Visualize the within-TLS radial-distance analysis — radial maps, composition by radial category, and expression along the radial axis. | `space2` (CPU) |
| `save_final_adata` | Drop reference-only obs columns and write the final processed AnnData object for the manuscript. | `space2` (CPU) |
| `05_cellchat/` | CellChat ligand–receptor analysis across TLS radial sub-regions (see below). | mixed |

### `05_cellchat/`
| Notebook | Purpose | Environment |
|---|---|---|
| `00_cellchat_tls-category` | Python prep — relabel activated CD4 T cells by TLS radial category and save the object for CellChat. | `space2` (CPU) |
| `01_cellchat_R_tls-category` | R — run CellChat (spatial, distance-constrained) and save the fitted object. | `cellchat` (R) |
| `02_cellchat_tls-category` | R — downstream CellChat plots, with custom adaptations restricting signaling heatmaps to the CD4 TLS categories. | `cellchat` (R) |

Each notebook pins its environment at the top; full specs are in `../conda_envs/`. Data and output
locations are configured once in `../config/paths.py` (and `../config/paths.R` for the R notebooks).
