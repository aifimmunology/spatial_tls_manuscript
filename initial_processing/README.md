# initial_processing

Builds the fully annotated Xenium AnnData object — predicting cell-type labels, refining them,
defining spatial zones, and identifying TLS structures. Notebooks are **run in numeric order**;
each consumes outputs written to disk by earlier ones. The final object is consumed by
`../downstream_analysis/`.

| Notebook | Purpose | Environment |
|---|---|---|
| `00_scRNAseq_reference` | Build the Hurskainen et al. scRNA-seq reference AnnData, subset to the Xenium gene panel. | `space2` (CPU) |
| `01_label_prediction_scanvi` | Co-embed Xenium + reference with scANVI and predict cell-type labels. | `cellcharter2` (GPU) |
| `02_label_refine` | Refine predicted labels by subclustering broad types and manually annotating (visual inspection of markers). | `space2` (CPU) |
| `03_cd4_cite_integration` | Integrate CITE-seq CD4 subsets with the Xenium CD4 cells (Seurat) to sub-label activated CD4 states. | `r_seurat` (R) |
| `04_cellcharter` | CellCharter spatial-neighborhood clustering (HDM + Flu samples) to define spatial niches. | `cellcharter2` (GPU) |
| `05_merge_cell_labels` | Merge refined labels + CITE CD4 sublabels; add the label hierarchy (fine/medium/coarse/coarser) and color palettes. | `space2` (CPU) |
| `06_label_zones` | Assign CellCharter clusters to interpretable zones (TLS, bronchi, vessels, adventitia, parenchyma, capsule). | `space2` (CPU) |
| `07_identify_structures` | Delineate individual TLS structures via alpha shapes and add per-cell structure IDs. | `space2` (CPU) |
| `08_plot_labels_zones` | Summary plots of cell labels and spatial zones for the manuscript figure. | `space2` (CPU) |

Each notebook pins its environment at the top; full specs are in `../conda_envs/`. Data and output
locations are configured once in `../config/paths.py` (and `../config/paths.R` for the R notebook).
