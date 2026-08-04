# Spatial Transcriptomics Analysis for Mouse Lung Allergy Model

Analysis for Xenium figure in Bangs et al., 2026:  
*Tertiary lymphoid structures support the development of allergen-specific TCF1+ progenitor CD4+ T cells*

---
![](images/xenium_image_crop.png)

## Analysis overview

This repo processes mouse lung 10x Xenium data (HDM/Flu) through two notebook pipelines, run in
numeric order. Each pipeline folder has its own `README.md` with a per-notebook table and the conda
environment used for each notebook (environment specs live in `conda_envs/`).

- **`initial_processing/`** — builds the fully annotated AnnData object: assemble the scRNA-seq
  reference, predict cell-type labels with scANVI, refine them (including CITE-seq CD4 integration),
  define spatial zones with CellCharter, and identify TLS structures.
- **`downstream_analysis/`** — hypothesis testing and manuscript figures: squidpy neighborhood
  interactions, differential expression, distance-to-TLS / distance-to-bronchi and within-TLS radial
  analyses, and CellChat ligand–receptor analysis (`05_cellchat/`).

Other folders: `inputs/` (small committed CSVs and color palettes), `config/` (data/output paths),
`utils/` (shared functions), `geo/` (GEO submission prep), and `hise_download/` (data download).

Analysis done by [@kathleenabadie](https://github.com/kathleenabadie).  
