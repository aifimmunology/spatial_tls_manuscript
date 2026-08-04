# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repo holds the spatial transcriptomics (10x Xenium) analysis for the Xenium figure in Bangs et al., 2025 (*Tertiary lymphoid structures support the development of allergen-specific TCF1+ progenitor CD4+ T cells*). It is a manuscript analysis repo, not a software package: the work product is a sequence of Jupyter notebooks that process mouse lung Xenium data (HDM/Flu, acute/memory timepoints), label cells, define spatial zones/structures, and run hypothesis-testing analyses.

Zenodo DOI: https://doi.org/10.5281/zenodo.15858877. Data is available by request only.

## Environments

Multiple conda environments are required; YAML specs live in `conda_envs/`. Notebooks are written for a specific environment and some require GPU:

- `space2` (CPU IDE) — most processing and downstream analysis notebooks (incl. `geo/` and the Python cellchat prep)
- `cellcharter2` (GPU IDE) — `initial_processing/01_label_prediction_scanvi.ipynb`, `initial_processing/04_cellcharter.ipynb`
- `r_seurat` — R/Seurat notebook: `initial_processing/03_cd4_cite_integration.ipynb`
- `cellchat` — CellChat R analysis: `downstream_analysis/05_cellchat/01_cellchat_R_tls-category.ipynb` and `02_cellchat_tls-category.ipynb`

Create an env, e.g.: `conda env create -n space2 -f conda_envs/space2_20250604.yml`. See `README.md` for the full notebook→environment mapping.

## Git
Never commit changes. I handle all commits myself.

## Pipeline architecture

The analysis is ordered by the numeric notebook prefixes; later notebooks consume outputs written to disk by earlier ones. Run them in order.

**`initial_processing/`** (build the annotated AnnData object):
- `00`–`02`: build scRNAseq reference, predict labels via scANVI, refine labels
- `03`: integrate CITE-seq CD4 subsets (R/Seurat)
- `04`: CellCharter spatial niche/zone clustering (GPU)
- `05`–`07`: merge cell labels, label spatial zones, identify structures (TLS, bronchi)
- `08`: plot labels and zones

**`downstream_analysis/`** (hypothesis testing + manuscript figures):
- `00_sq_interactions`: squidpy neighborhood/interaction analysis across zones
- `01_DE`: differential expression (CD4 T cells, between zones); DE functions imported from `utils/de_utils.py`
- `02_distance_tls_bronchi`: compute distance-to-TLS / distance-to-bronchi and within-TLS radial distance, and gate spatial regions
- `03_distance_tls_bronchi_vis`, `04_distance_tls_radial_vis`: visualization of those distance analyses
- `05_cellchat/`: CellChat ligand-receptor analysis across TLS radial sub-regions. The `00`→`01`→`02` chain is the active pipeline (`00` Python prep → `01` R CellChat → `02` R downstream plots). `05_cellchat/cellchat_variations/` holds abandoned parameterizations, and `02_cellchat_tls-category-backup.ipynb` is a backup — ignore both as non-canonical.

**`geo/`** (GEO submission prep): `00` extract/organize the required Xenium output files (adding sample prefixes), `01` prepare the final processed AnnData object for the manuscript, `02` consolidate + zip + upload.

**`hise_download/`**: `download-hise.ipynb` — download the input data from HISE.

Each analysis folder has a `README.md` summarizing its notebooks. (Superseded `prev_organization/` dirs and the abandoned cellchat `not_using/` set have been removed.)

## Paths and shared code

Notebooks do not hardcode data locations. They bootstrap shared config/utilities like this:

```python
sys.path.append(str(Path.cwd().resolve().parents[0]))   # repo root onto path
from config.paths import BASE_OUTDIR, INPUTS_DIR, FUNCTIONS_DIR
sys.path.append(str(FUNCTIONS_DIR))
import plotting_utils as pu
```

- **`config/paths.py`** / **`config/paths.R`** — the single place data/output locations are defined (`DATA_DIR`, `BASE_OUTDIR`, `REF_DIR`, `INPUTS_DIR`, `FUNCTIONS_DIR`). To run on a new machine, edit these (currently set to `/home/workspace/...`). Keep the `.py` and `.R` versions in sync. The bootstrap `parents[0]` reaches the repo root for notebooks one dir deep; notebooks two dirs deep (`downstream_analysis/05_cellchat/`) use `parents[1]`. R notebooks `source()` `config/paths.R` via `dirname(getwd())` (or `dirname(dirname(getwd()))` two dirs deep).
- **`utils/`** (= `FUNCTIONS_DIR`) — shared functions imported as modules: `plotting_utils.py` (`pu`; e.g. `plot_general`, `zone_composition`, `zone_composition_heatmap_2`) and `de_utils.py` (`DE_test`, `DE_volcano`, `run_zone_DE_analysis`, `filter_adata_expressed_in_n_cells`). Notebooks `importlib.reload` these during dev — edit the module, re-run the import cell.
- **`inputs/`** (= `INPUTS_DIR`) — committed small CSVs: cell-label hierarchy mapping (`label_fine`→`label_coarser`), CellChatDB interactions, KEGG cytokines. Downstream code joins on these.

Outputs are written under `BASE_OUTDIR` into subdirs like `downstream_analysis/distance`, `cell_labeling/...`, `cellcharter_analysis/...` — directory names are constructed inline in each notebook's setup cell.

## Conventions

- `.gitignore` excludes all data and large/binary artifacts (`*.h5ad`, `*.parquet`, `*.pt`, `*.joblib`, `*.pkl`, `*.geojson`, `*.pdf`, `*.gz`, etc.). Only notebooks, config, utils, small input CSVs, and images are tracked. Do not commit data outputs.
- Cell labels exist at multiple granularities (`label_fine`, `label_medium`, `label_coarse`, `label_coarser`) defined in the inputs mapping CSV; analyses pick a granularity column.

## Publication readiness

The repo is being prepared for a public GitHub release. `codereview.md` (repo root) holds a tiered (P0–P3) publication-readiness review. The `.claude/skills/notebook-cleanup` skill encodes the per-notebook cleanup procedure (config-path loading, import trimming, docstrings, leak/secret removal, markdown documentation) — read it before cleaning a notebook, and check whether a notebook has real cell `id`s first (some lack them and must be rebuilt programmatically rather than edited via `NotebookEdit`).
