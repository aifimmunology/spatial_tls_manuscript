# Spatial transcriptomics analysis for mouse lung allergy model
Analysis for Xenium figure in Bangs et al., 2025: Tertiary lymphoid structures support the development of allergen-specific TCF1+ progenitor CD4+ T cells

## Contents
1. inputs: required input csvs 
2. conda_envs: yml files to create conda environments require to run notebookds
3. initial_processing: notebooks for initial data processing, including cell labeling and spatial zone identification
4. downstream_analysis: notebooks for spatial hypothesis testing and plot generation

## Data download
Data is currently available for download by request only. 
The provided data includes the Xenium outputs for the HDM and Flu acute and memory time points, a processed adata object for the HDM Xenium data, and the processed CITE-seq seurat object. 

spatial_tls_manuscript_data
├── Xenium_outputs
│   ├── output-XETG00195__0036990__TIS08778-00-003__20241016__230427
│   ├── output-XETG00195__0037002__TIS08779-002-003__20241016__230427
│   ├── output-XETG00195__0036990__TIS08780-002-003__20241016__230427
│   └── output-XETG00195__0037002__TIS08781-003-003__20241016__230427
├── adata_processed
│   └── adata_distance_zones_structure.h5ad
└── CITEseq_processed 
    └── D3_D30_20250519.Rds


## Running the code:
1. Clone the repo
2. Create the required conda environments:

In cpu IDE:
- conda env create -n space2 --prefix /home/workspace/environment/space2 -f /home/workspace/spatial_tls_manuscript/conda_envs/space2_20250604.yml
- conda env create -n r_seurat --prefix /home/workspace/environment/space2 -f /home/workspace/spatial_tls_manuscript/conda_envs/r_seurat_20250604.yml
  
In gpu IDE: 
- conda env create -n cellcharter2 --prefix /home/workspace/environment/space2 -f /home/workspace/spatial_tls_manuscript/conda_envs/cellcharter2_20250604.yml

### IDE and conda environments
- [initial_processing/01_label_prediction_scanvi.ipynb](initial_processing/01_label_prediction_scanvi.ipynb) and [initial_processing/04_cellcharter.ipynb](initial_processing/04_cellcharter.ipynb) should be run in a gpu IDE with the cellcharter2 conda environment. All other notebooks can be run in a cpu IDE.
- [initial_processing/03_cd4_cite_integration.ipynb](initial_processing/03_cd4_cite_integration.ipynb) is in R and should be run in the r_seurat conda environment
- [initial_processing/00_scRNAseq_reference.ipynb](initial_processing/00_scRNAseq_reference.ipynb) and [initial_processing/02_label_refine.ipynb](initial_processing/02_label_refine.ipynb) can be run in conda environment space2
- Both notebooks in downstream_processing/ can be run in a cpu IDE with conda environment space2

