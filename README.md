# Spatial transcriptomics analysis for mouse lung allergy model
Analysis for Xenium figure in Bangs et al., 2025: Tertiary lymphoid structures support the development of allergen-specific TCF1+ progenitor CD4+ T cells

## Contents
1. inputs: required input csvs 
2. conda_envs: yml files to create conda environments require to run notebookds
3. initial_processing: notebooks for initial data processing, including cell labeling and spatial zone identification
4. downstream_analysis: notebooks for spatial hypothesis testing and plot generation


## Running the code:
1. Clone the repo
2. Create the required conda environments:

In cpu IDE:
- conda env create -n space2 --prefix /home/workspace/environment/space2 -f /home/workspace/spatial_tls_manuscript/conda_envs/space2_20250604.yml
- conda env create -n r_seurat --prefix /home/workspace/environment/space2 -f /home/workspace/spatial_tls_manuscript/conda_envs/r_seurat_20250604.yml
  
In gpu IDE: 
- conda env create -n cellcharter2 --prefix /home/workspace/environment/space2 -f /home/workspace/spatial_tls_manuscript/conda_envs/cellcharter2_20250604.yml

### IDE and conda environments
- [intial_processing/01_label_prediction_scanvi](intial_processing/01_label_prediction_scanvi) and [intial_processing/04_cellcharter](intial_processing/04_cellcharter) should be run in a gpu IDE with the cellcharter2 conda environment. All other notebooks can be run in a cpu IDE.
- [intial_processing/03_cd4_cite_integration](intial_processing/03_cd4_cite_integration) is in R and should be run in the r_seurat conda environment
- [intial_processing/00_scRNAseq_reference](intial_processing/00_scRNAseq_reference) and [intial_processing/02_label_refine](intial_processing/02_label_refine) can be run in conda environment space2
- Both notebooks in downstream_processing/ can be run in a cpu IDE with conda environment space2

