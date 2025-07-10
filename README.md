# Spatial Transcriptomics Analysis for Mouse Lung Allergy Model

Analysis for Xenium figure in Bangs et al., 2025:  
*Tertiary lymphoid structures support the development of allergen-specific TCF1+ progenitor CD4+ T cells*

---

## Contents

1. **inputs**: Required input CSV files  
2. **conda_envs**: YAML files to create conda environments required to run notebooks  
3. **initial_processing**: Notebooks for initial data processing, including cell labeling and spatial zone identification  
4. **downstream_analysis**: Notebooks for spatial hypothesis testing and plot generation  

---

## Data Download

> **Data access:** Currently available for download **by request only**.

The provided data includes:  
- Xenium outputs for HDM and Flu acute and memory time points  
- Processed adata object for the HDM Xenium data  
- Processed CITE-seq Seurat object  
```
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
```

## Running the code:
1. Clone the repo
2. Create the required conda environments:

In cpu IDE:
- conda env create -n space2 --prefix /home/workspace/environment/space2 -f /home/workspace/spatial_tls_manuscript/conda_envs/space2_20250604.yml
- conda env create -n r_seurat --prefix /home/workspace/environment/space2 -f /home/workspace/spatial_tls_manuscript/conda_envs/r_seurat_20250604.yml
  
In gpu IDE: 
- conda env create -n cellcharter2 --prefix /home/workspace/environment/space2 -f /home/workspace/spatial_tls_manuscript/conda_envs/cellcharter2_20250604.yml

### IDE and conda environments 
```
| Notebook(s)                                           | Recommended Environment   | Notes                                |
| ----------------------------------------------------| -------------------------| ------------------------------------|
| `initial_processing/01_label_prediction_scanvi.ipynb` | `cellcharter2` (GPU IDE) | Requires GPU environment             |
| `initial_processing/04_cellcharter.ipynb`             | `cellcharter2` (GPU IDE) | Requires GPU environment             |
| `initial_processing/03_cd4_cite_integration.ipynb`    | `r_seurat`               | R environment for Seurat processing  |
| `initial_processing/00_scRNAseq_reference.ipynb`      | `space2` (CPU IDE)       | CPU environment                      |
| `initial_processing/02_label_refine.ipynb`            | `space2` (CPU IDE)       | CPU environment                      |
| `downstream_analysis/00_downstream_analysis.ipynb`    | `space2` (CPU IDE)       | CPU environment                      |
| `downstream_analysis/01_downstream_analysis_hypothesis-testing.ipynb`    | `space2` (CPU IDE)       | CPU environment                      |
```


