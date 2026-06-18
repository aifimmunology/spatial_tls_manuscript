from pathlib import Path
import sys
import os
import time
import warnings
import re
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import anndata
import scanpy as sc
import squidpy as sq


import scipy
from scipy.spatial import distance_matrix
from scipy.spatial.distance import cdist
from scipy import stats
import scipy.ndimage as ndi
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter
from scipy.signal import find_peaks

from shapely.geometry import Point

from shapely.geometry import MultiPolygon
from alphashape import alphashape
from alphashape import optimizealpha
import geopandas as gpd

from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import PowerTransformer

import pickle
import joblib

import random
import matplotlib.colors as mcolors

import matplotlib as mpl
mpl.rcParams['axes.titlesize'] = 24
mpl.rcParams['pdf.fonttype'] = 42

import warnings
warnings.simplefilter(action='ignore', category=Warning)

sys.version_info

def filter_adata_expressed_in_n_cells(adata, fraction=0.05):
    """
    Filter the adata to only include genes expressed in more than a certain fraction of cells

    Parameters:
    - adata (anndata): The anndata object containing the gene expression
    - fraction (float): The fraction of cells a gene must be expressed in to be included

    Return:
    - adata (anndata): The filtered anndata object
    """
    adata = adata.copy()
    adata.layers["bin"] = adata.X > 0
    gene_expressed_in_fraction_cells = np.mean(adata.layers["bin"], axis=0)
    keep = gene_expressed_in_fraction_cells > fraction
    adata = adata[:, keep]
    return adata
    

def DE_test(adata, cell_column, cell_list=None, groupby=None, groups = 'all', reference= 'rest', stat_method = 'wilcoxon', use_raw = False, expression_threshold=None):
    """
    Performs differential expression (DE) testing on an AnnData object using Scanpy's `rank_genes_groups`.

    Parameters:
        adata (AnnData): The input AnnData object.
        cell_column (str): Column in `adata.obs` used to filter cells by type.
        cell_list (list, optional): List of cell types to subset before DE testing. If None, all cells are used.
        groupby (str): Column in `adata.obs` used to define groups for DE testing.
        groups (str or list, optional): Groups to test. Default is 'all'.
        reference (str, optional): Reference group for comparison. Default is 'rest'.
        stat_method (str, optional): Statistical test to use. 
        use_raw (bool, optional): Whether to use `.raw` for DE testing. Default is False.
        expression_threshold (float, optional): If set, filters genes to only those expressed in at least this fraction of cells.

    Returns:
        AnnData: A subsetted copy of the input `AnnData` object with DE results stored in `adata.uns['rank_genes_groups']`.
    """
    
    if cell_list:
        print('subsetting by celltype')
        adata_sub = adata[adata.obs[cell_column].isin(cell_list), :]
    else:
        print('not subsetting by cell type')
        adata_sub = adata.copy()
        
    # Optionally to include only genes that are expressed in 5% of the cells
    if expression_threshold is not None:
        print('removing genes by expression threshold: ', expression_threshold)
        adata_sub=filter_adata_expressed_in_n_cells(adata_sub,fraction=expression_threshold)
              
    sc.tl.rank_genes_groups(adata_sub, groupby, method=stat_method, corr_method='benjamini-hochberg', groups=groups, reference= reference, use_raw=use_raw)
    return(adata_sub)

from adjustText import adjust_text

def DE_volcano(
    df_DEG,
    pval_col="pvals_adj",
    gene_col="names",
    lfc_col="logfoldchanges",
    lfc_threshold=1.0,
    genes_to_label=None,
    pval_cutoff=0.05,
    axis_lims = None,
    down_color="darkblue",
    up_color="darkred",
    nonsig_color="gray",
    title=""
):
    """
    Generate a volcano plot from DE results.

    Parameters
    ----------
    df_DEG : pd.DataFrame
        DE result table with logFC and p-values.
    pval_col : str
        Name of column with adjusted p-values.
    gene_col : str
        Name of column with gene names.
    lfc_col : str
        Name of column with log2 fold changes.
    lfc_threshold : float
        Minimum |log2FC| to consider significant.
    genes_to_label : set or list, optional
        Genes to label regardless of significance.
    pval_cutoff : float
        Adjusted p-value threshold.
    down_color : str
        Color for significantly downregulated genes.
    up_color : str
        Color for significantly upregulated genes.
    nonsig_color : str
        Color for non-significant points.
    title : str
        Title for the plot.

    Returns
    -------
    matplotlib.figure.Figure
        The volcano plot figure.
    """

    df_DEG = df_DEG.copy()

    # Avoid log(0) issues
    df_DEG[pval_col] = df_DEG[pval_col].replace(0, 1e-300)
    df_DEG['-log10p'] = -np.log10(df_DEG[pval_col])

    # Default label set
    if genes_to_label is None:
        genes_to_label = set()

    # Assign color for each point
    def get_color(row):
        sig = row[pval_col] < pval_cutoff
        up = row[lfc_col] > lfc_threshold
        down = row[lfc_col] < -lfc_threshold
        if sig and up:
            genes_to_label.add(row[gene_col])
            return up_color
        elif sig and down:
            genes_to_label.add(row[gene_col])
            return down_color
        return nonsig_color

    df_DEG['color'] = df_DEG.apply(get_color, axis=1)

    # Define axis lims
    if axis_lims is None:
        lim = max(abs(min(df_DEG[lfc_col])), max(df_DEG[lfc_col]))
        axis_lims = [-lim, lim]
    # Clip axes
    df_DEG["lfc_clipped"] = df_DEG[lfc_col].clip(axis_lims[0], axis_lims[1])

    # Start plot
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(
        df_DEG["lfc_clipped"],
        df_DEG["-log10p"],
        c=df_DEG["color"],
        s=20,
        zorder=1,
        alpha=0.7
    )

    # Significance threshold lines
    ax.axhline(-np.log10(pval_cutoff), color="gray", linestyle="--", zorder=2)
    ax.axvline(-lfc_threshold, color="gray", linestyle="--", zorder=2)
    ax.axvline(lfc_threshold, color="gray", linestyle="--", zorder=2)

    # Label significant genes
    text_objects = []
    for _, row in df_DEG.iterrows():
        if row[gene_col] in genes_to_label:
            ax.scatter(row["lfc_clipped"], row["-log10p"], color=row["color"], s=20, zorder=3)
            text = ax.text(
                row["lfc_clipped"], row["-log10p"], row[gene_col],
                fontsize=10, zorder=4
            )
            text_objects.append(text)

    adjust_text(
        text_objects,
        expand_points=(1.5, 1.5),
        expand_text=(1.2, 1.2),
        force_text=(0.2, 0.5),
        force_points=(0.3, 0.5),
        lim=300,
        arrowprops=dict(arrowstyle="-", color="black", lw=0.5)
    )

    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(p-adj)")
    ax.set_title(title)
    ax.set_xlim(axis_lims[0], axis_lims[1])

    plt.tight_layout()
    return fig

def run_zone_DE_analysis(
    adata,
    celltype,
    zone_col,
    zone_order,
    celltype_col = None,
    reference_zone='rest',
    pval_cutoff=0.05,
    lfc_cutoff=1,
    expression_threshold = 0.01,
    plot=True, dot_save_path = None
):
    """
    Run differential expression within a cell type across two zones, return DE tables.

    Parameters
    ----------
    adata : AnnData
        Input AnnData object with raw counts or log-normalized data.
    celltype : str
        The cell type to isolate for DE analysis.
    celltype_col : str
        Column in adata.obs containing cell type labels.
    zone_col : str
        Column in adata.obs indicating zone (e.g. 'Tzone', 'Foll').
    zone_order : list of str
        Order of zones for plotting and sorting.
    reference_zone : str
        Zone to use as reference in DE.
    pval_cutoff : float
        Adjusted p-value threshold.
    lfc_cutoff : float
        Log fold change threshold.
    plot : bool
        Whether to generate volcano plots and dotplot.

    Returns
    -------
    df_DEG_all : pd.DataFrame
        Full DE table across zones.
    df_DEG_filt : pd.DataFrame
        Filtered DE genes.
    """

    # Subset to cell type and zones of interest
    if celltype_col is not None:
        adata_subset = adata[adata.obs[celltype_col] == celltype, :]
        adata_subset = adata_subset[adata_subset.obs[zone_col].isin(zone_order), :]
        print(f"Subset shape: {adata_subset.shape} for celltype: {celltype}")
    else:
        adata_subset = adata.copy()
        print(f"Subset shape: {adata_subset.shape} - no celltype subsetting")
        

    # Run DE (assumes DE_test adds results to .uns['rank_genes_groups'])
    adata_subset = DE_test(
        adata_subset,
        celltype_col,
        groupby=zone_col,
        groups='all',
        reference=reference_zone,
        stat_method='wilcoxon',
        use_raw=False,
        expression_threshold = expression_threshold
    )

    results = adata_subset.uns['rank_genes_groups']
    groups = results['names'].dtype.names

    # Build full DE table and volcano plots
    df_DEG_list = []
    for group in groups:
        df_DEG = pd.DataFrame({
            key: results[key][group]
            for key in ['names', 'logfoldchanges', 'pvals_adj']
        })
        df_DEG['group'] = group
        df_DEG_list.append(df_DEG)

        if plot:
            _ = DE_volcano(
                df_DEG,
                pval_col='pvals_adj',
                gene_col='names',
                lfc_col='logfoldchanges',
                lfc_threshold=lfc_cutoff,
                title=f'{celltype}: {group}'
            )

    df_DEG_all = pd.concat(df_DEG_list, ignore_index=True)

    # Get max-LFC per gene across zones
    df_max_lfc = df_DEG_all.loc[
        df_DEG_all.groupby("names")["logfoldchanges"].idxmax()
    ].copy()

    df_max_lfc["group_rank"] = df_max_lfc["group"].map(
        {z: i for i, z in enumerate(zone_order)}
    )
    df_max_lfc = df_max_lfc.sort_values(["group_rank", "logfoldchanges"])
    gene_order = df_max_lfc["names"].tolist()

    print(f"Unique DEGs before pval/lfc filter: {len(gene_order)}")

    # Apply filtering
    df_DEG_filt = df_DEG_all[
        (df_DEG_all["pvals_adj"] < pval_cutoff) &
        (df_DEG_all["logfoldchanges"].abs() > lfc_cutoff)
    ]
    gene_order_filt = [g for g in gene_order if g in set(df_DEG_filt["names"])]
    print(f"Unique DEGs after pval/lfc filter: {len(gene_order_filt)}")

    # Dotplot
    if plot and gene_order_filt:
        sc.pl.dotplot(
            adata_subset,
            var_names=gene_order_filt,
            groupby=zone_col,
            categories_order=zone_order,
            standard_scale='var',
            swap_axes=False,
            show=False
        )
    if dot_save_path is not None:
        plt.savefig(dot_save_path, dpi=300, bbox_inches="tight")
            

    return df_DEG_all, df_DEG_filt
