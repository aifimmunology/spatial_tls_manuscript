import math
import anndata as ad
import squidpy as sq
import pandas as pd
import scanpy as sc
import scvi
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_general(adata, x_col, y_col, color_col, palette, alpha=1, size=0.5, colors_sep=False, save_path=False, legend_outside=True, legend_remove = False, 
                 labels_on_plot=False, poly_list=False, poly_colors = ['black', 'blue', 'darkred']):
    """
    Generate a scatter plot from AnnData metadata with flexible coloring, labeling, and shape overlay.

    This function is useful for plotting 2D embeddings or spatial coordinates stored in `adata.obs`, 
    colored by any categorical or continuous column.
    """
    if adata.obs[color_col].isna().any():
        palette[np.nan] = 'black'
    
    if colors_sep:
        num_plots = len(adata.obs[color_col].unique().tolist())
        rows_plots = math.ceil(num_plots/2)
        print('rows_plots', rows_plots)
        fig, axs = plt.subplots(rows_plots, 2,
                                sharex = True, sharey = True,
                                figsize=(rows_plots*2, rows_plots*4)
                               )
        plot_num = 0
        for plot_num, id in enumerate(sorted(adata.obs[color_col].unique().tolist())):
            if plot_num < rows_plots:
                col_num = 0
                row_num = plot_num
            else:
                col_num = 1
                row_num = plot_num - rows_plots
    
            adata_zone = adata[adata.obs[color_col]==id, :]
            sns.scatterplot(data=adata_zone.obs, x=x_col, y=y_col, 
                    hue=color_col,
                    ax = axs[row_num, col_num], 
                    palette=[palette[plot_num]], 
                    alpha=alpha, s=size, 
                   )  
    else:
    
        fig, axs = plt.subplots(1, 1, figsize=(10, 10))
        sns.scatterplot(data=adata.obs, x=x_col, y=y_col, 
                hue=color_col,
                ax = axs, 
                palette=palette, 
                alpha=alpha, s=size, 
               )
        axs.spines['top'].set_visible(False)
        axs.spines['right'].set_visible(False)
        axs.spines['bottom'].set_visible(False)
        axs.spines['left'].set_visible(False)
        axs.get_xaxis().set_ticks([])
        axs.get_yaxis().set_ticks([])
        axs.set_xlabel('')
        axs.set_ylabel('')
    
         # add shape
        if poly_list:
            for i, poly in enumerate(poly_list):
                for p in poly.geoms:
                    axs.plot(*p.exterior.xy, c=poly_colors[i], linewidth=2)
                
         # Modify the legend
        if legend_outside:
            handles, labels = axs.get_legend_handles_labels()
            axs.legend(handles=handles, labels=labels, markerscale=5, bbox_to_anchor=(.8, 1.0), loc='upper left')  # Adjust the markerscale to increase the size of the markers
        elif labels_on_plot:
            # Annotate identities on the plot
            for i, txt in enumerate(adata.obs[color_col]):
                axs.text(adata.obs[x_col][i], adata.obs[y_col][i], txt, fontsize=8)
    
        elif legend_remove:
            axs.legend().remove()
        else:
            axs.legend(handles=handles, labels=labels, markerscale=5)  # Adjust the markerscale to increase the size of the markers
    
    if save_path:
        fig.savefig(save_path)



def zone_composition(
    adata,
    label_col,
    zone_col,
    zones_labels=None,
    remove_unlabeled=False,
    save_path=False,
    palette_map=None,
    plot_counts=False,
    zone_ncell_filter=None,
    horizontal=False,
    whitespace = 0.7
):
    """
    Plot stacked bar charts showing the composition of cell labels within each spatial zone.

    This function groups cells by a spatial zone column and a label column (e.g., cell type), 
    calculates either the count or fraction of each label in each zone, and visualizes the result 
    as a stacked bar chart (horizontal or vertical).
    """
    
    plt.rcParams.update({'font.size': 14})

    # Set palette
    if palette_map is not None:
        palette = palette_map
    else:
        palette = sns.color_palette(
            colorcet.glasbey_dark,
            n_colors=len(adata.obs[label_col].unique().tolist())
        )

    # Optionally remove unlabeled
    if remove_unlabeled:
        adata = adata[adata.obs[label_col] != 'unlabeled', :]

    # Summarize counts
    summary_df = adata.obs.groupby([zone_col, label_col]).size().reset_index(name='count')

    # Pivot table: zones x celltypes
    summary_counts = pd.pivot_table(
        summary_df, values='count', index=zone_col, columns=label_col, fill_value=0
    )

    # Filter by cell count if requested
    if zone_ncell_filter:
        summary_counts = summary_counts[summary_counts.sum(axis=1) > zone_ncell_filter]

    # Calculate fractions
    row_sums = summary_counts.sum(axis=1)
    summary_frac = summary_counts.div(row_sums, axis=0)

    if horizontal:
        summary_counts = summary_counts.iloc[::-1]
        summary_frac = summary_frac.iloc[::-1]

    if plot_counts:
        if horizontal:
            ax1 = summary_counts.plot.barh(stacked=True, color=palette, width=whitespace)
        else:
            ax1 = summary_counts.plot.bar(stacked=True, color=palette, rot=0, width=whitespace)
    else:
        if horizontal:
            ax1 = summary_frac.plot.barh(stacked=True, color=palette, width=whitespace)
        else:
            ax1 = summary_frac.plot.bar(stacked=True, color=palette, rot=0, width=whitespace)
    
    # Remove plot box (top/right)
    sns.despine(ax=ax1, top=True, right=True, left=False, bottom=False)

    # Titles and axis labels
    ax1.set_title('')
    if horizontal:
        ax1.set_ylabel('')
        ax1.set_xlabel('Count' if plot_counts else 'Fraction')
    else:
        ax1.set_xlabel('')
        ax1.set_ylabel('Count' if plot_counts else 'Fraction')

    # Handle custom zone labels
    if zones_labels is not None:
        if horizontal:
            ax1.set_yticklabels(zones_labels)
        else:
            ax1.set_xticklabels(zones_labels)

    # Legend styling
    ax1.legend(
        bbox_to_anchor=(1.8, 1),
        loc='upper right',
        fontsize=12,
        ncol=2,
        frameon=True,
        framealpha=0,
        edgecolor='black',
        handletextpad=0.5,
        columnspacing=0.8
    )

    return ax1


def zone_composition_heatmap_2(adata, label_col, zone_col, scale_by = 'rows', swap_axes = False, save_path=None):

    summary_df = adata.obs.groupby([zone_col, label_col]).size().reset_index(name='count')
    summary_counts = pd.pivot_table(summary_df, values='count', index=zone_col, columns=label_col, fill_value=0)

    # calculate fractions by row
    row_sums= summary_counts.sum(axis=1)
    summary_frac_row = summary_counts.div(row_sums, axis=0)  

    # calculate fractions by column
    col_sums= summary_counts.sum(axis=0)
    summary_frac_col = summary_counts.div(col_sums, axis=1)  

    if scale_by == 'rows':
        summary_frac_scaled = summary_frac_row
    elif scale_by == 'cols':
        summary_frac_scaled = summary_frac_col
        

    if swap_axes:
        print('swap')
        xticklabels=summary_frac_scaled.index
        yticklabels=summary_frac_scaled.columns
        summary_frac_scaled = summary_frac_scaled.T
        fig, ax = plt.subplots(1, 1, figsize=(2, 6))
    else:
        xticklabels=summary_frac_scaled.columns
        yticklabels=summary_frac_scaled.index
        fig, ax = plt.subplots(1, 1, figsize=(6, 2))
        
    
    
    # fig, ax = plt.subplots(1, 1, figsize=(10, 1))
    sns.heatmap(
        summary_frac_scaled,
        xticklabels=xticklabels,
        yticklabels=yticklabels,
        ax=ax,
        cmap='Blues',
        linecolor="white",
        linewidths=0.5,
    )
    ax.tick_params(axis='x', labelsize=14)

    if save_path is not None:
        fig.savefig(save_path, dpi=300, bbox_inches='tight', transparent=True)
    return(summary_frac_scaled)



def plot_fraction_stacked_bar_horizontal(
    adata,
    x_col: str,
    stack_col: str,
    palette: dict,
    sample_col: str = None,
    figsize=(10,6),
    show_legend: bool = True,
    gap: float = 0.05,  # fraction of bar_width reserved as gap
    return_df: bool = False
):
    """
    Create horizontal stacked bar plots of fractional composition of `stack_col` categories
    within each `x_col` group. If `sample_col` is provided, bars are split side-by-side by
    `sample_col`; otherwise a single bar is drawn per `x_col` group.

    Parameters
    ----------
    adata : AnnData
        AnnData object with obs metadata.
    x_col : str
        Column in adata.obs for category grouping (e.g., 'zone_consol').
    stack_col : str
        Column in adata.obs used for stacked bar coloring (e.g., 'label_fine').
    palette : dict
        Dictionary mapping `stack_col` categories to colors.
    sample_col : str, optional
        Column in adata.obs to split bars side-by-side (e.g., 'sample_label').
        If None (default), a single bar is drawn per `x_col` group.
    figsize : tuple
        Figure size for the matplotlib plot.
    gap : float
        Fraction of bar_width to leave as gap between sample_col bars.
    return_df : bool
        If True, return (ax, plot_df) where plot_df holds the fractional
        composition used for plotting. If False (default), return only ax.
    """

    # Step 1: dataframe
    group_cols = [x_col, stack_col] if sample_col is None else [sample_col, x_col, stack_col]
    df = adata.obs[group_cols].copy()

    # Step 2: counts
    frac_index = [x_col] if sample_col is None else [x_col, sample_col]
    counts = (
        df.groupby(frac_index + [stack_col])
          .size()
          .reset_index(name='count')
    )

    # Step 3: fractions
    counts['fraction'] = counts.groupby(frac_index)['count'].transform(lambda x: x / x.sum())

    # Step 4: pivot
    plot_df = counts.pivot_table(
        index=frac_index,
        columns=stack_col,
        values='fraction',
        fill_value=0
    )

    # Step 5: plot
    fig, ax = plt.subplots(figsize=figsize)

    zones = plot_df.index.get_level_values(x_col).unique()
    samples = (
        [None] if sample_col is None
        else plot_df.index.get_level_values(sample_col).unique()
    )

    bar_height = 0.8 / len(samples)
    effective_height = bar_height * (1 - gap)  # shrink to leave a gap

    y = np.arange(len(zones))

    for i, sample in enumerate(samples):
        left = [0]*len(zones)
        for lf in plot_df.columns:
            widths = []
            for z in zones:
                key = z if sample is None else (z, sample)
                widths.append(plot_df.loc[key, lf] if key in plot_df.index else 0)

            offset = 0 if sample is None else (i-0.5)*bar_height
            ax.barh(
                [pos + offset for pos in y],
                widths,
                effective_height,
                left=left,
                label=lf if i==0 else "_nolegend_",
                color=palette[lf],
                linewidth=0.25
            )

            left = [l+w for l, w in zip(left, widths)]

    # axis labels
    ax.set_yticks(y)
    ax.set_yticklabels(zones)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xlabel('')
    ax.set_ylabel('')

    
    # Reverse the order of categories
    ax.invert_yaxis()

    if show_legend:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title=stack_col)
    else:
        ax.legend().set_visible(False)

    sns.despine()
    plt.tight_layout()

    if return_df:
        return ax, plot_df
    return ax

    