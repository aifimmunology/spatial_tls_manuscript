"""Shared plotting utilities for the TLS/bronchi distance analyses.

Consolidated from the downstream_analysis distance notebooks (02-04) so the
plotting functions live in one place. Imported in those notebooks via
``from distance_utils import ...`` (FUNCTIONS_DIR is on sys.path).

Functions
---------
- _draw_gates / plot_celltype_density_2d : 2D distance-space cell density with gates
- scatter_with_gaussian_kde / distance_xy_kde : expression-weighted 2D KDE in distance space
- scvelo_heatmap : gene-expression-along-a-trajectory heatmap (scVelo)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import scanpy as sc
from scipy.stats import gaussian_kde
import scvelo as scv

from de_utils import filter_adata_expressed_in_n_cells


def _draw_gates(ax, gates, x_clip, y_clip, show_labels=True, xlim=None, ylim=None):
    """Overlay named rectangular gates on a 2D distance plot.

    Each gate is {'x_min','x_max','y_min','y_max'} (any subset of keys), in RAW
    distance units. Missing/None bounds default to the plot edges
    (0 .. x_clip / y_clip).

    Gate borders are drawn as four independent line segments. If `xlim` / `ylim`
    (the view limits) are given, each border is drawn only where it lies within
    the view: a border whose constant edge is outside the limits is omitted, and
    borders running past an edge are clipped to it. This lets a gate whose bound
    sits beyond the axis (e.g. y_min below the y-axis minimum) still draw its
    in-view borders without a stray line off the plot.
    """
    def _clip_seg(a, b, lo, hi):
        """Clip the 1D segment spanning (a, b) to [lo, hi]; None if fully outside."""
        s0, s1 = min(a, b), max(a, b)
        s0, s1 = max(s0, lo), min(s1, hi)
        return (s0, s1) if s0 <= s1 else None

    for name, b in gates.items():
        x0 = b.get('x_min') if b.get('x_min') is not None else 0
        x1 = b.get('x_max') if b.get('x_max') is not None else x_clip
        y0 = b.get('y_min') if b.get('y_min') is not None else 0
        y1 = b.get('y_max') if b.get('y_max') is not None else y_clip

        xlo, xhi = xlim if xlim is not None else (min(x0, x1), max(x0, x1))
        ylo, yhi = ylim if ylim is not None else (min(y0, y1), max(y0, y1))

        # Horizontal borders (bottom=y0, top=y1): draw only if the edge's y is
        # within the y-view, clipped to the x-view.
        for ty in (y0, y1):
            if ylo <= ty <= yhi:
                seg = _clip_seg(x0, x1, xlo, xhi)
                if seg is not None:
                    ax.plot([seg[0], seg[1]], [ty, ty],
                            color='red', linestyle='-', linewidth=2)
        # Vertical borders (left=x0, right=x1): draw only if the edge's x is
        # within the x-view, clipped to the y-view.
        for tx in (x0, x1):
            if xlo <= tx <= xhi:
                seg = _clip_seg(y0, y1, ylo, yhi)
                if seg is not None:
                    ax.plot([tx, tx], [seg[0], seg[1]],
                            color='red', linestyle='-', linewidth=2)

        if show_labels:
            seg_x = _clip_seg(x0, x1, xlo, xhi)
            seg_y = _clip_seg(y0, y1, ylo, yhi)
            if seg_x is not None and seg_y is not None:
                # center the label on the VISIBLE portion of the gate
                ax.text((seg_x[0] + seg_x[1]) / 2, (seg_y[0] + seg_y[1]) / 2,
                        name, ha='center', va='center', fontsize=8, color='black')


def plot_celltype_density_2d(
    adata: sc.AnnData,
    sample_label: list[str],
    celltype: str = None,
    x_axis: str = 'avg_distance_to_bronchi_zone',
    y_axis: str = 'distance_to_tls',
    celltype_col: str = 'label_fine',
    gates: dict = None,
    x_clip: float = 450,
    y_clip: float = 450,
    x_ticks: list = None,
    y_ticks: list = None,
    max_cells: int = None,
    random_state: int = 0,
    figsize: tuple = (3, 3),
    point_size: float = 2,
    show_gate_labels: bool = True,
):
    """
    2D distance plot of the *density* of a cell population.

    Like `distance_xy_kde`, but instead of weighting by gene expression this
    shows where cells concentrate in the (x_axis, y_axis) distance space. Each
    cell is a point colored by its local 2D Gaussian-KDE density, with KDE
    contours.

    By default uses ALL cells; pass `celltype` to restrict to a single subset
    from `adata.obs[celltype_col]`. Optionally overlays rectangular `gates`
    (see `_draw_gates`) so different region definitions can be eyeballed before
    applying `assign_spatial_region`. Gate borders that fall outside the axis
    limits are omitted / clipped to the view.

    Subsampling
    -----------
    gaussian_kde evaluates an O(n^2) pairwise density, so plotting 300k+ cells
    (e.g. when celltype=None) is very slow and memory-heavy. Set `max_cells` to
    randomly downsample to at most that many points before the KDE; the density
    pattern is stable as long as max_cells is reasonably large (~20k+). Use
    `random_state` to control / vary the subsample.

    Use `x_ticks` / `y_ticks` to label the axes at chosen distance values, e.g.
    `y_ticks=[0, 25, 50, 100, 450]`. If omitted, matplotlib's default ticks are
    used.

    Parameters
    ----------
    adata : AnnData
    sample_label : list[str]
        sample_label values to include.
    celltype : str, optional
        Value in adata.obs[celltype_col] whose density is plotted. None (default)
        uses all cell types.
    celltype_col : str
        obs column holding the cell-type labels (default 'label_fine').
    x_axis, y_axis : str
        Distance columns for the x and y axes.
    gates : dict, optional
        region name -> {'x_min','x_max','y_min','y_max'} bounds (raw units) to overlay.
    x_clip, y_clip : float
        Axis limits in RAW units; cells beyond these are dropped from the plot.
    x_ticks, y_ticks : list, optional
        Distance values at which to place labeled ticks.
    max_cells : int, optional
        Randomly subsample to at most this many cells before the KDE. None = use all.
    random_state : int
        Seed for the subsample (default 0).
    point_size : float
        Scatter point size.

    Returns
    -------
    plt : module
        The pyplot module (call .show() or .savefig()).
    """
    fig = plt.figure(figsize=figsize, dpi=200)

    # Subset by sample_label
    adata = adata[adata.obs['sample_label'].isin(sample_label)]

    # Clip to the plotting window in RAW units (matches distance_xy_kde behavior)
    if x_clip:
        adata = adata[adata.obs[x_axis] <= x_clip, :]
    if y_clip:
        adata = adata[adata.obs[y_axis] <= y_clip, :]

    # Optionally subset to a cell type (celltype=None -> use all cell types)
    if celltype is not None:
        print(f"Subset to '{celltype_col}' == '{celltype}'")
        sub_adata = adata[adata.obs[celltype_col] == celltype].copy()
    else:
        print('Using all cell types')
        sub_adata = adata.copy()
    print('cells in subset: ', sub_adata.n_obs)

    # Downsample before the KDE if there are too many cells. gaussian_kde
    # evaluates an O(n^2) pairwise density, so 300k+ cells is very slow and
    # memory-heavy; cap the number of points used for the plot.
    if max_cells is not None and sub_adata.n_obs > max_cells:
        print(f'Subsampling {sub_adata.n_obs} -> {max_cells} cells for the density')
        sub_adata = sc.pp.subsample(sub_adata, n_obs=max_cells, copy=True,
                                    random_state=random_state)

    x = sub_adata.obs[x_axis].values
    y = sub_adata.obs[y_axis].values

    ax = fig.add_subplot(1, 1, 1)

    # color each cell by its local 2D density
    xy = np.vstack([x, y])
    z = gaussian_kde(xy)(xy)
    ax.scatter(x, y, c=z, s=point_size, marker='.', cmap='viridis')

    # density contours
    sns.kdeplot(x=x, y=y, ax=ax, color='#444444', linewidths=1)

    # axis limits with a small relative pad
    xpad = 0.02 * x_clip
    ypad = 0.02 * y_clip
    xlim = (-xpad, x_clip + xpad)
    # ylim = (-ypad, y_clip + ypad)
    ylim = (-100, y_clip + ypad)

    # overlay gates; borders outside the view limits are omitted / clipped
    if gates:
        _draw_gates(ax, gates, x_clip, y_clip,
                    show_labels=show_gate_labels, xlim=xlim, ylim=ylim)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    # ticks: label chosen distances if provided, else matplotlib defaults
    if x_ticks is not None:
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([str(t) for t in x_ticks])
    if y_ticks is not None:
        ax.set_yticks(y_ticks)
        ax.set_yticklabels([str(t) for t in y_ticks])

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.grid(False)
    plt.title(celltype if celltype is not None else 'all cells', fontsize=12)
    fig.tight_layout()
    return plt


def scatter_with_gaussian_kde(ax, x, y, weights=None, s=2, **kwargs):
    xy = np.vstack([x, y])
    try:
        z = gaussian_kde(xy, weights=weights)(xy)
    except:
        z = gaussian_kde(xy)(xy)

    kwargs.setdefault('marker', '.')
    ax.scatter(x, y, c=z, s=s, **kwargs)

def distance_xy_kde(
    adata: sc.AnnData,
    sample_label: list[str],
    gene_name: str = None,
    signature_col: str = None,
    x_axis: str = 'avg_distance_to_bronchi_zone',
    y_axis: str = 'distance_to_tls',
    celltype: str = None,
    celltype_col: str = 'label_fine',
    gates: dict = None,
    x_clip: float = 450,
    y_clip: float = 450,
    x_ticks: list = None,
    y_ticks: list = None,
    max_cells: int = None,
    random_state: int = 0,
    figsize: tuple = (3, 3),
    point_size: float = 2,
    show_gate_labels: bool = True,
    vline: float = None,
    hline: float = None,
    palette=None,
):
    """
    2D distance plot of gene-expression (or signature) density.

    Companion to `plot_celltype_density_2d`: same subsetting / clipping / gating
    parameters, but the KDE is weighted by `gene_name` expression (or a
    `signature_col` in obs) instead of plain cell density.
    """
    fig = plt.figure(figsize=figsize, dpi=200)

    # Subset by sample_label
    adata = adata[adata.obs["sample_label"].isin(sample_label)]

    # Optionally subset to a cell type (celltype=None -> use all cells)
    if celltype is not None:
        print(f"Subset to '{celltype_col}' == '{celltype}'")
        adata = adata[adata.obs[celltype_col] == celltype]
    else:
        print("Using all cells")

    # Clip x and y axes
    if x_clip:
        print('Clipping x to', x_clip)
        adata = adata[adata.obs[x_axis] <= x_clip, :]
    if y_clip:
        print('Clipping y to', y_clip)
        adata = adata[adata.obs[y_axis] <= y_clip, :]

    sub_adata = adata.copy()

    # Downsample before the KDE if there are too many cells. gaussian_kde
    # evaluates an O(n^2) pairwise density, so 300k+ cells is very slow and
    # memory-heavy; cap the number of points used for the plot.
    if max_cells is not None and sub_adata.n_obs > max_cells:
        print(f'Subsampling {sub_adata.n_obs} -> {max_cells} cells for the density')
        sub_adata = sc.pp.subsample(sub_adata, n_obs=max_cells, copy=True,
                                    random_state=random_state)

    if signature_col:
        print('plotting signature')
        gene_expr = sub_adata.obs[signature_col]
        title = signature_col
    else:
        print('plotting gene expr')
        gene_expr = sub_adata[:, gene_name].X.toarray().flatten() if hasattr(sub_adata[:, gene_name].X, "toarray") else sub_adata[:, gene_name].X.flatten()
        gene_expr = np.nan_to_num(gene_expr)
        title = gene_name

    ax = fig.add_subplot(1, 1, 1)


    # KDE scatter plot
    scatter_with_gaussian_kde(
        ax=ax,
        x=sub_adata.obs[x_axis],
        y=sub_adata.obs[y_axis],
        s=point_size,
        weights=(gene_expr - np.min(gene_expr)) ** 2,
        cmap="viridis",
    )

    # KDE weighted by expression
    sns.kdeplot(
        data=sub_adata.obs,
        x=x_axis,
        y=y_axis,
        ax=ax,
        color="#444444",
        linewidths=1,
        weights=(gene_expr - np.min(gene_expr)) ** 2,
        cmap="viridis",
    )



    # Add vline and hline if specified
    if vline is not None:
        # ax.axvline(x=vline, color='black', linestyle='--', linewidth=2)
        ax.plot([vline, vline], [hline, y_clip], color='black', linestyle='--', linewidth=2)
    if hline is not None:
        ax.axhline(y=hline, color='black', linestyle='--', linewidth=2)

    # # axis limits
    # xlim = (-10, x_clip + 10)
    # ylim = (-10, y_clip + 10)
    # axis limits with a small relative pad
    xpad = 0.02 * x_clip
    ypad = 0.02 * y_clip
    xlim = (-xpad, x_clip + xpad)
    # ylim = (-ypad, y_clip + ypad)
    ylim = (-100, y_clip + ypad)

    # overlay gates (same drawing as plot_celltype_density_2d); borders outside
    # the view limits are omitted / clipped to the axis
    if gates:
        _draw_gates(ax, gates, x_clip, y_clip,
                    show_labels=show_gate_labels, xlim=xlim, ylim=ylim)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    # ticks: label chosen distances if provided, else blank ticks at the extremes
    if x_ticks is not None:
        ax.set_xticks(x_ticks)
        ax.set_xticklabels([str(t) for t in x_ticks])
    # else:
    #     ax.set_xticks([0, x_clip])
    #     ax.set_xticklabels([])
    if y_ticks is not None:
        ax.set_yticks(y_ticks)
        ax.set_yticklabels([str(t) for t in y_ticks])
    # else:
    #     ax.set_yticks([0, y_clip])
    #     ax.set_yticklabels([])

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.grid(False)

    fig.tight_layout()
    # plt.title(title, fontsize=20)
    plt.title(gene_name)

    return(plt)


def scvelo_heatmap(
    adata: sc.AnnData,
    sample_label: list[str],
    sortby: str,
    key_name: str,
    key_value: str = None,
    highlight: list[str] = None,
    n_bins: int = 5,
    col_color = None, 
    palette = None,
    x_clip = None,
    expression_threshold = None, 
    figsize: tuple = None 
):
    """
    Create a heatmap to visualize gene expression trends in single-cell RNA-seq data,
    with options for subsetting, sorting, and highlighting genes.

    Parameters:
    - adata (sc.AnnData): Annotated data object containing single-cell RNA-seq data.
    - sample_label (List[str]): List of batch identifiers to subset the data.
    - key_name (str): String representing the key in `adata.obs` to use for subsetting cells.
    - key_value (str): String representing the value of `key_name` to subset to.
    - sortby (str): Variable to sort the heatmap by (e.g., "crypt_villi_axis").
    - highlight (List[str]): List of labels to highlight on the heatmap.
    - n_bins (int, optional): Integer specifying the number of bins to use for convolution (default: 5).

    Returns:
    - s (seaborn.matrix.ClusterGrid): Matplotlib figure object representing the heatmap.

    This function subsets the input data based on specified sample_labels and key-value pairs,
    filters genes expressed in a minimum percentage of cells, and creates a heatmap
    to visualize gene expression trends along a specified variable. The function also allows
    highlighting specific labels on the y-axis.
    """
    print("Creating Heatmap for sample labels ", " + ".join(sample_label))

    # Subset sample_label
    adata = adata[adata.obs["sample_label"].isin(sample_label)]
    if key_value:
        print(f"Subset to '{key_name}'=='{key_value}'")
        # Subset to key
        adata = adata[adata.obs[key_name] == key_value]
    else:
        print("using all cells")
        
    # Filter to include only genes that are expressed in 5% of the cells
    if expression_threshold:
        print('removing genes by expression threshold: ', expression_threshold, ' fraction')
        adata=filter_adata_expressed_in_n_cells(adata,fraction=expression_threshold)
    else:
        adata = adata.copy()
    print('adata shape: ', adata.shape)

    # Clip x axis
    if x_clip:
        print('clipping x to ', x_clip)
        print('max value before clip ', np.max(adata.obs[sortby]))
        adata = adata[adata.obs[sortby] <= x_clip, :]
        # adata.obs.loc[adata.obs[sortby] > x_clip, sortby] = x_clip
    else: 
        adata = adata.copy()
        
    n_convolve = len(adata) // n_bins
    print(f"Setting `n_convolve` to {n_convolve} ({n_bins} bins, {len(adata)} cells) ")
    # Plot
    s = scv.pl.heatmap(
        adata,
        var_names=adata.var_names,
        sortby=sortby,
        n_convolve=n_convolve,
        show=False,
        yticklabels=True,
        rasterized=True,
        # color_map=colormap,
        color_map='viridis',
        col_color = col_color,
        palette = palette, 
        figsize=figsize
        
    )
    ax = s.ax_heatmap

    # Loop through the x-axis tick labels and show/hide based on the 'highlight' list
    if highlight:
        for i, label in enumerate(ax.get_yticklabels()):
            if label.get_text() not in highlight:
                label.set_visible(False)
                ax.get_yticklines()[2 * i + 1].set_visible(False)
            ax.get_yticklines()[2 * i].set_visible(False)
    else:
        # Make sure all labels and ticks are visible
        for label in ax.get_yticklabels():
            label.set_visible(True)
    
        for tickline in ax.get_yticklines():
            tickline.set_visible(False)  # optional: hide all tick lines
        

    ax.set_xlabel("")
    ax.set_title("")
    return s
