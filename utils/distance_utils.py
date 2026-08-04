"""Shared plotting utilities for the TLS/bronchi distance analyses.

Functions
---------
- _draw_gates : overlay named rectangular gates as shaded boxes on a 2D axis (internal helper)
- plot_celltype_density_2d : 2D cell scatter in distance space, colored by local KDE density, with optional gate overlays
- scatter_with_gaussian_kde : scatter plot colored by 2D Gaussian KDE density (internal helper)
- distance_xy_kde : 2D cell scatter in distance space, KDE-weighted by gene expression or an obs signature column
- scvelo_heatmap : gene expression heatmap sorted by a spatial distance axis (wraps scv.pl.heatmap)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import seaborn as sns
import scanpy as sc
from scipy.stats import gaussian_kde
import scvelo as scv

from de_utils import filter_adata_expressed_in_n_cells


def _draw_gates(ax, gates, x_clip, y_clip, palette=None, gate_alpha=0.25,
                xlim=None, ylim=None):
    """Overlay named rectangular gates as light shaded boxes on a 2D distance plot.

    Each gate is {'x_min','x_max','y_min','y_max'} (any subset of keys), in RAW
    distance units. Missing/None bounds default to the plot edges
    (0 .. x_clip / y_clip).

    Each gate is drawn as a single filled rectangle with `zorder=0`, so it sits
    in the background behind the cell scatter dots. If `xlim` / `ylim` (the view
    limits) are given, the rectangle is clipped to the visible window: a gate
    fully outside the view is omitted, and one running past an edge is clipped to
    it. This lets a gate whose bound sits beyond the axis (e.g. y_min below the
    y-axis minimum) still shade its in-view portion without spilling off the plot.

    Colors are resolved per gate (in insertion order, index ``i``):
    - `palette` is a dict {gate_name: color} -> look up by name, falling back to
      the default matplotlib color cycle for any missing name.
    - `palette` is a list/tuple -> cycle by index.
    - `palette` is a str -> that single color for all gates.
    - `palette` is None -> default matplotlib color cycle by index.
    """
    def _clip_seg(a, b, lo, hi):
        """Clip the 1D segment spanning (a, b) to [lo, hi]; None if fully outside."""
        s0, s1 = min(a, b), max(a, b)
        s0, s1 = max(s0, lo), min(s1, hi)
        return (s0, s1) if s0 <= s1 else None

    cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']

    for i, (name, b) in enumerate(gates.items()):
        x0 = b.get('x_min') if b.get('x_min') is not None else 0
        x1 = b.get('x_max') if b.get('x_max') is not None else x_clip
        y0 = b.get('y_min') if b.get('y_min') is not None else 0
        y1 = b.get('y_max') if b.get('y_max') is not None else y_clip

        xlo, xhi = xlim if xlim is not None else (min(x0, x1), max(x0, x1))
        ylo, yhi = ylim if ylim is not None else (min(y0, y1), max(y0, y1))

        # clip the gate to the visible window; skip if fully outside
        seg_x = _clip_seg(x0, x1, xlo, xhi)
        seg_y = _clip_seg(y0, y1, ylo, yhi)
        if seg_x is None or seg_y is None:
            continue

        # resolve the fill color for this gate
        if isinstance(palette, dict):
            color = palette.get(name, cycle[i % len(cycle)])
        elif isinstance(palette, (list, tuple)):
            color = palette[i % len(palette)]
        elif isinstance(palette, str):
            color = palette
        else:
            color = cycle[i % len(cycle)]

        ax.add_patch(mpatches.Rectangle(
            (seg_x[0], seg_y[0]), seg_x[1] - seg_x[0], seg_y[1] - seg_y[0],
            facecolor=color, edgecolor='none', alpha=gate_alpha, zorder=0))


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
    y_min: float = None,
    x_ticks: list = None,
    y_ticks: list = None,
    max_cells: int = None,
    random_state: int = 0,
    figsize: tuple = (3, 3),
    point_size: float = 2,
    gate_palette=None,
    gate_alpha: float = 0.25,
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
    `y_ticks=[0, 25, 50, 100, 450]`. If omitted, ticks are placed every 200 um.

    `distance_xy_kde` uses identical axis-limit / tick handling so the two plots
    are directly comparable. The axes fill the figure box, so both axes have the
    same physical length regardless of their distance range.

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
        region name -> {'x_min','x_max','y_min','y_max'} bounds (raw units) to
        overlay as light shaded background boxes (cell dots draw on top).
    gate_palette : dict | list | tuple | str, optional
        Colors for the gate fills. dict {gate_name: color} colors gates by name
        (default cycle for any missing name); list/tuple cycles by gate order;
        str applies one color to all gates; None (default) uses the matplotlib
        color cycle.
    gate_alpha : float
        Opacity of the shaded gate boxes (default 0.25).
    x_clip, y_clip : float, optional
        Axis limits in RAW units; cells beyond these are dropped from the plot.
        None = no clipping; the axis limits, ticks and gate boxes then fall back
        to the data range.
    y_min : float, optional
        Lower y-axis limit. None (default) uses the (padded) data minimum; set a
        fixed value to keep the lower edge constant across panels/samples.
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

    # axis limits with a small relative pad. When x_clip / y_clip is None (no
    # clipping requested) fall back to the data max so the limits, ticks and
    # gate boxes still have a finite upper bound.
    x_hi = x_clip if x_clip is not None else (float(np.nanmax(x)) if len(x) else 1.0)
    y_hi = y_clip if y_clip is not None else (float(np.nanmax(y)) if len(y) else 1.0)
    xpad = 0.02 * x_hi
    ypad = 0.02 * y_hi
    xlim = (-xpad, x_hi + xpad)
    # lower y limit: explicit y_min if given, else the padded data minimum
    if y_min is not None:
        y_lo = y_min
    else:
        y_lo = (float(np.nanmin(y)) if len(y) else 0.0) - ypad
    ylim = (y_lo, y_hi + ypad)

    # overlay gates; borders outside the view limits are omitted / clipped
    if gates:
        _draw_gates(ax, gates, x_hi, y_hi,
                    palette=gate_palette, gate_alpha=gate_alpha,
                    xlim=xlim, ylim=ylim)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    # ticks: the values in x_ticks / y_ticks if given, else every 200 um
    xt = x_ticks if x_ticks is not None else np.arange(0, int(x_hi) + 1, 200)
    yt = y_ticks if y_ticks is not None else np.arange(0, int(y_hi) + 1, 200)
    ax.set_xticks(xt)
    ax.set_xticklabels([str(t) for t in xt], fontsize=18)
    ax.set_yticks(yt)
    ax.set_yticklabels([str(t) for t in yt], fontsize=18)

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.grid(False)
    plt.title(celltype if celltype is not None else 'all cells', fontsize=18)
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
    y_min: float = None,
    x_ticks: list = None,
    y_ticks: list = None,
    max_cells: int = None,
    random_state: int = 0,
    figsize: tuple = (3, 3),
    point_size: float = 2,
    gate_alpha: float = 0.25,
    vline: float = None,
    hline: float = None,
    gate_palette=None,
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
        adata = adata[adata.obs[celltype_col] == celltype]

    # Clip x and y axes
    if x_clip:
        adata = adata[adata.obs[x_axis] <= x_clip, :]
    if y_clip:
        adata = adata[adata.obs[y_axis] <= y_clip, :]

    sub_adata = adata.copy()

    # Downsample before the KDE if there are too many cells. gaussian_kde
    # evaluates an O(n^2) pairwise density, so 300k+ cells is very slow and
    # memory-heavy; cap the number of points used for the plot.
    if max_cells is not None and sub_adata.n_obs > max_cells:
        print(f'Subsampling {sub_adata.n_obs} -> {max_cells} cells for the density')
        sub_adata = sc.pp.subsample(sub_adata, n_obs=max_cells, copy=True,
                                    random_state=random_state)

    # effective upper axis bounds: fall back to the data max when no clip is
    # requested (x_clip / y_clip = None) so the limits, ticks and gate boxes
    # still have a finite upper bound.
    xv = sub_adata.obs[x_axis].to_numpy(dtype=float)
    yv = sub_adata.obs[y_axis].to_numpy(dtype=float)
    x_hi = x_clip if x_clip is not None else (float(np.nanmax(xv)) if len(xv) else 1.0)
    y_hi = y_clip if y_clip is not None else (float(np.nanmax(yv)) if len(yv) else 1.0)

    if signature_col:
        gene_expr = sub_adata.obs[signature_col]
        title = signature_col
    else:
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
    )

    # Add vline and hline if specified
    if vline is not None:
        # ax.axvline(x=vline, color='black', linestyle='--', linewidth=2)
        ax.plot([vline, vline], [hline, y_hi], color='black', linestyle='--', linewidth=2)
    if hline is not None:
        ax.axhline(y=hline, color='black', linestyle='--', linewidth=2)

    # axis limits with a small relative pad
    xpad = 0.02 * x_hi
    ypad = 0.02 * y_hi
    xlim = (-xpad, x_hi + xpad)
    # lower y limit: explicit y_min if given, else the padded data minimum
    if y_min is not None:
        y_lo = y_min
    else:
        y_lo = (float(np.nanmin(yv)) if len(yv) else 0.0) - ypad
    ylim = (y_lo, y_hi + ypad)

    # overlay gates (same shaded-box drawing as plot_celltype_density_2d); gates
    # outside the view limits are omitted / clipped to the axis
    if gates:
        _draw_gates(ax, gates, x_hi, y_hi,
                    palette=gate_palette, gate_alpha=gate_alpha,
                    xlim=xlim, ylim=ylim)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    # ticks: the values in x_ticks / y_ticks if given, else every 200 um
    xt = x_ticks if x_ticks is not None else np.arange(0, int(x_hi) + 1, 200)
    yt = y_ticks if y_ticks is not None else np.arange(0, int(y_hi) + 1, 200)
    ax.set_xticks(xt)
    ax.set_xticklabels([str(t) for t in xt], fontsize=18)
    ax.set_yticks(yt)
    ax.set_yticklabels([str(t) for t in yt], fontsize=18)

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.grid(False)

    plt.title(gene_name, fontsize=18)
    fig.tight_layout()

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
