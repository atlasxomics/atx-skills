"""atx_figures — report figure generators built on atx_plots_lib + atx_project.

Each generator returns (figure, dataframe): the figure for the slide, the table behind it.
"""
from __future__ import annotations

from typing import List, Optional, Tuple

import os
import warnings

import numpy as np
import pandas as pd
import plotly.express as px

import atx_plots_lib as L
from atx_project import LATCH


def qc_violin(adata, metric: str, group_by: str = "SampleName",
              plot_type: str = "violin") -> Tuple[object, pd.DataFrame]:
    """The QC distribution panels on the 'Samples pass QCs for each assay' slide.

    Reproduces plots_epigenome/violin_plots/violin_plot.py exactly: same
    `create_violin_data` input, same `sort_group_categories` ordering, same
    `build_discrete_color_map` over DEFAULT_H5_CATEGORICAL_PALETTE, same px.violin
    construction and layout. Title is "Distribution of <metric> by <group_by>", matching
    the delivered Demo_002 deck (log10_nFrags and TSSEnrichment by SampleName).

    plot_type: "violin" (box=True inside, as delivered) or "box" (the cell's fast path).
    """
    df = L.create_violin_data(adata, group_by, metric, data_type="obs")
    categories = L.sort_group_categories(df["group"].unique().tolist())
    color_map = L.build_discrete_color_map(categories, L.DEFAULT_H5_CATEGORICAL_PALETTE)

    kwargs = dict(x="group", y="value", points=False, color="group",
                  color_discrete_map=cmap, category_orders={"group": categories})
    if plot_type == "violin":
        fig = px.violin(df, box=True, **kwargs)
    elif plot_type == "box":
        fig = px.box(df, **kwargs)
    else:
        raise ValueError("plot_type must be 'violin' or 'box'")

    fig.update_layout(
        title=f"Distribution of {metric} by {group_by}",
        xaxis_title=group_by,
        yaxis_title=metric,
        plot_bgcolor="rgba(0,0,0,0)",
        showlegend=False,
    )
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor="lightgrey")
    fig.update_xaxes(categoryorder="array", categoryarray=categories)
    return fig, df


def qc_violin_panels(project, metrics: Optional[List[str]] = None,
                     group_by: str = "SampleName") -> dict:
    """Both panels as delivered: log10_nFrags and TSSEnrichment, grouped by SampleName."""
    if metrics is None:
        metrics = ["log10_nFrags", "TSSEnrichment"]
    adata = project.cell_adata()
    out = {}
    for m in metrics:
        fig, df = qc_violin(adata, m, group_by=group_by)
        out[m] = {"figure": fig, "data": df,
                  "medians": df.groupby("group", observed=True)["value"].median().to_dict()}
    return out


# ---------------------------------------------------------------- cluster maps

def cluster_maps(adata, samples=None, color_by="cluster", pt_size=3,
                 coords="spatial", flipY=True, color_scheme="bright",
                 vmin=None, vmax=None, palette="h5", color_map=None,
                 panel_labels=None, legend_header=False, equal_aspect=True):
    """Spatial maps per sample (catalog #1) or spatial QC maps (catalog #10).

    Thin wrapper over plots_epigenome's `plot_umap_for_samples`. Pass color_by="cluster"
    for the cluster slide, or a QC obs column ("TSSEnrichment", "log10_nFrags") for the
    QC slide's spatial-map style, which is what Demo_001 p14 and Demo_003 p13 use.

    NOTE: the delivered panels were rendered by Latch's proprietary w_h5 widget (the
    `Mean: x` / `Genes: y` colorbar annotation). This reproduces the same data with Plotly
    cosmetics; it is not a pixel match and cannot be.
    """
    if samples is None:
        samples = L.sort_group_categories(adata.obs["sample"].unique().tolist())
    if coords != "spatial" and flipY:
        raise ValueError(
            f"cluster_maps(coords={coords!r}, flipY=True): the ported plot_umap_for_samples "
            "ignores `coords` when flipY is true -- it rebuilds spatial_flippedY from "
            "obsm['spatial'] and would silently plot the unaligned frame. Pass flipY=False "
            "with coordinates that are already screen-oriented, e.g. from "
            "align_spatial_to_image().")
    if color_map:
        cats = adata.obs[color_by].astype(str).unique().tolist()
        ctx = use_palette(ordered_palette(color_map, cats))
    elif palette == "h5":
        ctx = use_palette()
    elif palette:
        ctx = use_palette(palette)
    else:
        ctx = _contextlib.nullcontext()
    with ctx:
        fig = L.plot_umap_for_samples(
            adata, list(samples), color_by=color_by, pt_size=pt_size, coords=coords,
            flipY=flipY, color_scheme=color_scheme, vmin=vmin, vmax=vmax,
        )
    if equal_aspect:
        # A tissue section is a physical object, so x and y must carry the same number of
        # microns per pixel. Plotly stretches each subplot to fill its box by default,
        # which made the QC maps noticeably more elongated than the same sections on the
        # brightfield. Anchoring y to its own x restores the true outline.
        fig.for_each_yaxis(
            lambda ax: ax.update(scaleanchor=(ax.anchor if ax.anchor and
                                              ax.anchor != "free" else "x"),
                                 scaleratio=1))
    if panel_labels:
        relabel_panels(fig, {str(k): v for k, v in panel_labels.items()})
    if legend_header:
        fig.update_layout(legend=dict(
            title=dict(text=observation_legend_title(color_by), font=dict(size=12)),
            itemsizing="constant", font=dict(size=11),
            y=0.5, yanchor="middle"))
    return fig, list(samples)


def spatial_qc_maps(adata, metrics=None, **kw):
    """The Demo_001/Demo_003 QC slide panels: one spatial map per QC metric."""
    if metrics is None:
        metrics = ["TSSEnrichment", "log10_nFrags"]
    out = {}
    for m in metrics:
        fig, samples = cluster_maps(adata, color_by=m, **kw)
        out[m] = {"figure": fig, "samples": samples,
                  "mean": float(adata.obs[m].mean()),
                  "min": float(adata.obs[m].min()), "max": float(adata.obs[m].max())}
    return out


# ------------------------------------------------------------- proportion plot

def cluster_proportion(adata, group_by="sample", stack_by="cluster",
                       return_type="proportion", color_map=None, group_labels=None,
                       font_pt=None):
    """'Cluster proportion by sample' (catalog #4).

    Reproduces plots_epigenome/proportion_plot/proportion_plot.py exactly. Demo_001's deck
    uses return_type="proportion" (y 0-1); Demo_002's uses "counts" (raw cell counts).
    """
    df = L.create_proportion_dataframe(adata, group_by, stack_by, return_type=return_type)
    stack_categories = natural_sort(df["stack_by"].unique().tolist())
    if color_map:
        cmap = {c: color_map[c] for c in stack_categories}
    else:
        cmap = L.build_discrete_color_map(stack_categories,
                                          L.DEFAULT_H5_CATEGORICAL_PALETTE)
    if group_labels:
        df["group_by"] = df["group_by"].astype(str).map(lambda g: group_labels.get(g, g))

    fig = px.bar(
        df, x="group_by", y="value", color="stack_by", barmode="stack",
        color_discrete_map=cmap,
        category_orders={"stack_by": stack_categories},
        title=f"Distribution of {stack_by} by {group_by}",
    )
    fig.update_layout(
        xaxis_title=group_by,
        yaxis_title="Proportion" if return_type == "proportion" else "Count",
        plot_bgcolor="rgba(0,0,0,0)", showlegend=True, legend_title=stack_by,
    )
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor="lightgrey")
    fig.update_xaxes(categoryorder="array",
                     categoryarray=natural_sort(df["group_by"].unique().tolist()))
    if font_pt:
        # Plots' own defaults are sized for a notebook; on a slide the sample labels
        # under the bars are the thing a reader needs first.
        fig.update_layout(font=dict(size=font_pt),
                          title=dict(font=dict(size=font_pt + 3)),
                          legend=dict(font=dict(size=font_pt - 2),
                                      title=dict(font=dict(size=font_pt - 1))))
        fig.update_xaxes(tickfont=dict(size=font_pt + 2))
    return fig, df


# --------------------------------------------------- neighborhood enrichment

def neighborhood_enrichment(adata, sample_key="sample", cluster_key="cluster",
                            spatial_key="spatial", vmin=-100.0, vmax=100.0):
    """'Neighborhood analysis by sample' (catalog #3) — one z-score heatmap per sample.

    The delivered slides show a per-sample grid on a shared -100..100 colour range, so the
    enrichment is computed per sample rather than once over the pooled object.
    Uses `squidpy_analysis` then `plotly_heatmap`, both from the Plots code.
    """
    uns_key = f"{cluster_key}_nhood_enrichment"
    out = {}
    for s in L.sort_group_categories(adata.obs[sample_key].unique().tolist()):
        sub = adata[adata.obs[sample_key] == s].copy()
        sub.obs[cluster_key] = sub.obs[cluster_key].astype("category")
        L.squidpy_analysis(sub, cluster_key=cluster_key, sample_key=None,
                           spatial_key=spatial_key)
        fig, data = L.plotly_heatmap(sub, uns_key, key=cluster_key, title=str(s),
                                     mode="zscore", vmin=vmin, vmax=vmax)
        out[s] = {"figure": fig, "matrix": data,
                  # axis labels as plotly_heatmap ordered them, not the raw categories:
                  # it re-sorts the matrix, so reading them back off the trace is the
                  # only way to be sure labels and rows still correspond.
                  "labels": [str(v) for v in fig.data[0].x],
                  "clusters": list(sub.obs[cluster_key].cat.categories)}
    return out


def neighborhood_grid(nb, labels=None, cols=2, vmin=-100.0, vmax=100.0,
                      colorscale="RdBu_r", title_pt=36, tick_pt=11,
                      cell_px=340, colorbar_title="zscore"):
    """Compose the per-sample neighborhood heatmaps into ONE figure.

    Montaging the individual PNGs gives every panel its own colourbar and its own
    (identical) scale, which wastes a third of the width and invites the reader to
    compare panels that are only coincidentally comparable. Here the panels share a
    single `coloraxis`, so one bar serves all three and the shared range is enforced by
    the figure rather than by convention. Panel titles are sized for a slide, not for a
    notebook.

    Font sizes are in LAYOUT pixels and the canvas is deliberately small
    (`cell_px=340`), because a picture placed on a slide is scaled by
    slide_width / layout_width: a 1200px-wide figure dropped into a 6in column
    shrinks every label by a factor of three. Keep the canvas near its final size
    and set the type large, rather than rendering big and letting placement shrink it.

    `nb` is the dict from `neighborhood_enrichment`; `labels` maps sample -> panel title.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    samples = list(nb)
    rows = -(-len(samples) // cols)
    titles = [(labels or {}).get(s, str(s)) for s in samples]
    fig = make_subplots(rows=rows, cols=cols, subplot_titles=titles,
                        horizontal_spacing=0.10, vertical_spacing=0.12)
    for k, s in enumerate(samples):
        r, cc = divmod(k, cols)
        lab = nb[s].get("labels") or [str(v) for v in nb[s]["clusters"]]
        fig.add_trace(go.Heatmap(z=nb[s]["matrix"], x=lab, y=lab,
                                 coloraxis="coloraxis", hoverongaps=False),
                      row=r + 1, col=cc + 1)
        fig.update_xaxes(showgrid=False, tickfont=dict(size=tick_pt),
                         row=r + 1, col=cc + 1)
        fig.update_yaxes(showgrid=False, autorange="reversed",
                         tickfont=dict(size=tick_pt), row=r + 1, col=cc + 1)
    fig.update_layout(
        coloraxis=dict(colorscale=colorscale, cmin=vmin, cmax=vmax,
                       colorbar=dict(title=colorbar_title, tickformat=".2f", len=0.55,
                                     y=0.5, yanchor="middle", thickness=18,
                                     tickfont=dict(size=tick_pt + 1))),
        width=cols * cell_px + 170, height=rows * cell_px + 60,
        margin=dict(l=55, r=20, t=50, b=40), showlegend=False,
        plot_bgcolor="white", paper_bgcolor="white")
    for a in fig.layout.annotations:
        a.font.size = title_pt
    return fig


# -------------------------------------------------------------- peak barplots

#: Approximation of the ArchR peak-annotation palette read off the delivered slides.
#: The exact colours are not defined in any repository available here.
PEAK_ANNO_COLORS = {"Distal": "#4CAF50", "Exonic": "#64B5F6",
                    "Intronic": "#5E35B1", "Promoter": "#FFC107"}


def peaks_per_group(project, level="cluster"):
    """'Total number of peaks vary across clusters' (catalog #7).

    Built from the ArchR peak-call summary the workflow already delivers, so the bar
    heights and the `C1(n = ...)` tick labels are the pipeline's own numbers, not recomputed.
    UnionPeaks is placed last, as in the decks.
    """
    df = project.peak_call_summary(level)
    groups = [g for g in L.sort_group_categories(
        df.loc[~df.is_union, "Group"].unique().tolist())] + ["UnionPeaks"]
    fig = px.bar(
        df, x="Group", y="peaks_thousands", color="peak_anno",
        color_discrete_map=PEAK_ANNO_COLORS,
        category_orders={"Group": groups,
                         "peak_anno": ["Distal", "Exonic", "Intronic", "Promoter"]},
        title="Peaks per Cluster" if level == "cluster" else "Peaks per Sample",
    )
    fig.update_layout(
        xaxis_title="", yaxis_title="Number of Peaks (x10^3)",
        plot_bgcolor="rgba(0,0,0,0)", legend_title="PeakAnno",
    )
    fig.update_xaxes(showgrid=False, tickangle=90)
    fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor="lightgrey")
    return fig, df


# -------------------------------------------- already-delivered ArchR figures

#: Catalog class-A figures: the secondary-analysis workflow already wrote these, so the
#: right action is to locate the file, not to re-implement ArchR in Python.
WORKFLOW_FIGURES = {
    "marker_heatmaps": "heatmaps_all.pdf",              # catalog #5
    "neighborhoods_all": "all_neighborhoods.pdf",
    "peak_call_summary_cluster_pdf": None,              # Clusters_Peak-Call-Summary_Date-*
    "ripleys_L": "ripleys_L.pdf",
}


def workflow_figure(project, name):
    """Fetch a figure the workflow already produced, by catalog name.

    `volcano` figures are per-condition (volcano_plots_<condition>.pdf), so pass the exact
    filename instead for those.
    """
    fname = WORKFLOW_FIGURES.get(name, name)
    if fname is None:
        raise ValueError(f"{name} is date-stamped; pass the exact filename from figures/")
    import os
    return project._cp(f"{project.latch_root}/figures/{fname}",
                       os.path.join(project.local_root, "figures", fname))


# ------------------------------------------------------- morphology overlays

import re as _re


class MorphologyImageMissing(FileNotFoundError):
    """Raised when a run has no *_fix image. Report this to the scientist."""


def short_run_id(run_id):
    """RUN1_SEQ001 -> RUN1. Image folders drop the leading zeros (D02778 -> D2778)."""
    head = str(run_id).split("_")[0]
    m = _re.match(r"^([A-Za-z]+)0*(\d+)$", head)
    return f"{m.group(1)}{int(m.group(2))}" if m else head


#: Searched in order; the first that has <short>_fix.* wins.
IMAGE_ROOTS = ("latch://13502.account/spatials", "latch://atx-illumina.mount/Images")


def fetch_morphology_image(project, run_id, roots=IMAGE_ROOTS):
    """Fetch the registered brightfield image for a run.

    Only files containing `_fix` are used — the other TIFs in these folders (_PostB,
    _PostB_BSA) are pre-registration captures and must not be substituted silently.
    Raises MorphologyImageMissing when no _fix image exists, which is a note for the
    scientist rather than a failure to work around.
    """
    import os, shlex, subprocess
    short = short_run_id(run_id)
    for root in roots:
        folder = f"{root}/{short}"
        from atx_project import latch_ls
        try:
            names = [n for n in latch_ls(folder) if "_fix" in n and not n.endswith("/")]
        except FileNotFoundError:
            continue
        if not names:
            continue
        name = sorted(names)[0]
        local = os.path.join(project.local_root, "images", name)
        return project._cp(f"{folder}/{name}", local), f"{folder}/{name}"
    raise MorphologyImageMissing(
        f"No *_fix image for run {run_id} (searched {', '.join(r + '/' + short for r in roots)}). "
        "Ask the scientist whether a registered image exists; do not substitute _PostB."
    )


def fetch_tissue_positions(project, run_id):
    """tissue_positions_list.csv for a run: barcode, in_tissue, row, col, y_px, x_px."""
    short = short_run_id(run_id)
    local = os.path.join(project.local_root, "spatial", f"{short}_tissue_positions_list.csv")
    path = project._cp(
        f"latch://13502.account/spatials/{short}/spatial/tissue_positions_list.csv", local)
    df = pd.read_csv(path, header=None,
                     names=["barcode", "in_tissue", "row", "col", "y_px", "x_px"])
    return df


def fetch_registered_image(project, run_id):
    """The image tixel coordinates are actually registered to: spatial/tissue_hires_image.png.

    tissue_positions_list.csv is in the AtlasXBrowser crop space (metadata.json `crop_area`,
    7686x7686 for RUN1), and tissue_hires_image.png is that same space scaled by
    `tissue_hires_scalef`. Overlay geometry is therefore exact, with no guessing.
    """
    import json, os
    short = short_run_id(run_id)
    base = f"latch://13502.account/spatials/{short}/spatial"
    img = project._cp(f"{base}/tissue_hires_image.png",
                      os.path.join(project.local_root, "spatial", f"{short}_hires.png"))
    sf_path = project._cp(f"{base}/scalefactors_json.json",
                          os.path.join(project.local_root, "spatial", f"{short}_scalefactors.json"))
    return img, json.load(open(sf_path))


def check_fix_image(project, run_id):
    """Report whether a *_fix image exists and whether it is usable for overlay.

    The _fix TIFs are NOT in the tixel coordinate space: they differ in size per run and
    match neither the crop_area extent nor a uniform scaling of it. Registering onto one
    needs a transform that is not stored anywhere in the run output, so this returns a
    note for the scientist instead of silently misaligning the overlay.
    """
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    try:
        path, remote = fetch_morphology_image(project, run_id)
    except MorphologyImageMissing as e:
        return {"fix_image": None, "note": str(e)}
    w, h = Image.open(path).size
    return {"fix_image": remote, "size": [w, h],
            "note": (f"_fix image found ({w}x{h}) but it is not in the tixel coordinate space; "
                     "the overlay was drawn on tissue_hires_image.png instead. Ask the "
                     "scientist for the _fix registration transform if the _fix image is required.")}


def morphology_overlay(project, adata, run_id, cluster_key="cluster",
                       alpha=0.75, point_size=4, cmap_scheme="bright", figsize=None,
                       color_map=None, crop_pad=0.04):
    """'Clustering aligns to morphology' (catalog #2): overlay + side-by-side.

    Drawn on the registered hires image so the tixels land where they belong. Both panels
    are cropped to the tixel bounding box (plus `crop_pad` of its larger side) because the
    capture area is several times the tissue: an uncropped panel is mostly empty slide,
    which is what the delivered decks crop away. Pass ``crop_pad=None`` for the full frame.
    Returns (figure, info); info carries a `fix_note` when a *_fix image exists but is
    unusable, and when none exists at all, and `crop_box` when a crop was applied.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None

    img_path, sf = fetch_registered_image(project, run_id)
    scalef = float(sf["tissue_hires_scalef"])
    img = np.asarray(Image.open(img_path).convert("L"))
    pos = fetch_tissue_positions(project, run_id)

    sub = adata.obs[adata.obs["sample"].astype(str) == str(run_id)]
    bc = pd.Index(sub.index).str.split("#").str[-1].str.replace(r"-\d+$", "", regex=True)
    joined = (pd.DataFrame({"barcode": bc, cluster_key: sub[cluster_key].astype(str).values})
              .merge(pos, on="barcode", how="inner"))

    cats = natural_sort(joined[cluster_key].unique().tolist())
    if color_map:
        cmap = {cv: color_map[cv] for cv in cats}
    else:
        cmap = L.build_discrete_color_map(
            cats, L.generate_color_palette(len(cats), cmap_scheme))

    px = joined["x_px"].to_numpy() * scalef
    py = joined["y_px"].to_numpy() * scalef
    crop = None
    if crop_pad is not None and len(px):
        x0, x1 = float(px.min()), float(px.max())
        y0, y1 = float(py.min()), float(py.max())
        pad = crop_pad * max(x1 - x0, y1 - y0)
        crop = (max(x0 - pad, 0.0), min(x1 + pad, float(img.shape[1])),
                max(y0 - pad, 0.0), min(y1 + pad, float(img.shape[0])))

    if figsize is None:
        if crop is not None:
            ph = 7.0
            pw = ph * (crop[1] - crop[0]) / (crop[3] - crop[2])
            figsize = (2 * pw + 1.3, ph)   # +legend column, or the panels float apart
        else:
            figsize = (14, 7)
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    for ax in axes:
        ax.imshow(img, cmap="gray"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        if crop is not None:
            cx0, cx1, cy0, cy1 = crop
            ax.set_xlim(cx0, cx1)
            ax.set_ylim(cy1, cy0)          # image rows run downwards
            ax.set_box_aspect((cy1 - cy0) / (cx1 - cx0))
    axes[0].scatter(px, py, s=point_size,
                    c=[cmap[v] for v in joined[cluster_key]], alpha=alpha, linewidths=0)
    axes[0].set_title(f"Clusters on morphology — {run_id}")
    axes[1].set_title("Brightfield")
    axes[0].legend(handles=[Line2D([], [], marker="o", linestyle="", markersize=5,
                                   color=cmap[cv], label=cv) for cv in cats],
                   loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=7)
    fig.tight_layout()
    info = {"run_id": run_id, "short_id": short_run_id(run_id),
            "registered_image": img_path, "hires_scalef": scalef,
            "image_shape": [int(img.shape[0]), int(img.shape[1])],
            "tixels_plotted": int(len(joined)), "clusters": cats,
            "crop_box": None if crop is None else [round(v, 1) for v in crop]}
    info.update({f"fix_{k}": v for k, v in check_fix_image(project, run_id).items()})
    return fig, info


# ----------------------------------------------- gene-activity / motif panels


def _backed_column(a, feature):
    """Read one feature column out of a backed AnnData as a dense 1-D array.

    `a[:, f].X` on a backed object does not materialise the column (it yields a length-1
    object array), so the view is pulled into memory first — one column, not the matrix.
    """
    from scipy import sparse as _sp
    sub = a[:, [feature]]
    sub = sub.to_memory() if hasattr(sub, "to_memory") else sub
    x = sub.X
    if _sp.issparse(x):
        x = x.toarray()
    return np.asarray(x).ravel()


def add_features(project, adata_cells, features, kind="ge"):
    """Attach feature columns to the light cell AnnData, reading from the backed object.

    Only the named columns are pulled off disk, so this costs megabytes rather than the
    18 GB the full matrix would. Values are aligned on the barcode index, never by position.
    Reads combined_ge.h5ad / combined_motifs.h5ad — never the float16 *_sm_* objects.
    """
    a = project.feature_adata(kind, backed=True)
    idx = pd.Index(a.obs_names)
    missing = [f for f in features if f not in set(a.var_names)]
    for f in features:
        if f in missing:
            continue
        vals = _backed_column(a, f)
        if len(vals) != a.n_obs:
            raise ValueError(f"{f}: read {len(vals)} values for {a.n_obs} cells")
        adata_cells.obs[f] = pd.Series(vals, index=idx).reindex(adata_cells.obs_names).to_numpy()
    return adata_cells, missing


def feature_panels(adata_cells, features, samples=None, pt_size=3,
                   color_scheme="viridis", include_umap=True):
    """Gene-activity (#6) and motif (#8) panels.

    One spatial figure per feature across samples, plus the matching UMAP, mirroring the
    'Gene Activity Scores identify distinct biological clusters' slides. Colour limits are
    set per feature from its own range, as the H5 viewer does.
    """
    if samples is None:
        samples = L.sort_group_categories(adata_cells.obs["sample"].unique().tolist())
    out = {}
    for f in features:
        if f not in adata_cells.obs:
            out[f] = {"status": "missing"}
            continue
        v = adata_cells.obs[f].astype(float)
        vmin, vmax = float(np.nanmin(v)), float(np.nanmax(v))
        entry = {"status": "ok", "min": vmin, "max": vmax, "mean": float(np.nanmean(v))}
        entry["spatial"] = L.plot_umap_for_samples(
            adata_cells, list(samples), color_by=f, pt_size=pt_size, coords="spatial",
            flipY=True, color_scheme=color_scheme, vmin=vmin, vmax=vmax)
        if include_umap:
            entry["umap"] = L.plot_umap_for_samples(
                adata_cells, list(samples), color_by=f, pt_size=pt_size, coords="X_umap",
                flipY=False, color_scheme=color_scheme, vmin=vmin, vmax=vmax)
        out[f] = entry
    return out


# ------------------------------------------------- palette and UMAP overview

import contextlib as _contextlib

#: The delivered reports are coloured by the Plots H5 viewer, whose categorical palette is
#: DEFAULT_H5_CATEGORICAL_PALETTE. Legend swatches sampled from the Demo_001 deck confirm it
#: (#292F61 vs #282E66, #F6E76B vs #FAE64D, #723582 vs #7E2F8A, #AF423A vs #C33530).
#: `generate_color_palette` cannot produce it — it maps matplotlib colormap names — so the
#: ported functions are given the palette explicitly rather than being edited.
H5_PALETTE = L.DEFAULT_H5_CATEGORICAL_PALETTE


@_contextlib.contextmanager
def use_palette(colors=None):
    """Make the ported plotting functions use an explicit discrete palette.

    Patching the lookup for the duration of a call keeps atx_plots_lib byte-identical to
    its upstream source, which matters because that file is a verbatim extraction.
    """
    if colors is None:
        colors = H5_PALETTE
    original = L.generate_color_palette

    def _explicit(length, scheme=None, *a, **kw):
        if length <= len(colors):
            return list(colors[:length])
        reps = -(-length // len(colors))
        return (list(colors) * reps)[:length]

    L.generate_color_palette = _explicit
    try:
        yield
    finally:
        L.generate_color_palette = original


def umap_overview(adata, color_by="cluster", palette=None, point_size=2.2,
                  width=700, height=620, legend_title=None, color_map=None,
                  category_order=None):
    """The single combined UMAP in the right column of the clustering slide.

    plot_umap_for_samples draws one panel per sample; the report's right column instead
    shows every cell in one UMAP, coloured by cluster (top) and by sample (bottom).
    """
    if palette is None:
        palette = H5_PALETTE
    import plotly.graph_objects as go

    if color_by not in adata.obs:
        raise KeyError(f"{color_by} not in obs; available: {list(adata.obs.columns)[:12]}")
    xy = adata.obsm["X_umap"]
    cats = natural_sort(adata.obs[color_by].astype(str).unique().tolist())
    if category_order:
        cats = [c for c in category_order if c in set(cats)]
    if color_map:
        cmap = {c: color_map[c] for c in cats}
    else:
        pal = (list(palette) if len(cats) <= len(palette)
               else list(palette) * (-(-len(cats) // len(palette))))
        cmap = L.build_discrete_color_map(cats, pal[:len(cats)])
    fig = go.Figure()
    vals = adata.obs[color_by].astype(str).to_numpy()
    for cv in cats:
        m = vals == cv
        fig.add_trace(go.Scattergl(
            x=xy[m, 0], y=xy[m, 1], mode="markers", name=str(cv),
            marker=dict(size=point_size, color=cmap[cv]), showlegend=True,
        ))
    fig.update_layout(
        width=width, height=height, plot_bgcolor="white",
        legend=dict(title=dict(text=observation_legend_title(legend_title or color_by)),
                    itemsizing="constant",
                    font=dict(size=10), y=0.5, yanchor="middle"),
        margin=dict(l=10, r=10, t=10, b=10),
    )
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    return fig, cats


#: Per-cluster colours read off the delivered Demo_001 slide 2 legend and snapped to
#: the nearest DEFAULT_H5_CATEGORICAL_PALETTE entry (a clean bijection onto 15 of its 24
#: colours, so the deck does use this palette). The H5 viewer assigns colours in its own
#: internal category order, which is not stored anywhere in the project — so an exact
#: per-cluster match for a NEW project cannot be derived, only a palette-consistent one.
#: Supply a mapping like this when a delivered deck must be matched colour for colour.
DEMO_CLUSTER_COLORS = {
    "C1": "#96C56C", "C2": "#9584B9", "C3": "#E38180", "C4": "#E48341", "C5": "#C33530",
    "C6": "#B570A8", "C7": "#9FD3E2", "C8": "#8E9ECD", "C9": "#E0C3DA", "C10": "#D0A970",
    "C11": "#FAE64D", "C12": "#7E2F8A", "C13": "#63B9A8", "C14": "#282E66", "C15": "#43884A",
}


def ordered_palette(color_map, categories):
    """Order an explicit {category: colour} map the way the ported plotters consume it.

    `plot_umap_for_samples` does `sorted(unique())` then zips colours positionally, so a
    palette listed in that same lexicographic order (C1, C10, C11, ... C2, C3) reproduces
    any desired per-category assignment without touching the ported code.
    """
    cats = sorted(str(c) for c in categories)
    missing = [c for c in cats if c not in color_map]
    if missing:
        raise KeyError(f"no colour supplied for {missing}")
    return [color_map[c] for c in cats]


def natural_sort(values):
    """Sort labels so C1, C2, ... C10, C15 come out in numeric order.

    The ported `sort_group_categories` is lexicographic (C1, C10, C11, ... C2), which is
    the order `plot_umap_for_samples` assigns colours in. The delivered decks DISPLAY the
    legend numerically, so draw/legend order and colour-assignment order differ; keep the
    two concerns separate rather than changing the ported function.
    """
    def key(v):
        parts = _re.split(r"(\d+)", str(v))
        return [int(p) if p.isdigit() else p for p in parts]
    return sorted((str(v) for v in values), key=key)


def observation_legend_title(name):
    """The two-line legend header the delivered decks use: grey 'Observation:' over the
    bold obs column name."""
    return f'<span style="font-size:9px;color:#808080">Observation:</span><br><b>{name}</b>'


def relabel_panels(fig, labels):
    """Rewrite per-sample subplot titles, e.g. RUN1_SEQ001 -> 'RUN1|CondA'.

    `labels` maps the sample id used in obs to the label wanted on the panel.
    """
    anns = list(fig.layout.annotations or ())
    for a in anns:
        if a.text in labels:
            a.text = labels[a.text]
    fig.layout.annotations = anns
    return fig


def deck_panel_labels(adata, samples, condition_col="condition"):
    """Build {sample: 'RUN1|CondA'} from obs, using the short run id and its condition."""
    out = {}
    for s in samples:
        sub = adata.obs[adata.obs["sample"].astype(str) == str(s)]
        cond = ""
        if condition_col in sub.columns and len(sub):
            vals = [v for v in sub[condition_col].astype(str).unique() if v and v != "nan"]
            cond = vals[0] if vals else ""
        out[str(s)] = f"{short_run_id(s)}|{cond}" if cond else short_run_id(s)
    return out


def align_spatial_to_image(project, adata, samples=None, coords="spatial",
                           out_key="spatial_img"):
    """Re-orient each sample's spatial coordinates into its morphology-image frame.

    The ArchR spatial table and the AtlasXBrowser `tissue_positions_list.csv` are both
    barcode-indexed, but nothing guarantees they share axis order or sign -- so a spatial
    panel can come out mirrored or transposed relative to the registered brightfield
    image. The image is the ground truth the scientist checks alignment against, so the
    mapping is DERIVED per sample by correlating the two coordinate sets over shared
    barcodes, never assumed.

    Writes `adata.obsm[out_key]` already oriented for the screen -- x grows with image
    column, y is negated against image row -- so it must be drawn with **flipY=False**.
    That is not a style choice: the ported `plot_umap_for_samples` IGNORES its `coords`
    argument whenever flipY is true (it rebuilds `spatial_flippedY` from `obsm['spatial']`
    and overrides `coords`), so an aligned frame passed with flipY=True is silently
    discarded and the old orientation is plotted instead. `cluster_maps` raises on that
    combination rather than letting it through.

    Returns a per-sample report naming the axis mapping and signs actually used, so a
    surprising one is visible rather than silent.
    """
    if samples is None:
        samples = natural_sort(adata.obs["sample"].unique().tolist())
    xy = np.asarray(adata.obsm[coords], dtype=float)
    out = np.full_like(xy, np.nan)
    bc_all = (pd.Index(adata.obs.index).str.split("#").str[-1]
              .str.replace(r"-\d+$", "", regex=True))
    report = {}
    for s in samples:
        mask = (adata.obs["sample"].astype(str) == str(s)).to_numpy()
        pos = fetch_tissue_positions(project, s)[["barcode", "x_px", "y_px"]]
        sub = pd.DataFrame({"barcode": bc_all[mask],
                            "s1": xy[mask, 0], "s2": xy[mask, 1]})
        sub["_i"] = np.nonzero(mask)[0]
        j = sub.merge(pos, on="barcode", how="inner")
        if len(j) < 10:
            raise ValueError(f"{s}: only {len(j)} barcodes shared with tissue_positions")
        corr = {(a, b): float(np.corrcoef(j[a], j[b])[0, 1])
                for a in ("s1", "s2") for b in ("x_px", "y_px")}
        # assign each image axis the spatial axis it tracks most strongly
        if abs(corr[("s1", "x_px")]) + abs(corr[("s2", "y_px")]) >= \
           abs(corr[("s2", "x_px")]) + abs(corr[("s1", "y_px")]):
            ax_x, ax_y = "s1", "s2"
        else:
            ax_x, ax_y = "s2", "s1"
        # x grows with image column; y is NEGATED against image row so that the result is
        # screen-ready for plotly (whose y grows upward) and must be drawn with flipY=False
        sx = 1.0 if corr[(ax_x, "x_px")] >= 0 else -1.0
        sy = -1.0 if corr[(ax_y, "y_px")] >= 0 else 1.0
        idx = j["_i"].to_numpy()
        out[idx, 0] = sx * j[ax_x].to_numpy()
        out[idx, 1] = sy * j[ax_y].to_numpy()
        report[str(s)] = {
            "x_from": ax_x, "x_sign": sx, "y_from": ax_y, "y_sign": sy,
            "matched_barcodes": int(len(j)),
            "corr": {f"{a}~{b}": round(v, 3) for (a, b), v in corr.items()},
            "transposed": ax_x != "s1", "mirrored_x": sx < 0, "mirrored_y": sy < 0,
        }
    unmatched = int(np.isnan(out[:, 0]).sum())
    if unmatched:
        warnings.warn(f"align_spatial_to_image: {unmatched} cells had no tissue_positions "
                      f"row and will not be plotted from {out_key}")
    adata.obsm[out_key] = out
    return report


def assert_sample_coverage(samples, raise_on_gap=True, **groups):
    """Check that every per-sample figure family covers every sample.

    Call it with the enumerated sample list and one keyword per figure family, each
    holding whatever was produced -- a dict keyed by sample, a list of file paths, or a
    list of sample ids:

        assert_sample_coverage(samples,
                               neighborhood=nb,            # dict keyed by sample
                               morphology=morph_paths,     # paths containing the run id
                               spatial=panel_samples)

    Rendering one sample and shipping a deck that silently covers a third of the
    experiment is not something any figure will complain about, so this is the check that
    has to notice. Returns {family: [missing samples]}; raises by default when non-empty.
    """
    samples = [str(s) for s in samples]
    missing = {}
    for family, produced in groups.items():
        if produced is None:
            missing[family] = list(samples)
            continue
        if isinstance(produced, dict):
            have = {str(k) for k in produced}
        else:
            items = [str(p) for p in produced]
            have = {s for s in samples if any(s in it for it in items)}
        gap = [s for s in samples if s not in have]
        if gap:
            missing[family] = gap
    if missing and raise_on_gap:
        lines = [f"  {fam}: missing {gap}" for fam, gap in missing.items()]
        raise ValueError(
            f"per-sample figures do not cover all {len(samples)} samples:\n"
            + "\n".join(lines))
    return missing


# ------------------------------------------------- genome browser coverage tracks
#
# The delivered "Gene Activity Scores identify distinct biological clusters" slide shows
# per-cluster ATAC coverage over a gene locus. In the Plots app that panel is the Track
# Browser tab (plots_epigenome/track_browser/), which hands bigwigs to Latch's proprietary
# `w_igv` widget -- an igv.js wrapper. The widget cannot be driven headlessly, but the
# inputs are ordinary ArchR bigwigs, so the same figure is reproducible with pyBigWig and
# matplotlib. It is not a pixel match to an IGV screenshot; the data is the same.

COVERAGE_DIRS = {
    "cluster": "cluster_coverages",
    "sample": "sample_coverages",
    "condition": "condition_coverages",
}

#: ArchR writes coverage bigwigs as "<group>-TileSize-<n>-normMethod-<m>-ArchR.bw";
#: the group label is everything before the first hyphen.
def coverage_group_label(filename):
    return os.path.basename(filename).split("-")[0]


def fetch_coverages(project, group="cluster", members=None):
    """Download the ArchR coverage bigwigs for a grouping; returns {label: local path}.

    These are ~100 MB each (15 clusters is ~1.4 GB for Demo_001), so they are cached
    under the project cache and only fetched once. pyBigWig cannot read latch:// URLs, so
    they must be local -- there is no range-request shortcut.
    """
    import subprocess
    sub = COVERAGE_DIRS.get(group, group)
    dst = os.path.join(project.local_root, "coverages")
    os.makedirs(dst, exist_ok=True)
    listing = subprocess.run([LATCH, "ls", f"{project.latch_root}/{sub}"],
                             capture_output=True, text=True)
    if listing.returncode != 0:
        raise FileNotFoundError(f"no {sub} directory for {project.name}")
    names = [ln.split()[-1] for ln in listing.stdout.splitlines()
             if ln.strip().endswith(".bw")]
    if not names:
        raise FileNotFoundError(
            f"{sub} holds no .bw files. ArchR may have written bedgraphs (.gz) instead; "
            "the Track Browser accepts both but this reader needs bigwigs.")
    out = {}
    for n in names:
        label = coverage_group_label(n)
        if members is not None and label not in set(members):
            continue
        local = os.path.join(dst, n)
        if not os.path.exists(local):
            project._cp(f"{project.latch_root}/{sub}/{n}", local)
        out[label] = local
    return out


def coverage_matrix(bigwigs, chrom, start, end, nbins=1000):
    """Binned mean coverage per group over a locus -> (labels, array[n_groups, nbins]).

    Returns labels in natural order (C1, C2, ... C15) so track order matches the legend.
    """
    import pyBigWig
    labels = natural_sort(bigwigs)
    rows = []
    for lab in labels:
        bw = pyBigWig.open(bigwigs[lab])
        try:
            if chrom not in bw.chroms():
                alt = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
                if alt not in bw.chroms():
                    raise KeyError(f"{chrom} not in {os.path.basename(bigwigs[lab])}; "
                                   f"first chroms: {list(bw.chroms())[:4]}")
                chrom_use = alt
            else:
                chrom_use = chrom
            vals = bw.stats(chrom_use, int(start), int(end), type="mean", nBins=nbins)
        finally:
            bw.close()
        rows.append([0.0 if v is None else float(v) for v in vals])
    return labels, np.asarray(rows, dtype=float)


def coverage_tracks(bigwigs, locus, gene_model=None, color_map=None, palette=None,
                    nbins=1000, title=None, track_height=0.30, width=4.6,
                    share_scale=True, highlight=None):
    """Stacked per-group coverage tracks over a locus (the Track Browser panel).

    `locus` is (chrom, start, end); `gene_model` is the optional dict returned by a RefSeq
    lookup, with `exons` as [(start, end), ...], `strand` and `name`, drawn underneath as
    the delivered slide's "Refseq Curated" row.

    All tracks share one y scale by default, which is what makes the panel readable as a
    comparison -- the delivered slide prints that shared maximum at the top left of every
    track. `highlight` is an optional (first_label, last_label) pair boxed to mark the
    clusters being called out.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    chrom, start, end = locus
    labels, mat = coverage_matrix(bigwigs, chrom, start, end, nbins=nbins)
    if color_map:
        colors = [color_map[l] for l in labels]
    else:
        pal = list(palette or H5_PALETTE)
        colors = [pal[i % len(pal)] for i in range(len(labels))]

    n = len(labels)
    model_h = 0.55 if gene_model else 0.12
    fig, axes = plt.subplots(
        n + 1, 1, figsize=(width, n * track_height + model_h + (0.3 if title else 0.05)),
        gridspec_kw={"height_ratios": [1] * n + [model_h / track_height]})
    vmax = float(mat.max()) if share_scale else None
    x = np.linspace(start, end, mat.shape[1])

    for i, (lab, col) in enumerate(zip(labels, colors)):
        ax = axes[i]
        top = vmax if share_scale else float(mat[i].max())
        top = top or 1.0
        ax.fill_between(x, 0, mat[i], color=col, linewidth=0)
        ax.set_ylim(0, top)
        ax.set_xlim(start, end)
        ax.set_yticks([0, top])
        ax.set_yticklabels([0, f"{top:.2f}"], fontsize=4)
        ax.set_xticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0, pad=1)
        ax.text(0.012, 0.5, lab, transform=ax.transAxes, fontsize=4.5, va="center",
                ha="left", bbox=dict(boxstyle="square,pad=0.18", fc="white",
                                     ec="#444444", lw=0.4))

    axm = axes[-1]
    axm.set_xlim(start, end)
    axm.set_ylim(0, 1)
    axm.axis("off")
    if gene_model:
        y = 0.62
        axm.plot([gene_model["txStart"], gene_model["txEnd"]], [y, y],
                 color="#17418C", lw=0.7, zorder=1)
        for (es, ee) in gene_model["exons"]:
            axm.add_patch(Rectangle((es, y - 0.12), max(ee - es, (end - start) / 600),
                                    0.24, color="#17418C", zorder=2))
        axm.text((gene_model["txStart"] + gene_model["txEnd"]) / 2, y - 0.42,
                 gene_model.get("name", ""), ha="center", va="top", fontsize=4.5,
                 color="#17418C", style="italic")
        axm.text(start, y + 0.26, "Refseq Curated", fontsize=4, va="bottom",
                 bbox=dict(boxstyle="square,pad=0.18", fc="white", ec="#444444", lw=0.4))

    if title:
        fig.suptitle(title, fontsize=9, y=0.995)
    # a small hspace is required, not cosmetic: with touching axes the "0" of one track
    # lands on the shared maximum printed at the top of the one below it
    fig.subplots_adjust(hspace=0.45, left=0.17, right=0.995,
                        top=0.955 if title else 0.995, bottom=0.005)

    if highlight:
        # positions must be read AFTER subplots_adjust, or the box lands on whichever
        # tracks happened to occupy those coordinates in the pre-layout figure
        try:
            i0, i1 = labels.index(highlight[0]), labels.index(highlight[1])
        except ValueError:
            warnings.warn(f"highlight {highlight} not among tracks {labels}")
        else:
            lo, hi = min(i0, i1), max(i0, i1)
            x0, x1 = (highlight[2], highlight[3]) if len(highlight) > 3 else (start, end)
            fig.canvas.draw()
            pb = axes[hi].get_position()
            pt = axes[lo].get_position()
            fx0 = pb.x0 + pb.width * (x0 - start) / (end - start)
            fx1 = pb.x0 + pb.width * (x1 - start) / (end - start)
            fig.patches.append(Rectangle(
                (fx0, pb.y0), fx1 - fx0, pt.y1 - pb.y0, transform=fig.transFigure,
                fill=False, ec="#17418C", lw=1.1, zorder=10))
    info = {"locus": f"{chrom}:{start}-{end}", "tracks": labels,
            "shared_max": float(mat.max()), "nbins": nbins,
            "gene": (gene_model or {}).get("name")}
    return fig, info


UCSC_API = "https://api.genome.ucsc.edu"


def refseq_gene(gene, genome="mm39", track="ncbiRefSeqCurated", pad=5000, timeout=90):
    """Locus and exon structure for a gene, from the UCSC track the decks cite.

    The delivered track panels label their gene row "Refseq Curated", which is UCSC's
    `ncbiRefSeqCurated`. Two calls: `/search` to find the locus, then `/getData/track`
    for the transcripts in it. The longest transcript whose `name2` matches the gene is
    used, so a gene with several isoforms gets the full span rather than an arbitrary one.

    ALWAYS pass the genome build the data was aligned to -- read it from a coverage
    bigwig header, do not assume. Demo_001 is mm39, and the same query against mm10
    returns a plausible-looking locus that is simply the wrong place.

    Returns {"name", "chrom", "txStart", "txEnd", "strand", "exons", "locus"} where
    `locus` is the padded (chrom, start, end) tuple to plot.
    """
    import json as _json
    import urllib.parse
    import urllib.request

    def _get(path):
        with urllib.request.urlopen(UCSC_API + path, timeout=timeout) as fh:
            return _json.load(fh)

    q = urllib.parse.quote(str(gene))
    hits = _get(f"/search?search={q}&genome={genome}").get("positionMatches", [])
    pos = None
    for group in hits:
        for m in group.get("matches", []):
            label = (m.get("posName") or "").split(" ")[0]
            if label.lower() == str(gene).lower():
                pos = m["position"]
                break
        if pos:
            break
    if pos is None:
        raise LookupError(f"UCSC has no {genome} match named {gene!r}")
    chrom, span = pos.split(":")
    lo, hi = (int(v.replace(",", "")) for v in span.split("-"))

    data = _get(f"/getData/track?genome={genome};track={track};chrom={chrom};"
                f"start={max(0, lo - pad)};end={hi + pad}")
    rows = next((v for v in data.values() if isinstance(v, list)), [])
    cand = [r for r in rows if str(r.get("name2", "")).lower() == str(gene).lower()]
    if not cand:
        raise LookupError(f"{track} has no transcript named {gene!r} at {pos}")
    best = max(cand, key=lambda r: r["txEnd"] - r["txStart"])
    starts = [int(v) for v in str(best["exonStarts"]).rstrip(",").split(",") if v != ""]
    ends = [int(v) for v in str(best["exonEnds"]).rstrip(",").split(",") if v != ""]
    return {
        "name": best.get("name2", gene), "transcript": best.get("name"),
        "chrom": best["chrom"], "txStart": best["txStart"], "txEnd": best["txEnd"],
        "strand": best.get("strand", "+"), "exons": list(zip(starts, ends)),
        "genome": genome,
        "locus": (best["chrom"], max(0, best["txStart"] - pad), best["txEnd"] + pad),
    }


# --------------------------------------------- paired spatial + UMAP gene blocks
#
# The delivered gene-activity slide does NOT show spatial panels alone: each marker gene
# appears as its spatial maps across samples AND a UMAP coloured by that same gene, under
# one shared colour scale. The pair is what makes the slide readable -- the spatial maps
# show where the signal sits in the tissue, the UMAP shows which cluster carries it.
# Build this pairing by default for marker genes when composing a deck from scratch.

#: Continuous scale used by the delivered gene-activity panels: blue (low) through pale
#: yellow to red (high). Set explicitly rather than relying on a library default.
FEATURE_COLORSCALE = "RdYlBu_r"


def feature_umap(adata, feature, colorscale=FEATURE_COLORSCALE, point_size=2.0,
                 vmin=None, vmax=None, width=430, height=430, coloraxis=None):
    """One combined UMAP of every cell, coloured by a continuous feature.

    `feature_panels(include_umap=True)` returns a per-SAMPLE grid of UMAPs; the delivered
    slide shows a single UMAP over all cells, which is this.
    """
    import plotly.graph_objects as go
    if feature not in adata.obs:
        raise KeyError(f"{feature!r} not in obs; call add_features() first")
    xy = np.asarray(adata.obsm["X_umap"], dtype=float)
    v = adata.obs[feature].astype(float).to_numpy()
    fig = go.Figure(go.Scattergl(
        x=xy[:, 0], y=xy[:, 1], mode="markers",
        marker=dict(size=point_size, color=v, colorscale=colorscale,
                    cmin=float(np.nanmin(v)) if vmin is None else vmin,
                    cmax=float(np.nanmax(v)) if vmax is None else vmax,
                    showscale=coloraxis is None),
        hoverinfo="skip", showlegend=False))
    fig.update_layout(width=width, height=height, plot_bgcolor="white",
                      margin=dict(l=4, r=4, t=18, b=4))
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    return fig


def gene_activity_block(adata, feature, samples=None, sample_labels=None,
                        row_label=None, colorscale=FEATURE_COLORSCALE, pt_size=2.6,
                        umap_point_size=1.8, panel_w=240, height=300):
    """One marker gene as spatial panels per sample PLUS a UMAP, on a shared colour scale.

    This is the unit the delivered gene-activity slide is built from; four of them tile
    the slide. `row_label` is written down the left edge ("Cd3e: T cell") the way the
    delivered slide labels each row, and `sample_labels` lets the panels be titled by
    condition (CondA, CondB, CondC) rather than by run id, which is what that slide uses.

    Spatial coordinates come from `obsm["spatial_img"]` when present, so the panels keep
    the brightfield orientation used everywhere else in the deck.
    """
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    if feature not in adata.obs:
        raise KeyError(f"{feature!r} not in obs; call add_features() first")
    if samples is None:
        samples = natural_sort(adata.obs["sample"].unique().tolist())
    samples = [str(s) for s in samples]
    labels = [(sample_labels or {}).get(s, s) for s in samples] + ["UMAP"]

    coords = "spatial_img" if "spatial_img" in adata.obsm else "spatial"
    flip = coords == "spatial"
    v_all = adata.obs[feature].astype(float).to_numpy()
    vmin, vmax = float(np.nanmin(v_all)), float(np.nanmax(v_all))

    fig = make_subplots(rows=1, cols=len(samples) + 1, subplot_titles=labels,
                        horizontal_spacing=0.012)
    smp = adata.obs["sample"].astype(str).to_numpy()
    xy_sp = np.asarray(adata.obsm[coords], dtype=float)
    for i, s in enumerate(samples, start=1):
        m = smp == s
        y = -xy_sp[m, 1] if flip else xy_sp[m, 1]
        fig.add_trace(go.Scattergl(
            x=xy_sp[m, 0], y=y, mode="markers",
            marker=dict(size=pt_size, color=v_all[m], coloraxis="coloraxis"),
            hoverinfo="skip", showlegend=False), row=1, col=i)
    xy_um = np.asarray(adata.obsm["X_umap"], dtype=float)
    fig.add_trace(go.Scattergl(
        x=xy_um[:, 0], y=xy_um[:, 1], mode="markers",
        marker=dict(size=umap_point_size, color=v_all, coloraxis="coloraxis"),
        hoverinfo="skip", showlegend=False), row=1, col=len(samples) + 1)

    fig.update_layout(
        width=panel_w * (len(samples) + 1) + 90, height=height, plot_bgcolor="white",
        margin=dict(l=70 if row_label else 10, r=10, t=26, b=8),
        coloraxis=dict(colorscale=colorscale, cmin=vmin, cmax=vmax,
                       colorbar=dict(thickness=9, len=0.78, x=1.005, tickfont=dict(size=7),
                                     title=dict(text=f"Mean: {np.nanmean(v_all):.2f}<br>"
                                                     f"Genes:<br>{feature}",
                                                font=dict(size=7), side="bottom"))))
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False, scaleanchor=None)
    for a in fig.layout.annotations:
        a.font.size = 11
    if row_label:
        fig.add_annotation(text=row_label, textangle=-90, xref="paper", yref="paper",
                           x=-0.055, y=0.5, showarrow=False, font=dict(size=12))
    info = {"feature": feature, "min": vmin, "max": vmax,
            "mean": float(np.nanmean(v_all)), "samples": samples, "coords": coords}
    return fig, info


def top_clusters_for_features(adata, features, cluster_key="cluster", n=1):
    """Cluster(s) with the highest mean activity for each feature, and their UMAP centroid.

    Used to place the "Cd19+" / "Cd3e+" / "Hmox1+" callouts on the cluster UMAP from the
    data instead of by eye. The delivered slide marks a second Hmox1+ region by hand, so
    raise `n` or pass explicit positions when a gene is high in more than one place.
    """
    xy = np.asarray(adata.obsm["X_umap"], dtype=float)
    cl = adata.obs[cluster_key].astype(str).to_numpy()
    out = {}
    for f in features:
        if f not in adata.obs:
            continue
        v = adata.obs[f].astype(float).to_numpy()
        means = {c: float(np.nanmean(v[cl == c])) for c in np.unique(cl)}
        top = sorted(means, key=means.get, reverse=True)[:n]
        out[f] = [{"cluster": c, "mean": means[c],
                   "x": float(np.median(xy[cl == c, 0])),
                   "y": float(np.median(xy[cl == c, 1]))} for c in top]
    return out


def annotate_umap(fig, callouts, suffix="+", font_size=13, halo=True):
    """Add "<gene>+" text callouts at given UMAP positions (see top_clusters_for_features).

    `halo` puts a translucent white plate behind each label. Without it a callout placed
    at a cluster centroid can be invisible: dark text over a dense dark cluster (Cd19+
    over the navy C14 blob) simply disappears, and the label looks like it was never
    drawn.
    """
    for feature, spots in callouts.items():
        for sp in spots:
            fig.add_annotation(
                x=sp["x"], y=sp["y"], text=f"{feature}{suffix}", showarrow=False,
                font=dict(size=font_size, color="#111111"),
                xanchor="center", yanchor="middle",
                bgcolor="rgba(255,255,255,0.78)" if halo else None,
                borderpad=2 if halo else 0)
    return fig


# ------------------------------------------------- motif logos and Compare results
#
# The differential-motif slide is the one figure set on a Data Report whose NUMBERS are
# not computable here. MeanDiff and FDR come from ArchR::getMarkerFeatures, run as a
# Latch "Compare Workflow" launched from the Plots app, which writes
# latch:///compare_outs/<run>/motif_results/all_motifs.csv. Recomputing them in Python
# (a scanpy t-test, say) produces different numbers from the delivered deck, so FETCH the
# workflow output -- never re-derive it.

#: CisBP position-probability matrices ship as JSON with the Plots repo, one per genome.
#: There is NO mm39 file: the repo has hg38, mm10 and rn6 only. The CisBP mouse motif set
#: is the same for mm10 and mm39 (motif identity does not depend on the assembly), so an
#: mm39 project uses the mm10 file -- that substitution is deliberate and is reported.
SEQLOGO_GENOMES = {"hg38": "seqlogo_hg38.json", "mm10": "seqlogo_mm10.json",
                   "mm39": "seqlogo_mm10.json", "rn6": "seqlogo_rn6.json"}
SEQLOGO_BASES = ["A", "C", "G", "T"]

#: Commit the ported Plots code and the CisBP PPM files are pinned to.
PLOTS_EPIGENOME_COMMIT = "3d05a74f7683"


def motif_key(name):
    """Normalise a motif name to the CisBP JSON convention.

    The h5ad var_names use a hyphen before the index ("Ascl1-780") while the seqlogo JSON
    and the Compare results use an underscore ("Ascl1_780"). Every one of the 884 mouse
    entries uses the underscore form, so the mapping is a straight substitution of the
    LAST hyphen.
    """
    s = str(name)
    return s[::-1].replace("-", "_", 1)[::-1] if "-" in s else s


def load_seqlogos(genome="mm10", assets_dir=None):
    """Load the CisBP position-probability matrices for a genome -> {motif: DataFrame}."""
    import json as _json
    fname = SEQLOGO_GENOMES.get(genome)
    if fname is None:
        raise KeyError(f"no seqlogo set for {genome!r}; have {sorted(SEQLOGO_GENOMES)}")
    if genome == "mm39":
        warnings.warn("no mm39 CisBP seqlogo file exists; using the mm10 set "
                      "(same mouse motif collection)")
    root = assets_dir or os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets")
    os.makedirs(root, exist_ok=True)
    path = os.path.join(root, fname)
    if not os.path.exists(path):
        # These PPM files are 0.5-1.4 MB, past the per-file limit for skill assets, so
        # they are fetched once from the pinned public plots_epigenome commit and cached
        # rather than vendored. Pinning the commit keeps the motif set reproducible.
        import urllib.request
        url = (f"https://raw.githubusercontent.com/atlasxomics/plots_epigenome/"
               f"{PLOTS_EPIGENOME_COMMIT}/motif_logo/{fname}")
        try:
            with urllib.request.urlopen(url, timeout=120) as fh, open(path, "wb") as out:
                out.write(fh.read())
        except Exception as exc:
            raise FileNotFoundError(
                f"{path} not found and fetching {url} failed ({exc}). Copy {fname} from "
                f"plots_epigenome/motif_logo/ into {root} by hand.") from exc
    with open(path) as fh:
        raw = _json.load(fh)
    return {k: pd.DataFrame(v)[SEQLOGO_BASES].astype(float) for k, v in raw.items()}


def motif_logo(name, genome="mm10", logos=None, assets_dir=None,
               information_content=True, title=None):
    """Sequence logo for a CisBP motif, as shown beside each motif's spatial panels."""
    logos = logos if logos is not None else load_seqlogos(genome, assets_dir)
    key = motif_key(name)
    if key not in logos:
        near = [k for k in logos if k.split("_")[0].lower() == key.split("_")[0].lower()]
        raise KeyError(f"{name!r} (as {key!r}) not in the {genome} CisBP set"
                       + (f"; same TF: {near[:6]}" if near else ""))
    df = logos[key]
    fig = L.plot_motif_logo(df, title=title if title is not None else
                            f"{key} (len={len(df)})",
                            information_content=information_content)
    return fig, {"motif": key, "length": int(len(df)), "genome": genome}


def fetch_compare_results(project, run, feature="motif"):
    """Download an ArchR Compare-workflow result table from latch:///compare_outs/<run>.

    `feature` is "motif", "gene" or "peak". Returns the DataFrame of all features with
    columns group, group_name, seqnames, idx, name, FDR, MeanDiff, Pval -- the numbers
    the delivered slide plots.
    """
    sub = {"motif": "motif_results", "gene": "gene_results", "peak": "peak_results"}[feature]
    fname = f"all_{feature}s.csv"
    local = os.path.join(project.local_root, "compare", run, fname)
    os.makedirs(os.path.dirname(local), exist_ok=True)
    path = project._cp(f"latch:///compare_outs/{run}/{sub}/{fname}", local)
    return pd.read_csv(path)


def fetch_compare_annotation(project, run, column="CD19_high"):
    """The per-barcode group assignment the Compare run actually used.

    Reading this is better than re-deriving the selection from a threshold: it is the
    exact set of cells the delivered statistics were computed on, so the figure and the
    numbers cannot drift apart.
    """
    local = os.path.join(project.local_root, "compare", run, "UpdateClustName_by_barcode.csv")
    os.makedirs(os.path.dirname(local), exist_ok=True)
    path = project._cp(f"latch:///compare_outs/{run}/UpdateClustName_by_barcode.csv", local)
    df = pd.read_csv(path)
    return df


def motif_rank_plot(df, group_name="GroupA", y_col="MeanDiff", color_col="FDR",
                    n_labels=4, colorscale="PuBu_r", marker_size=6,
                    width=760, height=620, title=""):
    """Rank plot of differential motifs, matching the delivered slide's settings.

    These arguments are not arbitrary: plots_epigenome/compare/motif_results.py calls
    plot_ranked_feature_plotly with exactly rank-by MeanDiff, colour-by FDR, PuBu_r,
    4 labels and marker size 6.
    """
    sub = df[df["group_name"] == group_name] if "group_name" in df else df
    fig = L.plot_ranked_feature_plotly(
        sub, y_col=y_col, x_col=None, n_labels=n_labels, label_col="name",
        color_col=color_col, colorscale=colorscale, marker_size=marker_size,
        title=title, y_label=y_col)
    fig.update_layout(width=width, height=height)
    return fig, {"group": group_name, "n_motifs": int(len(sub)),
                 "y": y_col, "color": color_col}


def feature_violin(adata, feature, group_by, order=None, colors=None,
                   width=520, height=430, title=None, y_label=None):
    """Violin of a continuous obs feature split by a categorical column.

    Used for "<motif> by Cd19+ region", comparing the motif deviation between conditions
    within the selected population.
    """
    import plotly.graph_objects as go
    if feature not in adata.obs:
        raise KeyError(f"{feature!r} not in obs; call add_features() first")
    df = L.create_violin_data(adata, group_by, feature, data_type="obs").dropna()
    groups = order or natural_sort(df["group"].astype(str).unique().tolist())
    fig = go.Figure()
    for g in groups:
        vals = df.loc[df["group"].astype(str) == g, "value"].astype(float)
        fig.add_trace(go.Violin(y=vals, name=str(g), box_visible=True, meanline_visible=False,
                                line_color=(colors or {}).get(g), points=False,
                                showlegend=False))
    fig.update_layout(width=width, height=height, plot_bgcolor="white",
                      title=dict(text=title or f"{feature} by {group_by}", x=0.5,
                                 font=dict(size=15)),
                      yaxis_title=y_label or feature,
                      margin=dict(l=50, r=10, t=44, b=30))
    fig.update_yaxes(gridcolor="#EEEEEE", zeroline=False)
    stats = {g: {"n": int((df["group"].astype(str) == g).sum()),
                 "median": float(df.loc[df["group"].astype(str) == g, "value"].median())}
             for g in groups}
    return fig, {"feature": feature, "group_by": group_by, "groups": groups, "stats": stats}


#: Colours the delivered CD19_high map uses: the selected conditions in the H5 palette,
#: unselected cells in a neutral tint so they read as background rather than a category.
COMPARE_NULL_COLOR = "#E7DAD5"


def apply_compare_annotation(adata, ann, column="CD19_high", mapping=None,
                             null_label="NULL", barcode_col="cell_barcode",
                             value_col="UpdateClustName"):
    """Join a Compare run's per-barcode group assignment onto obs as a categorical column.

    `mapping` renames the workflow's GroupA/GroupB to meaningful labels, e.g.
    {"GroupA": "CondA", "GroupB": "CondB"} -- check which is which against the run name and
    the cell counts in the compare_config filename (`..._A<nA>_B<nB>_...`) rather than
    assuming.

    Cells outside the selection get `null_label`, so the whole tissue can still be drawn
    with the unselected area as background, which is what the delivered slide shows.
    Returns a counts dict.
    """
    m = dict(zip(ann[barcode_col].astype(str), ann[value_col].astype(str)))
    if mapping:
        m = {k: mapping.get(v, v) for k, v in m.items()}
    vals = pd.Index(adata.obs.index.astype(str)).map(lambda b: m.get(b, null_label))
    adata.obs[column] = pd.Categorical(vals)
    counts = adata.obs[column].value_counts().to_dict()
    matched = sum(v for k, v in counts.items() if k != null_label)
    if matched == 0:
        raise ValueError(
            f"no barcode in {barcode_col} matched the obs index. Compare writes the full "
            f"'<run_id>#<barcode>-1' form, e.g. {ann[barcode_col].iloc[0]!r}; obs index "
            f"looks like {adata.obs.index[0]!r}")
    return {str(k): int(v) for k, v in counts.items()}


def subset_to_groups(adata, column, groups):
    """Cells in the named groups of a categorical column, as a copy safe to plot."""
    keep = adata.obs[column].astype(str).isin([str(g) for g in groups]).to_numpy()
    if keep.sum() == 0:
        raise ValueError(f"no cells in {groups} of {column!r}")
    return adata[keep].copy()
