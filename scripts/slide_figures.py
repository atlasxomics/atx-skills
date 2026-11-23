"""Render every figure a Data Report deck needs, each at the exact inch size of its slide
slot, so labels land on the slide at the point size written here (12-18 pt).

    python slide_figures.py <project> [--spec deck_spec.json] [--one-shot] [--out .]

Runs on the box (it reads the h5ad objects, bigwigs and images). Writes `<project>__*.png`,
`samples_<project>.csv`, `qc_table_<project>.csv`, `marker_heatmap_p*.png` and
`slide_figures_report.json` into --out; `dense_deck.py` (local) turns that folder and the
same spec into the deck. Every figure runs under a guard: a failure is recorded in the
report and the next figure is drawn.

The spec is per project and holds every scientific choice. It accepts render_figures.py's
CONFIG keys too (genes, gene_rows, motifs, colors), so a one-shot config works unchanged.
All keys are optional:

    palette        {cluster: hex} for every cluster          default: H5 in viewer order
    sample_labels  {run_id: "RUN1|CondA"}                    default: deck_panel_labels
    sample_short   {run_id: "CondA"}                         default: text after "|"
    sample_colors  {run_id: hex}
    cell_types     {cluster: "putative label"}               adds the cell-type figures
    rows           {feature: {"kind": "ge"|"motif", "label": "...", "sublabel": "..."}}
                   one figure per feature; built from genes/gene_rows/motifs if absent
    track_genes    [3 genes] for the track-browser slide     default: first 3 rows genes
    track_sets     {"name": [genes]} extra 3-column track slides (no UMAP)
    he_images      {run_id: "file.png" | {"pdf": "file.pdf", "page": 1}}
    he_exclude     {run_id: [[x0, y0, x1, y1], ...]}  page fractions to ignore (titles,
                   scale bars, thumbnails on a slide-scanner export)
    genome_heatmap true -> median one-vs-rest log2FC per 5 Mb bin, per sample and cluster

The rules this follows are in references/deck.md ("Dense deck"). The short version: draw
at slot size, never let a label fall below ~10 pt, one palette everywhere, every spatial
panel at one physical scale (only when the chips match), and the H&E cropped to the same
aspect as the brightfield it sits beside.
"""
import matplotlib
matplotlib.use("Agg")
import argparse, json, os, re, shutil, traceback, urllib.request, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
from PIL import Image
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

import atx_figures as F
from atx_project import Project
from atx_slots import SLOTS

Image.MAX_IMAGE_PIXELS = None
HERE = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.size": 13, "font.family": "DejaVu Sans", "axes.titlesize": 14,
                     "axes.labelsize": 13, "xtick.labelsize": 12, "ytick.labelsize": 12,
                     "legend.fontsize": 13, "axes.spines.top": False, "axes.spines.right": False})

PEAK_COLORS = {"Distal": "#4CAF50", "Exonic": "#64B5F6", "Intronic": "#5E35B1", "Promoter": "#FFC107"}
SAMPLE_COLORS = ["#d62728", "#1f77b4", "#98df8a", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22"]
CYTO = {"gneg": "#ffffff", "gpos25": "#c8c8c8", "gpos50": "#969696", "gpos75": "#5a5a5a",
        "gpos100": "#141414", "acen": "#c0392b", "gvar": "#dcdcdc", "stalk": "#aaaaaa"}


# =========================================================================== geometry
def ax_in(fig, x, y, w, h):
    """An axes at (x, y, w, h) in INCHES from the figure's bottom-left."""
    FW, FH = fig.get_size_inches()
    return fig.add_axes([x / FW, y / FH, w / FW, h / FH])


def ftext(fig, x, y, s, **kw):
    FW, FH = fig.get_size_inches()
    return fig.text(x / FW, y / FH, s, **kw)


def save(fig, out_dir, proj, name):
    # never bbox_inches="tight": it changes the size, and the size is the contract with the slot
    fig.savefig(os.path.join(out_dir, f"{proj}__{name}.png"), dpi=300)
    plt.close(fig)


# =========================================================================== H&E
def he_image(spec_dir, src):
    """Load an H&E as an RGB array; a PDF page is rendered at 4x (pypdfium2)."""
    if isinstance(src, dict):
        import pypdfium2 as pdfium
        doc = pdfium.PdfDocument(os.path.join(spec_dir, src["pdf"]))
        return np.asarray(doc[int(src.get("page", 1)) - 1].render(scale=4).to_pil().convert("RGB"))
    return np.asarray(Image.open(os.path.join(spec_dir, src)).convert("RGB"))


def crop_he(A, aspect, exclude=(), pad=0.04):
    """Crop an H&E to its tissue, padded, then widen to `aspect` (w/h) on a white canvas.

    The brightfield panel beside it is cropped to the tixels; giving the H&E the same
    aspect makes the three morphology panels the same size. Slide-scanner exports carry
    large black margins and edge lines: only LARGE dark regions (>=1% of the image) are
    whitened, because small dark specks inside tissue are tissue, and whitening them punches
    holes. `exclude` (page fractions) removes titles, scale bars and overview thumbnails.
    """
    A = A.astype(float)
    H, W, _ = A.shape
    k = max(3, int(0.003 * max(H, W)))
    dark = ndi.binary_opening(A.max(2) < 35, iterations=3)
    lab, n = ndi.label(dark)
    if n:
        sz = ndi.sum(dark, lab, range(1, n + 1))
        dark = np.isin(lab, 1 + np.flatnonzero(sz >= 0.01 * H * W))
    kill = ndi.binary_dilation(dark, iterations=k)
    for x0, y0, x1, y1 in exclude:
        kill[int(y0 * H):int(y1 * H), int(x0 * W):int(x1 * W)] = True
    B = np.where(kill[..., None], 255.0, A)
    rgb = A / 255
    mx, mn = rgb.max(2), rgb.min(2)
    tissue = ((mx - mn) / (mx + 1e-6) > 0.25) & (mx < 0.95) & ~kill
    big = ndi.binary_closing(tissue, iterations=max(5, int(0.006 * max(H, W))))
    lab, n = ndi.label(big)
    if not n:
        raise ValueError("no tissue found in the H&E image")
    sz = ndi.sum(big, lab, range(1, n + 1))
    ys, xs = np.nonzero(np.isin(lab, 1 + np.flatnonzero(sz >= 0.02 * sz.max())))
    x0, x1, y0, y1 = xs.min(), xs.max(), ys.min(), ys.max()
    p = pad * max(x1 - x0, y1 - y0)
    x0, x1, y0, y1 = x0 - p, x1 + p, y0 - p, y1 + p
    cx, cy, w, h = (x0 + x1) / 2, (y0 + y1) / 2, x1 - x0, y1 - y0
    if w / h < aspect:
        w = h * aspect
    else:
        h = w / aspect
    box = [int(round(v)) for v in (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)]
    canvas = Image.new("RGB", (box[2] - box[0], box[3] - box[1]), "white")
    sb = (max(box[0], 0), max(box[1], 0), min(box[2], W), min(box[3], H))
    canvas.paste(Image.fromarray(B.astype(np.uint8)).crop(sb), (sb[0] - box[0], sb[1] - box[1]))
    return np.asarray(canvas)


# =========================================================================== context
class Ctx:
    """Everything the figure functions share: cells, coordinates, labels, palette."""

    def __init__(self, project, spec, spec_dir, out_dir, one_shot=False):
        self.proj, self.spec, self.spec_dir, self.out = project, spec, spec_dir, out_dir
        self.report = {"project": project, "figures": {}, "notes": {}, "assumptions": list(spec.get("assumptions", []))}
        self.p = Project(project, cache="/data/analysis_cache")
        c = self.cells = self.p.cell_adata()
        self.samples = F.natural_sort(c.obs["sample"].astype(str).unique().tolist())
        if one_shot and not (spec.get("genes") or spec.get("rows")):
            cfg, assumptions = F.auto_feature_config(self.p, c)
            spec.setdefault("genes", cfg["genes"]); spec.setdefault("gene_rows", cfg["gene_rows"])
            self.report["assumptions"] += assumptions
        # labels
        dl = F.deck_panel_labels(c, self.samples)
        self.lab = {s: spec.get("sample_labels", {}).get(s, dl.get(s, s)) for s in self.samples}
        self.short = {s: spec.get("sample_short", {}).get(s, self.lab[s].split("|")[-1]) for s in self.samples}
        self.scol = {s: spec.get("sample_colors", {}).get(s, SAMPLE_COLORS[i % len(SAMPLE_COLORS)])
                     for i, s in enumerate(self.samples)}
        # palette: spec palette > render_figures "colors" > H5 default in the viewer's order
        clusters = c.obs["cluster"].astype(str).unique()
        pal = spec.get("palette") or spec.get("colors")
        if isinstance(pal, str):
            pal = getattr(F, pal)
        self.cm = dict(pal) if pal else F.default_cluster_colors(clusters)
        missing = sorted(set(clusters) - set(self.cm))
        if missing:
            raise KeyError(f"palette has no colour for clusters {missing}")
        self.cl = F.natural_sort(clusters)
        self.cell_types = spec.get("cell_types", {})
        # geometry
        F.align_spatial_to_image(self.p, c, self.samples)
        self.xy = np.asarray(c.obsm["spatial_img"], float)
        self.um = np.asarray(c.obsm["X_umap"], float)
        self.smp = c.obs["sample"].astype(str).to_numpy()
        self.clu = c.obs["cluster"].astype(str).to_numpy()
        self.ok = np.isfinite(self.xy).all(1)
        self.ext = {}
        for s in self.samples:
            b = self.xy[(self.smp == s) & self.ok]
            self.ext[s] = (b[:, 0].min(), b[:, 0].max(), b[:, 1].min(), b[:, 1].max())
        big = max(self.samples, key=lambda s: ((self.smp == s) & self.ok).sum())
        b = self.xy[(self.smp == big) & self.ok]
        self.pitch = float(np.median(cKDTree(b).query(b, k=2)[0][:, 1]))
        self.shared, geoms = F.check_shared_scale(self.p, self.samples)
        self.report["notes"].update({"shared_scale": self.shared, "chips": geoms, "pitch": self.pitch,
                                     "palette": {k: self.cm[k] for k in self.cl},
                                     "sample_labels": self.lab})
        if not self.shared and len(self.samples) > 1:
            self.report["assumptions"].append(
                "Runs are on different chip types, so spatial panels are each drawn at their own "
                "scale; relative tissue size is not comparable between panels.")
        self.clu_col = np.array([self.cm[x] for x in self.clu])
        # features
        self.rows = self._rows()
        feats_ge = [f for f, r in self.rows.items() if r["kind"] == "ge"]
        feats_m = [f for f, r in self.rows.items() if r["kind"] != "ge"]
        tg = spec.get("track_genes") or feats_ge[:3]
        extra = [g for gs in spec.get("track_sets", {}).values() for g in gs]
        ge_needed = list(dict.fromkeys(feats_ge + tg + extra))
        if ge_needed:
            c, miss = F.add_features(self.p, c, ge_needed, kind="ge")
            self.report["notes"]["missing_genes"] = miss
        if feats_m:
            c, miss = F.add_features(self.p, c, feats_m, kind="motifs")
            self.report["notes"]["missing_motifs"] = miss
        self.track_genes = [g for g in tg if g in c.obs]
        c.obs["ln(nFrags)"] = c.obs["log10_nFrags"].astype(float)
        # the pipeline column called log10_nFrags is ln(nFrags); prove it before labelling it so
        if not np.allclose(np.exp(c.obs["ln(nFrags)"]), c.obs["n_fragment"].astype(float), rtol=1e-6):
            c.obs["ln(nFrags)"] = np.log(c.obs["n_fragment"].astype(float))
            self.report["notes"]["ln_nFrags"] = "recomputed from n_fragment"

    def _rows(self):
        sp = self.spec
        if sp.get("rows"):
            return {f: {"kind": r.get("kind", "ge"), "label": r.get("label", f),
                        "sublabel": r.get("sublabel", "")} for f, r in sp["rows"].items()}
        rows = {}
        for g in sp.get("genes", []):
            lab = sp.get("gene_rows", {}).get(g, g).replace("<br>", " ")
            head, _, tail = lab.partition(":")
            rows[g] = {"kind": "ge", "label": g, "sublabel": tail.strip()}
        for m in sp.get("motifs", []):
            rows[m] = {"kind": "motif", "label": f"{m.rsplit('-', 1)[0]} motif", "sublabel": ""}
        return rows

    def zscores(self, feat):
        m = self.cells.obs.groupby(self.clu)[feat].mean()
        return ((m - m.mean()) / m.std(ddof=0)).reindex(self.cl)

    def top_clusters(self, feat, zmin=1.0):
        z = self.zscores(feat)
        hi = [c for c in self.cl if z[c] >= zmin]
        return hi or [z.idxmax()]


# =========================================================================== spatial helpers
def spatial_row(C, fig, x0, yb, w_budget, h, title_h=0.32, min_w=1.15, gap=0.10, title_pt=14,
                center=False, labels=None):
    """One panel per sample, top-aligned, all at ONE scale when the chips match.

    Panels keep a minimum width so a tiny section still gets a readable title; the scale
    shrinks until everything fits. Returns ({sample: axes}, {sample: scale in/unit}, width).
    """
    labels = labels or C.lab
    S = C.samples
    xs = {s: C.ext[s][1] - C.ext[s][0] for s in S}
    ys = {s: C.ext[s][3] - C.ext[s][2] for s in S}
    hh = h - (title_h if labels is not False else 0)
    padu = 0.03 * max(ys.values())
    n = len(S)
    if C.shared or n == 1:
        sc0 = min(hh / (max(ys.values()) + 2 * padu),
                  (w_budget - gap * (n - 1)) / sum(x + 2 * padu for x in xs.values()))
        scales = {s: sc0 for s in S}
    else:
        scales = {s: hh / (ys[s] + 2 * padu) for s in S}
    for _ in range(500):
        widths = {s: max((xs[s] + 2 * padu) * scales[s], min_w) for s in S}
        if sum(widths.values()) + gap * (n - 1) <= w_budget + 1e-9:
            break
        scales = {s: v * 0.99 for s, v in scales.items()}
    total = sum(widths.values()) + gap * (n - 1)
    x = x0 + ((w_budget - total) / 2 if center else 0)
    FW = fig.get_size_inches()[0]
    axes = {}
    for s in S:
        aw, sc = widths[s], scales[s]
        ax = ax_in(fig, x, yb, aw, hh)
        cx = (C.ext[s][0] + C.ext[s][1]) / 2
        top = C.ext[s][3] + padu
        ax.set_xlim(cx - aw / sc / 2, cx + aw / sc / 2)
        ax.set_ylim(top - hh / sc, top)
        ax.set_axis_off()
        if labels is not False:
            tx = min(max(x + aw / 2, x0 + 0.06 * title_pt), FW - 0.06 * title_pt)
            ftext(fig, tx, yb + hh + 0.05, labels[s], ha="center", va="bottom", fontsize=title_pt)
        axes[s] = ax
        x += aw + gap
    return axes, scales, total


def spatial_points(C, ax, s, scale, c, **kw):
    """Square markers one tixel wide at this panel's scale, so tissue reads as solid."""
    m = (C.smp == s) & C.ok
    size = (C.pitch * scale * 72 * 1.08) ** 2
    cc = c[m] if isinstance(c, np.ndarray) and len(c) == len(C.smp) else c
    if isinstance(cc, np.ndarray) and cc.dtype.kind == "f":
        o = np.argsort(cc)
        return ax.scatter(C.xy[m, 0][o], C.xy[m, 1][o], c=cc[o], s=size, marker="s",
                          linewidths=0, rasterized=True, **kw)
    return ax.scatter(C.xy[m, 0], C.xy[m, 1], c=cc, s=size, marker="s", linewidths=0,
                      rasterized=True, **kw)


def legend_cols(n, per_row=10):
    return min(n, per_row)


# =========================================================================== figures
def fig_cluster(C):
    w, h = SLOTS["cluster_spatial"]
    fig = plt.figure(figsize=(w, h))
    ftext(fig, w / 2, h - 0.03, "Spatial UMAP", ha="center", va="top", fontsize=18)
    nrow = -(-len(C.cl) // 10)
    leg_h = 0.40 + 0.33 * (nrow - 1)
    axes, sc, _ = spatial_row(C, fig, 0.1, leg_h + 0.15, w - 0.2, h - leg_h - 0.15 - 0.42,
                              title_h=0.40, title_pt=16, center=True)
    for s, ax in axes.items():
        spatial_points(C, ax, s, sc[s], C.clu_col)
    hd = [Line2D([], [], marker="o", ls="", ms=10, color=C.cm[c], label=c) for c in C.cl]
    fig.legend(handles=hd, loc="lower center", ncol=legend_cols(len(C.cl)), frameon=False,
               fontsize=14, handletextpad=0.15, columnspacing=1.1, bbox_to_anchor=(0.5, 0.02 / h))
    save(fig, C.out, C.proj, "cluster_spatial")


def umap_fig(C, name, colors, items, header, slot, legend="right", fontsize=13):
    w, h = SLOTS[slot]
    fig = plt.figure(figsize=(w, h))
    n = len(items)
    ncol = 1 if (legend != "right" or n <= 9) else 2
    leg_w = 1.45 if ncol == 1 else 2.2
    ftext(fig, w / 2 if legend == "below" else (w - leg_w) / 2 + 0.05, h - 0.03, header,
          ha="center", va="top", fontsize=17)
    if legend == "right":
        ax = ax_in(fig, 0.05, 0.05, w - leg_w - 0.1, h - 0.45)
    else:
        side = min(w - 0.2, h - 0.45 - 0.27 * n - 0.1)
        ax = ax_in(fig, (w - side) / 2, h - 0.42 - side, side, side)
    o = np.random.default_rng(0).permutation(len(C.um))
    ax.scatter(C.um[o, 0], C.um[o, 1], c=colors[o], s=0.5, linewidths=0, rasterized=True)
    ax.set_axis_off()
    hd = [Line2D([], [], marker="o", ls="", ms=9, color=c, label=l) for l, c in items]
    if legend == "right":
        fig.legend(handles=hd, loc="center left", bbox_to_anchor=((w - leg_w - 0.03) / w, 0.45),
                   ncol=ncol, frameon=False, fontsize=fontsize if n <= 9 else 11,
                   handletextpad=0.2, labelspacing=0.35, columnspacing=0.6)
    else:
        fig.legend(handles=hd, loc="lower left", bbox_to_anchor=(0.12 / w, 0.03),
                   frameon=False, fontsize=fontsize, handletextpad=0.25, labelspacing=0.3)
    save(fig, C.out, C.proj, name)


def fig_umaps(C):
    umap_fig(C, "umap_cluster", C.clu_col, [(c, C.cm[c]) for c in C.cl], "UMAP by cluster", "umap_cluster")
    sc = np.array([C.scol[s] for s in C.smp])
    umap_fig(C, "umap_sample", sc, [(C.short[s], C.scol[s]) for s in C.samples], "UMAP by sample", "umap_sample")


def fig_celltypes(C):
    w, h = SLOTS["celltype_spatial"]
    fig = plt.figure(figsize=(w, h))
    axes, sc, _ = spatial_row(C, fig, 0.05, 0.05, w - 0.1, h - 0.1, title_h=0.40, title_pt=16, center=True)
    for s, ax in axes.items():
        spatial_points(C, ax, s, sc[s], C.clu_col)
    save(fig, C.out, C.proj, "celltype_spatial")
    items = [(f"{c} {C.cell_types.get(c, '')}".strip(), C.cm[c]) for c in C.cl]
    umap_fig(C, "celltype_umap", C.clu_col, items, "Putative cell types", "celltype_umap",
             legend="below", fontsize=12.5)


def fig_neighborhood(C):
    nb = F.neighborhood_enrichment(C.cells, spatial_key="spatial_img")
    w, h = SLOTS["neighborhood"]
    fig = plt.figure(figsize=(w, h))
    ftext(fig, w / 2, h - 0.03, "Neighborhood analysis by sample", ha="center", va="top", fontsize=17)
    idx = {c: i for i, c in enumerate(C.cl)}
    n = len(C.samples)
    cols = 1 if n == 1 else 2
    rows = -(-n // cols)
    cb_cell = n % 2 == 1 and n > 1        # odd count: the colourbar takes the empty cell
    cb_w = 0 if (cb_cell or n == 1) else 1.6
    if n == 1:
        cb_w = 1.6
    cellw, cellh = (w - 0.1 - cb_w) / cols, (h - 0.45) / rows
    side = min(cellh - 0.62, cellw - 0.75)
    tick = max(7, min(11, side * 72 / len(C.cl) * 0.75))
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#dddddd")
    im = None
    for k, s in enumerate(C.samples):
        tr = nb[s]["figure"].data[0]
        z = np.array(tr.z, float)
        M = np.full((len(C.cl), len(C.cl)), np.nan)
        for i, yv in enumerate([str(v) for v in tr.y]):
            for j, xv in enumerate([str(v) for v in tr.x]):
                M[idx[yv], idx[xv]] = z[i, j]
        r, cc = divmod(k, cols)
        x0 = 0.05 + cc * cellw + (cellw - side) / 2 + 0.25
        y0 = (rows - 1 - r) * cellh + 0.42
        ax = ax_in(fig, x0, y0, side, side)
        im = ax.imshow(M, cmap=cmap, vmin=-100, vmax=100, interpolation="nearest")
        ax.set_xticks(range(len(C.cl))); ax.set_xticklabels(C.cl, rotation=90, fontsize=tick)
        ax.set_yticks(range(len(C.cl))); ax.set_yticklabels(C.cl, fontsize=tick)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_title(C.lab[s], fontsize=15, pad=5)
    if cb_cell:
        cx, cy, ch = 0.05 + cellw + cellw * 0.30, 0.55, cellh - 0.6
    else:
        cx, cy, ch = w - cb_w + 0.1, 0.9, h - 1.8
    cax = ax_in(fig, cx, cy, 0.22, ch)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("Neighborhood enrichment\n(z-score)", fontsize=13)
    cb.ax.tick_params(labelsize=12)
    if cb_cell:
        ftext(fig, cx + 1.55, cy + ch / 2, "Grey: cluster\nabsent from\nthat sample",
              ha="left", va="center", fontsize=11, color="0.35")
    save(fig, C.out, C.proj, "neighborhood")


def fig_proportion(C):
    w, h = SLOTS["proportion"]
    fig = plt.figure(figsize=(w, h))
    ftext(fig, (w - 0.9) / 2 + 0.6, h - 0.03, "Cluster proportion by sample", ha="center", va="top", fontsize=17)
    ax = ax_in(fig, 0.85, 0.95, w - 0.85 - 1.05, h - 0.95 - 0.55)
    tab = pd.crosstab(C.smp, C.clu).reindex(index=C.samples, columns=C.cl).fillna(0)
    frac = tab.div(tab.sum(1), axis=0)
    bottom = np.zeros(len(C.samples))
    for c in C.cl:
        ax.bar(range(len(C.samples)), frac[c].values, bottom=bottom, color=C.cm[c], width=0.7)
        bottom += frac[c].values
    ax.set_xticks(range(len(C.samples)))
    ax.set_xticklabels([f"{C.lab[s].split('|')[0]}\n{C.short[s]}\n(n = {int(tab.loc[s].sum()):,})"
                        for s in C.samples], fontsize=12 if len(C.samples) <= 3 else 10)
    ax.set_ylabel("Proportion of spots", fontsize=14)
    ax.set_ylim(0, 1)
    ax.margins(x=0.04)
    hd = [Line2D([], [], marker="s", ls="", ms=11, color=C.cm[c], label=c) for c in C.cl[::-1]]
    ax.legend(handles=hd, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False,
              fontsize=13 if len(C.cl) <= 12 else 11, handletextpad=0.2)
    save(fig, C.out, C.proj, "proportion")


def fig_morph(C, s):
    """Clusters on brightfield | brightfield | H&E: one panel size, the brightfield crop's aspect."""
    p = C.p
    img_path, sf = F.fetch_registered_image(p, s)
    scalef = float(sf["tissue_hires_scalef"])
    img = np.asarray(Image.open(img_path).convert("L"))
    pos = F.fetch_tissue_positions(p, s)
    sub = C.cells.obs[C.smp == s]
    bc = pd.Index(sub.index).str.split("#").str[-1].str.replace(r"-\d+$", "", regex=True)
    j = (pd.DataFrame({"barcode": bc, "cluster": sub["cluster"].astype(str).values})
         .merge(pos, on="barcode", how="inner"))
    px, py = j["x_px"].to_numpy() * scalef, j["y_px"].to_numpy() * scalef
    x0, x1, y0, y1 = px.min(), px.max(), py.min(), py.max()
    pad = 0.04 * max(x1 - x0, y1 - y0)
    cx0, cx1 = max(x0 - pad, 0.0), min(x1 + pad, float(img.shape[1]))
    cy0, cy1 = max(y0 - pad, 0.0), min(y1 + pad, float(img.shape[0]))
    a = (cx1 - cx0) / (cy1 - cy0)
    he_src = C.spec.get("he_images", {}).get(s)
    he = None
    if he_src:
        he = crop_he(he_image(C.spec_dir, he_src), a, C.spec.get("he_exclude", {}).get(s, ()))
    npan = 3 if he is not None else 2
    w, h = SLOTS["morph"]
    fig = plt.figure(figsize=(w, h))
    gap, leg_h, tit = 0.18, 0.45, 0.42
    ph = min(h - leg_h - tit - 0.05, (w - 0.1 - (npan - 1) * gap) / (npan * a))
    pw = a * ph
    xl = (w - npan * pw - (npan - 1) * gap) / 2
    yb = leg_h + (h - leg_h - tit - ph)
    pitch = float(np.median(cKDTree(np.c_[px, py]).query(np.c_[px, py], k=2)[0][:, 1]))
    pt_per_px = pw * 72 / (cx1 - cx0)
    titles = ["Clusters on brightfield", "Brightfield", "H&E (adjacent section)"][:npan]
    for k in range(npan):
        ax = ax_in(fig, xl + k * (pw + gap), yb, pw, ph)
        if k < 2:
            ax.imshow(img, cmap="gray")
            ax.set_xlim(cx0, cx1); ax.set_ylim(cy1, cy0)
        else:
            ax.imshow(he)
        if k == 0:
            ax.scatter(px, py, c=[C.cm[c] for c in j["cluster"]], s=(pitch * pt_per_px * 1.05) ** 2,
                       marker="s", linewidths=0, alpha=0.85, rasterized=True)
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_visible(False)
        ftext(fig, xl + k * (pw + gap) + pw / 2, yb + ph + 0.07, titles[k], ha="center", va="bottom", fontsize=16)
    present = [c for c in C.cl if c in set(j["cluster"])]
    hd = [Line2D([], [], marker="o", ls="", ms=10, color=C.cm[c], label=c) for c in present]
    fig.legend(handles=hd, loc="lower center", ncol=legend_cols(len(hd), 12), frameon=False,
               fontsize=14, handletextpad=0.15, columnspacing=1.1, bbox_to_anchor=(0.5, 0.02))
    C.report["notes"][f"morph_{s}"] = {"crop": [cx0, cx1, cy0, cy1], "tixels": int(len(j)),
                                       "bf_aspect": a, "he": bool(he is not None),
                                       "he_aspect": (he.shape[1] / he.shape[0]) if he is not None else None}
    save(fig, C.out, C.proj, f"morph_{s}")


def fig_row(C, feat, slot):
    """One feature: label | spatial panels | UMAP | colourbar | violin by cluster."""
    r = C.rows.get(feat, {"kind": "ge", "label": feat, "sublabel": ""})
    w, h = SLOTS[slot]
    big = h > 2.0
    fig = plt.figure(figsize=(w, h))
    v = C.cells.obs[feat].astype(float).to_numpy()
    if r["kind"] == "ge":
        vmin, vmax = np.nanpercentile(v, 1), np.nanpercentile(v, 99)
        clab = "Gene score"
    else:
        lim = np.nanpercentile(np.abs(v), 99)
        vmin, vmax = -lim, lim
        clab = "chromVAR z"
    t_pt = 13 if big else 11.5
    ftext(fig, 0.10, h / 2, r["label"], rotation=90, ha="left", va="center", fontsize=16 if big else 14,
          fontweight="bold", style="italic" if r["kind"] == "ge" else "normal")
    if r.get("sublabel"):
        ftext(fig, 0.36, h / 2, r["sublabel"], rotation=90, ha="left", va="center", fontsize=11.5 if big else 10)
    lab_w, yb = 0.62, 0.04
    inner = h - 0.08
    title_h = 0.30 if big else 0.27
    axes, sc, tot = spatial_row(C, fig, lab_w, yb, 5.3, inner, title_h=title_h, title_pt=t_pt,
                                min_w=1.12 if big else 1.0, gap=0.08,
                                labels=None if big else C.short)
    im = None
    for s, ax in axes.items():
        im = spatial_points(C, ax, s, sc[s], v, cmap="RdYlBu_r", vmin=vmin, vmax=vmax)
    x = lab_w + tot + 0.12
    side = inner - title_h
    axu = ax_in(fig, x, yb, side, side)
    o = np.argsort(v)
    axu.scatter(C.um[o, 0], C.um[o, 1], c=v[o], cmap="RdYlBu_r", vmin=vmin, vmax=vmax, s=0.35,
                linewidths=0, rasterized=True)
    axu.set_axis_off()
    ftext(fig, x + side / 2, yb + side + 0.05, "UMAP", ha="center", va="bottom", fontsize=t_pt)
    x += side + 0.08
    cax = ax_in(fig, x, yb + side * 0.08, 0.13, side * 0.84)
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=11 if big else 10)
    cb.locator = MaxNLocator(4)
    cb.update_ticks()
    x += 0.13 + 0.70 + 0.62
    vb = yb + 0.30
    axv = ax_in(fig, x, vb, w - x - 0.08, yb + side - vb)
    lo, hi = np.nanpercentile(v, 0.5), np.nanpercentile(v, 99.5)
    data = [np.clip(v[(C.clu == c) & np.isfinite(v)], lo, hi) for c in C.cl]
    parts = axv.violinplot(data, positions=range(len(C.cl)), widths=0.8, showextrema=False, showmedians=True)
    for b_, c in zip(parts["bodies"], C.cl):
        b_.set_facecolor(C.cm[c]); b_.set_edgecolor("0.3"); b_.set_linewidth(0.5); b_.set_alpha(0.9)
    parts["cmedians"].set_color("black")
    parts["cmedians"].set_linewidth(1.2)
    if r["kind"] != "ge":
        axv.axhline(0, color="0.5", lw=0.8, zorder=0)
    axv.set_xticks(range(len(C.cl)))
    axv.set_xticklabels(C.cl, fontsize=(11 if big else 10) if len(C.cl) <= 10 else 9,
                        rotation=0 if len(C.cl) <= 10 else 90)
    axv.tick_params(axis="y", labelsize=11 if big else 10)
    axv.set_ylabel(clab, fontsize=12 if big else 10.5)
    axv.set_xlim(-0.6, len(C.cl) - 0.4)
    ftext(fig, x + (w - x - 0.08) / 2, yb + side + 0.05, "By cluster", ha="center", va="bottom", fontsize=t_pt)
    save(fig, C.out, C.proj, f"row_{feat}_{slot}")


_BW = {}


def cytobands(genome, chrom):
    url = f"{F.UCSC_API}/getData/track?genome={genome};track=cytoBandIdeo;chrom={chrom}"
    with urllib.request.urlopen(url, timeout=60) as r:
        d = json.load(r)
    return d.get("cytoBandIdeo") or d.get(chrom) or []


def fig_tracks(C, gene, slot, name, genome, nbins=360):
    """Track-browser column: ideogram, kb ruler, one track per cluster on a shared y-axis,
    RefSeq model; boxes on clusters whose gene score is >= 1 SD above the cluster mean."""
    if not _BW:
        _BW.update(F.fetch_coverages(C.p, "cluster"))
    w, h = SLOTS[slot]
    fig = plt.figure(figsize=(w, h))
    model = F.refseq_gene(gene, genome=genome)
    chrom, start, end = model["locus"]
    labels, mat = F.coverage_matrix(_BW, chrom, start, end, nbins=nbins)
    labels = [str(l) for l in labels]
    vmax = float(mat.max()) or 1.0
    highlight = C.top_clusters(gene) if gene in C.cells.obs else []
    L_, R_ = 0.62, 0.06
    tw = w - L_ - R_
    y = h - 0.03
    ftext(fig, w / 2, y, gene, ha="center", va="top", fontsize=18, fontweight="bold", style="italic")
    y -= 0.42
    try:
        bands = cytobands(genome, chrom)
        if bands:
            axi = ax_in(fig, L_, y - 0.20, tw, 0.18)
            clen = max(b["chromEnd"] for b in bands)
            for b in bands:
                axi.add_patch(Rectangle((b["chromStart"], 0), b["chromEnd"] - b["chromStart"], 1,
                                        color=CYTO.get(b["gieStain"], "#cccccc"), lw=0))
            axi.add_patch(Rectangle((0, 0), clen, 1, fill=False, ec="0.3", lw=0.8))
            axi.axvline((start + end) / 2, color="#e31a1c", lw=2.5)
            axi.set_xlim(0, clen); axi.set_ylim(0, 1); axi.set_axis_off()
            ftext(fig, L_ - 0.05, y - 0.11, chrom.replace("chr", "chr "), ha="right", va="center", fontsize=10.5)
            y -= 0.30
    except Exception as e:                      # the ideogram is context, not data
        C.report["notes"][f"cytoband_{gene}"] = str(e)
    y -= 0.30
    axk = ax_in(fig, L_, y, tw, 0.001)
    axk.set_xlim(start, end); axk.set_yticks([])
    for sp in ("left", "right", "bottom"):
        axk.spines[sp].set_visible(False)
    axk.xaxis.set_ticks_position("top")
    axk.tick_params(axis="x", pad=1, length=3)
    ticks = [t for t in axk.get_xticks() if start <= t <= end]
    ticks = ticks[::max(1, -(-len(ticks) // 3))][:3]
    axk.set_xticks(ticks)
    axk.set_xticklabels([f"{t / 1e3:,.0f} kb" for t in ticks], fontsize=10.5)
    y -= 0.12
    gm_h = 0.52
    th = (y - gm_h - 0.06) / len(labels)
    x = np.linspace(start, end, mat.shape[1])
    boxes = {}
    for i, lab in enumerate(labels):
        yb = y - (i + 1) * th
        ax = ax_in(fig, L_, yb + 0.05, tw, th - 0.10)
        ax.fill_between(x, 0, mat[i], color=C.cm.get(lab, "0.4"), lw=0, step="mid")
        ax.set_xlim(start, end); ax.set_ylim(0, vmax)
        ax.set_yticks([0, vmax]); ax.set_yticklabels(["0", f"{vmax:.2f}"], fontsize=9.5)
        ax.set_xticks([]); ax.spines["bottom"].set_visible(False)
        ax.tick_params(axis="y", length=2, pad=1)
        ax.text(0.012, 0.95, lab, transform=ax.transAxes, fontsize=11 if len(labels) <= 10 else 9,
                va="top", ha="left", bbox=dict(boxstyle="square,pad=0.15", fc="white", ec="0.4", lw=0.6))
        boxes[lab] = (L_, yb + 0.02, tw, th - 0.04)
    axm = ax_in(fig, L_, 0.06, tw, gm_h)
    axm.set_xlim(start, end); axm.set_ylim(0, 1); axm.set_axis_off()
    yy = 0.70
    axm.plot([model["txStart"], model["txEnd"]], [yy, yy], color="#17418C", lw=1.0)
    for es, ee in model["exons"]:
        axm.add_patch(Rectangle((es, yy - 0.14), max(ee - es, (end - start) / 400), 0.28, color="#17418C"))
    axm.text((model["txStart"] + model["txEnd"]) / 2, yy - 0.30, gene, ha="center", va="top",
             fontsize=12, color="#17418C", style="italic")
    ftext(fig, L_ - 0.05, 0.06 + gm_h * yy, "RefSeq", ha="right", va="center", fontsize=10)
    FW, FH = fig.get_size_inches()
    for lab in highlight:
        if lab in boxes:
            bx, by, bw_, bh = boxes[lab]
            fig.patches.append(Rectangle((bx / FW, by / FH), bw_ / FW, bh / FH, transform=fig.transFigure,
                                         fill=False, ec="#17418C", lw=2.2, zorder=10))
    C.report["notes"][f"tracks_{gene}"] = {"locus": f"{chrom}:{start}-{end}", "shared_max": vmax,
                                           "highlight": highlight}
    save(fig, C.out, C.proj, name)


def fig_track_umap(C, genes):
    w, h = SLOTS["track_umap"]
    fig = plt.figure(figsize=(w, h))
    ax = ax_in(fig, 0.04, 0.04, w - 0.08, h - 0.40)
    o = np.random.default_rng(0).permutation(len(C.um))
    ax.scatter(C.um[o, 0], C.um[o, 1], c=C.clu_col[o], s=0.25, linewidths=0, rasterized=True)
    ax.set_axis_off()
    ftext(fig, w / 2, h - 0.03, "UMAP by cluster", ha="center", va="top", fontsize=13)
    by_c = {}
    for g in genes:
        by_c.setdefault(C.zscores(g).idxmax(), []).append(f"{g}+")
    x0_, x1_ = ax.get_xlim()
    for c, gl in by_c.items():
        m = C.clu == c
        cx, cy = np.median(C.um[m, 0]), np.median(C.um[m, 1])
        frac = (cx - x0_) / (x1_ - x0_)
        ha = "left" if frac < 0.3 else ("right" if frac > 0.7 else "center")
        ax.text(cx, cy, "\n".join(gl), fontsize=11, fontweight="bold", ha=ha, va="center",
                path_effects=[pe.withStroke(linewidth=3, foreground="white")], clip_on=False)
    C.report["notes"]["track_callouts"] = by_c
    save(fig, C.out, C.proj, "track_umap")


def fig_peaks(C, level):
    w, h = SLOTS["peaks"]
    fig = plt.figure(figsize=(w, h))
    df = C.p.peak_call_summary(level)
    groups = list(df.loc[~df.is_union, "Group"].unique())
    order = C.cl if level == "cluster" else C.samples
    groups = sorted(groups, key=lambda g: order.index(g.split("(")[0]) if g.split("(")[0] in order else len(order))
    groups += ["UnionPeaks"]

    def lab(g):
        if g == "UnionPeaks":
            return "Union peaks"
        n = int(re.search(r"n = (\d+)", g).group(1))
        k = g.split("(")[0]
        return f"{C.lab.get(k, k)}\n(n = {n:,})"
    ftext(fig, w / 2, h - 0.03, "Peaks per cluster" if level == "cluster" else "Peaks per sample",
          ha="center", va="top", fontsize=17)
    bot = 1.45 if level == "cluster" else 0.85
    ax = ax_in(fig, 0.85, bot, w - 0.85 - 0.15, h - bot - 0.95)
    bottom = np.zeros(len(groups))
    for anno in ["Distal", "Exonic", "Intronic", "Promoter"]:
        vals = [float(df[(df.Group == g) & (df.peak_anno == anno)]["peaks_thousands"].sum()) for g in groups]
        ax.bar(range(len(groups)), vals, bottom=bottom, color=PEAK_COLORS[anno], width=0.72, label=anno)
        bottom += np.array(vals)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([lab(g) for g in groups], fontsize=11 if level == "cluster" else 12,
                       rotation=90 if level == "cluster" or len(groups) > 4 else 0)
    ax.set_ylabel("Number of peaks (×10³)", fontsize=14)
    ax.margins(x=0.02)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, frameon=False, fontsize=13,
              handlelength=1.0, columnspacing=1.0)
    save(fig, C.out, C.proj, f"peaks_{level}")


def fig_qc(C, metric, header, name):
    w, h = SLOTS["qc_map"]
    fig = plt.figure(figsize=(w, h))
    v = C.cells.obs[metric].astype(float).to_numpy()
    vmin, vmax = np.nanpercentile(v, 1), np.nanpercentile(v, 99)
    ftext(fig, (w - 1.0) / 2, h - 0.03, header, ha="center", va="top", fontsize=16)
    axes, sc, _ = spatial_row(C, fig, 0.05, 0.05, w - 1.15, h - 0.45, title_h=0.30, title_pt=12,
                              min_w=1.0, center=True)
    im = None
    for s, ax in axes.items():
        im = spatial_points(C, ax, s, sc[s], v, cmap="viridis", vmin=vmin, vmax=vmax)
    cax = ax_in(fig, w - 0.95, 0.35, 0.14, h - 1.0)
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=11)
    save(fig, C.out, C.proj, name)


def fig_genome(C):
    """Median one-vs-rest gene-score log2FC per 5 Mb bin, per sample and per cluster.

    Reads the pipeline's delivered ranked_genes tables; no statistic is recomputed. A block
    of bins shifted together is what a copy-number change looks like in ATAC gene scores.
    """
    from matplotlib.colors import TwoSlopeNorm
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    BIN = 5_000_000
    mats = {}
    for level in ("sample", "cluster"):
        df = pd.read_csv(C.p.table(f"ranked_genes_per_{level}.csv"))
        df = df[df["seqnames"].isin(chroms)].copy()
        df["ci"] = df["seqnames"].map({c: i for i, c in enumerate(chroms)})
        df["bin"] = (df[["start", "end"]].min(axis=1) // BIN).astype(int)
        g = df.groupby(["group_name", "ci", "bin"]).agg(m=("Log2FC", "median"), n=("Log2FC", "size")).reset_index()
        mats[level] = g[g.n >= 5]
    cols = (pd.concat(mats.values())[["ci", "bin"]].drop_duplicates().sort_values(["ci", "bin"])
            .reset_index(drop=True))
    cols["x"] = np.arange(len(cols))
    key = {(a, b): x for a, b, x in cols[["ci", "bin", "x"]].itertuples(index=False)}

    def mat(df, rows):
        M = np.full((len(rows), len(cols)), np.nan)
        for i, r in enumerate(rows):
            d = df[df.group_name.astype(str) == r]
            M[i, [key[(a, b)] for a, b in zip(d.ci, d.bin)]] = d.m.values
        return M
    rows_s = [s for s in C.samples if s in set(mats["sample"].group_name.astype(str))]
    Ms, Mc = mat(mats["sample"], rows_s), mat(mats["cluster"], C.cl)
    w, h = SLOTS["genome"]
    fig = plt.figure(figsize=(w, h))
    norm = TwoSlopeNorm(vcenter=0, vmin=-1.5, vmax=1.5)
    L0 = 1.75
    W0 = w - L0 - 1.85
    hs = 0.33 * max(len(rows_s), 1)
    hc = min(3.0, 0.34 * len(C.cl))
    yb = 0.72
    bounds = cols.groupby("ci").x.agg(["min", "max"])
    im = None
    for M, names, y0, hh, title in [(Ms, [C.lab[s] for s in rows_s], yb + hc + 0.5, hs,
                                     "Per sample: each sample vs the others"),
                                    (Mc, C.cl, yb, hc, "Per cluster: each cluster vs all other spots")]:
        ax = ax_in(fig, L0, y0, W0, hh)
        im = ax.imshow(M, aspect="auto", cmap="RdBu_r", norm=norm, interpolation="nearest")
        ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=13)
        for v in bounds["min"].values[1:]:
            ax.axvline(v - 0.5, color="0.35", lw=0.5)
        for s_ in ax.spines.values():
            s_.set_visible(False)
        ax.tick_params(length=0); ax.set_xticks([])
        ax.set_title(title, loc="left", fontsize=14, pad=4)
    ticks = ((bounds["min"] + bounds["max"]) / 2).values
    names = [chroms[i].replace("chr", "") for i in bounds.index]
    ax.set_xticks(ticks)
    ax.set_xticklabels([("\n" + n) if n in ("19", "21") else n for n in names], fontsize=12)
    ax.set_xlabel("Chromosome (5 Mb bins)", fontsize=14, labelpad=4)
    cax = ax_in(fig, L0 + W0 + 0.3, 0.9, 0.2, 3.5)
    cb = fig.colorbar(im, cax=cax)
    cb.set_ticks([-1.5, -0.75, 0, 0.75, 1.5]); cb.ax.tick_params(labelsize=12)
    cb.set_label("Median gene-score log2FC\n(ArchR, one-vs-rest)", fontsize=13)
    save(fig, C.out, C.proj, "genome")


def fig_heatmaps(C):
    """ArchR's own marker heatmaps (heatmaps_all.pdf), one PNG per page, white border cropped."""
    import pypdfium2 as pdfium
    doc = pdfium.PdfDocument(F.workflow_figure(C.p, "marker_heatmaps"))
    for i in range(len(doc)):
        im = doc[i].render(scale=3).to_pil().convert("RGB")
        a = np.asarray(im)
        nz = np.where((a < 245).any(2))
        box = (max(nz[1].min() - 12, 0), max(nz[0].min() - 12, 0),
               min(nz[1].max() + 12, a.shape[1]), min(nz[0].max() + 12, a.shape[0]))
        im.crop(box).save(os.path.join(C.out, f"{C.proj}__heatmap_p{i + 1}.png"))


# =========================================================================== driver
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("project")
    ap.add_argument("--spec", help="deck spec JSON (also accepts a render_figures config)")
    ap.add_argument("--one-shot", action="store_true", help="choose marker genes from the data")
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    spec = json.load(open(a.spec)) if a.spec else {}
    spec_dir = os.path.dirname(os.path.abspath(a.spec)) if a.spec else "."
    os.makedirs(a.out, exist_ok=True)
    C = Ctx(a.project, spec, spec_dir, a.out, one_shot=a.one_shot)
    genome = str(C.p.input_parameters().get("genome", "hg38"))

    def guard(name, fn):
        try:
            fn()
            C.report["figures"][name] = "ok"
            print("OK  ", name, flush=True)
        except Exception as e:
            C.report["figures"][name] = f"FAIL {type(e).__name__}: {e}"
            print("FAIL", name, type(e).__name__, e, traceback.format_exc()[-900:], flush=True)

    def tables():
        C.p.sample_table().to_csv(os.path.join(a.out, f"samples_{a.project}.csv"), index=False)
        C.p.qc_table().to_csv(os.path.join(a.out, f"qc_table_{a.project}.csv"), index=False)
    guard("tables", tables)
    guard("cluster_spatial", lambda: fig_cluster(C))
    guard("umaps", lambda: fig_umaps(C))
    if C.cell_types:
        guard("celltypes", lambda: fig_celltypes(C))
    guard("neighborhood", lambda: fig_neighborhood(C))
    if len(C.samples) > 1:
        guard("proportion", lambda: fig_proportion(C))
    for s in C.samples:
        guard(f"morph_{s}", lambda s=s: fig_morph(C, s))
    guard("heatmaps", lambda: fig_heatmaps(C))
    for f in C.rows:
        if f not in C.cells.obs:
            C.report["figures"][f"row_{f}"] = "SKIP feature absent from the h5ad"
            continue
        for slot in ("row2", "row3"):
            guard(f"row_{f}_{slot}", lambda f=f, slot=slot: fig_row(C, f, slot))
    if len(C.track_genes) >= 1:
        for g in C.track_genes[:3]:
            guard(f"track_{g}", lambda g=g: fig_tracks(C, g, "track_col", f"track_{g}", genome))
        guard("track_umap", lambda: fig_track_umap(C, C.track_genes[:3]))
    for name, genes in C.spec.get("track_sets", {}).items():
        for g in genes[:3]:
            guard(f"trackset_{name}_{g}", lambda g=g, name=name:
                  fig_tracks(C, g, "track_col_wide", f"trackset_{name}_{g}", genome))
    guard("peaks_cluster", lambda: fig_peaks(C, "cluster"))
    if len(C.samples) > 1:
        guard("peaks_sample", lambda: fig_peaks(C, "sample"))
    guard("qc_tss", lambda: fig_qc(C, "TSSEnrichment", "TSS enrichment", "qc_TSS"))
    guard("qc_frags", lambda: fig_qc(C, "ln(nFrags)", "ln(nFrags)", "qc_frags"))
    if C.spec.get("genome_heatmap"):
        guard("genome", lambda: fig_genome(C))
    wf = os.path.join(HERE, "..", "assets", "bioinformatics_workflow.png")
    if os.path.exists(wf):
        shutil.copyfile(wf, os.path.join(a.out, f"{a.project}__bioinformatics_workflow.png"))
    C.report["rows"] = C.rows
    C.report["track_genes"] = C.track_genes[:3]
    C.report["samples"] = C.samples
    C.report["sample_short"] = C.short
    C.report["clusters"] = C.cl
    C.report["zscores"] = {g: C.zscores(g).round(2).to_dict() for g in C.track_genes[:3]}
    tab = pd.crosstab(C.clu, C.smp).reindex(index=C.cl, columns=C.samples).fillna(0)
    C.report["cluster_sizes"] = {c: int(tab.loc[c].sum()) for c in C.cl}
    main_s = {}
    for c in C.cl:
        frac = tab.loc[c] / tab.loc[c].sum()
        top = frac.idxmax()
        main_s[c] = (f"{C.short[top]} ({frac[top]:.0%})" if frac[top] >= 0.75
                     else f"mixed ({int((frac > 0.05).sum())} samples)")
    C.report["cluster_main_sample"] = main_s
    json.dump(C.report, open(os.path.join(a.out, "slide_figures_report.json"), "w"), indent=1, default=str)
    bad = [k for k, v in C.report["figures"].items() if v.startswith("FAIL")]
    print(f"\n{len(C.report['figures']) - len(bad)} ok/skipped, {len(bad)} failed: {bad}")


if __name__ == "__main__":
    main()
