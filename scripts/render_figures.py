"""Render every figure the Data Report deck needs, from a Latch project.

This is the missing half of the pipeline: `build_demo_deck.py` assembles a deck from
figures that already exist, and this produces them. Run it on the analysis box (it needs
`/data/analysis_cache` and the Latch CLI), then run the deck builder.

    python render_figures.py Demo_001
    python render_figures.py <project> --one-shot             # pick features from the data
    python render_figures.py <project> --config features.json  # features chosen elsewhere

Every figure is produced under `guard()`, so one failure is reported in
`render_report.json` rather than taking the run down with it -- a 25-figure render that
dies on figure 3 wastes the twenty minutes already spent.

THE CONFIG BLOCK IS NOT BOILERPLATE. The marker genes, the motifs and the Compare run are
per-project scientific choices. Confirm them with whoever ordered the report before
rendering; see "Choosing features for gene/motif panels" (Part 3) in SKILL.md. The
exception is a one-shot report, where the scientist has asked for best assumptions
instead of confirmations: `--config` takes a JSON file with the same keys (what you
chose for the tissue), and `--one-shot` without it picks marker genes from the data
(`atx_figures.auto_feature_config`). Either way the choices are written to
`render_report.json` under "assumptions" for the deck to print.
"""
import argparse
import json
import os
import sys
import traceback

import matplotlib
matplotlib.use("Agg")

import atx_figures as F
from atx_project import Project

_ap = argparse.ArgumentParser()
_ap.add_argument("project", nargs="?", default="Demo_001")
_ap.add_argument("--one-shot", action="store_true",
                 help="no CONFIG entry needed: choose marker genes from the data")
_ap.add_argument("--config", help="JSON file of feature choices (CONFIG keys)")
ARGS = _ap.parse_args()
PROJ = ARGS.project

CONFIG = {
    "Demo_001": {
        "genome": "mm39",
        "logo_genome": "mm10",          # no CisBP seqlogo file ships for mm39
        "genes": ["Cd3e", "Cd19", "Hmox1"],
        "gene_rows": {"Cd3e": "Cd3e: T cell", "Cd19": "Cd19: B cell",
                      "Hmox1": "Hmox1:<br>red pulp macrophage"},
        "motifs": ["Ascl1-780", "Runx1-712"],
        "compare_run": "Demo_CD19_high_CondA_vs_CondB",
        "compare_column": "CD19_high",
        # Which of the workflow's GroupA/GroupB is which: confirm against the cell
        # counts in the compare_config filename (`..._A<nA>_B<nB>_...`), do not assume
        # from the run name. apply_compare_annotation() leaves them as GroupA/GroupB
        # without this mapping.
        "compare_mapping": {"GroupA": "CondA", "GroupB": "CondB"},
        "compare_groups": ["CondA", "CondB"],
        "colors": "DEMO_CLUSTER_COLORS",
    },
}

#: Used for any project with no CONFIG entry: the standard figures only. Empty selections
#: SKIP their panels rather than failing, so a first pass can run before the scientist has
#: chosen marker genes, motifs or a comparison. The guard that matters is still in place —
#: nothing here guesses a feature; those panels simply are not drawn.
#:
#: DO NOT add real project, customer or run names to this file. It is published.
STANDARD_ONLY = {
    "genome": None,                 # taken from the project's input_parameters.csv
    "logo_genome": None,            # from the genome via atx_figures.LOGO_GENOME
    "genes": [],
    "gene_rows": {},
    "motifs": [],
    "compare_run": None,
    "compare_column": None,
    "compare_mapping": None,
    "compare_groups": [],
    "colors": None,                 # pinned below: H5 default in the viewer's order
}
p = Project(PROJ, cache="/data/analysis_cache")
assumptions = []
if ARGS.config:
    with open(ARGS.config) as fh:
        cfg = {**STANDARD_ONLY, **json.load(fh)}
    assumptions = list(cfg.pop("assumptions", []))
else:
    cfg = CONFIG.get(PROJ)
cells = p.cell_adata()
if cfg is None and ARGS.one_shot:
    cfg, assumptions = F.auto_feature_config(p, cells)
    print(f"One-shot: marker genes chosen from the data: {cfg['genes']}")
elif cfg is None:
    cfg = dict(STANDARD_ONLY)
    print(f"No feature configuration for {PROJ!r}: rendering the standard figures only. "
          "Gene-activity, coverage-track and differential-motif panels are skipped until "
          "marker genes, motifs and a comparison are agreed with the scientist.")
if not cfg.get("genome"):
    cfg["genome"] = str(p.input_parameters()["genome"])
if cfg.get("logo_genome") is None:
    cfg["logo_genome"] = F.LOGO_GENOME.get(cfg["genome"], "mm10")
# One palette for every figure. A named palette (e.g. "DEMO_CLUSTER_COLORS") or an
# explicit {cluster: hex} dict from a JSON config wins; otherwise pin the H5 default in
# the viewer's order. Never leave this None: the generators' own defaults disagree with
# each other, so an unpinned run colours the same cluster differently on different slides.
_colors = cfg.get("colors")
if isinstance(_colors, str):
    CM = getattr(F, _colors)
elif isinstance(_colors, dict):
    CM = dict(_colors)
else:
    CM = F.default_cluster_colors(cells.obs["cluster"].astype(str).unique())
_missing = sorted(set(cells.obs["cluster"].astype(str)) - set(CM))
if _missing:
    raise KeyError(f"palette has no colour for clusters {_missing}")
samples = F.natural_sort(cells.obs["sample"].unique().tolist())
labels = F.deck_panel_labels(cells, samples)

report = {"project": PROJ, "n_obs": int(cells.n_obs), "samples": samples,
          "panel_labels": labels, "figures": {}, "notes": [],
          "config": {k: v for k, v in cfg.items() if k != "assumptions"},
          "palette": {c: CM[c] for c in F.natural_sort(CM)},
          "assumptions": assumptions}


def guard(name, fn):
    try:
        out = fn()
        png = f"{PROJ}__{name}.png" if out is None else out
        size = os.path.getsize(png) if os.path.exists(png) else None
        report["figures"][name] = {"png": png, "bytes": size}
        print(f"OK   {name}  {size}")
    except Exception as e:
        report["figures"][name] = {"status": "FAIL",
                                   "error": f"{type(e).__name__}: {e}",
                                   "tb": traceback.format_exc()[-600:]}
        print(f"FAIL {name}: {type(e).__name__}: {str(e)[:200]}")


def skip(name, reason):
    """Record a panel that is not applicable, distinct from one that failed."""
    report["figures"][name] = {"status": "SKIPPED", "reason": reason}
    print(f"SKIP {name}: {reason}")


def emit(fig, name, w, h, scale=2):
    fig.update_layout(width=w, height=h)
    fig.write_image(f"{PROJ}__{name}.png", scale=scale)


def emit_spatial(fig, name, w, h, axes=None):
    """emit() for a per-sample spatial figure, recording the rendered scale of each panel.

    With SHARED_SCALE on, every spatial panel must come out at the same pixels per unit;
    a spread over 1% means the shared range did not survive layout, and the figure is
    failed rather than shipped with one section silently enlarged.
    """
    emit(fig, name, w, h)
    if not (SHARED_SCALE and len(samples) > 1):
        return
    sc = F.panel_scales(fig, axes=axes)
    report.setdefault("panel_scales", {})[name] = sc
    xs = [v[0] for v in sc.values()]
    if len(xs) > 1 and (max(xs) - min(xs)) > 0.01 * max(xs):
        raise RuntimeError(f"{name}: spatial panels rendered at different scales {sc}")


# --- the sample table: the deck's first slide, and what the scientist checks first -----
def _samples():
    t = p.sample_table()
    t.to_csv(f"samples_{PROJ}.csv", index=False)
    report["figures"]["sample_table"] = {"csv": f"samples_{PROJ}.csv", "rows": int(len(t))}
    print(f"OK   sample_table  {len(t)} runs")
guard("samples", _samples)   # records itself under "sample_table"

# --- orientation first: every spatial panel below is drawn in the image frame ----------
report["alignment"] = F.align_spatial_to_image(p, cells, samples)

# --- one physical scale: a small section must look small ------------------------------
# Every spatial panel shares one span so relative tissue size is true. Only valid when the
# runs are on the same chip type (the coordinates are crop-space pixels, ~7,800 px across
# any chip); otherwise each panel keeps its own scale and the report says so.
SHARED_SCALE, report["chip_geometry"] = F.check_shared_scale(p, samples)
report["shared_scale"] = SHARED_SCALE
if len(samples) > 1 and not SHARED_SCALE:
    report["notes"].append("runs are on different chip types: spatial panels are drawn "
                           "at per-panel scale, so relative tissue size is NOT comparable")
    print("NOTE chip geometry differs between runs: per-panel spatial scale")

# --- slide 2: spatial clusters + two UMAPs --------------------------------------------
guard("cluster_spatial_h5", lambda: emit_spatial(
    F.cluster_maps(cells, samples=samples, color_by="cluster", color_map=CM,
                   coords="spatial_img", flipY=False, panel_labels=labels,
                   legend_header=True, shared_scale=SHARED_SCALE)[0],
    "cluster_spatial_h5", 1150, 640))

guard("umap_cluster", lambda: emit(
    F.umap_overview(cells, color_by="cluster", color_map=CM,
                    category_order=F.natural_sort(
                        cells.obs["cluster"].astype(str).unique().tolist()))[0],
    "umap_cluster", 640, 460))

if len(samples) > 1:
    guard("umap_sample", lambda: emit(
        F.umap_overview(cells, color_by="SampleName")[0], "umap_sample", 640, 460))
else:
    skip("umap_sample", "single-sample project: one colour")

# --- slide 3: neighborhood grid + cluster proportion ----------------------------------
def _nb():
    nb = F.neighborhood_enrichment(cells, spatial_key="spatial_img")
    report["notes"].append(f"neighborhood panels: {list(nb)}")
    F.neighborhood_grid(nb, labels=labels).write_image(
        f"{PROJ}__neighborhood_grid.png", scale=2)
guard("neighborhood_grid", _nb)

# A proportion bar compares samples; with one sample it is a single bar that says nothing.
if len(samples) > 1:
    guard("cluster_proportion", lambda: emit(
        F.cluster_proportion(cells, color_map=CM, group_labels=labels, font_pt=16)[0],
        "cluster_proportion", 820, 640))
else:
    report["figures"]["cluster_proportion"] = {
        "status": "SKIPPED", "reason": "single-sample project: nothing to compare"}
    print("SKIP cluster_proportion: single-sample project")

# --- slides 4-6: one morphology overlay per run (cropped to the tixels) ---------------
def _morph():
    for s in samples:
        fig, info = F.morphology_overlay(p, cells, s, color_map=CM)
        fig.savefig(f"{PROJ}__morphology_{s}.png", dpi=150, bbox_inches="tight")
        report["figures"][f"morphology_{s}"] = {
            "png": f"{PROJ}__morphology_{s}.png", "crop_box": info["crop_box"],
            "tixels": info["tixels_plotted"]}
        print(f"OK   morphology_{s}  crop {info['crop_box']}")
guard("morphology", _morph)
F.assert_sample_coverage(samples, morphology=[f"morphology_{s}" for s in samples])

# --- slide 7: marker heatmaps ArchR already wrote -------------------------------------
def _heatmaps():
    """Fetch heatmaps_all.pdf and split it into the three page PNGs the deck places.

    The deck shows marker genes, marker peaks and marker motifs side by side; they are
    three pages of one ArchR PDF, so the deck cannot use the PDF directly.
    """
    import pypdfium2 as pdfium
    pdf_path = F.workflow_figure(p, "marker_heatmaps")
    doc = pdfium.PdfDocument(pdf_path)
    report["notes"].append(f"heatmaps_all.pdf pages: {len(doc)}")
    for i in range(len(doc)):
        out = f"marker_heatmap_p{i + 1}.png"
        doc[i].render(scale=3).to_pil().save(out)
        report["figures"][f"marker_heatmap_p{i + 1}"] = {
            "png": out, "bytes": os.path.getsize(out)}
        print(f"OK   marker_heatmap_p{i + 1}")
guard("marker_heatmaps", _heatmaps)

# --- slide 8: genome-browser coverage tracks ------------------------------------------
def _tracks():
    bw = F.fetch_coverages(p, "cluster")
    for gene in cfg["genes"]:
        model = F.refseq_gene(gene, genome=cfg["genome"])
        fig, info = F.coverage_tracks(bw, model["locus"], gene_model=model,
                                      color_map=CM, track_height=0.33, width=3.5,
                                      title=gene)
        # No bbox_inches="tight" here: the track stack is sized deliberately
        # (width=3.5in at 300 dpi = 1050 px) and trimming it changes the width each
        # gene ends up with, so the three panels stop matching each other on the slide.
        fig.savefig(f"{PROJ}__tracks_{gene}.png", dpi=300)
        report["figures"][f"tracks_{gene}"] = {
            "png": f"{PROJ}__tracks_{gene}.png", "locus": model["locus"],
            "shared_max": info.get("shared_max")}
        print(f"OK   tracks_{gene}  {model['locus']}  max {info.get('shared_max')}")
if cfg["genes"]:
    guard("tracks", _tracks)
else:
    skip("tracks", "no marker genes chosen yet (coverage bigwigs not fetched)")

# --- slide 9: gene-activity blocks + annotated UMAP -----------------------------------
def _genes():
    cells2, missing = F.add_features(p, cells, cfg["genes"], kind="ge")
    if missing:
        report["notes"].append(f"gene features absent from the h5ad: {missing}")
    for g in cfg["genes"]:
        if g in missing:
            continue
        fig, info = F.gene_activity_block(cells2, g, samples=samples,
                                          sample_labels=labels,
                                          row_label=cfg["gene_rows"].get(g, g),
                                          panel_w=250, height=420,
                                          shared_scale=SHARED_SCALE)
        fig.write_image(f"{PROJ}__block_{g}.png", scale=2)
        report["figures"][f"block_{g}"] = {"png": f"{PROJ}__block_{g}.png",
                                           "mean": info.get("mean")}
        print(f"OK   block_{g}  mean {info.get('mean')}")
    tops = F.top_clusters_for_features(cells2, cfg["genes"])
    report["callouts"] = {k: [sp["cluster"] for sp in v] for k, v in tops.items()}
    fig, _ = F.umap_overview(cells2, color_by="cluster", color_map=CM,
                             category_order=F.natural_sort(
                                 cells2.obs["cluster"].astype(str).unique().tolist()))
    F.annotate_umap(fig, tops)
    emit(fig, "umap_cluster_annotated", 720, 600)
if cfg["genes"]:
    guard("gene_blocks", _genes)
else:
    skip("gene_blocks", "no marker genes chosen yet")

# --- slide 10: peak counts ------------------------------------------------------------
levels = ("cluster", "sample") if len(samples) > 1 else ("cluster",)
for lvl in levels:
    guard(f"peaks_per_{lvl}",
          lambda lvl=lvl: emit(F.peaks_per_group(p, lvl)[0],
                               f"peaks_per_{lvl}", 1000, 650))
if "sample" not in levels:
    skip("peaks_per_sample", "single-sample project")

# --- slide 11: differential motifs, from the Compare workflow's own statistics ---------
def _motifs():
    ann = F.fetch_compare_annotation(p, cfg["compare_run"], column=cfg["compare_column"])
    cells3 = cells.copy()
    counts = F.apply_compare_annotation(cells3, ann, column=cfg["compare_column"],
                                        mapping=cfg["compare_mapping"])
    report["notes"].append(f"compare groups: {counts}")

    emit_spatial(F.cluster_maps(cells3, samples=samples, color_by=cfg["compare_column"],
                                coords="spatial_img", flipY=False, panel_labels=labels,
                                color_map=None, palette=None,
                                shared_scale=SHARED_SCALE)[0],
                 "cd19_high_spatial", 720, 520)

    sub = F.subset_to_groups(cells3, cfg["compare_column"], cfg["compare_groups"])
    sub, missing = F.add_features(p, sub, cfg["motifs"], kind="motifs")
    if missing:
        report["notes"].append(f"motifs absent from the h5ad: {missing}")
    for m in cfg["motifs"]:
        if m in missing:
            continue
        emit_spatial(F.cluster_maps(sub, samples=samples, color_by=m, coords="spatial_img",
                                    flipY=False, palette=None, panel_labels=labels,
                                    shared_scale=SHARED_SCALE)[0],
                     f"motif_sel_{m}", 620, 470)
        emit(F.feature_violin(sub, m, cfg["compare_column"],
                              order=cfg["compare_groups"][::-1])[0],
             f"motif_violin_{m}", 520, 430)
        # motif_logo returns a PLOTLY figure (not matplotlib), so it saves the same way
        # as every other plotly panel here.
        lg, linfo = F.motif_logo(m, genome=cfg["logo_genome"])
        emit(lg, f"logo_{m}", 700, 320)
        print(f"OK   logo_{m}  len {linfo['length']}")

    df = F.fetch_compare_results(p, cfg["compare_run"], feature="motif")
    df.to_csv(f"{PROJ}__all_motifs.csv", index=False)
    report["notes"].append(f"compare motif rows: {len(df)}")
    emit(F.motif_rank_plot(df)[0], "motif_rank", 760, 620)
if cfg["compare_run"] and cfg["motifs"]:
    guard("motifs", _motifs)
else:
    skip("motifs", "no compare run / motifs chosen yet")

# --- appendix: the QC table and its spatial maps --------------------------------------
def _qc():
    qc = p.qc_table()
    qc.to_csv(f"qc_table_{PROJ}.csv", index=False)
    report["figures"]["qc_table"] = {"csv": f"qc_table_{PROJ}.csv", "rows": int(len(qc))}
    for m in ("TSSEnrichment", "log10_nFrags"):
        emit_spatial(F.cluster_maps(cells, samples=samples, color_by=m, coords="spatial_img",
                                    flipY=False, palette=None, panel_labels=labels,
                                    shared_scale=SHARED_SCALE)[0],
                     f"qc_spatial_{m}", 1150, 620)
        print(f"OK   qc_spatial_{m}")
guard("qc", _qc)

# --- the closing workflow diagram is a reusable asset, not a per-project render -------
def _workflow():
    import shutil
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "assets", "bioinformatics_workflow.png")
    dst = f"{PROJ}__bioinformatics_workflow.png"
    if not os.path.exists(src):
        raise FileNotFoundError(f"missing skill asset {src}")
    shutil.copyfile(src, dst)
    return dst
guard("bioinformatics_workflow", _workflow)

with open("render_report.json", "w") as f:
    json.dump(report, f, indent=1, default=str)

failed = [k for k, v in report["figures"].items() if v.get("status") == "FAIL"]
skipped = [k for k, v in report["figures"].items() if v.get("status") == "SKIPPED"]
print(f"\n{len(report['figures']) - len(failed) - len(skipped)} ok, {len(failed)} failed"
      + (f": {failed}" if failed else "")
      + (f", {len(skipped)} skipped: {skipped}" if skipped else ""))
