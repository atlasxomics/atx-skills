"""Render every figure the Data Report deck needs, from a Latch project.

This is the missing half of the pipeline: `build_demo_deck.py` assembles a deck from
figures that already exist, and this produces them. Run it on the analysis box (it needs
`/data/analysis_cache` and the Latch CLI), then run the deck builder.

    python render_figures.py Demo_001

Every figure is produced under `guard()`, so one failure is reported in
`render_report.json` rather than taking the run down with it -- a 25-figure render that
dies on figure 3 wastes the twenty minutes already spent.

THE CONFIG BLOCK IS NOT BOILERPLATE. The marker genes, the motifs and the Compare run are
per-project scientific choices. Confirm them with whoever ordered the report before
rendering; see "Choosing features for gene/motif panels" (Part 3) in SKILL.md.
"""
import json
import os
import sys
import traceback

import matplotlib
matplotlib.use("Agg")

import atx_figures as F
from atx_project import Project

PROJ = sys.argv[1] if len(sys.argv) > 1 else "Demo_001"

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
cfg = CONFIG.get(PROJ)
if cfg is None:
    raise SystemExit(f"No feature configuration for {PROJ!r}. Add one to CONFIG after "
                     f"agreeing the marker genes and motifs with the requester.")

p = Project(PROJ, cache="/data/analysis_cache")
cells = p.cell_adata()
CM = getattr(F, cfg["colors"])
samples = F.natural_sort(cells.obs["sample"].unique().tolist())
labels = F.deck_panel_labels(cells, samples)

report = {"project": PROJ, "n_obs": int(cells.n_obs), "samples": samples,
          "panel_labels": labels, "figures": {}, "notes": []}


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


def emit(fig, name, w, h, scale=2):
    fig.update_layout(width=w, height=h)
    fig.write_image(f"{PROJ}__{name}.png", scale=scale)


# --- orientation first: every spatial panel below is drawn in the image frame ----------
report["alignment"] = F.align_spatial_to_image(p, cells, samples)

# --- slide 2: spatial clusters + two UMAPs --------------------------------------------
guard("cluster_spatial_h5", lambda: emit(
    F.cluster_maps(cells, samples=samples, color_by="cluster", color_map=CM,
                   coords="spatial_img", flipY=False, panel_labels=labels,
                   legend_header=True)[0], "cluster_spatial_h5", 1150, 640))

guard("umap_cluster", lambda: emit(
    F.umap_overview(cells, color_by="cluster", color_map=CM,
                    category_order=F.natural_sort(
                        cells.obs["cluster"].astype(str).unique().tolist()))[0],
    "umap_cluster", 640, 460))

guard("umap_sample", lambda: emit(
    F.umap_overview(cells, color_by="SampleName")[0], "umap_sample", 640, 460))

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
guard("tracks", _tracks)

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
                                          panel_w=250, height=420)
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
guard("gene_blocks", _genes)

# --- slide 10: peak counts ------------------------------------------------------------
for lvl in ("cluster", "sample"):
    guard(f"peaks_per_{lvl}",
          lambda lvl=lvl: emit(F.peaks_per_group(p, lvl)[0],
                               f"peaks_per_{lvl}", 1000, 650))

# --- slide 11: differential motifs, from the Compare workflow's own statistics ---------
def _motifs():
    ann = F.fetch_compare_annotation(p, cfg["compare_run"], column=cfg["compare_column"])
    cells3 = cells.copy()
    counts = F.apply_compare_annotation(cells3, ann, column=cfg["compare_column"],
                                        mapping=cfg["compare_mapping"])
    report["notes"].append(f"compare groups: {counts}")

    emit(F.cluster_maps(cells3, samples=samples, color_by=cfg["compare_column"],
                        coords="spatial_img", flipY=False, panel_labels=labels,
                        color_map=None, palette=None)[0],
         "cd19_high_spatial", 720, 520)

    sub = F.subset_to_groups(cells3, cfg["compare_column"], cfg["compare_groups"])
    sub, missing = F.add_features(p, sub, cfg["motifs"], kind="motifs")
    if missing:
        report["notes"].append(f"motifs absent from the h5ad: {missing}")
    for m in cfg["motifs"]:
        if m in missing:
            continue
        emit(F.cluster_maps(sub, samples=samples, color_by=m, coords="spatial_img",
                            flipY=False, palette=None,
                            panel_labels=labels)[0],
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
guard("motifs", _motifs)

# --- appendix: the QC table and its spatial maps --------------------------------------
def _qc():
    qc = p.qc_table()
    qc.to_csv(f"qc_table_{PROJ}.csv", index=False)
    report["figures"]["qc_table"] = {"csv": f"qc_table_{PROJ}.csv", "rows": int(len(qc))}
    for m in ("TSSEnrichment", "log10_nFrags"):
        emit(F.cluster_maps(cells, samples=samples, color_by=m, coords="spatial_img",
                            flipY=False, palette=None, panel_labels=labels)[0],
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
