"""Assemble a Data Report deck for any project from a render_figures.py output folder.

    python build_report_deck.py <render_dir> [--out deck.pptx] [--draft]

`render_dir` holds what `render_figures.py` wrote: the `<project>__*.png` panels,
`marker_heatmap_p*.png`, `qc_table_<project>.csv`, `samples_<project>.csv` and
`render_report.json`. Slides are built only from figures that exist and did not fail, so a
skipped or failed panel drops its slide rather than leaving a hole; what was left out is
listed in the speaker notes of the last slide.

The slide structures follow the delivered decks (references/deck.md). Titles are
descriptive, not findings: a one-shot draft has not been read by a scientist, so it does
not claim biology. `--draft` marks the title slide and appends the assumptions slide.

`build_demo_deck.py` is the hand-tuned Demo_001 deck; this is the general case.
"""
import argparse
import json
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atx_deck  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATE = os.path.join(HERE, "..", "assets", "atx_template.pptx")

#: Columns of the sample table that fit on one slide, in display order.
SAMPLE_COLUMNS = ("run_id", "sample_name", "condition", "n_cells", "Species", "Organ",
                  "Tissue type", "Disease", "Workflow")


def report_slides(d: str, report: dict):
    """Slide specs for build_deck, plus the names of slides left out."""
    proj = report["project"]
    figs = report.get("figures", {})
    cfg = report.get("config", {})
    labels = report.get("panel_labels", {})
    samples = report.get("samples", [])
    left_out = []

    def ok(name):
        return figs.get(name, {}).get("status") not in ("FAIL", "SKIPPED")

    def f(stem, prefixed=True):
        path = os.path.join(d, f"{proj}__{stem}.png" if prefixed else stem)
        return path if os.path.exists(path) else None

    slides = []

    st = os.path.join(d, f"samples_{proj}.csv")
    if os.path.exists(st):
        df = pd.read_csv(st)
        cols = [c for c in SAMPLE_COLUMNS if c in df.columns]
        slides.append({"title": "Samples in this report", "table": df[cols]})
    else:
        left_out.append("sample table")

    sp, uc, us = f("cluster_spatial_h5"), f("umap_cluster"), f("umap_sample")
    if sp and uc:
        slides.append({"title": "Spatial clustering",
                       "columns": [[sp], [x for x in (uc, us) if x]],
                       "column_titles": ["Spatial UMAP", "UMAP"]})
    else:
        left_out.append("clustering")

    nb, prop = f("neighborhood_grid"), f("cluster_proportion")
    if nb and prop:
        slides.append({"title": "Neighborhood enrichment and cluster proportions",
                       "subtitle": ("Neighborhood analysis identifies how often each cluster "
                                    "is adjacent to another cluster (higher number is more "
                                    "related)"),
                       "columns": [[nb], [prop]],
                       "column_titles": ["Neighborhood analysis by sample",
                                         "Cluster proportion by sample"]})
    elif nb:
        slides.append({"title": "Neighborhood enrichment",
                       "subtitle": ("How often each cluster is adjacent to another cluster "
                                    "(higher number is more related)"),
                       "images": [nb]})
    else:
        left_out.append("neighborhood")

    for s in samples:
        m = f(f"morphology_{s}")
        if m:
            slides.append({"title": f"Clustering on morphology: {labels.get(s, s)}",
                           "images": [m]})
        else:
            left_out.append(f"morphology {s}")

    hm = [os.path.join(d, f"marker_heatmap_p{i}.png") for i in (1, 2, 3)]
    if all(os.path.exists(h) for h in hm):
        slides.append({"title": "Differential genes, peaks and motifs between clusters",
                       "subtitle": "P-value <=0.05 & Log2FC >0.1 for all plots.",
                       "columns": [[h] for h in hm],
                       "column_titles": ["Gene activity scores", "Peak enrichment",
                                         "Motif enrichment"]})
    else:
        left_out.append("marker heatmaps")

    genes = cfg.get("genes", [])
    tracks = [(g, f(f"tracks_{g}")) for g in genes]
    tracks = [(g, t) for g, t in tracks if t]
    if tracks and uc:
        slides.append({"title": "Coverage at marker gene loci",
                       "columns": [[t] for _, t in tracks] + [[uc]],
                       "column_titles": [g for g, _ in tracks] + [""],
                       "column_weights": [1] * len(tracks) + [0.85]})
    elif genes:
        left_out.append("coverage tracks")

    blocks = [f(f"block_{g}") for g in genes]
    blocks = [b for b in blocks if b]
    ann = f("umap_cluster_annotated")
    cells = blocks + ([ann] if ann else [])
    if blocks:
        rows = [cells[i:i + 2] for i in range(0, len(cells), 2)]  # equal cells, 2 across
        slides.append({"title": "Gene activity of marker genes", "grid": rows,
                       "notes": "Marker genes: " + ", ".join(
                           cfg.get("gene_rows", {}).get(g, g) for g in genes)})
    elif genes:
        left_out.append("gene activity blocks")

    pc, ps = f("peaks_per_cluster"), f("peaks_per_sample")
    if pc:
        slides.append({"title": "Total number of peaks per cluster"
                                + (" and sample" if ps else ""),
                       "subtitle": "N= number of datapoints per cluster",
                       "images": [x for x in (pc, ps) if x]})
    else:
        left_out.append("peaks")

    motifs = cfg.get("motifs", [])
    mrows = cfg.get("motif_rows", {})
    have = [(m, f(f"motif_block_{m}"), f(f"logo_{m}")) for m in motifs]
    have = [(m, b, lg) for m, b, lg in have if b]
    for i in range(0, len(have), 4):  # four motif rows fill a slide
        part = have[i:i + 4]
        spec = {"title": "Motif activity of selected motifs",
                "subtitle": "chromVAR deviation z-score per spot",
                "notes": "Motifs: " + ", ".join(mrows.get(m, m) for m, _, _ in part)}
        if all(lg for _, _, lg in part):
            spec.update({"columns": [[b for _, b, _ in part], [lg for _, _, lg in part]],
                         "column_titles": ["", "Motif (CisBP)"],
                         "column_weights": [3, 1]})
        else:
            blocks_m = [b for _, b, _ in part]
            spec["grid"] = [blocks_m[j:j + 2] for j in range(0, len(blocks_m), 2)]
        slides.append(spec)
    if motifs and not have:
        left_out.append("motif activity")

    if motifs and ok("motifs") and cfg.get("compare_run"):
        grid = []
        for m in motifs:
            sel, logo, vio = f(f"motif_sel_{m}"), f(f"logo_{m}"), f(f"motif_violin_{m}")
            if sel and vio:
                grid.append([[x for x in (sel, logo) if x], vio])
        extra = [x for x in (f("cd19_high_spatial"), f("motif_rank")) if x]
        for row, e in zip(grid, extra):
            row.append(e)
        if grid:
            slides.append({"title": "Differential motifs between the compared groups",
                           "grid": grid,
                           "notes": (f"MeanDiff / FDR from ArchR getMarkerFeatures via the "
                                     f"Compare workflow run {cfg.get('compare_run')}; "
                                     f"groups as assigned by that run.")})
        else:
            left_out.append("differential motifs")

    slides.append({"title": "Appendix", "section": True})

    qc = os.path.join(d, f"qc_table_{proj}.csv")
    tss, frags = f("qc_spatial_TSSEnrichment"), f("qc_spatial_log10_nFrags")
    if os.path.exists(qc):
        spec = {"title": "Samples pass QCs for each assay", "table": pd.read_csv(qc)}
        if tss and frags:
            spec.update({"columns": [[tss], [frags]],
                         "column_titles": ["TSS Enrichment", "log10nFrags"]})
        slides.append(spec)
    else:
        left_out.append("QC table")

    wf = f("bioinformatics_workflow")
    if wf:
        slides.append({"title": "Bioinformatics Workflow",
                       "subtitle": "For more information go to our website",
                       "images": [wf]})
    return slides, left_out


def build(render_dir: str, out: str = None, draft: bool = False, template: str = TEMPLATE):
    with open(os.path.join(render_dir, "render_report.json")) as fh:
        report = json.load(fh)
    proj = report["project"]
    slides, left_out = report_slides(render_dir, report)
    if draft:
        rows = list(report.get("assumptions", []))
        if left_out:
            rows.append("Left out (failed, skipped or not applicable): "
                        + ", ".join(left_out) + ".")
        rows.append("Not reviewed by a scientist. Check every slide before it goes to "
                    "the customer, and remove this slide.")
        slides.append({"title": "Assumptions in this draft",
                       "table": pd.DataFrame({"Assumption": rows})})
    elif left_out:
        slides[-1]["notes"] = (slides[-1].get("notes", "") + "\nLeft out: "
                               + ", ".join(left_out)).strip()
    out = out or os.path.join(render_dir, f"{proj}_Data_Report{'_draft' if draft else ''}.pptx")
    # The delivered decks open on the brand line over "Data Report" (references/deck.md).
    atx_deck.build_deck(slides, out, template=template,
                        deck_title="Chart the unexplored depths of spatial biology",
                        subtitle=f"Data Report: {proj}" + (" (draft)" if draft else ""))
    return out, left_out


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("render_dir")
    ap.add_argument("--out")
    ap.add_argument("--draft", action="store_true")
    a = ap.parse_args()
    path, missing = build(a.render_dir, a.out, a.draft)
    print(f"wrote {path}" + (f"; left out: {', '.join(missing)}" if missing else ""))
