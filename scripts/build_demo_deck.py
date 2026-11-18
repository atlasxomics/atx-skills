"""Assemble the Demo_001 Data Report deck from saved figure artifacts.

Kept as a script rather than inline cells so the slide list survives a kernel reset and
can be re-run or edited without reconstructing it from memory. Figures are resolved by
FILENAME against the artifact store, so re-rendering a figure and saving a new version is
picked up automatically.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import atx_deck

TEMPLATE = "atx_template.pptx"
OUT = "Demo_001_Data_Report_branded.pptx"
QC_TABLE = "qc_table_Demo_001.csv"


def resolve(host):
    """Map figure filename -> local path via the artifact store (latest version)."""
    items = host.artifacts()["artifacts"]
    latest = {}
    for a in items:
        latest.setdefault(a["filename"], a.get("version_id") or a.get("latest_version_id"))
    return latest


def build(host, out_path=OUT):
    vid = resolve(host)
    A = host.artifact_path

    def F(name):
        if name not in vid:
            raise KeyError(f"no artifact named {name!r}; render it first")
        return A(vid[name])

    # Rendered as ONE plotly figure with a shared coloraxis by
    # atx_figures.neighborhood_grid() -- not montaged from the per-sample PNGs, which
    # would give each panel its own duplicate colourbar.
    nb_grid = F("Demo_001__neighborhood_grid.png")
    qc = pd.read_csv(F(QC_TABLE))

    slides = [
        {"title": "Clustering shows distinct spatial features across tissues",
         "columns": [[F("Demo_001__cluster_spatial_h5.png")],
                     [F("Demo_001__umap_cluster.png"), F("Demo_001__umap_sample.png")]],
         "column_titles": ["Spatial UMAP", "UMAP"]},
        {"title": "Spatial distribution of clusters differs between samples",
         "subtitle": ("Neighborhood analysis identifies how often each cluster is adjacent "
                      "to another cluster (higher number is more related)"),
         "columns": [[nb_grid], [F("Demo_001__cluster_proportion.png")]],
         "column_titles": ["Neighborhood analysis by sample", "Cluster proportion by sample"]},
        {"title": "Clustering aligns to morphology: CondA|RUN1",
         "images": [F("Demo_001__morphology_RUN1_SEQ001.png")]},
        {"title": "Clustering aligns to morphology: CondB|RUN2",
         "images": [F("Demo_001__morphology_RUN2_SEQ002.png")]},
        {"title": "Clustering aligns to morphology: CondC|RUN3",
         "images": [F("Demo_001__morphology_RUN3_SEQ003.png")]},
        {"title": "Numerous differential genes, motifs and peaks between clusters",
         "subtitle": "P-value <=0.05 & Log2FC >0.1 for all plots.",
         "columns": [[F("marker_heatmap_p1.png")], [F("marker_heatmap_p2.png")],
                     [F("marker_heatmap_p3.png")]],
         "column_titles": ["Gene activity scores", "Peak enrichment", "Motif enrichment"]},
        {"title": "Gene Activity Scores identify distinct biological clusters",
         "subtitle": "Alternative clustering and annotation could provide further insights into biology",
         "columns": [[F("Demo_001__tracks_Cd19.png")], [F("Demo_001__tracks_Cd3e.png")],
                     [F("Demo_001__tracks_Hmox1.png")], [F("Demo_001__umap_cluster.png")]],
         "column_titles": ["Cd19", "Cd3e", "Hmox1", ""],
         "column_weights": [1, 1, 1, 0.85]},
        # Each marker gene is shown as its spatial panels AND a UMAP coloured by the same
        # gene, two blocks per column, with the annotated cluster UMAP filling the fourth
        # cell. The spatial/UMAP pairing is the point of the slide -- see SKILL.md.
        {"title": "Gene Activity Scores identify distinct biological clusters",
         "subtitle": "Gene markers from https://panglaodb.se/markers.html",
         "grid": [[F("Demo_001__block_Cd3e.png"), F("Demo_001__block_Hmox1.png")],
                  [F("Demo_001__block_Cd19.png"),
                   F("Demo_001__umap_cluster_annotated.png")]]},
        {"title": "Total number of peaks vary across clusters",
         "subtitle": "N= number of datapoints per cluster",
         "images": [F("Demo_001__peaks_per_cluster.png"), F("Demo_001__peaks_per_sample.png")]},
        # MeanDiff/FDR come from the ArchR Compare workflow output, never recomputed here.
        # Each motif pairs its spatial panels (selected cells only) with its CisBP logo.
        {"title": "Distinct differential motifs found between Cd19 populations",
         "subtitle": ("Selected datapoints with Cd19 gene activity score higher than 0.6 and "
                      "compared differential motifs. No Cd19 high was found in CondC"),
         # The spatial panel and its sequence logo are a LIST, so they are stacked as two
         # separate pictures -- the scientist can resize the logo without dragging the
         # panel above it, which a montaged single image would not allow.
         "grid": [[[F("Demo_001__motif_sel_Ascl1-780.png"),
                    F("Demo_001__logo_Ascl1-780.png")],
                   F("Demo_001__motif_violin_Ascl1-780.png"),
                   F("Demo_001__cd19_high_spatial.png")],
                  [[F("Demo_001__motif_sel_Runx1-712.png"),
                    F("Demo_001__logo_Runx1-712.png")],
                   F("Demo_001__motif_violin_Runx1-712.png"),
                   F("Demo_001__motif_rank.png")]],
         # Slide notes carry the provenance of the statistics: which workflow produced
         # them, which groups were compared, and the literature the biology rests on.
         "notes": ("MeanDiff / FDR from ArchR getMarkerFeatures via the Compare "
                   "workflow; GroupA vs GroupB as assigned by that run. Add the "
                   "citations supporting the biological claim here.")},
        {"title": "Appendix", "section": True},
        # Both QC metrics are shown, each as its own titled column of per-sample maps.
        {"title": "Samples pass QCs for each assay", "table": qc,
         "columns": [[F("Demo_001__qc_spatial_TSSEnrichment.png")],
                     [F("Demo_001__qc_spatial_log10_nFrags.png")]],
         "column_titles": ["TSS Enrichment", "log10nFrags"]},
        {"title": "Bioinformatics Workflow",
         "subtitle": "For more information go to our website",
         "images": [F("Demo_001__bioinformatics_workflow.png")]},
    ]
    return atx_deck.build_deck(
        slides, out_path, template=TEMPLATE,
        deck_title="Chart the unexplored depths of spatial biology",
        subtitle="Data Report")
