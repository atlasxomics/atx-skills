# AtlasXomics Data Report — Figure Catalog

Derived from three delivered decks (Demo_003, Demo_002, Demo_001)
and from source inspection of atlasxomics/plots_epigenome, plots_wt, plots_copro,
archrproject_latch and ATX_snap. Every "generator" entry below was confirmed by reading the
code; every "absent" entry was confirmed by searching all five repositories.

## Provenance classes

- **A — workflow-generated**: already written as a PDF/PNG by the secondary-analysis workflow.
  Reproducing means locating the file, not re-running R.
- **B — portable**: produced by a GUI-free Python function in a Plots repo. Verified to contain
  no `lplots`/widget references, so it runs headless unchanged.
- **C — H5 Viewer**: rendered by the proprietary `w_h5` lplots widget in the Plots GUI
  (recognisable in the decks by the `Mean: <x>` / `Genes: <feature>` colorbar annotation).
  Not reproducible from the repositories; a class-B function gives the equivalent plot with
  different cosmetics.
- **D — hand-assembled**: no generator exists in any repository; the scientist builds it from a
  delivered table or image. These are the highest-value targets for automation.

## Catalog

| # | Figure | Deck pages | Class | Generator | Data source |
|---|--------|-----------|-------|-----------|-------------|
| 1 | Cluster spatial maps + UMAP | Demo_003 2, Demo_002 2, Demo_001 2 | C (equiv. B) | `plot_umap_for_samples(color_by='cluster')` — plots_epigenome/welcome/init.py:1468 | `combined_ge.h5ad`: `obs['cluster']`, `obsm['spatial']`, `obsm['X_umap']` |
| 2 | Clustering vs brightfield/H&E | Demo_001 4-6, Demo_002 4-7, Demo_003 3 | D | none | `spatials/<run>/` tissue images + `tissue_positions_list.csv` |
| 3 | Neighborhood enrichment by sample | Demo_001 3, Demo_002 3, Demo_003 4 | B | `plot_neighborhood_groups`, `plotly_heatmap` (init.py:1266, 1863); fallback `squidpy_analysis` (init.py:2176) | `uns['cluster_nhood_enrichment']` (precomputed per sample/condition) |
| 4 | Cluster proportion by sample | Demo_001 3, Demo_002 3 | B | `create_proportion_dataframe` (init.py:208) + `proportion_plot/proportion_plot.py` | `obs['cluster']` x `obs['sample']` |
| 5 | Marker gene / peak / motif heatmaps | Demo_002 8, Demo_003 5, Demo_001 7 | A (equiv. B) | ArchR `plotMarkerHeatmap` -> `figures/heatmaps_all.pdf` (archr_objs.R:524-1234). Portable: `build_archr_like_heatmap_df`, `build_motif_archr_like_heatmap_from_wide`, `select_archr_like_heatmap_features` (init.py:966-1222) | `tables/genes_per_cluster_hm.csv`, `motif_per_cluster_hm.csv`, `marker_peaks_per_cluster.csv` |
| 6 | Gene-activity feature panels (gene x sample + UMAP) | Demo_002 9-10, Demo_001 8-9, Demo_003 6-8 | C (equiv. B) | `plot_umap_for_samples(color_by=<gene>)` + `convert_feature_expression` (init.py:198) | `combined_ge.h5ad` feature matrix |
| 7 | Peaks per cluster / per sample | Demo_003 9, Demo_002 11, Demo_001 10 | A (equiv. B) | ArchR peak-call summary. `PeakAnno` occurs in archrproject_latch/wf/archr_objs.R; `UnionPeaks`/`Number of Peaks` were not in the files searched (all .py in the three Plots repos plus archr_objs.R, archr.R, figure_gallery.py, plotting.py, peaks.py — **not** every R file in archrproject_latch/ATX_snap). The figure is nonetheless workflow-generated: the table and PDF are both delivered in the project output | `tables/Clusters_Peak-Call-Summary.csv`, `Sample_Peak-Call-Summary.csv` — columns `Group,Var1,Freq`, Group already carrying the `C1(n = ...)` label and Freq already in thousands. Figures also delivered as `figures/*_Peak-Call-Summary_Date-*` |
| 8 | Differential motif panels | Demo_002 12, Demo_001 11 | B | `plot_motif_logo` (init.py:641), `plotly_heatmap`, `violin_plots/violin_plot.py`, `plot_ranked_feature_plotly` (init.py:1718) | `combined_motifs.h5ad`, `tables/enrichedMotifs_cluster.csv`, `motif_per_cluster_hm.csv`, `seurat_objects/seqlogo.rds` |
| 9 | Appendix QC table | Demo_003 13, Demo_002 14, Demo_001 14 | D | none | `tables/medians.csv` + `run_metrics/pycistopic_metrics/summary.csv` (read count, duplicate fraction) |
| 10 | Spatial QC maps (TSS enrichment, log10 nFrags) | Demo_001 14 | C (equiv. B) | `plot_umap_for_samples(color_by='TSSEnrichment'\|'log10_nFrags')` | `obs` QC columns in `combined_ge.h5ad` |
| 12 | QC distribution violins | Demo_002 14 | B | `qc_violin` in atx_figures.py, reproducing plots_epigenome/violin_plots/violin_plot.py (`create_violin_data` + `sort_group_categories` + `build_discrete_color_map` + `px.violin(box=True, points=False)`) | `tables/obs.csv`: `log10_nFrags`, `TSSEnrichment`, grouped by `SampleName` |
| 11 | Volcano plots (not used in these three decks) | — | A (equiv. B) | `figures/volcano_plots_<cond>.pdf` (archr_objs.R:926); portable `plot_volcano` (init.py:1628) + `make_volcano_df` (init.py:1223) | `tables/volcanoMarkers_genes_<n>.csv`, `volcanoMarkers_motifs_<n>.csv` |

Static slides (AtlasXPlots blurb, "Additional Data Outputs", Bioinformatics Workflow) carry no data
and are out of scope.

## Conventions observed in the delivered decks

- Significance threshold stated on every differential slide: **P-value <= 0.05 & Log2FC > 0.1**.
- Cluster ordering is lexicographic on the label string, so the axis reads
  C1, C10, C11, C12, C13, C14, C15, C2, C3, ... — matching `sort_group_categories` (init.py:2155).
  Faithful reproduction must preserve this, not sort numerically.
- Neighborhood heatmaps share a fixed colour range of -100 to 100, labelled `zscore`.
- Peak bar charts annotate each x tick with the cluster's cell count, `C1(n = 924)`, and append a
  `UnionPeaks` bar; the y axis is "Number of Peaks (x10^3)" stacked by PeakAnno category.
- Cluster proportion is plotted as a fraction in one deck (Demo_001) and as a raw count in another
  (Demo_002); `create_proportion_dataframe` supports both via `return_type`.

## Correctness rule

Computation must read the full `combined_ge.h5ad` / `combined_motifs.h5ad` objects. The reduced
`*_sm_*.h5ad` objects that Plots loads cast `.X` to float16 for the ATAC modality, so any value
derived from them is imprecise.

## Verified against Demo_001 (2026-11-11)

All catalog inputs exist in `latch:///epi_analysis_archr/Demo_001/`. Checks made against the
delivered deck:

- **QC table rounding is reproducible exactly.** `tables/medians.csv` holds
  `median_TSS, median_fragments, median_FRIP` per run; the deck rounds TSS and FRiP to two
  decimals and median fragments to the nearest hundred (two decimals for TSS and FRiP,
  nearest hundred for median fragments). Reproduced exactly on a delivered deck.
- **The QC table needs a join across two workflows.** `Number of reads` and
  `Fraction Duplicates` are not in `medians.csv`; they come from
  `fastq2frags/<run_id>/run_metrics/pycistopic_metrics/summary.csv`. The generator must resolve
  the preprocessing directory for each run, not just the analysis project.
- **`condition` in medians.csv is the deck's sample label** (CondA / CondB / CondC).
- **Object sizes force a table-first design.** In Demo_001 `combined_ge.h5ad` is 19 GB and
  `combined_motifs.h5ad` 685 MB, against 64 GB RAM and a 100 GB `/data`. Anything obtainable
  from `tables/obs.csv`, `spatial.csv` or `UMAPHarmony.csv` (clusters, QC metrics, embeddings,
  proportions, neighborhood input) must be read from those, not from the h5ad. Gene- or
  motif-level feature plots are the only figures that require the matrices, and those should use
  backed reads — `read_backed_h5ad` and `tests/test_backed_plotting.py` in plots_copro are the
  team's existing precedent.

## The QC slide differs by project

"Samples pass QCs for each assay" is not one layout. The table at the top is common, but the
panels beneath it differ:

- **Demo_002 (p14)**: two violin plots — "Distribution of log10_nFrags by SampleName" and
  "Distribution of TSSEnrichment by SampleName", four samples each, x ordered
  AP_4TFF_DM, AP_55PK_DR, AP_JTZH_XA, AP_S3DN_2B.
- **Demo_001 (p14)** and **Demo_003 (p13)**: spatial tissue maps coloured by
  TSSEnrichment and log10_nFrags, with the H5-viewer colorbar annotation (`Mean: 5.34`).

So there are two panel styles, not three, and the spatial-map style is the more common one.
A generator must take the style as a parameter: `qc_violin` (catalog #12) or
`plot_umap_for_samples` on a QC obs column (catalog #10).

The table itself also varies by project, which the formatting rules above do NOT cover:

- Demo_003 adds an **Assay** column ("ATAC + RNA") absent from the other two decks.
- Demo_003 prints FRiP to **three** decimals (0.211) where Demo_001 and Demo_002 use two. The
  rounding rules verified above hold for Demo_001 and Demo_002 only; treat column set
  and precision as per-deck settings rather than constants.
- Demo_003 `Sample` is an external specimen label (`NL39-6wk-1069349-S6-R4`), matching neither
  `medians.csv:condition` nor `obs.csv:SampleName` patterns seen elsewhere.

### Violin generator verification (Demo_002, 2026-11-11)

Generated from 159,513 cells. Figure titles and x-axis order match the delivered slide exactly.
Per-sample medians of `TSSEnrichment` from the violin data reproduce `medians.csv:median_TSS` for
all four samples to three decimals (5.649, 4.746, 3.989, 5.876), confirming the violins and the
appendix table are computed from the same per-cell values.

## Morphology overlays (catalog #2) — measured geometry, 2026-11-16

Images live under `latch://13502.account/spatials/<SHORT_ID>/`, where SHORT_ID drops the
leading zeros from the run (RUN1_SEQ001 -> RUN1). `latch://atx-illumina.mount/Images` is
the documented alternative root; its folders are not named by run id, so the spatials root is
the one that resolves.

**The `_fix` TIFs are not in the tixel coordinate space.** All three Demo_001 runs have one
(RUN1 3261x6456, RUN2 4607x8436, RUN3 3258x5447), but `tissue_positions_list.csv` is in the
AtlasXBrowser crop space defined by `metadata.json:crop_area` — [745, 572, 8431, 8258], i.e.
7686x7686 for RUN1, with observed coordinates x 224-7394, y 153-7373. The `_fix` images match
neither that extent nor a uniform scaling of it, and they differ in size per run. No registration
transform is stored anywhere in the run output, and `spatial/figure/` holds only postB.tif and
postB_BSA.tif, not an overlay.

The overlay is therefore drawn on `spatial/tissue_hires_image.png`, which IS that crop space
scaled by `scalefactors_json.json:tissue_hires_scalef` (0.2602 for RUN1, 0.2527 RUN2,
0.2540 RUN3) — exact geometry, no fitting. `check_fix_image()` still reports whether a `_fix`
image exists and returns a note when it does, so the scientist can supply the transform or
confirm the hires image is acceptable.

**H&E images** are not tracked in a standard location; the generator cannot find them and the
scientist must supply a path or upload.

## Open discrepancy: the H5 viewer `Mean:` annotation

Spatial QC maps reproduce the delivered colorbar RANGE exactly for both metrics, but not
the `Mean:` figure printed under it, which comes out slightly different each time. Matching
endpoints with a drifting mean indicates the same underlying data and a different
denominator in the viewer — displayed cells, or one sample rather than the pooled object.
Not resolved; deferred.

## Feature panels: how they were verified

Gene and motif columns are read from `combined_ge.h5ad` and `combined_motifs.h5ad`, opened
backed, one column at a time. For each feature, check the value range and mean against the
delivered panel's colourbar endpoints and its `Mean:` annotation — both reproduce.

`a[:, feature].X` on a backed AnnData does NOT materialise the column (it yields a length-1
object array); the view must be pulled with `.to_memory()` first. Cross-check the per-run
tixel counts in the morphology overlays against `Sample_Peak-Call-Summary.csv`; they match
exactly, which is a cheap way to confirm the barcode join is right.
