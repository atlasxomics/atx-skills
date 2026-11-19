# AtlasXomics data processing pipeline

Distilled from the AtlasXomics docs (<https://docs.atlasxomics.com/en/latest/>, source in
`atlasxomics-docs/`). Read this when you need to know **where a number or file came from**,
what produced it, or what the parameters were. The figure generators do not depend on this
file; it explains the data they read.

Focus is the epigenomics path, which is what Data Reports are built from. Whole
transcriptome and co-profiling are summarised at the end.

## Terms

- **Run**: one tissue section assayed by DBiT-seq. `run_id` is `Dxxxxx_NGxxxxx`. `D…` is the
  slide, `NG…` the sequencing library.
- **Project**: several Runs merged for analysis. The `<project>` under `epi_analysis_*`.
- **Tixel**: one spatial barcode position on the grid (50×50 or 96×96). AnnData rows are
  tixels, not cells, even where the code says "cell".
- **Fragments file**: `fragments.tsv.gz`, bgzipped BED-like (chrom, start, end, barcode,
  count), tabix-indexed. The handoff between preprocessing and everything downstream.
- **Spatial folder**: AtlasXBrowser output per Run. `tissue_positions_list.csv` (barcode,
  on/off-tissue flag, row, col, x, y), `tissue_lowres_image.png`, `tissue_highres_image.png`,
  `metadata.json`.
- **condition**: optional per-Run label set at launch. Drives the `condition_*` tables and
  volcano outputs, and is the sample label in the deck (CondA / CondB …).

## Epigenomics (spatial ATAC-seq and CUT&Tag)

```text
 read-2 FASTQ ─► BarcodeQC (optional) ─► cleanTable.csv / xtalkTable.csv
                                               │ ds_table
 FASTQ (R1,R2) ─────────────────────────► ATX epigenomic preprocessing ─► fastq2frags/<run_id>/
 (or aligned CRAM ─► cram2frags ─► cram2frags/<project>/)                      │ fragments
                                                                               ▼
 tissue images ─► AtlasXBrowser ─► spatials/<run_id>/ ──────────► optimize archr | optimize_snap
                                                                   (parameter sweep, first spatial maps)
                                                                               │
                                   missing lanes ─► ATX impute ─┐              │
                                   out-of-bounds ─► fragment fixer ─► re-optimize
                                                                               ▼
                                                    create ArchRProject | ATX_snap
                                                    epi_analysis_archr/  epi_analysis_snap/
                                                                               │
                         ┌─────────────────────────────┬───────────────────────┤
                         ▼                             ▼                       ▼
                 Epigenomics Plots ──launch──► compare clusters        atx_glue (co-profiling)
                 (reads *_sm_*.h5ad)           compare_outs/
```

The ArchR and SnapATAC2 branches are alternatives, chosen per project: `optimize archr` →
`create ArchRProject`, or `optimize_snap` → `ATX_snap`. Both emit the same `tables/` set, so
report figures are generated the same way from either.

### 0. BarcodeQC (optional, before preprocessing)

Repo `barcodeqc_wf`. Reads read 2 only, matches barcodes to the whitelist and makes an HTML
report of match rates and their spatial distribution. It also **detects** two microfluidic
artifacts and writes the correction tables that preprocessing applies:

- **CLEAN** (Channel-Linked Efficiency Artifact Normalization): rows/columns whose summed
  counts deviate from a LOWESS trend by more than ±0.5 normalised difference are masked;
  each masked tixel gets a downsampling fraction from its neighbours. `cleanTable.csv`
  (`mer`, `finalDownFract`). A regional variant works in 1/16 sub-regions.
- **Cross-talk**: bleed-through between adjacent channels. Positions above
  `fold_above_background` are flagged in `xtalkTable.csv`; in preprocessing their fragments
  are replaced by those of the nearest clean neighbour (±1, then ±2, no diagonals).

If BarcodeQC runs first, correction happens inline and there is no post-optimisation repair loop.

### 1. ATX epigenomic preprocessing → `latch:///fastq2frags/<run_id>/`

Repo `atlasXomics_pipeline_wf`. Read 1 is genomic; read 2 is
`linker1 | barcodeA | linker2 | barcodeB | genomic`.

| Task | What it does |
|---|---|
| `filtering` | bbduk on read-2 linkers; reads with >3 mismatches in linker 1 or 2 are dropped (`skip1`/`skip2` disable each). Applies CLEAN / cross-talk if `ds_table` given. Bulk and no-ligation-bulk modes substitute random barcodes. |
| `alignment` | Chromap against the `species` reference, barcodes from read 2. Barcodes get a `-1` suffix; out-of-bounds fragments removed using chrom sizes; optional blacklist filtering. |
| `metrics_task` | `atac_analysis_summary_prep.py` (GENCODE annotation, `--min-fragments 100`): per-tixel fragments, TSS enrichment, FRiP, fragment sizes, duplicate and mito rates, a coverage bigWig, a peak list, QC histograms, a lane-QC plot. `fragReport.py` builds the HTML report and writes the run-level `run_metrics/summary.csv` (the row pushed to the Registry and SLIMS). |
| `statistics` | **Off by default** (`skipExtendedMetrics=True`). When on: pycisTopic per-cell QC, optional MACS consensus peaks (`peak_calling`), and `singlecellsummary.py` writing its own `summary.csv` and `singlecell.csv` to `run_metrics/pycistopic_metrics/`. |

Key outputs:

```text
fastq2frags/<run_id>/
├── fragments.tsv.gz, fragment_analysis_report.html, fragmentscell.bw   # top-level copies
├── filtered_fastqs/      <run_id>_linker2_R{1,2}.fastq.gz, <run_id>_l{1,2}_stats.txt
├── chromap_output/       fragments.tsv.gz(.tbi), aln.bed, chromap_log.txt
│                         (fragments.final.tsv.gz when cross-talk correction ran)
└── run_metrics/          summary.csv, fragments_adata_obs.csv (per tixel), lane_qc.png,
                          fragments_*_hist.png, fragments_peakList.bed.gz, …
    └── pycistopic_metrics/   only if the statistics task ran: singlecell.csv, cistopic_cell_data.csv, …
```

`run_metrics/summary.csv` is where the report QC table gets **Number of reads** and
**Fraction Duplicates**, and is the path `atx_project.py` reads. Most runs have no
`pycistopic_metrics/` folder; if one exists, its `summary.csv` is a different, pycisTopic-based
file. Raw FASTQs are not delivered to customers; `filtered_fastqs/` stands in for them.

`cram2frags` (repo `cram2frag_wf`) is the alternative entry for aligned single-end CRAM
(e.g. Ultima). It writes `cram2frags/<project>/fragments.sort.bed.gz` plus `frag_metrics/`, and
the docs list no `summary.csv` for it, so read count / duplicate fraction must come from elsewhere for those runs.

### 2. AtlasXBrowser → `latch:///spatials/<run_id>/`

Image app (Latch Pod or local), not a workflow. Registers the tissue image to the barcode grid
and marks tixels on/off tissue. Its output is the Spatial folder that every downstream step
takes as `spatial_dir`. Separate docs: <https://atlasxbrowser-docs.readthedocs.io/>.

### 3. Optimisation (parameter sweep)

All Runs are merged once, then each parameter combination is evaluated in parallel. Output is
galleries for choosing parameters, not report figures. This is also the first point where
spatial maps exist, which is when failed lanes or out-of-bounds coordinates show up.

| | optimize archr (`archr_latch`) | optimize_snap (`optimize_snap`) |
|---|---|---|
| Output | `epi_optimize_archr/<project>/` | `epi_optimize_snap/<project>/` |
| Cell filter | `min_TSS` 2.0, `min_frags` 0 | `min_tss` 2.0, `min_frags` 10 |
| Features | 5 kb tiles; IterativeLSI (`lsi_varfeatures` 25000, `lsi_resolution` 0.5, 2 iterations, 30 dims), optional Harmony | 5 kb tiles; `n_features` 25000, spectral embedding `n_comps` 30 |
| Clustering | Seurat via `addClusters`, `clustering_resolution` 1.0, ≤25 clusters, ≥20 cells | Leiden, `resolution` 1.0, ≥20 cells |
| Swept | `lsi_resolution`, `lsi_varfeatures`, `clustering_resolution`, `maximum_dims` | `n_features`, `n_comps`, `resolution`, `varfeat_iters` |
| Key files | `medians.csv`, `figures/{umap,spatialdim,qc}_plots_<NNN>.png`, `cluster_csvs/` | `combined.h5ad`, `medians.csv`, `all_umaps.html`, `all_spatialdim.html`, `spatial_qc.html` |

### 4. Repair helpers (only when the maps show a problem)

- **ATX impute** (`ATX_impute`) → `impute/<dir>/`: fills missing rows/columns (a failed
  channel) from neighbouring tixels, guided by cluster labels.
- **fragment fixer** (`fragment_fixer`) → `ff_outs/<dir>/`: strips fragments outside the
  genome's chromosome bounds. New preprocessing does this already; needed only for old runs or
  external fragment files.

After either, re-optimise. A project whose fragments came from `impute/` or `ff_outs/`
rather than `fastq2frags/` has had its counts altered; note this if QC numbers look off.

### 5. Secondary analysis → the objects Data Reports are built from

**create ArchRProject** (repo `archrproject_latch`, script `archr_objs.R`) →
`latch:///epi_analysis_archr/<project>/`. One `archr_task` does everything in ArchR + Seurat:
Arrow files and cell filtering, IterativeLSI, clustering, UMAP, spatial alignment, then

- **gene accessibility**: `GeneScoreMatrix`. Tn5 signal over the gene body and a regulatory
  window, exponentially distance-weighted, bounded by neighbouring genes, imputation-smoothed.
  A per-tixel proxy for expression.
- **peaks**: MACS2 per group, then reproducible peak sets.
- **motifs**: chromVAR deviations (`addDeviationsMatrix`, per-tixel z-scores against
  GC/accessibility-matched background peaks) and hypergeometric enrichment in marker peaks
  (`peakAnnoEnrichment`).
- **differential testing**: `getMarkerFeatures` per cluster / sample / condition. Genes use a
  t-test against background matched on TSS and log2(fragments); peaks use Wilcoxon.

**ATX_snap** (repo `ATX_snap`) → `latch:///epi_analysis_snap/<project>/`. SnapATAC2 does the
tiles, spectral embedding and **Leiden clustering**, and writes the labels to `tables/obs.csv`.
ArchR then rebuilds a project from the fragments, **adopts those clusters**, and computes gene
scores, peaks, chromVAR motifs and marker stats exactly as above. Squidpy adds neighbourhood
enrichment and Moran's I spatially variable genes. Heavy stages checkpoint to `checkpoints/`,
which is deleted only after every branch succeeds. If it is present, the run failed.

Defaults common to both: `tile_size` 5000, `min_TSS` 2.0, ≤25 clusters (ArchR),
≥20 cells per cluster, chrY excluded unless `include_y_chromosome`.
`tables/input_parameters.csv` records what was actually used for a given project.

Output layout (both branches; ArchR names shown, SnapATAC differences in brackets):

```text
epi_analysis_{archr,snap}/<project>/
├── <project>_ArchRProject/        [snap: absent; checkpoint deleted after success]
├── anndata/
│   ├── combined_ge.h5ad           full gene accessibility (Demo_001: ~19 GB)
│   ├── combined_motifs.h5ad       full motif deviations
│   ├── combined_sm_ge.h5ad, combined_sm_motifs.h5ad   Plots copies, float16, NOT for computation
│   ├── <run_id>_{g,m}_converted.h5ad
│   └── [snap: combined.h5ad — tile matrix, spectral, Leiden]
├── seurat_objects/                <run_id>_SeuratObj(Motif).rds, combined(_m).rds, seqlogo.rds,
│                                  <run_id>_BP/ on very large projects (part of the object)
├── {cluster,sample,condition}_coverages/   bigWig   [snap: *_coverage/, bedgraph]
├── {cluster,sample,condition}_peak_beds/          [snap: peaks.bed, filtering_summary.csv]
├── figures/                       heatmaps_all.pdf, volcano_plots_<cond>.pdf, *_Peak-Call-Summary_*, …
├── tables/
└── Launch_Plots/artifact.json
```

`tables/`, the files report figures are read from:

| File | Content |
|---|---|
| `obs.csv`, `spatial.csv`, `UMAPHarmony.csv` [snap: `X_umap.csv`, `spectral.csv`] | per-tixel metadata and clusters, spatial coordinates, embedding |
| `medians.csv` | per-Run median fragments, FRiP, TSS, with `condition` |
| `input_parameters.csv` | parameters this project ran with |
| `ranked_genes_per_{cluster,sample,condition}.csv` | all genes: log2FC, mean diff, FDR |
| `genes_per_{cluster,sample}_hm.csv` | marker genes, `Pval ≤ 0.05` & `Log2FC ≥ 0.10`, clipped to ±2: the gene heatmap |
| `marker_peaks_per_{cluster,sample,condition}.csv` | Wilcoxon marker peaks, `Pval < 0.05` & `Log2FC ≥ 0.1` |
| `complete_peak_list_{cluster,sample}.csv` | unfiltered peaks + annotation, closest gene. Group index = cluster + 1 (R vs Python) |
| `enrichedMotifs_{cluster,sample,condition}.csv` | hypergeometric motif enrichment (`mlog10Padj`, …) |
| `motif_per_{cluster,sample}_hm.csv` | `-log(Padj)` scaled 0–100 per group: the motif heatmap |
| `volcanoMarkers_{genes,motifs}_<n>.csv` | condition vs rest, only when conditions were set |
| `svg_genes.csv`, `svg_motifs.csv` | spatially variable features (Moran's I) |
| `Clusters_Peak-Call-Summary.csv`, `Sample_Peak-Call-Summary.csv` | ArchR peak-call summary: peaks per cluster/sample |

### 6. Plots and compare clusters

**Epigenomics Plots** (repo `plots_epigenome`) is the interactive Latch Plots notebook,
launched from `Launch_Plots/artifact.json`. It loads only the `_sm` objects and the coverage
folders. Tabs: H5 Viewer, Heatmap, Violin, Volcano, Proportion, Track Browser, Motif Logo,
Neighborhood, Compare, Bulk Score Genes. `atx_plots_lib.py` in this skill is lifted from it.

**compare clusters** (repo `combined_cluster_wf`) → `latch:///compare_outs/<name>/`. Pairwise
`getMarkerFeatures` between GroupA and GroupB (by cluster / condition / sample, barcode
lists, or a file) on an ArchRProject. At most 500 cells per group by default; `closest` is
True for a direct launch but False from Plots. Writes `gene_results/`, `peak_results/`,
`motif_results/`, `coverages/` and `UpdateClustName_by_barcode.csv`. Method, defects and how
to read it: `references/differential_analysis.md`. It is normally launched from the Plots
Compare tab. The skill can launch it too, but only after the scientist confirms. Either way,
a deck's comparison slide is fetched from a finished run, never recomputed.

### Other epigenomic helpers

- **atx_convert** (`atx_convert`) → `converted/<dir>/`: AnnData ↔ Seurat ↔ BPCells, optional `_sm`.

## Whole transcriptome (spatial RNA-seq)

```text
FASTQ ─► (BarcodeQC) ─► RNAQC ─► rna_preprocessing/<sample>_star/
images ─► AtlasXBrowser ──────► optimize_wt ─► rna_analysis/<project>/ ─► WholeTxn Plots
```

- **RNAQC** (`rnaqc_wf`): read-2 barcode processing, optional CLEAN / cross-talk, STARsolo
  alignment, MultiQC, and a FastQ-Screen contamination screen of read 1. Optional quick
  Scanpy sweep.
- **optimize_wt**: filter, normalise, HVGs, optional Harmony, cluster with Scanpy Leiden or
  STAGATE, marker genes, squidpy neighbourhood enrichment, spatially variable genes. This
  **is** the secondary analysis; there is no separate object-building step.
- Objects live in **per-parameter-set folders**, e.g.
  `rna_analysis/<project>/set1_backend-scanpy_cr1-0-nc30-nn15-md0-5-sp1-0/combined.h5ad`,
  each with its own `Launch_Plots/`. The top level has only `medians.csv` (UMIs, genes,
  % mito), `metadata.csv`, `svg_genes.csv` and galleries.
- The WT `_sm` object is float32 with exact counts in `layers["counts"]`. It is still
  incomplete, so it is not for computation either.

## Co-profiling (epigenome × transcriptome)

```text
create ArchRProject | ATX_snap ─┐
                                ├─► atx_glue ─► copro_integration_analysis/<project>/ ─► CoPro Plots
optimize_wt ────────────────────┘
```

- **atx_glue** (`atx_SpatialGlue`): aligns the transcriptome AnnData and the gene-accessibility
  AnnData to a shared grid, trains SpatialGlue and sweeps Leiden resolutions for **joint
  clusters** (`sg_clusters`, min cluster size 200). Optionally adds per-cluster coverage
  bigWigs, RNA–ATAC Spearman correlations and ArchR Peak2Gene links.
- Outputs `rna_copro.h5ad`, `atac_gs_copro.h5ad` (+ `_sm`), `SpatialGlue_model.pickle`, and
  `tables/`, `figures/`, `coverages/`, `peak2gene/`.
- Both modality paths must finish first. Co-profiling report panels are out of scope for this
  skill today (see README).

## Workspace map

| Directory | Produced by | Unit |
|---|---|---|
| `spatials/` | AtlasXBrowser | run |
| `fastq2frags/` | ATX epigenomic preprocessing | run |
| `cram2frags/` | cram2frags | project |
| `epi_optimize_archr/`, `epi_optimize_snap/` | optimize archr / optimize_snap | project |
| `epi_analysis_archr/`, `epi_analysis_snap/` | create ArchRProject / ATX_snap | project |
| `compare_outs/` | compare clusters | comparison |
| `impute/`, `ff_outs/`, `converted/` | ATX impute / fragment fixer / atx_convert | run or dir |
| `rna_preprocessing/`, `rna_analysis/` | RNAQC / optimize_wt | sample / project |
| `copro_integration_analysis/` | atx_glue | project |

Older workspaces may instead use `Raw_Data/`, `Processed_Data/`, `Optimized_Data/`.
