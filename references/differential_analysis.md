# Differential analysis (compare clusters)

How AtlasXomics computes differential genes, peaks and motifs between two groups of cells,
and how to read the result without misreporting it. Read this before quoting any
GroupA-vs-GroupB number or building a differential slide. `scripts/atx_compare.py`
implements everything below.

**Not every differential question needs compare clusters.** Each cluster, sample and
condition against all other cells is already in the project's `tables/` (archrproject and
ATX_snap `archr_objs.R`):

| Question | Table | Reader |
|---|---|---|
| Marker genes per group | `ranked_genes_per_{cluster,sample,condition}.csv` | `Project.top_ranked_genes` |
| Motifs enriched in each group's marker peaks (hypergeometric) | `enrichedMotifs_{cluster,sample,condition_<n>}.csv` | `Project.top_motifs` |
| Motif activity (chromVAR deviation) higher in one condition, all cells and per cluster (Wilcoxon) | `volcanoMarkers_motifs_<n>_<condition>.csv` | `Project.top_condition_motifs` |

The volcano tables are filtered before they are written (rows with adjusted p == 0 or
|MeanDiff| >= 1.2 are dropped), and their `Significance` column uses raw p < 0.01. Compare
clusters is for two groups someone defines, such as a region, a set of clusters, or one
cluster's cells in two conditions.

Sources (read in full, not summarised from the docs):

- workflow `atlasxomics/combined_cluster_wf` v0.10.8 (`8cdf2d1`): `wf/compare_task.py`,
  `wf/compare_clusters.R`
- Plots `atlasxomics/plots_epigenome` `compare/` (last change `70a0b032824c`, main at
  `3d05a74f7683`): `select_cells.py`, `check_selection.py`, `launch_workflow.py`,
  `fetch_results.py`, `gene_results.py`, `motif_results.py`, `track_browser.py`
- docs page <https://docs.atlasxomics.com/en/latest/epigenomics/compare-clusters/>. Its
  output layout is out of date; the layout below is what the workflow writes.
- checked against a real run, `latch:///compare_outs/jm_test/` (828 vs 1575 cells)

## The rule

**Fetch the statistics, never recompute them.** The numbers come from
`ArchR::getMarkerFeatures` with bias matching and a cell cap. A scanpy or scipy test on the
same cells gives different numbers, and the deck will not match anything the customer
can reproduce. If no run exists for the comparison, launch one, but only after the
scientist approves it (see "Launching a run" below).

## End to end

```text
Plots Compare tab                                   compare clusters workflow
─────────────────                                   ─────────────────────────
pick a categorical obs column (< 30 values)
pick GroupA value, GroupB value
  └─ check: both non-empty, no overlap
  └─ upload {"groupA":[barcodes], "groupB":[...]}   groupings = that file (File mode)
     to <project>/compare_config_<dataset>_         archrproject = <project>_ArchRProject
        A<nA>_B<nB>_<YYYYmmdd_HHMMSS>.json   ──────► compare_task.py → compare_clusters.R
"Launch Workflow" (pins workflow 0.10.6)                 │
"Fetch Workflow Results" reads all_*.csv  ◄──────── latch:///compare_outs/<run>/
volcano + rank plot per feature, IGV of the 2 bigwigs
```

Plots offers only the categorical columns of the loaded `_sm` gene AnnData, including custom
annotations made with the H5 Viewer lasso. A continuous selection (e.g. "Cd19 gene score
> 0.6") has to be thresholded into a categorical column first.

Compare needs an **ArchRProject**, so it runs only on `epi_analysis_archr` projects.
ATX_snap deletes its ArchR checkpoint after a successful run, and the Plots tab refuses
without one.

## How groups are defined (workflow inputs)

| Mode | How cells are chosen |
|---|---|
| **File** (what Plots sends) | JSON `{"groupA": [...], "groupB": [...]}` of `<run_id>#<barcode>-1` names; treated as Barcodes |
| **Barcodes** | comma-separated barcode strings for each group |
| **Groupings** | `clusterA/conditionA/sampleA` and `…B`; comma lists, ranges like `C5-C7` expanded to C5,C6,C7; cluster, condition and sample filters combine with AND; an empty field means no filter. `multipleA/B` matches every condition containing the value as a whole word (case-insensitive). **Up to v0.10.8, `_` was not a word break (`WT` missed `WT_young`), and when nothing matched the condition filter silently became empty, putting every condition in the group.** From v0.10.9 any non-alphanumeric run separates words and no match stops the run. |

Cells in neither group are dropped. A cell in both is labelled GroupA (Plots blocks overlap
before this point). `UpdateClustName_by_barcode.csv` records the exact assignment. It is
the ground truth for which cells were compared, so use it rather than re-deriving a
selection.

GroupA/GroupB are **positional**. The only record of the intended sizes is the
`_A<nA>_B<nB>_` part of the compare_config file name in the project root. The run name is
not recorded there, so match on counts and time: `find_compare_configs()` and
`check_group_mapping()`.

## Cell cap and background: three parameters that change the numbers

| Parameter | Workflow default | Plots default | Effect |
|---|---|---|---|
| `max_cells` | 500 | 500 | `maxCells` for every test: each group is capped at this many cells |
| `use_max_possible_cells` | — | False | when True, `maxCells = min(nA, nB)` and `max_cells` is ignored |
| `closest` | **True** | **False** | ArchR `closest`: background cells matched as nearest neighbours in bias space vs. sampled |

**`closest` differs between a Plots launch and a direct workflow launch, and the run
output does not record which was used** (only the Latch execution's inputs do). Two runs
on the same cells can disagree for this reason alone. Ask, or check the execution,
before comparing runs.

## The tests

All four use `groupBy = "UpdateClustName"`. The bias-matched tests use
`bias = c("TSSEnrichment", "log10(nFrags)")`, so each group's background is matched on
TSS enrichment and depth.

| Feature | Call | Matrix | Test | Effect size | Output |
|---|---|---|---|---|---|
| Genes | `getMarkerFeatures` | GeneScoreMatrix (gene accessibility) | t-test | `Log2FC` | `gene_results/all_genes.csv` |
| Peaks | `getMarkerFeatures` | PeakMatrix | Wilcoxon | `Log2FC` | `peak_results/all_peaks.csv` |
| Motif activity | `getMarkerFeatures`, `useSeqnames = "z"`, `normBy = "none"` | MotifMatrix (chromVAR deviation z-scores) | Wilcoxon | `MeanDiff` of z | `motif_results/all_motifs.csv` |
| Motif enrichment | `peakAnnoEnrichment` on the peak test | motif annotation of differential peaks | hypergeometric | `mlog10Padj` | `motif_results/{up,down}Regulated_motifs.csv` |

`FDR` is ArchR's Benjamini–Hochberg value. The workflow adds raw `Pval` as an extra column.

**There are two different "differential motif" answers.** `all_motifs.csv` asks whether a TF's
accessibility (chromVAR deviation) differs between the groups. It is per cell, Wilcoxon,
and what the Plots tab and the delivered motif slide show. The up/down tables ask whether a
motif is over-represented in the peaks that are more open in one group. That is a property
of peak sets with a much looser input cut (`Pval <= 0.1`). Never mix the two on one slide,
and say which one a statement rests on.

Coverage tracks: `getGroupBW(normMethod = "ReadsInTSS", tileSize = 100, maxCells = 1000,
ceiling = 4)` writes one bigwig per group to `coverages/`. They are capped at 1000 cells,
not `max_cells`, so they can include cells the tests did not.

Peak and motif steps are optional. If one fails, the workflow logs it to
`analysis_errors.log` and still succeeds. A run with gene results but no motif results is
normal, so report the gap rather than treating it as "no differential motifs".

## Reading the output

```text
compare_outs/<run>/
├── UpdateClustName_by_barcode.csv   cell_barcode, UpdateClustName (GroupA | GroupB)
├── gene_results/   all_genes.csv, marker_genes.csv, volcano_gene.pdf
├── peak_results/   all_peaks.csv, MA_peaks.pdf
├── motif_results/  all_motifs.csv, marker_motifs.csv, volcano_motif.pdf,
│                   up/downRegulated_motifs.csv, up/downRegulated_motif_enrichment.pdf
└── coverages/      GroupA-TileSize-100-normMethod-ReadsInTSS-ArchR.bw, GroupB-…bw
```

`all_*.csv` columns: `group, group_name, name, Log2FC, FDR, MeanDiff, Pval` (+ `seqnames,
start, end, strand, idx` for genes; motifs have `seqnames = "z"`). `all_peaks.csv` adds the
peak annotation (`nearestGene, peakType, distToTSS, GC, …`).

**GroupA rows and GroupB rows are two separate tests, not one test shown twice.** ArchR
tests each group against the rest, which here is the other group, and samples its own
bias-matched background each time. On `jm_test` the two sides' gene Log2FC correlate at only
−0.64, and FDRs for the same gene differ by up to 0.98. Pick one side, normally GroupA
(what Plots shows first and what the decks use), and name it: `one_direction(df, "GroupA")`.
A positive effect then means higher in GroupA.

### Three different significance rules are in play

| Where | Rule |
|---|---|
| Plots Compare tab volcano (and `atx_compare.significant` default) | `FDR <= 0.05` and \|effect\| `>= 0.01` (both editable in the UI) |
| workflow `marker_genes.csv` / `marker_motifs.csv` and the volcano PDFs | raw `p < 0.05` and \|effect\| `>= 0.4` |
| workflow `MA_peaks.pdf` | `FDR <= 0.1` and \|Log2FC\| `>= 0.4` |

On `jm_test`, the Plots rule gives 741 genes higher in GroupA and 1038 higher in GroupB.
The workflow's marker file calls 2330 and 2565. Every count quoted needs its rule beside it.

## Known defects in the workflow output (v0.10.8 and earlier)

Both are fixed in `combined_cluster_wf` v0.10.9 (`wf/compare_clusters.R`). Every run made
before that, which is every run in `compare_outs/` today, still has them. Plots also
launches a pinned version (0.10.6), so a Plots launch keeps producing them until that pin
moves. Keep the guards on for any run you did not see launched at 0.10.9 or later.

1. **`marker_motifs.csv` `Significance` is wrong for direction.** The R code labels
   direction with `mpvalue > 0`, which is always true, so every significant motif is
   called GroupA whatever the sign of `mmean`. On `jm_test`, 72 of the 100 "GroupA" motifs
   are actually higher in GroupB. `volcano_motif.pdf` colours are affected the same way.
   Use `fixed_marker_motifs()`, or better, read `all_motifs.csv`.
2. **`all_peaks.csv` can attach another chromosome's annotation.** The peak annotation is
   merged on `start` and `end` only. The peak ranges carry no chromosome, so two peaks on
   different chromosomes with the same coordinates are cross-joined: each appears twice,
   once with the wrong `nearestGene` and distances. The statistics are unaffected. On
   `jm_test` this affected 32 peaks (64 rows). `flag_ambiguous_peaks()` marks them; don't
   quote `nearestGene` for a flagged row.

The guards are harmless on fixed runs: `fixed_marker_motifs()` relabels nothing and
`flag_ambiguous_peaks()` flags nothing.

## Using it in the skill

```python
import atx_compare as C
C.list_compare_runs()                          # does a run already exist?
C.find_compare_configs(project)                # intended group sizes, by launch time
run = C.load_run(project, "<run>")             # tables + assignment; run.missing lists gaps
C.check_group_mapping(run, n_groupA, n_groupB) # raises if A and B are swapped
s = C.summarize(run, group="GroupA")           # counts + top features, with test and rule
fig, info = C.compare_volcano(run.tables["motif"], "motif")
fig, info = C.compare_rank_plot(run.tables["motif"], "motif")   # = atx_figures.motif_rank_plot
```

To draw the selected cells on tissue, `atx_figures.fetch_compare_annotation()` and
`apply_compare_annotation()` join `UpdateClustName_by_barcode.csv` onto obs.

When reporting, state: the run name, the side tested (GroupA vs GroupB, with what each
group is), the cell counts, the test, the threshold, and `closest` if known. Example:
"Motif activity, GroupA (<what GroupA is>, 828 cells) vs GroupB (<what GroupB is>, 1575
cells), ArchR Wilcoxon on chromVAR deviations, FDR ≤ 0.05 and |MeanDiff| ≥ 0.01: 92 motifs
higher in GroupA" (the `jm_test` numbers).

## Launching a run

Only after the scientist says yes to the exact launch. A launch costs compute in their
workspace and writes a new `compare_outs/<run>/`.

```python
version = C.latest_release()      # newest non-wip build; C.registered_versions() lists all
info = C.write_groupings_file(obs, "condition", "WT", "KO", out_dir=scratch, dataset=project.name)
params = C.launch_parameters(project, "<run>", f"{project.latch_root}/{info['file']}",
                             closest=False)               # False = what Plots sends
C.check_launch(params, version, info)   # version registered, inputs convert, folder free
print(C.describe_launch(params, version, info))           # show this; ask; wait for yes
# --- only after a clear yes, for exactly these params and version ---
C.upload_groupings(project, info)                          # project root, as Plots does
rec = C.launch_compare(params, version, confirmed=True)   # keep rec: id, version, closest
C.wait_for_run(rec["execution_id"])                        # run in the background
```

- **Version.** `C.PLOTS_PINNED_VERSION` (`0.10.6-c5165d-8acb83`) is what Plots launches
  and has the defects above. Prefer `latest_release()` unless the scientist wants to
  match an earlier Plots run. Builds with `-wip-` in the version were registered from
  uncommitted code; don't launch those for a customer. Today 0.10.6 is also the latest
  release (0.10.7 and 0.10.8 exist only as `wip` builds). The fixes reach a launch once
  0.10.9 is registered. `check_launch` raises if a version is not registered.
- **The run name must be new.** `check_launch` refuses an existing `compare_outs/<run>/`.
- **`closest`.** Default False, as Plots sends. `describe_launch` shows it, and the
  launch record keeps it, because the output folder does not.
- **Failure.** A run can succeed with peak or motif steps missing
  (`analysis_errors.log`), so check `load_run(...).missing` before reporting.
