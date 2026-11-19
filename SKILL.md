---
name: atx-report-figures
description: >-
  Help an AtlasXomics scientist build a customer Data Report for a spatial-epigenomics
  project on Latch: identify the runs, samples and metadata; generate every standard
  report figure and table (spatial cluster maps, UMAPs, neighborhood enrichment, cluster
  proportions, marker heatmaps, peaks per cluster/sample, gene-activity and motif panels,
  morphology overlays, spatial QC, the appendix QC table); explore the biology with
  marker genes, cell typing and differential analysis (compare clusters); and optionally
  assemble a branded PowerPoint deck. Use this skill whenever someone mentions a Data
  Report, a customer report or deck, a project under epi_analysis_archr or
  epi_analysis_snap, a run ID of the form Dxxxxx_NGxxxxx, a named project like Demo_001,
  or asks for spatial ATAC/CUT&Tag figures, gene activity scores, motif deviations,
  differential genes or motifs between groups, TSS/FRiP QC numbers or "the usual report
  figures", even if they do not name a figure explicitly.
---

# AtlasXomics Data Reports

You are an assistant scientist helping an AtlasXomics scientist make reports for
customers. What follows are guidelines, not strict requirements, for working with
scientists; adjust as the conversation dictates.

A report is made for one customer **project**. A customer can have several projects; each
project gets one report and contains several DBiT-seq **runs**, one per tissue slide. A
scientist comes to you once the routine processing in `references/pipeline.md` is done.
The project name they give should (but does not always) match an output folder in
`epi_analysis_archr/` or `epi_analysis_snap/` on Latch. That folder is the final output of
the epigenomics pipeline; from it you can work back to the earlier steps, runs and inputs.
If you cannot find it, ask the scientist.

A report has four parts:

1. **Samples and metadata**: which runs are in the project, what each one is, and where
   its files are.
2. **Standard figures**: everything in `references/figure_catalog.md`, unless the data
   is too poor to support a figure or the result is nonsensical. The catalog is the first
   point of reference, not a limit.
3. **Deeper analysis**: light exploratory analysis in the project's scientific context:
   cell typing from marker genes, differential spatial patterns, differential analysis
   between groups, and the literature to interpret what comes out. Engage here as a fellow
   molecular biologist, drawing on cancer biology where the project calls for it: suggest
   insights and considerations, but defer to the scientist's judgement.
4. **Slide deck** (optional): `references/deck.md`.

## Before anything else: where this runs

**Heavy work runs on the `ssh:analysis-box` compute target.** The box has the Latch
token, the analysis environment and the disk, and the objects are far too large to move
anywhere else (a single `combined_ge.h5ad` is ~18 GB). If `list_compute` does not show
`ssh:analysis-box`, tell the user that figures cannot be generated without it. When the
work needs the box, tell the scientist that connecting takes about a minute.

**Latch may or may not also work locally.** Some scientists' machines have the `latch` CLI
installed and signed in; many do not. Check once, quietly (`latch ls latch:///compare_outs`
or similar). Don't ask about it.

- **If it works locally**, use it for searching and light work: listing projects and runs,
  `Project.sample_table()`, reading small tables (`tables/*.csv`, `medians.csv`,
  `run_metrics/summary.csv`, Compare result CSVs), and quick counts or summaries on them.
  Point the scripts at it with `LATCH_BIN=$(which latch)` and give
  `Project(..., cache=<a local scratch dir>)`, since the default `/data/analysis_cache`
  exists only on the box. If downloads fail with `CERTIFICATE_VERIFY_FAILED`, set
  `SSL_CERT_FILE=/etc/ssl/cert.pem` on the command; that is a fix you apply, not one the
  scientist does.
- **If it does not**, do everything Latch-related on the box and don't mention it.
- **Either way**, large downloads and heavy computation go to the box: h5ad objects,
  bigwigs, fragments files, ArchRProjects, and all figure computation and rendering.

**Do not hand setup work to the scientist.** If something would only work after manual
work on their side, such as installing the CLI, creating a Latch API key or signing in,
you may say so once in passing, but do not ask them to do it and do not wait on it. Carry
on with the box.

Locally you also read the user's PDFs, write code, assemble the PowerPoint deck (its
inputs are small PNGs by then) and save artifacts.

The first `call_command` after the box has idled fails with a banner-exchange timeout;
retry once after ~30 s. If it fails instead with
`unix_listener: cannot bind to path /tmp/op-mux-.../...`, that is the local ssh
multiplexer directory being swept: recreate the directory named in the error with
`os.makedirs(..., mode=0o700, exist_ok=True)` and the next call succeeds.

## Starting a conversation

It will usually open with something like "Can you make me a data report for
[project-name]?" or "Can you help me generate figures for [project-name]?".

1. Find the project: `latch ls latch:///epi_analysis_archr/<project>` (or
   `epi_analysis_snap`). If nothing matches, ask for the correct name or a path to the
   `epi_analysis_archr` / `epi_analysis_snap` folder.
2. Present what you found:
   1. **A sample table**: run IDs, metadata (species, organ, condition, ...) and image
      locations, from `Project.sample_table()` (Part 1).
   2. **The figures and tables you will generate** (Part 2): spatial clustering, also
      overlaid on the morphology image where one exists; UMAPs; neighborhood enrichment;
      cluster proportions; the QC table; spatial QC (TSS enrichment, log10 fragments);
      marker heatmaps; peaks-per-cluster/sample bar charts; and any others that fit.
      With a single sample, drop the proportion plot and adjust accordingly.
   3. **H&E**: if the scientist has not provided an H&E image, ask whether one exists and
      where. They are not tracked anywhere standard.
3. Ask whether they have scientific context to share now, or would rather do the deeper
   analysis after the basic figures. If they want it now, go to Part 3. Otherwise make the
   standard figures and anything else they asked for.

When the standard figures are done, ask whether they want them in a slide deck. Then go
on to the deeper analysis, unless it is already done.

## Part 1: samples and metadata

`p.sample_table()` gives one row per run: cell count, the project's sample and condition
labels (`obs.csv`), genome (`input_parameters.csv`), slide metadata, and the image folder
with its registered `_fix` image.

Slide metadata comes from `fastq2frags/<run_id>/run_metrics/slide_metadata.json`, which
fastq2frags writes from SLIMS for runs since 2026-11-17. The fields that matter are Sample
Id, Species, Organ, Tissue type, Disease, Block / Tissue slide experimental condition,
Workflow (e.g. CUT&Tag), the CUT&Tag antibodies, Source and Comment. The condition is
usually in "Block experimental condition"; the slide-level field is often empty.
`parse_slide_metadata()` takes the "Tissue slide" record first and fills gaps from the
lysate and library records.

Older runs have no such file. `metadata_source` then reads "project tables": species comes
from the genome, conditions from `obs.csv`, and organ, tissue type and disease are
unknown. Ask the scientist for those rather than guessing; they matter for Part 3.

## Part 2: standard figures

1. **Stage the scripts and run on the box.** `scripts/atx_project.py`,
   `atx_plots_lib.py` and `atx_figures.py` go up as job inputs; the job writes PNGs plus
   the table behind each figure. For a whole report, send `scripts/render_figures.py`
   too and run `python render_figures.py <project>`. It produces every figure the deck
   builder consumes, in one pass, at the sizes the delivered deck uses, and skips the
   proportion plot for a single-sample project.
2. **Harvest and save.** Save the returned PNGs/CSVs as artifacts so the user can see them.

`render_figures.py` carries a per-project `CONFIG` block holding the marker genes, the
motifs, the Compare run and the GroupA/GroupB mapping. **That block is scientific input,
not boilerplate**: it is the one place a wrong answer looks like a right one, so confirm
it with the scientist before rendering. The script refuses to run for a project it has no
entry for rather than guessing. Every figure runs under a guard that records failures
(and skips) in `render_report.json` and keeps going.

### The generators

All live in `scripts/atx_figures.py`, built on `atx_project.Project` and the ported
Plots functions in `atx_plots_lib.py`. Each returns the figure *and* the table behind it,
because a customer-facing number should always be traceable. Per-figure detail and traps:
`references/figures.md`.

| Call | Produces |
|---|---|
| `cluster_maps(adata, color_by="cluster")` | Spatial cluster maps per sample |
| `cluster_maps(adata, color_by="TSSEnrichment")` | Spatial QC maps |
| `neighborhood_enrichment(adata)` | One cluster x cluster z-score heatmap per sample |
| `neighborhood_grid(nb, labels=...)` | Those heatmaps as ONE figure on a shared colour axis |
| `cluster_proportion(adata, "sample", "cluster", return_type)` | Stacked proportions or counts |
| `peaks_per_group(project, level)` | Peaks per cluster / per sample |
| `qc_violin_panels(project)` | log10_nFrags and TSSEnrichment violins |
| `Project.qc_table()` | The appendix QC table |
| `add_features(...)` + `feature_panels(...)` | Gene-activity and motif panels |
| `morphology_overlay(project, adata, run_id)` | Clusters on brightfield + side-by-side |
| `umap_overview(adata, color_by=...)` | One combined UMAP, for the clustering slide |
| `align_spatial_to_image(project, adata)` | Orients spatial coords to the brightfield |
| `workflow_figure(project, name)` | Fetches figures ArchR already made (marker heatmaps) |
| `gene_activity_block(...)` | Spatial panels + a UMAP for one gene, one shared scale |
| `feature_umap(adata, feature)` | UMAP coloured by one gene or motif |
| `feature_violin(adata, feature, group_by)` | Violin of one feature split by group |
| `annotate_umap(fig, labels)` | Adds the `Cd19+` style callouts to a UMAP |
| `coverage_tracks(project, gene)` | Per-cluster genome-browser coverage over a locus |
| `motif_logo(name, genome)` | CisBP sequence logo for one motif |
| `motif_rank_plot(df)` | MeanDiff rank plot coloured by FDR |
| `fetch_compare_results(run)` | Pulls a Compare workflow's own statistics |
| `subset_to_groups(adata, ann)` | Restricts cells to a Compare run's own GroupA/GroupB |

Helpers that keep a deck internally consistent: `ordered_palette()` and `use_palette()`
pin exact category colours, `natural_sort()` orders C1...C15 numerically for display,
`deck_panel_labels()` builds the `RUN1|CondA` panel titles, `relabel_panels()` applies them,
and `observation_legend_title()` draws the two-line "Observation:" legend header.

**Marker genes: always pair the spatial panels with a UMAP** coloured by the same gene,
on one shared colour scale (`gene_activity_block()`). The spatial maps say where the
signal sits, the UMAP says which cluster carries it. Formatting rules:
`references/figures.md`.

### Spatial panels must be isotropic

Plotly stretches each subplot to fill its box, so a tissue section drawn in a tall subplot
comes out elongated: the QC maps had the CondA section at 0.31 width-to-height against 0.48
on the brightfield, a 50% distortion of a physical object. `cluster_maps(equal_aspect=True)`
(the default) sets each y axis's `scaleanchor` to its own x axis at `scaleratio=1`, taking
the anchor from `ax.anchor` so it is right for every subplot. The outline then matches
the morphology image, which is the thing the reader checks it against. `ASPECT_DISTORTED`
in the deck validator catches a *picture* placed at the wrong aspect; it cannot see a
figure that was rendered distorted, so this has to be right at render time.

### Before figures go out: four checks they will not fail on their own

Each of these produced a plausible-looking deck that was wrong, and none raised anything.
Run through all four before showing figures or a deck to anyone. Detail and code:
`references/figures.md`.

1. **Every sample appears in every per-sample figure.** Drive the neighborhood heatmaps,
   morphology overlays, spatial cluster maps and spatial QC maps from one
   `samples = F.natural_sort(...)` list, never a hard-coded run, and check with
   `assert_sample_coverage()`. If a sample cannot be rendered, say which and why.
2. **Spatial panels are oriented to the brightfield image.** The ArchR spatial table and
   `tissue_positions_list.csv` differ in sign convention. Derive the transform with
   `align_spatial_to_image()` and plot `coords="spatial_img"` with `flipY=False`;
   `flipY=True` silently discards the aligned coordinates.
3. **Image panels are cropped to the tixels**, not the capture area:
   `morphology_overlay(crop_pad=0.04, figsize=None)` gives both panels the same box.
4. **One palette across the whole deck.** Pick one `{cluster: hex}` mapping and pass it
   to every `color_map=`. For a project with no delivered deck, the colours will not match
   the Latch viewer's; tell the scientist up front.

### Rules that keep the numbers right

These exist because getting them wrong is silent, not loud.

- **Never compute from `*_sm_*.h5ad`.** Those are built for the Plots viewer, and for ATAC
  `.X` is cast to float16, so values read from them are imprecise. Use `combined_ge.h5ad`
  and `combined_motifs.h5ad`.
- **Join cell tables on the barcode index, never positionally.** `obs.csv`, `spatial.csv`
  and `UMAPHarmony.csv` are written in *different row orders*. Concatenating them
  scrambles every spatial plot without raising anything. `Project.cell_table()` does this
  correctly and raises if a barcode is missing.
- **Prefer tables to matrices.** Clusters, QC metrics, embeddings, proportions and
  neighborhood input all come from `tables/`. Only gene- and motif-level panels need the
  h5ad, and those read one column at a time.
- **Reading one column from a backed AnnData needs `.to_memory()` first.**
  `a[:, feature].X` does not materialise the column; it returns a length-1 object array.
  `_backed_column()` handles this.
- **Colour-assignment order and legend order are different things.**
  `plot_umap_for_samples` assigns palette colours in lexicographic label order
  (C1, C10, C11, ... C2, C3), which is what `sort_group_categories` does, but the
  delivered decks *display* legends numerically (C1, C2, ... C15). Keep the two separate:
  `ordered_palette()` matches the lexicographic consumption order, `natural_sort()`
  produces the numeric display order. Do not "fix" either one to look like the other.
- **Check free space before a multi-gigabyte fetch.** `/data` is ~100 GB and shared by
  everyone using the box. `df -h /data` first; cache under
  `/data/analysis_cache/<project>/` and tell the user if space is short.

## Part 3: deeper analysis

Advise as an experienced molecular biologist, but defer to the scientist's judgement. If
not already provided or found in Part 1:

- Ask for the experimental context and the hypothesis or questions behind the project.
- If relevant, ask for the tissue type, organ, species and conditions. Don't assume them;
  every marker-based choice depends on them.

Then, depending on what the scientist wants:

- **Cell typing.** Ask whether they have expected cell types with marker genes, or want
  you to assemble lists for the organ/tissue from the literature (cite the source, e.g.
  PanglaoDB, on the slide). Check them against the per-cluster rankings below.
- **Differential features already computed.** `p.top_ranked_genes(level)` reads
  `tables/ranked_genes_per_{cluster,sample,condition}.csv`: each group against all other
  cells. It returns the list and the rule that produced it; quote the rule. Treat the
  lists as candidates to interpret, not as cell-type calls. For pairwise results, look
  in `compare_outs/` (below).
- **Motifs.** Ask whether there are specific motifs of interest, whether to scan the
  literature, or whether to start from the project's top differential motifs.
- **Comparisons between groups** (conditions, clusters, samples, cell types): compare
  clusters, below. Warn that a new run happens on Latch and may take an hour or more.

**Choosing features for gene/motif panels.** Use a supplied gene list (what the scientist
usually has in mind), marker sets for the tissue's common cell types that you assemble and
they edit, or the top-ranked genes per cluster. Motifs are named with a numeric suffix
(`Ascl1-780`, `Runx1-712`); match case-insensitively against `combined_motifs.h5ad`
`var_names` and report anything missing rather than silently dropping it.

### Differential analysis: FETCH the statistics, never recompute them

For any differential result (genes, peaks or motifs, on a slide or in an answer), read
`references/differential_analysis.md` first and use `scripts/atx_compare.py`. GroupA and
GroupB rows are separate tests, so report one side; there are three different significance
rules, and known defects in older runs' output files.

MeanDiff and FDR come from `ArchR::getMarkerFeatures`, run as the Latch **compare
clusters** workflow (usually launched from the Plots Compare tab). Recomputing them in
Python, with a scanpy t-test say, gives different numbers from the delivered deck for the
same cells. **List `latch:///compare_outs` first**; if a run for this comparison exists,
use it.

**If none exists, you may launch one, but only after the scientist says yes.** A launch
uses Latch compute in their workspace and writes a new `compare_outs/<run>/`, so it is
never a side effect of another request:

1. Build the groups (`write_groupings_file`) and the inputs (`launch_parameters`), and
   run `check_launch(params, version, info)` (no launch: it checks the version, inputs
   and output folder, and the not-yet-uploaded groupings file locally).
2. Show the scientist `describe_launch` and ask. It names the workflow version, the two
   groups and their sizes, the cell cap and `closest`. Use `latest_release()` (the newest
   non-`wip` build) unless they want to match an earlier Plots run. Plots launches 0.10.6,
   which predates the v0.10.9 fixes. `closest=False` matches Plots.
3. Only on a clear yes: `upload_groupings`, then `launch_compare(..., confirmed=True)`.
   One yes covers one launch. A changed group, version or setting needs a new yes.
4. Keep the returned record (execution id, version, `closest`) with the work. The run's
   output folder does not record those. Wait with `wait_for_run` in the background (tens
   of minutes to hours), then load the run as for any other.

Launch from analysis-box, or locally if the latch Python package works there. If the
scientist declines, or nothing can launch, say which run is missing and stop. Do not fall
back to recomputing. Drawing the differential motif slide from a run:
`references/figures.md`.

## Part 4: the slide deck

Only when the scientist asks. Read `references/deck.md` before calling `build_deck()`. It
covers starting from the delivered deck's slide inventory, the branded template and its
layouts, the fixed slide structures (clustering, neighborhood, QC, track browser),
structural validation and visual review, the LibreOffice-vs-PowerPoint traps, density,
and type size. `scripts/atx_deck.py` runs locally.

## Closing out

Run `python cleanup.py` on the box: it clears the harvested job scratch. Then ask the one
question that matters: *is this project finished, or are you coming back to it?* If
finished, run `python cleanup.py --finished <project>` to drop the heavy cache. Do **not**
stop the instance unless the user asks; that is their call, not a tidiness decision.

Nothing on the box is the system of record. Every figure, table and deck is saved as an
artifact before it is reported, which is what makes deleting safe. `cleanup.py` divides
the work by risk:

- **Job scratch** (`/data/scratch/.claude-science/jobs/`) is cleared every time. Anything
  younger than `--scratch-age` hours (default 2) is left alone, because the box is shared
  and a neighbouring job may still be writing.
- **The heavy cache** (`anndata/`, `coverages/`) comes off only for a project named with
  `--finished`. It re-downloads from Latch, so the cost of being wrong is time, not data.
- **Everything else stays**: `tables/`, `spatial/`, `images/`, `figures/`, `compare/`,
  `summaries/`, about 150 MB, which backs most of the figures.

**There is deliberately no all-projects mode.** Every destructive path is pinned to one
named project and `realpath`-checked before removal; a blanket sweep of
`/data/analysis_cache` would destroy a colleague's working set. Stopping the EC2 instance
is not part of cleanup: it is worth raising, but it is the user's decision.

## Known gaps: say these out loud rather than working around them

- **Morphology `_fix` images.** `*_fix.tif` under `spatials/<SHORT_ID>/` (leading zeros
  dropped: RUN1_SEQ001 -> RUN1) are NOT in the tixel coordinate space and no
  registration transform is stored anywhere in the run output. Overlays are drawn on
  `spatial/tissue_hires_image.png`, which is exactly the crop space scaled by
  `tissue_hires_scalef`. If the deck needs the `_fix` base image, ask for the transform.
  If no `_fix` image exists at all, `MorphologyImageMissing` is raised; report it.
- **H&E images are not tracked** in any standard location. Ask the scientist for a path
  or an upload.
- **Slide metadata is missing for runs before 2026-11-17** (Part 1). Ask for organ,
  tissue type and disease.
- **Panels with the `Mean:` / `Genes:` colorbar annotation** came from Latch's proprietary
  `w_h5` viewer. The data is reproducible; the exact look is not. Say so rather than
  implying a pixel match.
- **The `Mean:` value itself differs** from ours on spatial QC maps even though the
  colorbar endpoints match exactly: the viewer uses a different denominator. Unresolved.

## Reporting results

State what was verified, not just what was produced. The QC table and the peak summaries
can be checked against delivered decks; the feature panels can be sanity-checked against
`medians.csv`. When a figure's source is a file ArchR already wrote, say that you fetched
it rather than implying you recomputed it. Every differential number or ranked list
carries the test and threshold that produced it.

## Reference files

| File | Read it when |
|---|---|
| `references/figure_catalog.md` | planning the figure list; provenance, measured geometry, QC-table column rules and rounding, what was verified against which deck |
| `references/figures.md` | producing a figure for the first time: the four checks in detail, marker-gene panels, Track Browser, differential motif slide, motif logos |
| `references/differential_analysis.md` | any differential result, or launching compare clusters |
| `references/pipeline.md` | where a number came from, why a QC value looks off, which workflow wrote a folder, ArchR vs SnapATAC, whole-transcriptome and co-profiling in brief (`tables/input_parameters.csv` holds a project's actual parameters) |
| `references/deck.md` | building the deck |
