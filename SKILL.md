---
name: atx-report-figures
description: >-
  Generate any figure or table from an AtlasXomics spatial-epigenomics Data Report
  directly from a Latch project — spatial cluster maps, neighborhood enrichment, cluster
  proportions, peaks per cluster/sample, gene-activity and motif panels, morphology
  overlays, QC violins and the appendix QC table — and optionally assemble them into a
  PowerPoint deck. Use this skill whenever someone mentions a Data Report, a customer
  report or deck, an ArchR project under epi_analysis_archr or epi_analysis_snap, a run
  ID of the form Dxxxxx_NGxxxxx, a named project like Demo_001 or Demo_002, or asks
  for spatial ATAC/CUT&Tag figures, gene activity scores, motif deviations, TSS/FRiP QC
  numbers or "the usual report figures" — even if they do not name a figure explicitly.
---

# AtlasXomics report figures

Regenerate the figures and tables that appear in AtlasXomics Data Reports, from the
project outputs on Latch.

## Before anything else: where this runs

**Every Latch operation must run on the `ssh:analysis-box` compute target.** The Latch
token lives only on that box, latch.bio is not reachable from the local sandbox, and the
objects are far too large to move (a single `combined_ge.h5ad` is ~18 GB). If
`list_compute` does not show `ssh:analysis-box`, stop and tell the user — nothing here
works without it.

The split that matters:

- **On the box**: every `latch ls` / `latch cp`, all table and h5ad reads, all figure
  computation and PNG rendering.
- **Locally**: reading the user's PDFs, writing code, assembling the PowerPoint deck
  (inputs are small harvested PNGs by then), saving artifacts.

The first `call_command` after the box has idled fails with a banner-exchange timeout;
retry once after ~30 s. If it fails instead with
`unix_listener: cannot bind to path /tmp/op-mux-.../...`, that is the local ssh
multiplexer directory being swept — recreate the directory named in the error with
`os.makedirs(..., mode=0o700, exist_ok=True)` and the next call succeeds.

## Workflow

1. **Establish the project.** Ask for the project name if not given, then confirm it
   exists: `latch ls latch:///epi_analysis_archr/<project>`. Projects are also under
   `epi_analysis_snap` for SnapATAC runs.
2. **Ask what to plot.** If the user named figures, use those. Otherwise show the catalog
   (`references/figure_catalog.md`) and let them choose.
3. **For gene-activity or motif panels, ask two questions up front** — see
   "Choosing features" below. Do not guess genes.
4. **Stage the scripts and run on the box.** `scripts/atx_project.py`,
   `atx_plots_lib.py` and `atx_figures.py` go up as job inputs; the job writes PNGs plus
   the table behind each figure. For a whole report, send `scripts/render_figures.py`
   too and run `python render_figures.py <project>` — it produces every figure the deck
   builder consumes, in one pass, with the sizes the delivered deck uses.
5. **Harvest and save.** Save the returned PNGs/CSVs as artifacts so the user can see them.
6. **Assemble a deck only if asked.** Read `references/deck.md` first;
   `scripts/atx_deck.py` runs locally.
7. **Close out the session.** Run `python cleanup.py` on the box — it clears the harvested
   job scratch unconditionally. Then ask the one question that matters: *is this project
   finished, or are you coming back to it?* If finished, run
   `python cleanup.py --finished <project>` to drop the heavy cache. Do **not** stop the
   instance unless the user asks; that is their call, not a tidiness decision.

`render_figures.py` carries a per-project `CONFIG` block holding the marker genes, the
motifs, the Compare run and the GroupA/GroupB mapping. **That block is scientific input,
not boilerplate** — it is the one place a wrong answer looks like a right one, so confirm
it with whoever ordered the report before rendering. The script refuses to run for a
project it has no entry for rather than guessing. Every figure is produced under a guard
that records failures in `render_report.json` and keeps going: a 25-figure render that
aborts on figure 3 wastes the twenty minutes already spent.

## Cleanup: the artifact store is the record, the box is a cache

Nothing on the box is the system of record. Every figure, table and deck is saved and
versioned as an artifact before it is reported, which is what makes deleting safe — and
is the reason not to hoard a cache "just in case."

`scripts/cleanup.py` divides the work by risk, not by a level the user has to choose:

- **Job scratch** (`/data/scratch/.claude-science/jobs/`) is cleared every time, no
  question. Those directories are copies of outputs already harvested. Anything younger
  than `--scratch-age` hours (default 2) is left alone, because the box is shared and a
  neighbouring job may still be writing.
- **The heavy cache** (`anndata/`, `coverages/`) comes off only for a project named with
  `--finished`. For Demo_001 that is 19.6 GB against 95 MB kept — and both re-download
  from Latch, so the cost of being wrong is time, not data.
- **Everything else stays**: `tables/`, `spatial/`, `images/`, `figures/`, `compare/`,
  `summaries/`. About 150 MB, and it backs most of the figures.

The trigger is a question to the user — *finished with this project, or coming back to
it?* — not a disk threshold. Disk pressure is not what makes a cache stale.

**There is deliberately no all-projects mode.** Every destructive path is pinned to one
named project and `realpath`-checked before removal. A blanket sweep of
`/data/analysis_cache` would destroy a colleague's working set; at the time of writing
`Demo_002` sits beside `Demo_001` there.

Stopping the EC2 instance is not part of cleanup. It dominates the cost and restarting is
free, so it is worth raising — but it is the user's decision, and `cleanup.py` only
reminds them it is still running.

## Building the deck

Deck assembly is in `references/deck.md`. Read it before calling `build_deck()`. It
covers starting from the delivered deck's slide inventory, the branded template and its
layouts, the fixed slide structures (clustering, neighborhood, QC, track browser),
structural validation and visual review, the LibreOffice-vs-PowerPoint traps, density,
and type size on slides. The figure rules below apply whether or not a deck is built.

## Spatial panels must be isotropic

Plotly stretches each subplot to fill its box, so a tissue section drawn in a tall subplot
comes out elongated — the QC maps had the CondA section at 0.31 width-to-height against 0.48
on the brightfield, a 50% distortion of a physical object. `cluster_maps(equal_aspect=True)`
(the default) sets each y axis's `scaleanchor` to its own x axis at `scaleratio=1`, taking
the anchor from `ax.anchor` so it is right for every subplot. The outline then matches
the morphology image, which is the thing the reader checks it against. `ASPECT_DISTORTED`
in the deck validator catches a *picture* placed at the wrong aspect; it cannot see a
figure that was rendered distorted, so this has to be right at render time.

## The generators

All live in `scripts/atx_figures.py`, built on `atx_project.Project` and the ported
Plots functions in `atx_plots_lib.py`. Each returns the figure *and* the table behind it,
because a customer-facing number should always be traceable.

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
| `workflow_figure(project, name)` | Fetches figures ArchR already made |
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

## Differential motifs: FETCH the statistics, never recompute them

This is the one slide whose NUMBERS cannot be produced here. MeanDiff and FDR come from
`ArchR::getMarkerFeatures`, run as a Latch **Compare Workflow** launched from the Plots
app (`plots_epigenome/compare/`). Recomputing them in Python -- a scanpy t-test, say --
gives different numbers from the delivered deck for the same cells. Fetch the workflow
output instead:

```python
mot = F.fetch_compare_results(project, "Demo_CD19_high_CondA_vs_CondB", "motif")
ann = F.fetch_compare_annotation(project, "Demo_CD19_high_CondA_vs_CondB")
F.apply_compare_annotation(cells, ann, "CD19_high",
                           mapping={"GroupA": "CondA", "GroupB": "CondB"})
```

Runs live at `latch:///compare_outs/<run>/` with `motif_results/all_motifs.csv`,
`gene_results/`, `peak_results/`, per-group `coverages/`, and
`UpdateClustName_by_barcode.csv`. **List `latch:///compare_outs` first** -- if a run for
this comparison already exists, use it; if none does, the scientist must launch the
Compare workflow from the Plots app, which this skill cannot do.

- **Use the run's own cell assignment**, not a re-derived threshold. The delivered
  subtitle says "Cd19 gene activity score higher than 0.6", but re-thresholding can select
  a slightly different set than the statistics were computed on. `UpdateClustName_by_barcode.csv`
  is the exact selection; its barcodes are already in the `<run_id>#<barcode>-1` obs form.
- **Check which group is which.** GroupA/GroupB are positional. Confirm against the run
  name and the cell counts, which are embedded in the `compare_config_*` filename in the
  project root (`..._A<nA>_B<nB>_...` gives each group's cell count). Cells outside the
  selection are left unassigned, so the whole tissue can still be drawn as background.
- **The rank plot settings are fixed by the Plots code**, not chosen:
  `motif_rank_plot()` mirrors `compare/motif_results.py` -- rank by MeanDiff, colour by
  FDR, `PuBu_r`, 4 labels, marker size 6.

### Motif sequence logos

`motif_logo()` draws a CisBP logo from position-probability matrices that live in
`plots_epigenome/motif_logo/seqlogo_<genome>.json`. At 0.5-1.4 MB these are too large to
ship as skill assets, so `load_seqlogos()` fetches the file once from the pinned public
commit (`PLOTS_EPIGENOME_COMMIT`) and caches it next to the scripts; pass `assets_dir=` to
use a copy you already have. Two traps:

- **There is no mm39 file.** The repo ships hg38, mm10 and rn6 only. CisBP motif identity
  does not depend on the assembly, so an mm39 project uses the mm10 set; `load_seqlogos`
  makes that substitution explicitly and warns rather than failing.
- **Naming differs between sources.** The h5ad `var_names` use a hyphen before the index
  (`Ascl1-780`); the seqlogo JSON and the Compare results use an underscore
  (`Ascl1_780`). All 884 mouse entries use the underscore form. `motif_key()` converts,
  and raises listing same-TF alternatives when a motif is absent.

Validated end to end against a delivered slide: the fetched MeanDiff and FDR values match
it, and the logos render the expected consensus sequences (CAGCTG for Ascl1, CCACA for
Runx1 -- these are the published CisBP motifs, not project-specific results).

## Marker genes: always pair the spatial panels with a UMAP

**When building a deck from scratch, show each marker gene as its spatial panels AND a
UMAP coloured by that same gene.** This pairing is the recommended default, not a styling
choice: the spatial maps say where in the tissue the signal sits, the UMAP says which
cluster carries it, and neither answers the other's question. Spatial panels alone are an
incomplete slide.

`gene_activity_block()` builds one gene's pair as a single figure -- spatial panels per
sample plus the UMAP, all on one shared colour scale, with the gene label written down
the left edge:

```python
fig, info = F.gene_activity_block(cells, "Cd19", samples=samples,
                                  sample_labels=cond,          # CondA / CondB / CondC
                                  row_label="Cd19: B cell")
```

Formatting rules this slide follows, all taken from the delivered deck:

- **Title panels by CONDITION, not run id.** This slide uses `CondA`, `CondB`, `CondC`,
  unlike the clustering slide's `RUN1|CondA`. Pass `sample_labels`.
- **Label each row with gene and cell type** -- "Cd3e: T cell", "Cd19: B cell",
  "Hmox1: red pulp macrophage" -- down the left edge via `row_label`.
- **One shared colour scale per gene**, spanning its spatial panels and its UMAP, with the
  colourbar between them. Scale per gene, never across genes.
- **Subtitle cites the marker source**, e.g. "Gene markers from
  https://panglaodb.se/markers.html".
- **Lay the blocks out with `grid`, not `columns`.** `_place_columns` scales a whole
  column as one stack, so a square UMAP sharing a column with a wide gene block shrinks
  the block to match; `{"grid": [[a, b], [c, d]]}` gives every cell the same size, which
  is what the delivered slide does.
- **Callouts come from the data.** `top_clusters_for_features()` returns the highest-mean
  cluster per gene and its UMAP centroid, and `annotate_umap()` writes "Cd19+" there. On
  Demo_001 this reproduces the delivered callouts exactly (Cd3e->C7, Cd19->C14,
  Hmox1->C10). Keep the label halo: dark text at the centroid of a dark cluster is
  invisible, and the callout looks like it was never drawn.

The per-gene means this produces match the delivered panels exactly, so the `Mean:`
annotation on gene-activity panels is reproducible even though the one on the spatial QC
maps is not.

## Genome browser tracks (the Track Browser slide)

The "Gene Activity Scores identify distinct biological clusters" slide with per-cluster
coverage over a gene locus is the **Track Browser** tab, and it is easy to miss: the
ported `atx_plots_lib` contains none of it, because it comes from
`plots_epigenome/track_browser/`, not `welcome/init.py`. That module hands bigwigs to
Latch's proprietary `w_igv` widget (an igv.js wrapper), which cannot run headlessly --
but its inputs are ordinary ArchR bigwigs, so the figure is reproducible with pyBigWig
and matplotlib. It is not a pixel match to an IGV screenshot; the data is identical.

```python
bws = F.fetch_coverages(project, "cluster")        # also "sample", "condition"
fig, info = F.coverage_tracks(bws, ("chr7", 126002623, 126019067),
                              color_map=CM, title="Cd19",
                              gene_model=model,    # optional RefSeq row
                              highlight=("C14", "C15", 126009000, 126014500))
```

Things that matter here:

- **Check the genome build; do not assume it.** Demo_001 is **mm39**, not mm10. Read it
  from the bigwig header rather than the pipeline docs: mm39 `chr1` is 195,154,279,
  mm10 is 195,471,971, hg38 is 248,956,422. The Track Browser offers hg38/mm10/mm39/rn6
  and the wrong choice silently plots the wrong locus.
- **The bigwigs must be local.** pyBigWig cannot read `latch://` URLs, so there is no
  range-request shortcut: `fetch_coverages` downloads the set once (~100 MB each, 1.37 GB
  for 15 Demo_001 clusters) and caches it under the project cache.
- **All tracks share one y scale.** That shared maximum is what makes the panel a
  comparison rather than fifteen unrelated plots, and the delivered slides print it at the
  top left of every track, so it is directly checkable against a delivered slide.
- **Leave a gap between tracks.** With touching axes the `0` of one track prints on top of
  the shared maximum of the one below it.
- **Read axis positions after `subplots_adjust`, not before**, or a highlight box lands on
  whichever tracks happened to sit at those coordinates in the pre-layout figure.

Gene loci and the RefSeq Curated gene model come from the UCSC API
(`api.genome.ucsc.edu`) via `refseq_gene(gene, genome)`, which runs
`/search?search=<gene>&genome=<build>` for the locus and then
`/getData/track?track=ncbiRefSeqCurated` for the exons, taking the longest matching
transcript so a multi-isoform gene gets its full span.

**The view window is a human choice and is not recorded anywhere.** `refseq_gene` pads
the transcript by 5 kb. In validation that reproduced one delivered panel's shared maximum
exactly and missed the other two, because those IGV views had been scrolled to different
windows. A shared maximum only matches when the window does. Report the computed maximum
and the locus you used rather than tuning the window until the number agrees -- and if a
panel must match a delivered one exactly, ask the scientist for the locus they viewed.

(Placing this figure on its slide: see `references/deck.md`.)

## Before figures go out: four checks they will not fail on their own

Each of these produced a plausible-looking deck that was wrong, and none of them raised
anything. Run through all four before showing figures or a deck to anyone.

### 1. Every sample appears in every per-sample figure

`neighborhood_enrichment()` returns a dict keyed by sample and `morphology_overlay()`
takes one run at a time — so it is easy to render one sample, put it on the slide, and
ship a deck that silently covers a third of the experiment. The delivered reports show
**one panel per sample** for the neighborhood heatmaps, the morphology overlays, the
spatial cluster maps and the spatial QC maps.

Enumerate samples once and drive every per-sample figure from that list:

```python
samples = F.natural_sort(cells.obs["sample"].unique().tolist())
for s in samples:                      # never a single hard-coded run id
    fig, info = F.morphology_overlay(p, cells, s, color_map=CM)
```

Then assert the count before building the deck: the number of neighborhood PNGs, of
morphology PNGs and of spatial panels must each equal `len(samples)`. If a sample cannot
be rendered — a missing `_fix` image, no tissue positions — say which sample and why
rather than quietly producing a shorter deck.

### 2. Spatial panels must be oriented to the brightfield image

The ArchR spatial table and the AtlasXBrowser `tissue_positions_list.csv` do not share a
sign convention, so the default spatial panels come out **mirrored** against the
registered image. The scientist checks alignment by eye against the morphology slide, so
the image is the ground truth and the spatial panels move to match it — never the reverse.

Do not assume a transform; derive it. Both tables are barcode-indexed, so correlating
them over shared barcodes recovers the axis mapping and signs exactly:

```python
report = F.align_spatial_to_image(p, cells, samples)   # writes obsm["spatial_img"]
fig, _ = F.cluster_maps(cells, coords="spatial_img", flipY=False, ...)
```

For all three Demo_001 runs the relation is `x = +Spatial_1`, `y = -Spatial_2` at
|r| = 1.00, but check the returned report rather than hard-coding that — it names the
mapping and signs used per sample so a surprising one is visible.

**The trap:** `plot_umap_for_samples` IGNORES its `coords` argument whenever
`flipY=True`. It rebuilds `spatial_flippedY` from `obsm['spatial']` and overwrites
`coords`, so an aligned frame passed with `flipY=True` is discarded and the old
orientation is plotted — producing a figure that looks fine and is unchanged.
`cluster_maps` now raises on that combination; pass `flipY=False` with coordinates that
are already screen-oriented.

### 3. Crop image panels to the tixels, not to the capture area

The brightfield capture is several times the size of the tissue, so a panel drawn at the
full image extent is mostly empty slide — the tissue sat in the middle of a grey field at
roughly a third of the width it has in the delivered deck. The bench scientists crop by
hand before pasting; `morphology_overlay()` now does it with `crop_pad=0.04`.

The crop box is the bounding box of the **plotted tixels** (in hires pixel space), padded
by 4% of its larger side, and **both panels get the same box** so the overlay and the
brightfield stay registered and comparable. Cropping each panel to its own content would
silently put them at different scales. The box is reported as `info["crop_box"]`; pass
`crop_pad=None` to keep the full frame.

Sizing the figure matters as much as the crop: with `set_box_aspect()` holding the panel
to the crop's aspect, a fixed `figsize=(14, 7)` leaves the two panels floating apart in a
wide white figure, and `bbox_inches="tight"` cannot remove space that is *between* axes.
`figsize=None` (the default) computes the width from the crop aspect — two panels plus a
1.3in legend column — and the panels end up adjacent.

The same reasoning applies to any panel drawn on a microscopy image. Before placing one,
ask what fraction of it is tissue; if it is well under half, crop.

### 4. One palette across the whole deck

A cluster must be the same colour on the clustering slide, the morphology overlays and
the proportion bar. The reports use `DEFAULT_H5_CATEGORICAL_PALETTE` from the Plots code,
but the Latch H5 viewer assigns it in its own internal category order, which **is not
stored anywhere in the project** — so the colours cannot be re-derived from the data.

Pick one `{cluster: hex}` mapping at the top of the run and pass it to every generator
that takes `color_map=` (`cluster_maps`, `umap_overview`, `morphology_overlay`,
`cluster_proportion`). `DEMO_CLUSTER_COLORS` is the mapping for that project,
recovered by sampling the delivered slide's legend and snapping to the palette.

For a project with no delivered deck to sample, the figures will be palette-consistent
but assigned differently from whatever the viewer would have chosen. Tell the scientist
that up front — it is not something to discover during a customer review.

## Rules that keep the numbers right

These exist because getting them wrong is silent, not loud.

- **Never compute from `*_sm_*.h5ad`.** Those are built for the Plots viewer, and for ATAC
  `.X` is cast to float16 — values read from them are imprecise. Use `combined_ge.h5ad`
  and `combined_motifs.h5ad`.
- **Join cell tables on the barcode index, never positionally.** `obs.csv`, `spatial.csv`
  and `UMAPHarmony.csv` are written in *different row orders*. Concatenating them
  scrambles every spatial plot without raising anything. `Project.cell_table()` does this
  correctly and raises if a barcode is missing.
- **Prefer tables to matrices.** Clusters, QC metrics, embeddings, proportions and
  neighborhood input all come from `tables/`. Only gene- and motif-level panels need the
  h5ad, and those read one column at a time.
- **Reading one column from a backed AnnData needs `.to_memory()` first.**
  `a[:, feature].X` does not materialise the column — it returns a length-1 object array.
  `_backed_column()` handles this.
- **Colour-assignment order and legend order are different things.**
  `plot_umap_for_samples` assigns palette colours in lexicographic label order
  (C1, C10, C11, ... C2, C3), which is what `sort_group_categories` does — but the
  delivered decks *display* legends numerically (C1, C2, ... C15). Keep the two separate:
  `ordered_palette()` matches the lexicographic consumption order, `natural_sort()`
  produces the numeric display order. Do not "fix" either one to look like the other.
- **Check free space before a multi-gigabyte fetch.** `/data` is ~100 GB and shared by
  everyone using the box. `df -h /data` first; cache under
  `/data/analysis_cache/<project>/` and tell the user if space is short.

## Choosing features for gene/motif panels

Ask these before plotting, in this order:

1. **What tissue and species is this?** Needed for any marker-based choice; don't assume.
2. **Which selection route?**
   - a supplied gene list (what the scientist usually has in mind),
   - marker sets for the tissue's common cell types, which you assemble and they edit,
   - or the top markers per cluster from `tables/ranked_genes_per_cluster.csv`.

Motifs are named with a numeric suffix (`Ascl1-780`, `Runx1-712`); match case-insensitively
against `combined_motifs.h5ad` `var_names` and report anything missing rather than
silently dropping it.

## Known gaps — say these out loud rather than working around them

- **Morphology `_fix` images.** `*_fix.tif` under `spatials/<SHORT_ID>/` (leading zeros
  dropped: RUN1_SEQ001 -> RUN1) are NOT in the tixel coordinate space and no
  registration transform is stored anywhere in the run output. Overlays are drawn on
  `spatial/tissue_hires_image.png`, which is exactly the crop space scaled by
  `tissue_hires_scalef`. If the deck needs the `_fix` base image, ask for the transform.
  If no `_fix` image exists at all, `MorphologyImageMissing` is raised — report it.
- **H&E images are not tracked** in any standard location. Ask the scientist for a path
  or an upload.
- **Panels with the `Mean:` / `Genes:` colorbar annotation** came from Latch's proprietary
  `w_h5` viewer. The data is reproducible; the exact look is not. Say so rather than
  implying a pixel match.
- **The `Mean:` value itself differs** from ours on spatial QC maps even though the
  colorbar endpoints match exactly — the viewer uses a different denominator. Unresolved.

## Reporting results

State what was verified, not just what was produced. The QC table and the peak summaries
can be checked against delivered decks; the feature panels can be sanity-checked against
`medians.csv`. When a figure's source is a file ArchR already wrote, say that you fetched
it rather than implying you recomputed it.

`references/figure_catalog.md` has the full provenance table, the measured geometry, the
QC-table column rules and rounding, and what was verified against which deck.

`references/pipeline.md` describes the processing pipeline that produced the project:
which workflow wrote each Latch directory, the default parameters and filters, and how the
ArchR and SnapATAC branches differ. It also briefly covers the whole-transcriptome and
co-profiling paths. Read it when a user asks where a number came from, why a QC value looks
off, or which parameters a project ran with (`tables/input_parameters.csv` holds the actual
values).
