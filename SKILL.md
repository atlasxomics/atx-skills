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
6. **Assemble a deck only if asked.** `scripts/atx_deck.py`, run locally.
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

## Branding the deck

`build_deck(..., template="<branded>.pptx")` is how output looks like an AtlasXomics
report rather than a generic one. A template is any delivered deck with its slides
stripped (masters, layouts, theme, logo and footer bar all survive; the file drops from
~47 MB to ~1 MB). Content slides are built on the **Title Only** layout, which 14 of the
16 slides in the delivered Demo_001 deck use, and the title and subtitle go into that
layout's own placeholders so fonts, sizes and positions are inherited rather than
restated in code.

Two details that are easy to get wrong:

- Placeholder types must be matched with `PP_PLACEHOLDER`, not raw integers — `13` is
  `SLIDE_NUMBER`, not a title. Treating it as one deletes the page numbers.
- `add_slide` does not clone slide-number, footer or date placeholders from the layout,
  so `_inherit_chrome()` copies them across explicitly; without it a branded template
  still yields unnumbered slides.

### The opening slide is a different layout

Slide 1 is **not** `Title Only`. The blue branded artwork is a *background picture on the
`6_Title Slide` layout*, not part of the theme — so it appears only if the slide is built
on that layout, and a deck whose first slide uses any other layout comes out plain white
no matter how good the template is. `_title_layout()` prefers `6_Title Slide`, and the
deck title and subtitle are written into that layout's own text placeholder. If a scientist
reports "the branded opening slide is missing", this is almost always why.

### Match the delivered slide structures, not just the figure set

Two slides have a fixed arrangement that the generic "row or column" placement will not
produce. Both are built with `_place_columns()`, which takes explicit column groups:

- **Clustering slide** — left column: the per-sample spatial maps under the heading
  "Spatial UMAP". Right column: one combined UMAP coloured by cluster stacked over one
  coloured by `SampleName`, under the heading "UMAP". A single wide spatial panel alone
  leaves half the slide empty and is missing half the content.
- **Neighborhood slide** — left column: the per-sample heatmaps as one figure on a shared
  colour axis, from `atx_figures.neighborhood_grid()`. Right column: the cluster-proportion
  bar. Stacking the heatmaps in a deck column instead shrinks each one to unreadable.

### Visual review is not optional

Geometry checks pass on slides that look wrong. Rendering the deck is what caught a
figure sitting on top of the QC table and wide panels leaving half a slide empty.
LibreOffice is installed on analysis-box, so render and look before calling a deck done:

```
soffice --headless --convert-to pdf --outdir out deck.pptx
python -c "import pypdfium2 as p; d=p.PdfDocument('out/deck.pdf'); [d[i].render(scale=120/72).to_pil().save(f'slide_{i+1:02d}.png') for i in range(len(d))]"
```

Three layout rules came out of that loop and are built into `atx_deck`: images are
arranged as a row or a column by whichever renders them larger (wide panels stack);
table columns are sized from the measured text extent of their widest value — relative
"units" with a floor were not enough, since a 9-unit floor still gave the Sample column
0.93in while "CondC" at 11pt needs ~0.91in plus padding, so it wrapped anyway and grew
the table past its estimated height; and a wrapped title pushes the subtitle and the
content below it down.

Check a claim about appearance against a fresh render, not the last one you looked at.
Reporting a layout as fixed while showing an image rendered before the fix is an easy
mistake to make and an obvious one to the reader.

### Every build is validated structurally

`build_deck()` runs `validate_deck()` on the saved file and raises unless it is clean, so
a broken deck cannot be handed over by accident. Pass `validate=False` only to inspect a
deliberately broken one. The checks, and the bug each one caught:

| Code | Catches |
|---|---|
| `OFF_WITHOUT_EXT` | a placeholder with an offset and no extent — zero width in PowerPoint |
| `CROSSES_FOOTER` | content running over the master's footer rule |
| `OFF_SLIDE` | a shape past a slide edge |
| `SHAPE_OVERLAP` | content printed on top of other content (a figure over the QC table) |
| `MISSING_AUTOFIT` | a title placeholder without `<a:normAutofit/>`, so a long title wraps |
| `ASPECT_DISTORTED` | a picture placed at an aspect ratio its source does not have |
| `MISSING_SLIDE_NUMBER` | chrome not inherited, stripping page numbers |
| `EMPTY_PLACEHOLDER` | a leftover placeholder rendering "Click to add" |

`scripts/validate_selftest.py` plants one fault per check and asserts each one fires —
run it after touching `atx_deck`, because a validation pass nobody has seen fail is not
evidence of anything. It needs the template as its first argument and any PNG as its
second: `python validate_selftest.py assets/atx_template.pptx some.png` prints
`8/8 checks fired as expected`.

For figure coverage rather than slide geometry, `atx_figures.assert_sample_coverage()`
takes the enumerated sample list and one keyword per per-sample figure family, and raises
naming the family and the samples it is missing.

**When a check fires, fix the measurement, not the threshold.** The first clean run
reported six `SHAPE_OVERLAP` findings, all of them real: the subtitle placeholder's
declared box was taller than the one line of text in it, so figures were being placed
inside it. Sizing the box to its content removed the finding; widening the tolerance
would have hidden it.

### ...but LibreOffice is not PowerPoint: two bugs it will not show you

A clean soffice render is necessary and **not sufficient**. These two reached a customer
deck looking perfect in every render, and both must be checked in the XML instead.

- **Set all four of left/top/width/height, never just one.** python-pptx emits `<a:off>`
  as soon as `left` or `top` is assigned, but emits `<a:ext>` only when a size is assigned
  too. A placeholder that inherits its size from the layout and then has only `top` moved
  ends up with an offset and no extent. PowerPoint reads that width as **zero** and
  collapses the text to one character per line spilling off the slide; LibreOffice quietly
  falls back to the layout extent and looks right. Use `_set_box()`, and assert after
  building that no placeholder has an `off` without an `ext`:

  ```python
  tags = re.findall(r'<a:(off|ext)[^/]*/>', ph._element.xml)
  assert not ('off' in tags and 'ext' not in tags)
  ```

- **Content must stop above the footer furniture, which is measured, not assumed.** The
  AtlasXomics master draws a thin cyan rule at 6.83in, a dot at 6.805in and the logo at
  6.97in on a 7.5in slide. A fixed 0.5in bottom margin lets figures run to 7.0in, straight
  over the rule. `_footer_top()` finds the highest shape in the bottom quarter of the
  master and the layouts and returns it; placement stops 0.1in above that (6.705in for
  this template). Never hard-code the margin — a different template moves the rule.

### Density: never estimate where the title ends — turn autofit on and read the box

A bench-made deck is visibly denser than a naively generated one: title and subtitle sit
close together, figures start just under the subtitle and run nearly to the slide edges.
Measuring a delivered deck against a generated one gave the numbers: delivered subtitles
all at **0.99in** (the template's own position), generated at **1.42in**; delivered
picture lefts down to **0.09in**, generated at 0.45in. Every slide inherited the error, so
the deck looked airy throughout even though no single slide was obviously wrong.

The cause was not a margin constant. Delivered title placeholders carry
`<a:bodyPr><a:normAutofit/></a:bodyPr>`; python-pptx writes a bare `<a:bodyPr/>`. With
normAutofit a long title **shrinks to fit one line**; without it the same string **wraps
to two**. The generated builder had been compensating with a character-width heuristic
that predicted a two-line wrap for 57-62 character titles, and pushed the subtitle down
to clear a second line that PowerPoint never draws — then dragged all content with it.

Three rules follow, and the first is the one that matters:

1. Call `_set_autofit()` on both title and subtitle (`word_wrap = True`,
   `MSO_AUTO_SIZE.TEXT_TO_FIT_SHAPE`) so the title shrinks instead of wrapping. Content
   then starts at `title_ph.top + title_ph.height` and the subtitle stays wherever the
   template put it. **Delete any wrap-estimating code**; a character-width guess is a
   second, disagreeing layout engine, and it will drift from the renderer on every title
   it has not been calibrated against.
2. `SIDE_MARGIN = 0.20in` and `PANEL_GAP = 0.12in`, not the half-inch margins a general
   layout routine defaults to. The delivered decks run figures out to the slide edge.
3. Render a figure at roughly the aspect of the cell it will land in. A grid cell fits an
   image by `min(cw/w, ch/h)`, so a 3.5:1 strip in a 2.6:1 cell is width-limited and
   leaves the bottom third of the slide blank no matter how the margins are set. The
   gene-activity blocks were re-rendered at `panel_w=250, height=420` (2.6:1) for exactly
   this reason.

Diagnose this class of problem by rendering the **delivered** deck through the same
soffice path as the generated one. If the delivered title stays on one line and yours
wraps, it is not a missing font or a renderer difference — it is a property of your XML,
and diffing the two `<a:bodyPr>` elements will find it in a minute.

### The deck is an editable document, so do not fuse figures that belong apart

`montage()` composes several PNGs into one picture. That is right when the pieces are one
figure (the neighborhood panels were, until they became a real shared-axis figure), and
wrong when they are not: on the motif slide the spatial panel and its sequence logo were
montaged, so the scientist who opened the deck could not resize the logo without dragging
the panel above it. A deck is handed to someone who will adjust it.

A `grid` cell may therefore be a **list** of paths instead of a single path. The images
are stacked inside the cell as separate pictures, each with its own frame, with heights
allocated from each image's own aspect so the stack fills the cell without distortion:

```python
"grid": [[[spatial_png, logo_png], violin_png, other_png], ...]
```

Reach for `montage()` only when the result is genuinely one object.

### Spatial panels must be isotropic

Plotly stretches each subplot to fill its box, so a tissue section drawn in a tall subplot
comes out elongated — the QC maps had the CondA section at 0.31 width-to-height against 0.48
on the brightfield, a 50% distortion of a physical object. `cluster_maps(equal_aspect=True)`
(the default) sets each y axis's `scaleanchor` to its own x axis at `scaleratio=1`, taking
the anchor from `ax.anchor` so it is right for every subplot. The outline then matches
the morphology image, which is the thing the reader checks it against. `ASPECT_DISTORTED`
in the deck validator catches a *picture* placed at the wrong aspect; it cannot see a
figure that was rendered distorted, so this has to be right at render time.

### Type size is set by the canvas, not by the font number

Plotly font sizes are in layout pixels, and a picture placed on a slide is scaled by
`slide_width / layout_width`. A figure rendered 1200px wide and dropped into a 6in column
is reduced to about a third, so a `size=20` title lands at roughly 7pt — unreadable, and
no amount of raising the number inside a too-large canvas fixes the ratio for the panels
around it. Render near the final size and set the type large: `neighborhood_grid()`
defaults to `cell_px=340` with `title_pt=36`, which lands near 19pt on the slide.
`cluster_proportion(font_pt=16)` does the same for the bar chart's sample labels.

Two panels drawn on the same quantity share one colourbar. Montaging three separately
rendered heatmaps gave three identical bars, wasting a third of the width and implying
the scales were merely coincidentally equal; `neighborhood_grid()` puts all panels on one
`coloraxis` so the shared range is enforced by the figure. Read axis labels back off the
trace (`fig.data[0].x`) rather than from the categorical dtype — `plotly_heatmap()`
re-sorts the matrix, so the raw category order is not guaranteed to match the rows.

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

## Start from the delivered deck, not from the PDF

If a previous deck for this customer (or any recent Data Report) exists as a `.pptx`,
**read its slide inventory before choosing what to generate**:

```python
for s in atx_deck.deck_outline("Demo_001 AtlasXomics Data Report the reference deck.pptx"):
    print(s["index"], s["pictures"], s["title"])
```

and check the built deck against it:

```python
gaps = atx_deck.compare_outlines(delivered_pptx, built_pptx)
gaps["missing"]      # reference slides with no counterpart -- the list that matters
gaps["picture_gap"]  # matched slides reproduced at less depth
```

The report PDFs are a poor substitute for this. A slide list built from the PDF catalog
omitted slide 7, "Numerous differential genes, motifs and peaks between clusters" --
three ArchR marker heatmaps that were already sitting in `figures/heatmaps_all.pdf`.
Nothing failed and no panel was empty; the slide simply was not there, and no check
based on the deck itself can notice a slide that was never specified. Only the reference
inventory can.

`heatmaps_all.pdf` has three pages in order: gene activity scores (GeneScoreMatrix),
peak enrichment (PeakMatrix), motif enrichment. Render them with pypdfium2 at scale
200/72 and place them as three columns; `workflow_figure(project, "marker_heatmaps")`
fetches the PDF.

Two slides in the delivered decks **cannot be generated** and must come from the
scientist: "Visualize and explore with AtlasXPlots" and "Additional Data Outputs
available on LatchBio Workspace" are both screenshots of the Latch web UI. Ask for them
rather than quietly shipping a deck without them. "Bioinformatics Workflow" is a generic
diagram that can be lifted from any delivered deck and reused; the "Appendix" divider is
a title-only slide on the plain `Title Slide` layout (`{"title": "Appendix",
"section": True}`).

Slide wording is worth copying exactly. The opening slide reads "Chart the unexplored
depths of spatial biology" over "Data Report", and the Cd19 motif slide defines its
population as **Cd19 gene activity score higher than 0.6** -- a threshold, not a
hand-drawn lasso, so that slide is reproducible from the data.

## The QC slide carries BOTH spatial metrics

"Samples pass QCs for each assay" is the metrics table plus **two** titled groups of
per-sample spatial maps: **TSS Enrichment** and **log10nFrags**. Rendering only TSS
enrichment leaves out half the slide, and nothing complains -- the slide looks complete.

```python
{"title": "Samples pass QCs for each assay", "table": qc,
 "columns": [[F("..._qc_spatial_TSSEnrichment.png")],
             [F("..._qc_spatial_log10_nFrags.png")]],
 "column_titles": ["TSS Enrichment", "log10nFrags"]}
```

`build_deck` accepts `columns` (or `grid`) after a `table`, so the groups sit under the
metrics table with their headings. Both metrics are plain obs columns, so
`spatial_qc_maps()` or `cluster_maps(color_by=...)` draws them -- in the brightfield frame
like every other spatial panel, and titled by **condition** (CondA, CondB, CondC) rather
than run id, which is what this slide uses.

Validated against a delivered deck: the colourbar endpoints for both metrics reproduce
exactly. (The viewer's `Mean:` annotation still differs slightly on these two panels -- a
known, deferred discrepancy that does NOT affect the gene-activity panels, where the means
do reproduce.)

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

On the deck side this slide needs `column_weights` (four columns: three track panels and
a UMAP). Without it `_place_columns` fits every column to one common height, so the wide
UMAP drags the tall track panels down and a third of the slide is left empty; with
weights each column is fitted independently against its share of the width.

## Before a deck goes out: four checks the figures will not fail on their own

Each of these produced a plausible-looking deck that was wrong, and none of them raised
anything. Run through all three before showing a deck to anyone.

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
