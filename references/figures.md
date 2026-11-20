# Figure detail

How the individual report figures are built and the traps behind each one. `SKILL.md`
lists the generators and the rules every figure follows; read the matching section here
before producing that figure for the first time in a session.

Contents: the four pre-release checks in detail; marker-gene panels; the Track Browser
slide; the differential motif slide; motif sequence logos.

## The four pre-release checks in detail

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

## Spatial panels share one physical scale

Not one of the four checks, because the figure rarely looks wrong: each sample's panel is
fitted to its own tissue, so every section fills its box. On a project whose sections
differed tenfold in area (0.47, 4.4 and 5.7 mm² by H&E) the smallest was drawn the same size
as the largest, and the zoom spread its tixels into what looked like a sparse grid across a
full-size section.

`shared_spatial_scale(fig, axes=None, pad=0.04)` fixes the ranges after the figure is
built: one x span and one y span for every spatial subplot (the largest extent, padded),
each centred on its own tissue. With equal subplot domains, that means one
pixels-per-unit scale; `panel_scales(fig)` reads plotly's final layout and returns it per
panel, so the claim is checked rather than assumed. Pass `axes=` when the figure also holds
a UMAP, as `gene_activity_block` does.

Only for runs on the same chip type. The spatial coordinates are crop-space pixels and the
crop covers the whole capture area (~7,800 px) for any chip, so a 50-channel and a
220-channel run would be drawn at the same pixel span for very different physical sizes.
`check_shared_scale()` compares `numChannels`, `barcode_filename` and `chip_resolution`
from each run's `spatials/<SHORT_ID>/spatial/metadata.json`, and capture widths within 5%.
`chip_resolution` is often null, so no absolute scale bar is drawn; the relative scale is
exact.

Sanity check against H&E when sizes are disputed: measure each section against its own
scale bar (the H&E images are frequently at different magnifications) and compare the
area ratios with the on-tissue tixel counts. In the case above the two larger sections
agreed to 2% and the smallest, a lacy fragment, to about 20%.

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

Placing this figure on its slide: see `references/deck.md`.

## The differential motif slide

The statistics come from a compare clusters run (`references/differential_analysis.md`);
this section is about drawing the slide from one.

```python
mot = F.fetch_compare_results(project, "Demo_CD19_high_CondA_vs_CondB", "motif")
ann = F.fetch_compare_annotation(project, "Demo_CD19_high_CondA_vs_CondB")
F.apply_compare_annotation(cells, ann, "CD19_high",
                           mapping={"GroupA": "CondA", "GroupB": "CondB"})
```

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

## Motif sequence logos

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
