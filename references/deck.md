# Building the AtlasXomics deck

Read this before calling `atx_deck.build_deck()`. The figure rules in `SKILL.md`
(orientation, cropping, every sample, one palette, isotropic panels) still apply; this
file covers turning those figures into slides. `atx_deck.py` runs locally, not on the box.

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
  A single-sample project has no proportion bar (`render_figures.py` skips it): give the
  one heatmap the whole slide.

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

## Type size is set by the canvas, not by the font number

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

## Track Browser slide placement

On the deck side this slide needs `column_weights` (four columns: three track panels and
a UMAP). Without it `_place_columns` fits every column to one common height, so the wide
UMAP drags the tall track panels down and a third of the slide is left empty; with
weights each column is fitted independently against its share of the width.
