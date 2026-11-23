"""Slide geometry shared by slide_figures.py (renders on the box) and dense_deck.py (builds
locally), so a figure is drawn at exactly the inch size of the slot it will fill.

Why this exists: a picture placed on a slide is scaled to its slot, and every label in it
is scaled with it. A plotly canvas rendered 1150 px wide and dropped into a 6 in column
lands its 14 px labels at ~7 pt. Drawing each figure at its slot size in inches and
placing it at 100% makes the font size in the code the point size on the slide.

The numbers are measured from assets/atx_template.pptx (13.333 x 7.5 in, Title Only
layout): the title placeholder ends at 0.97 in, a one-line subtitle at ~1.31 in, and the
footer rule sits at 6.805 in (atx_deck._footer_top). A different template needs these
re-measured, not guessed.
"""

SLIDE_W = 13.333
X0 = 0.22            # left/right content margin
Y_SUB = 1.43         # content top on a slide with a one-line subtitle
Y_NOSUB = 1.03       # content top on a slide with a title only
BOTTOM = 6.70        # content bottom: 0.1 in above the footer rule
GAP = 0.12           # gap between side-by-side panels

CONTENT_W = SLIDE_W - 2 * X0                     # 12.89
H_SUB = BOTTOM - Y_SUB                           # 5.27
H_NOSUB = BOTTOM - Y_NOSUB                       # 5.67

#: (width, height) in inches of every figure slot. The slide builders in dense_deck.py
#: place each figure at exactly this size.
SLOTS = {
    # clustering slide (title only): spatial maps left, two stacked UMAPs right
    "cluster_spatial": (8.55, 5.58), "umap_cluster": (4.22, 2.73), "umap_sample": (4.22, 2.73),
    # putative cell types: spatial maps left, labelled UMAP right
    "celltype_spatial": (8.30, H_SUB), "celltype_umap": (4.47, H_SUB),
    # neighborhood grid left, cluster proportion right
    "neighborhood": (7.55, 5.58), "proportion": (5.22, 5.58),
    # morphology: clusters on brightfield | brightfield | H&E, one slide per run
    "morph": (CONTENT_W, H_SUB),
    # one feature per row: spatial panels + UMAP + colourbar + violin by cluster
    "row2": (CONTENT_W, 2.575), "row3": (CONTENT_W, 1.677),
    # track browser: three gene columns + a small annotated UMAP
    "track_col": (3.52, H_SUB), "track_umap": (1.97, 2.55),
    # track set without a UMAP (three columns across the full width)
    "track_col_wide": (4.21, H_SUB),
    "peaks": (6.385, H_SUB), "qc_map": (6.385, 3.25), "genome": (CONTENT_W, H_SUB),
}

#: Type sizes in points. Labels a reader must read on the slide are 12-18 pt; dense
#: secondary text (track y-limits, colourbar ticks on short rows) never below 9.5 pt.
PT = {"header": 17, "panel_title": 15, "label": 13, "tick": 12, "legend": 13, "small": 10.5}
