"""atx_deck — assemble generated figures and tables into a PowerPoint deck.

Runs LOCALLY, not on analysis-box: by this stage the inputs are small PNGs and CSVs that
have already been harvested, so assembly needs no Latch access and no project data. Keeping
it local also lets the box idle-stop while slide layout is iterated.

Pass `template=` a .pptx to inherit AtlasXomics branding (master, fonts, footer, logo).
Without one the deck is clean but unbranded — python-pptx cannot invent a theme.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional, Sequence

from PIL import Image
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN
from pptx.util import Emu, Inches, Pt

WIDESCREEN = (Inches(13.333), Inches(7.5))
TITLE_COLOR = RGBColor(0x1A, 0x1A, 0x1A)
SUB_COLOR = RGBColor(0x44, 0x44, 0x44)
HEADER_FILL = RGBColor(0x29, 0xB6, 0xE8)   # the cyan header band used in the decks

#: Width of one character as a fraction of the font size, for the template's bold title
#: face. Calibrated from the delivered deck (see _fill_placeholders); raising it causes
#: phantom title wraps.
TITLE_EM = 0.47

#: The delivered decks run figures close to the slide edge -- measured picture lefts go
#: down to 0.09in -- so content is not inset the way a generic layout would be.
SIDE_MARGIN = Inches(0.20)
PANEL_GAP = Inches(0.12)


#: The AtlasXomics decks put 14 of 16 slides on "Title Only" — it carries the title,
#: the subtitle line, the slide-number placeholder, and inherits the footer bar and logo
#: from the master. Prefer it so branding comes from the template rather than being redrawn.
PREFERRED_LAYOUTS = ("Title Only",)


def _content_layout(prs):
    """The layout to build content slides on."""
    by_name = {lay.name: lay for lay in prs.slide_layouts}
    for want in PREFERRED_LAYOUTS:
        if want in by_name:
            return by_name[want]
    best, best_n = prs.slide_layouts[0], 10 ** 6
    for lay in prs.slide_layouts:
        n = len(lay.placeholders)
        if n < best_n:
            best, best_n = lay, n
    return best


#: The branded opening slide (blue graphic + logo) is the template's "6_Title Slide"
#: layout, which carries the artwork as a layout BACKGROUND PICTURE. The stock
#: "Title Slide" layout is plain white — using it is why the opening slide looked unbranded.
TITLE_LAYOUTS = ("6_Title Slide", "Title Slide")

#: The delivered decks put section dividers ("Appendix") on the PLAIN `Title Slide`
#: layout, not the branded `6_Title Slide` the opening slide uses.
SECTION_LAYOUTS = ("Title Slide", "Section Header", "Title Only")


def _title_layout(prs):
    by_name = {lay.name: lay for lay in prs.slide_layouts}
    for want in TITLE_LAYOUTS:
        if want in by_name:
            return by_name[want]
    return _content_layout(prs)


def _section_layout(prs):
    """Layout for a divider slide, which is not the branded opening layout."""
    by_name = {l.name: l for l in prs.slide_layouts}
    for n in SECTION_LAYOUTS:
        if n in by_name:
            return by_name[n]
    return _content_layout(prs)


def _inherit_chrome(slide, layout):
    """Copy the layout's slide-number / footer / date placeholders onto the slide.

    python-pptx's add_slide clones only content placeholders, so a template whose master
    draws page numbers still produces slides without them. Copying the elements verbatim
    keeps the master's position and formatting, and the <a:fld> field keeps auto-updating.
    """
    import copy as _copy
    from pptx.enum.shapes import PP_PLACEHOLDER as PP

    CHROME = (PP.SLIDE_NUMBER, PP.FOOTER, PP.DATE)
    have = {ph.placeholder_format.type for ph in slide.placeholders}
    spTree = slide.shapes._spTree
    for ph in layout.placeholders:
        t_ = ph.placeholder_format.type
        if t_ in CHROME and t_ not in have:
            spTree.append(_copy.deepcopy(ph._element))


def _fill_placeholders(slide, title, subtitle):
    """Populate the layout's own title/subtitle placeholders.

    Using the placeholders rather than new textboxes is what makes the output inherit the
    template's fonts, sizes and position, so a branded template needs no styling code here.

    Placeholder types are matched by enum, not by number: PP_PLACEHOLDER.TITLE is 1 while
    13 is SLIDE_NUMBER, and confusing them silently strips the slide numbers the master
    supplies. Footer, date and slide-number placeholders are always left alone.
    Returns the y offset content should start below, or None if the layout has no title.
    """
    from pptx.enum.shapes import PP_PLACEHOLDER as PP

    TITLES = (PP.TITLE, PP.CENTER_TITLE)
    SUBS = (PP.SUBTITLE, PP.BODY)
    CHROME = (PP.SLIDE_NUMBER, PP.FOOTER, PP.DATE, PP.HEADER)

    title_ph = sub_ph = None
    for ph in slide.placeholders:
        t_ = ph.placeholder_format.type
        if title_ph is None and t_ in TITLES:
            title_ph = ph
        elif sub_ph is None and t_ in SUBS:
            sub_ph = ph
    if title_ph is None:
        return None

    title_ph.text_frame.text = title
    _set_autofit(title_ph)
    t_pt = 28.0
    try:
        t_pt = float(title_ph.text_frame.paragraphs[0].font.size.pt)
    except Exception:
        pass
    # bold title faces run wider than body text, and the text frame's own insets shrink
    # the usable width -- ignoring both under-counts the wrap and drops the subtitle on
    # top of the title's second line
    # No line-count estimate and no repositioning. The delivered decks set <a:normAutofit/>
    # on the title, so PowerPoint shrinks an over-long title to fit its 0.57in box instead
    # of wrapping; the subtitle therefore never needs to move and stays at the template's
    # own 0.99in. Estimating the wrap instead was what pushed the subtitle -- and with it
    # every figure -- a third of an inch down the slide on EVERY slide.
    bottom = title_ph.top + title_ph.height
    if subtitle is not None and sub_ph is not None:
        sub_ph.text_frame.text = subtitle
        _set_autofit(sub_ph)
        s_pt = 13.0
        try:
            s_pt = float(sub_ph.text_frame.paragraphs[0].font.size.pt)
        except Exception:
            pass
        s_lines = _est_lines(subtitle, sub_ph.width / 914400.0, s_pt)
        # resize the box to the text it actually holds. The layout's subtitle placeholder
        # is taller than one line of 14pt, so trusting its declared height either wastes a
        # strip of slide or -- if content is placed from an estimate instead -- starts the
        # figures inside the still-declared box. Size it, then use its real bottom.
        _set_box(sub_ph, sub_ph.left, sub_ph.top, sub_ph.width,
                 Inches(max(0.30, s_lines * s_pt * 1.35 / 72.0)))
        bottom = max(bottom, sub_ph.top + sub_ph.height)
    elif sub_ph is not None:
        sub_ph._element.getparent().remove(sub_ph._element)

    # drop only empty, non-chrome placeholders so the slide has no "Click to add" prompts
    for ph in list(slide.placeholders):
        if ph is title_ph or ph is sub_ph:
            continue
        if ph.placeholder_format.type in CHROME:
            continue
        if not ph.text_frame.text.strip():
            ph._element.getparent().remove(ph._element)
    return bottom + Inches(0.12)


def _add_title(slide, title, subtitle=None, width=Inches(12.5)):
    tb = slide.shapes.add_textbox(Inches(0.45), Inches(0.28), width, Inches(0.8))
    p = tb.text_frame.paragraphs[0]
    p.text = title
    p.font.size = Pt(28); p.font.bold = True; p.font.color.rgb = TITLE_COLOR
    if subtitle:
        sb = slide.shapes.add_textbox(Inches(0.45), Inches(1.03), width, Inches(0.45))
        sp = sb.text_frame.paragraphs[0]
        sp.text = subtitle
        sp.font.size = Pt(13); sp.font.color.rgb = SUB_COLOR
    return Inches(1.55) if subtitle else Inches(1.25)


def _est_lines(text, width_in, pt, em=0.52, pad_in=0.0):
    """Rough wrapped-line count. python-pptx cannot measure rendered text, so estimate
    from an average glyph width. Table cells need a wider `em` and the cell padding
    subtracted — under-estimating here is what lets a figure land on top of a table."""
    if not text:
        return 0
    usable = max(0.2, width_in - pad_in)
    chars_per_line = max(1, int(usable * 72.0 / (em * pt)))
    return max(1, -(-len(str(text)) // chars_per_line))


def _best_arrangement(sizes, avail_w, avail_h, gap):
    """Choose a row or a column, whichever renders the images larger.

    A row of very wide panels (the gene-activity figures are ~2.5:1) is width-bound and
    leaves most of the slide empty; stacking them uses the height instead. Picking by
    total rendered area keeps both cases right without a per-slide flag.
    """
    n = len(sizes)
    ars = [w / h for w, h in sizes]

    row_h = min((avail_w - gap * (n - 1)) / sum(ars), avail_h)
    row_area = sum((row_h * ar) * row_h for ar in ars)

    col_each = (avail_h - gap * (n - 1)) / n
    col_hs = [min(col_each, avail_w / ar) for ar in ars]
    col_area = sum((h * ar) * h for h, ar in zip(col_hs, ars))

    return ("row", row_h) if row_area >= col_area else ("column", col_each)


def _place_images(slide, images, top, prs, gap=PANEL_GAP,
                  bottom_margin=None, side_margin=SIDE_MARGIN):
    """Lay images out, preserving aspect, choosing a row or a column by which is larger."""
    images = [i for i in images if i and os.path.exists(i)]
    if not images:
        return
    floor = (float(prs.slide_height - bottom_margin) if bottom_margin is not None
             else _footer_top(prs) - float(Inches(0.10)))
    avail_w = float(prs.slide_width - 2 * side_margin)
    avail_h = float(floor - top)
    if avail_h <= 0:
        return
    sizes = [Image.open(i).size for i in images]
    mode, h = _best_arrangement(sizes, avail_w, avail_h, float(gap))

    if mode == "row":
        widths = [h * (w / hh) for w, hh in sizes]
        x = side_margin + (avail_w - sum(widths) - float(gap) * (len(images) - 1)) / 2
        for img, wpx in zip(images, widths):
            slide.shapes.add_picture(img, Emu(int(x)), Emu(int(top)),
                                     width=Emu(int(wpx)), height=Emu(int(h)))
            x += wpx + float(gap)
    else:
        # h is the per-image height budget; cap each by the width it would need
        placed = []
        for (w, hh) in sizes:
            ar = w / hh
            ih = min(h, avail_w / ar)
            placed.append((ih * ar, ih))
        used = sum(ih for _, ih in placed) + float(gap) * (len(placed) - 1)
        y = float(top) + max(0.0, (avail_h - used) / 2)
        for img, (iw, ih) in zip(images, placed):
            slide.shapes.add_picture(img, Emu(int(side_margin + (avail_w - iw) / 2)),
                                     Emu(int(y)), width=Emu(int(iw)), height=Emu(int(ih)))
            y += ih + float(gap)





def _set_box(shape, left, top, width, height):
    """Position a shape by writing ALL FOUR of left/top/width/height.

    python-pptx writes `<a:off>` as soon as one of left/top is set, but only emits
    `<a:ext>` when a size is set too. A placeholder that inherits its size from the layout
    and then has only `top` assigned ends up with an offset and no extent: PowerPoint reads
    the width as zero and collapses the text to one character per line running off the
    slide. LibreOffice falls back to the layout extent instead, so the bug is INVISIBLE in
    a soffice render -- it only shows up in PowerPoint. Always set the full box.
    """
    shape.left, shape.top, shape.width, shape.height = (
        Emu(int(left)), Emu(int(top)), Emu(int(width)), Emu(int(height)))


def _set_autofit(shape):
    """Mark a text frame as shrink-to-fit, the way the delivered decks' titles are.

    Without <a:normAutofit/> a long title WRAPS and overflows its placeholder; with it the
    renderer scales the text down and keeps one line. Matching the template here removes
    the need to predict wrapping at all.
    """
    try:
        from pptx.enum.text import MSO_AUTO_SIZE
        shape.text_frame.word_wrap = True
        shape.text_frame.auto_size = MSO_AUTO_SIZE.TEXT_TO_FIT_SHAPE
    except Exception:
        pass


def _footer_top(prs, default_margin=Inches(0.5)):
    """Y of the highest piece of footer furniture, so content can stop above it.

    The AtlasXomics master puts a thin cyan rule, a dot and the logo along the bottom
    (rule at 6.83in, dot at 6.80in, logo at 6.97in on a 7.5in slide). A fixed bottom
    margin of 0.5in leaves content running to 7.0in, straight over the rule. Measure the
    furniture instead of guessing: anything in the bottom quarter of the master or the
    content layouts counts.
    """
    limit = float(prs.slide_height)
    zone = limit * 0.75
    sources = [prs.slide_master.shapes]
    sources += [l.shapes for l in prs.slide_layouts]
    for shapes in sources:
        for sh in shapes:
            try:
                top = float(sh.top)
            except (TypeError, ValueError):
                continue
            if zone <= top < limit:
                limit = min(limit, top)
    if limit >= float(prs.slide_height):
        limit = float(prs.slide_height - default_margin)
    return limit


def montage(paths, out_path, cols=2, pad=14, bg=(255, 255, 255)):
    """Compose several PNGs into one grid image, left to right then top to bottom.

    The delivered slides put the three per-sample neighborhood heatmaps in a 2-up grid
    beside the proportion bar. Stacking them in a deck column would make each one tiny,
    so they are composed into a single picture first and placed as one image.

    Cells are sized to the widest and tallest input and each panel is centred in its cell,
    so panels of slightly different pixel sizes stay aligned and undistorted.
    """
    paths = [p for p in paths if p and os.path.exists(p)]
    if not paths:
        raise ValueError("montage: no existing input images")
    ims = [Image.open(p).convert("RGB") for p in paths]
    cw = max(i.width for i in ims)
    ch = max(i.height for i in ims)
    rows = -(-len(ims) // cols)
    W = cols * cw + (cols + 1) * pad
    H = rows * ch + (rows + 1) * pad
    canvas = Image.new("RGB", (W, H), bg)
    for k, im in enumerate(ims):
        r, cidx = divmod(k, cols)
        x = pad + cidx * (cw + pad) + (cw - im.width) // 2
        y = pad + r * (ch + pad) + (ch - im.height) // 2
        canvas.paste(im, (x, y))
    canvas.save(out_path)
    return out_path

def _place_columns(slide, columns, top, prs, column_titles=None, gap=PANEL_GAP,
                   bottom_margin=None, side_margin=SIDE_MARGIN,
                   title_pt=14, column_weights=None):
    """Place images in explicit side-by-side columns, each stacked vertically.

    Used for the clustering slide, whose delivered layout is a wide spatial panel on the
    left and two UMAPs stacked on the right -- an arrangement `_best_arrangement` would
    never pick, since it only chooses between one row and one column.

    Each column is scaled to a common height; its width follows from the aspect ratio of
    its stack, so no image is distorted.
    """
    cols = [[i for i in c if i and os.path.exists(i)] for c in columns]
    cols = [c for c in cols if c]
    if not cols:
        return float(top)
    g = float(gap)
    head_h = Inches(0.34) if column_titles else 0
    floor = (float(prs.slide_height - bottom_margin) if bottom_margin is not None
             else _footer_top(prs) - float(Inches(0.10)))
    avail_w = float(prs.slide_width - 2 * side_margin)
    avail_h = float(floor - top - head_h)
    if avail_h <= 0:
        return float(top)

    # aspect of each column's stack at a common column width
    aspects = []
    for c in cols:
        inv = sum(h / w for w, h in (Image.open(i).size for i in c))
        aspects.append(1.0 / inv)          # width / height of the stack, gaps aside
    n_gaps_v = [len(c) - 1 for c in cols]
    usable_h = avail_h - g * max(n_gaps_v)

    if column_weights:
        # Each column gets a share of the width and is fitted INDEPENDENTLY, capped by the
        # available height. A single common height is wrong when one column is much wider
        # than the others -- a square UMAP beside three tall track panels forces every
        # column down to the width-limited height and leaves a third of the slide empty.
        wts = [float(w) for w in column_weights][:len(cols)]
        wts += [1.0] * (len(cols) - len(wts))
        total_w = avail_w - g * (len(cols) - 1)
        shares = [total_w * w / sum(wts) for w in wts]
        widths, heights = [], []
        for a, s in zip(aspects, shares):
            h = min(usable_h, s / a)
            heights.append(h)
            widths.append(h * a)
    else:
        total_a = sum(aspects)
        h_fit_w = (avail_w - g * (len(cols) - 1)) / total_a if total_a else avail_h
        height = min(usable_h, h_fit_w)
        widths = [a * height for a in aspects]
        heights = [height] * len(cols)

    x = side_margin + (avail_w - sum(widths) - g * (len(cols) - 1)) / 2
    bottom = float(top)
    for ci, (col, cw) in enumerate(zip(cols, widths)):
        y = float(top)
        if column_titles and ci < len(column_titles) and column_titles[ci]:
            tb = slide.shapes.add_textbox(Emu(int(x)), Emu(int(y)), Emu(int(cw)), head_h)
            tf = tb.text_frame
            tf.word_wrap = True
            tf.paragraphs[0].text = column_titles[ci]
            tf.paragraphs[0].alignment = PP_ALIGN.CENTER
            tf.paragraphs[0].runs[0].font.size = Pt(title_pt)
            y += float(head_h)
        for img in col:
            w, h = Image.open(img).size
            ih = cw * h / w
            slide.shapes.add_picture(img, Emu(int(x)), Emu(int(y)),
                                     width=Emu(int(cw)), height=Emu(int(ih)))
            y += ih + g
        bottom = max(bottom, y)
        x += cw + g
    return bottom



def _place_grid(slide, rows, top, prs, gap=PANEL_GAP, bottom_margin=None,
                side_margin=SIDE_MARGIN):
    """Place images in a true grid of EQUAL cells, each image fitted inside its cell.

    Different from `_place_columns`, which scales a whole column as one stack: there, a
    square figure sharing a column with a wide one shrinks the wide one to match. A grid
    keeps every cell the same size, which is what the delivered gene-activity slide does
    with its three gene blocks and the annotated UMAP.

    A cell may be a LIST of paths instead of a single path. The images are then stacked
    inside that cell as SEPARATE pictures, each with its own frame, rather than being
    composed into one PNG with `montage()`. That distinction matters to whoever opens the
    deck: a montage is one object, so nudging the motif logo without also moving the
    spatial panel above it is impossible. Heights within the cell are allocated from each
    image's own aspect, so the stack fills the cell without distorting anything.
    """
    def _exists(c):
        if isinstance(c, (list, tuple)):
            return [i for i in c if i and os.path.exists(i)]
        return c if (c and os.path.exists(c)) else None

    rows = [[_exists(i) for i in r] for r in rows]
    rows = [[i for i in r if i] for r in rows]
    rows = [r for r in rows if r]
    if not rows:
        return float(top)
    g = float(gap)
    floor = (float(prs.slide_height - bottom_margin) if bottom_margin is not None
             else _footer_top(prs) - float(Inches(0.10)))
    avail_w = float(prs.slide_width - 2 * side_margin)
    avail_h = float(floor - top)
    if avail_h <= 0:
        return float(top)
    ncols = max(len(r) for r in rows)
    nrows = len(rows)
    cw = (avail_w - g * (ncols - 1)) / ncols
    ch = (avail_h - g * (nrows - 1)) / nrows

    for ri, row in enumerate(rows):
        y0 = float(top) + ri * (ch + g)
        for ci, cell in enumerate(row):
            cx = float(side_margin) + ci * (cw + g)
            if not isinstance(cell, (list, tuple)):
                cell = [cell]
            sizes = [Image.open(i).size for i in cell]
            n = len(cell)
            inner = g * 0.5
            # Height each image would take at full cell width, then shrunk to fit.
            nat = [cw * h / w for (w, h) in sizes]
            k = min(1.0, (ch - inner * (n - 1)) / sum(nat)) if sum(nat) else 1.0
            used = sum(nat) * k + inner * (n - 1)
            yy = y0 + (ch - used) / 2
            for img, (w, h) in zip(cell, sizes):
                ih = (cw * h / w) * k
                iw = cw * k
                slide.shapes.add_picture(img, Emu(int(cx + (cw - iw) / 2)), Emu(int(yy)),
                                         width=Emu(int(iw)), height=Emu(int(ih)))
                yy += ih + inner
    return float(top) + nrows * ch + (nrows - 1) * g


def _add_table(slide, df, top, prs, side_margin=Inches(0.7), font_pt=11):
    """Add the table and return its estimated bottom edge.

    Columns are sized to their widest content rather than split evenly: an evenly split
    table wraps `RUN1_SEQ001` onto two lines, which grows the row and pushes the real
    table past any fixed-height estimate — that is what put a figure over the QC table.
    """
    rows, cols = df.shape[0] + 1, df.shape[1]
    width = prs.slide_width - 2 * side_margin
    total_in = width / 914400.0

    # Width each column actually needs, in inches, from the widest string it holds.
    # Relative "units" were not enough: a 9-unit floor still gave Sample 0.93in, and
    # "CondC" at 11pt needs ~0.91in plus padding, so it wrapped anyway.
    pad_in = 0.18
    em = 0.62
    def _needed(s):
        return len(str(s)) * em * font_pt / 72.0 + pad_in
    required = [max([_needed(col)] + [_needed(v) for v in df[col]]) for col in df.columns]

    if sum(required) <= total_in:
        slack = total_in - sum(required)
        col_in = [r + slack * (r / sum(required)) for r in required]
    else:
        # genuinely too wide for the slide: shrink proportionally and accept wrapping
        col_in = [r * total_in / sum(required) for r in required]

    line_in = font_pt * 1.30 / 72.0
    def _row_lines(values):
        return max(_est_lines(v, col_in[j], font_pt, em=0.62, pad_in=pad_in)
                   for j, v in enumerate(values))
    heights = [_row_lines(list(df.columns)) * line_in + 0.14]
    heights += [_row_lines(list(r)) * line_in + 0.12 for _, r in df.iterrows()]
    total_h = Inches(sum(heights))

    shape = slide.shapes.add_table(rows, cols, side_margin, Emu(int(top)), width, total_h)
    table = shape.table
    for j, w_in in enumerate(col_in):
        table.columns[j].width = Inches(w_in)
    for j, name in enumerate(df.columns):
        cell = table.cell(0, j)
        cell.text = str(name)
        para = cell.text_frame.paragraphs[0]
        para.font.size = Pt(font_pt); para.font.bold = True
        cell.fill.solid(); cell.fill.fore_color.rgb = HEADER_FILL
    for i, (_, row) in enumerate(df.iterrows(), start=1):
        for j, val in enumerate(row):
            cell = table.cell(i, j)
            cell.text = "" if val is None else str(val)
            cell.text_frame.paragraphs[0].font.size = Pt(font_pt)
    for i, h_in in enumerate(heights):
        table.rows[i].height = Inches(h_in)
    return Emu(int(top)) + total_h


def build_deck(slides: Sequence[Dict], out_path: str, template: Optional[str] = None,
               deck_title: Optional[str] = None, subtitle: Optional[str] = None, validate: bool = True) -> str:
    """Assemble slides into a .pptx.

    Each slide dict takes: `title`, optional `subtitle`, optional `images` (list of paths,
    laid out in one row), optional `table` (a DataFrame rendered as a native, editable
    PowerPoint table), optional `notes` (speaker notes).
    """
    prs = Presentation(template) if template else Presentation()
    if not template:
        prs.slide_width, prs.slide_height = WIDESCREEN
    layout = _content_layout(prs)

    if deck_title:
        lay = _title_layout(prs) if template else layout
        s = prs.slides.add_slide(lay)
        # prefer the layout's own text placeholder so the branded artwork is preserved
        ph = next((p for p in s.placeholders), None)
        if template and ph is not None:
            tf = ph.text_frame
            tf.text = deck_title
            if subtitle:
                sp = tf.add_paragraph()
                sp.text = subtitle
        else:
            tb = s.shapes.add_textbox(Inches(0.8), Inches(2.6),
                                      prs.slide_width - Inches(1.6), Inches(1.4))
            p = tb.text_frame.paragraphs[0]
            p.text = deck_title; p.font.size = Pt(40); p.font.bold = True
            p.font.color.rgb = TITLE_COLOR
            if subtitle:
                sp = tb.text_frame.add_paragraph()
                sp.text = subtitle; sp.font.size = Pt(18); sp.font.color.rgb = SUB_COLOR

    for spec in slides:
        if spec.get("section"):
            sec = _section_layout(prs)
            s = prs.slides.add_slide(sec)
            _inherit_chrome(s, sec)
            _fill_placeholders(s, spec["title"], spec.get("subtitle"))
            if spec.get("notes"):
                s.notes_slide.notes_text_frame.text = spec["notes"]
            continue
        s = prs.slides.add_slide(layout)
        _inherit_chrome(s, layout)
        top = _fill_placeholders(s, spec["title"], spec.get("subtitle"))
        if top is None:
            top = _add_title(s, spec["title"], spec.get("subtitle"))
        if spec.get("table") is not None:
            top = Emu(int(top)) + Inches(0.1)
            bottom = _add_table(s, spec["table"], top, prs)
            below = Emu(int(bottom)) + Inches(0.3)
            # a table can be followed by titled column groups, which is how the QC slide
            # puts "TSS Enrichment" and "log10nFrags" side by side under the metrics table
            if spec.get("columns"):
                _place_columns(s, spec["columns"], below, prs,
                               column_titles=spec.get("column_titles"),
                               column_weights=spec.get("column_weights"))
            elif spec.get("grid"):
                _place_grid(s, spec["grid"], below, prs)
            elif spec.get("images"):
                _place_images(s, spec["images"], below, prs)
        elif spec.get("grid"):
            _place_grid(s, spec["grid"], top, prs)
        elif spec.get("columns"):
            _place_columns(s, spec["columns"], top, prs,
                           column_titles=spec.get("column_titles"),
                           column_weights=spec.get("column_weights"))
        elif spec.get("images"):
            _place_images(s, spec["images"], top, prs)
        if spec.get("notes"):
            s.notes_slide.notes_text_frame.text = spec["notes"]

    prs.save(out_path)
    if validate:
        # structural checks run on every build: they catch the PowerPoint-only failures
        # that a LibreOffice render shows as perfect. Pass validate=False only to inspect
        # a deliberately broken deck.
        assert_deck_valid(out_path)
    return out_path


# --------------------------------------------------------------- validation
#
# Every check below exists because that exact failure reached a built deck. Three of them
# (OFF_WITHOUT_EXT, CROSSES_FOOTER, SHAPE_OVERLAP) were invisible in a LibreOffice render,
# so none of this can be replaced by looking at a picture.

#: Shapes that belong to the template's furniture and are exempt from content rules.
_CHROME_NAMES = ("slide number", "footer", "date")


def _is_chrome(shape):
    try:
        from pptx.enum.shapes import PP_PLACEHOLDER as PP
        if shape.is_placeholder and shape.placeholder_format.type in (
                PP.SLIDE_NUMBER, PP.FOOTER, PP.DATE, PP.HEADER):
            return True
    except Exception:
        pass
    return any(k in (shape.name or "").lower() for k in _CHROME_NAMES)


def _box(shape):
    """(left, top, right, bottom) in EMU, or None when the shape has no geometry."""
    try:
        l, t_, w, h = shape.left, shape.top, shape.width, shape.height
        if None in (l, t_, w, h):
            return None
        return (float(l), float(t_), float(l) + float(w), float(t_) + float(h))
    except (TypeError, ValueError):
        return None


def _overlap_area(a, b):
    dx = min(a[2], b[2]) - max(a[0], b[0])
    dy = min(a[3], b[3]) - max(a[1], b[1])
    return dx * dy if dx > 0 and dy > 0 else 0.0


def validate_deck(deck, floor=None, aspect_tol=0.02, overlap_frac=0.02):
    """Structural checks on an assembled deck. Returns a list of findings.

    `deck` is a path or a Presentation. Each finding is a dict with `code`, `slide`,
    `shape` and `detail`. An empty list means every check passed.

    Checks, and the bug each one caught:

    - `OFF_WITHOUT_EXT`  a placeholder carries `<a:off>` with no `<a:ext>`. python-pptx
      writes the offset as soon as left/top is set but the extent only when a size is set,
      so a repositioned inherited placeholder ends up sizeless. PowerPoint reads width 0
      and runs the text down the page one character per line; LibreOffice falls back to
      the layout extent and looks perfect. Use `_set_box()` to avoid it.
    - `CROSSES_FOOTER`   content runs over the master's footer rule, because a fixed
      bottom margin was used instead of the measured furniture position.
    - `OFF_SLIDE`        a shape extends past a slide edge.
    - `SHAPE_OVERLAP`    two pieces of content sit on top of each other -- this is how a
      figure ended up printed over the QC table.
    - `MISSING_AUTOFIT` a title without `<a:normAutofit/>`, which wraps and overflows
      onto the content below instead of shrinking to fit its placeholder.
    - `ASPECT_DISTORTED` a picture was placed at an aspect ratio its source does not have.
    - `MISSING_SLIDE_NUMBER` the slide-number placeholder was not inherited, which strips
      page numbers from an otherwise branded deck.
    - `EMPTY_PLACEHOLDER`    an empty placeholder left behind renders as "Click to add".
    """
    import re as _re
    prs = deck if isinstance(deck, Presentation().__class__) else Presentation(deck)
    if floor is None:
        floor = _footer_top(prs)
    out = []

    def add(code, slide_i, shape, detail):
        out.append({"code": code, "slide": slide_i,
                    "shape": getattr(shape, "name", None), "detail": detail})

    for i, slide in enumerate(prs.slides, 1):
        is_title_slide = "title slide" in slide.slide_layout.name.lower()
        content = []

        for ph in slide.placeholders:
            tags = _re.findall(r"<a:(off|ext)[^/]*/>", ph._element.xml)
            if "off" in tags and "ext" not in tags:
                add("OFF_WITHOUT_EXT", i, ph,
                    "offset written without an extent: PowerPoint will read width 0")
            if not _is_chrome(ph) and not ph.text_frame.text.strip():
                add("EMPTY_PLACEHOLDER", i, ph, "empty placeholder renders a prompt")

        # A title must carry <a:normAutofit/>, as the delivered decks' do. Without it an
        # over-long title wraps and overflows onto whatever sits below instead of being
        # scaled down -- the defect this replaces a wrap ESTIMATE with a structural check.
        try:
            from pptx.enum.shapes import PP_PLACEHOLDER as PP
            t_ph = next((p for p in slide.placeholders
                         if p.placeholder_format.type in (PP.TITLE, PP.CENTER_TITLE)), None)
            if t_ph is not None and t_ph.text_frame.text.strip():
                if "normAutofit" not in t_ph._element.xml:
                    add("MISSING_AUTOFIT", i, t_ph,
                        "title has no <a:normAutofit/>: a long title will wrap and "
                        "overflow instead of shrinking to fit")
        except Exception:
            pass

        if not is_title_slide and not any(_is_chrome(sh) for sh in slide.shapes):
            add("MISSING_SLIDE_NUMBER", i, None,
                "no slide-number/footer placeholder was inherited from the layout")

        for sh in slide.shapes:
            b = _box(sh)
            if b is None or _is_chrome(sh):
                continue
            if b[0] < -1 or b[1] < -1 or b[2] > float(prs.slide_width) + 1 \
                    or b[3] > float(prs.slide_height) + 1:
                add("OFF_SLIDE", i, sh,
                    f"box {[round(v / 914400.0, 2) for v in b]} in leaves the slide")
            if b[3] > floor + 1:
                add("CROSSES_FOOTER", i, sh,
                    f"bottom {b[3] / 914400.0:.3f} in is below the footer line "
                    f"at {floor / 914400.0:.3f} in")
            if sh.shape_type == 13:          # PICTURE
                try:
                    pw, ph_ = sh.image.size
                    placed = float(sh.width) / float(sh.height)
                    native = float(pw) / float(ph_)
                    if abs(placed - native) / native > aspect_tol:
                        add("ASPECT_DISTORTED", i, sh,
                            f"placed {placed:.3f} vs native {native:.3f}")
                except Exception:
                    pass
            content.append((sh, b))

        for a_i in range(len(content)):
            for b_i in range(a_i + 1, len(content)):
                sa, ba = content[a_i]
                sb, bb = content[b_i]
                ov = _overlap_area(ba, bb)
                if ov <= 0:
                    continue
                smaller = min((ba[2] - ba[0]) * (ba[3] - ba[1]),
                              (bb[2] - bb[0]) * (bb[3] - bb[1]))
                if smaller > 0 and ov / smaller > overlap_frac:
                    add("SHAPE_OVERLAP", i, sa,
                        f"{sa.name!r} overlaps {sb.name!r} over "
                        f"{100 * ov / smaller:.0f}% of the smaller shape")
    return out


def assert_deck_valid(deck, **kw):
    """Raise with every finding listed, or return the (empty) finding list."""
    findings = validate_deck(deck, **kw)
    if findings:
        lines = [f"  slide {f['slide']}: {f['code']} {f['shape'] or ''} -- {f['detail']}"
                 for f in findings]
        raise ValueError("deck validation failed:\n" + "\n".join(lines))
    return findings


# ------------------------------------------------- comparing against a delivered deck

def deck_outline(deck):
    """Slide-by-slide inventory of a .pptx: index, title, layout, picture and table counts.

    Read the DELIVERED deck with this before choosing which figures to generate. Building
    the slide list from the report PDFs instead is how the marker-heatmap slide
    ("Numerous differential genes, motifs and peaks between clusters") was missed
    entirely: it is a class-A figure ArchR had already written to figures/heatmaps_all.pdf,
    so nothing failed and nothing was empty -- the slide simply was not there.
    """
    prs = deck if isinstance(deck, Presentation().__class__) else Presentation(deck)
    out = []
    for i, s in enumerate(prs.slides, 1):
        texts = [sh.text_frame.text.strip() for sh in s.shapes
                 if sh.has_text_frame and sh.text_frame.text.strip()]
        out.append({
            "index": i,
            "title": texts[0] if texts else "",
            "layout": s.slide_layout.name,
            "pictures": sum(1 for sh in s.shapes if sh.shape_type == 13),
            "tables": sum(1 for sh in s.shapes if getattr(sh, "has_table", False)),
        })
    return out


def compare_outlines(reference, built, cutoff=0.6):
    """Match a built deck's slides against a reference deck's by title similarity.

    Returns {"missing": [...], "extra": [...], "matched": [(ref, built), ...],
    "picture_gap": [...]}. `missing` is the list that matters: reference slides with no
    counterpart. `picture_gap` flags matched slides carrying fewer pictures than the
    reference, which is the signature of a slide that was reproduced at less depth (for
    example one gene-activity panel where the delivered deck shows five).
    """
    import difflib
    ref = deck_outline(reference) if not isinstance(reference, list) else reference
    got = deck_outline(built) if not isinstance(built, list) else built
    got_titles = [g["title"] for g in got]
    used, matched, missing, gap = set(), [], [], []
    for r in ref:
        hit = difflib.get_close_matches(r["title"], got_titles, n=1, cutoff=cutoff)
        if not hit:
            missing.append(r)
            continue
        g = next(x for x in got if x["title"] == hit[0])
        used.add(g["index"])
        matched.append((r, g))
        if g["pictures"] < r["pictures"]:
            gap.append({"title": r["title"], "reference_pictures": r["pictures"],
                        "built_pictures": g["pictures"]})
    extra = [g for g in got if g["index"] not in used]
    return {"missing": missing, "extra": extra, "matched": matched, "picture_gap": gap}
