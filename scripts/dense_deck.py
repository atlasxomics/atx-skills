"""Build a dense Data Report deck from a slide_figures.py folder: every figure placed at
100% in the slot it was drawn for, so labels read at 12-18 pt and slides carry little
white space.

    python dense_deck.py <figure_dir> [--spec deck_spec.json] [--out deck.pptx] [--draft]

Runs locally (python-pptx). `figure_dir` is what slide_figures.py wrote: the PNGs,
samples_/qc_table_ CSVs and slide_figures_report.json. The spec is the same JSON passed to
slide_figures.py; its optional "deck" block controls the slides:

    "deck": {
      "subtitle": "Data Report: <project>",
      "titles": {"clustering": ["Title", "Subtitle"], ...},   # override default titles
      "slides": [ ... ]                                         # explicit plan (else default)
    }

Slide types (each takes optional "title", "subtitle", "notes"):
    samples | clustering | celltype_table | celltype_map | neighborhood | morphology
    heatmaps | rows {"features": [..]} | tracks | trackset {"name"} | peaks | qc | genome
    table {"csv", "columns": {src: shown}, "col_w", "pt", "row_h", "bullets"}
    image {"path"} | section | workflow
`morphology` emits one slide per run ("subtitles": {run_id: text}). `rows` takes 2 or 3
features and picks the matching slot. Paths in the spec are relative to the spec file.

With no "slides", the default plan is the standard report: samples, clustering,
neighborhood, one morphology slide per run, marker heatmaps, gene rows two per slide,
track browser, motif rows three per slide, peaks; then Appendix with QC, the putative
cell-type table and map (if the spec names cell types), appendix rows and the workflow. Default titles are descriptive, not
findings: a scientist writes the finding titles (deck.titles) after reading the figures.

Every build is checked by atx_deck.validate_deck and refused if any finding remains.
Rules and the reasons for them: references/deck.md, "Dense deck".
"""
import argparse, json, os, sys
import pandas as pd
from PIL import Image
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import atx_deck as D  # noqa: E402
from atx_slots import SLOTS, SLIDE_W, X0, Y_SUB, Y_NOSUB, BOTTOM, GAP, CONTENT_W  # noqa: E402

TEMPLATE = os.path.join(HERE, "..", "assets", "atx_template.pptx")
DEFAULT_TITLES = {
    "samples": ("Samples in this report", None),
    "clustering": ("Spatial clustering", None),
    "celltype_table": ("Putative cell types assigned from marker gene activity",
                       "Putative: inferred from gene activity, not validated by IHC or molecular profiling"),
    "celltype_map": ("Putative cell types in tissue",
                     "Putative: not validated by IHC or molecular profiling; colors match the clusters throughout"),
    "neighborhood": ("Neighborhood enrichment and cluster proportions",
                     "Neighborhood enrichment: how often each cluster is adjacent to another (higher is more related)"),
    "morphology": ("Clustering on tissue morphology: {label}", None),
    "heatmaps": ("Differential genes, peaks and motifs between clusters",
                 "ArchR marker features, each cluster vs all others; Pval <= 0.05 & Log2FC >= 0.1"),
    "rows": ("Gene activity of marker genes", "Spatial panels share one color scale per feature"),
    "rows_motif": ("Transcription-factor motif activity",
                   "chromVAR deviation z-scores per spot; a motif is shared within its TF family"),
    "tracks": ("Coverage at marker gene loci",
               "Coverage per cluster on one shared y-axis per gene; boxes mark clusters whose gene "
               "score is at least 1 SD above the cluster mean"),
    "trackset": ("Coverage at {name} loci", "One shared y-axis per gene"),
    "peaks": ("Total number of peaks per cluster and sample", "n = number of spots per cluster or sample"),
    "qc": ("Samples pass QCs for each assay", None),
    "genome": ("Genome-wide gene activity by 5 Mb bin",
               "Median ArchR log2FC per 5 Mb bin, one-vs-rest; relative, not absolute, copy number"),
    "workflow": ("Bioinformatics Workflow", "For more information go to our website"),
}
SAMPLE_COLUMNS = {"run_id": "Run", "sample_name": "Sample", "n_cells": "Spots", "Species": "Species",
                  "Organ": "Organ", "Tissue type": "Tissue type", "Workflow": "Assay"}


class Deck:
    def __init__(self, figdir, spec, spec_dir, template=TEMPLATE):
        self.figdir, self.spec, self.spec_dir = figdir, spec, spec_dir
        self.rep = json.load(open(os.path.join(figdir, "slide_figures_report.json")))
        self.proj = self.rep["project"]
        self.prs = Presentation(template)
        self.layout = D._content_layout(self.prs)
        self.titles = {k: tuple(v) for k, v in spec.get("deck", {}).get("titles", {}).items()}
        self.left_out = []

    # ------------------------------------------------------------ primitives
    def fig(self, stem):
        p = os.path.join(self.figdir, f"{self.proj}__{stem}.png")
        return p if os.path.exists(p) else None

    def need(self, stem, what):
        p = self.fig(stem)
        if p is None:
            self.left_out.append(what)
        return p

    def title_for(self, key, spec, **fmt):
        t, s = self.titles.get(key, DEFAULT_TITLES.get(key, (key, None)))
        t = spec.get("title", t)
        s = spec.get("subtitle", s)
        return t.format(**fmt), (s.format(**fmt) if s else s)

    def slide(self, title, subtitle=None, notes=None):
        s = self.prs.slides.add_slide(self.layout)
        D._inherit_chrome(s, self.layout)
        D._fill_placeholders(s, title, subtitle)
        if notes:
            s.notes_slide.notes_text_frame.text = notes
        return s, (Y_SUB if subtitle else Y_NOSUB)

    def section(self, title):
        lay = D._section_layout(self.prs)
        s = self.prs.slides.add_slide(lay)
        D._inherit_chrome(s, lay)
        D._fill_placeholders(s, title, None)

    @staticmethod
    def pic(s, path, x, y, w=None, h=None):
        iw, ih = Image.open(path).size
        if w is None:
            w = h * iw / ih
        if h is None:
            h = w * ih / iw
        return s.shapes.add_picture(path, Inches(x), Inches(y), Inches(w), Inches(h))

    def fit_pic(self, s, path, x, y, w, h):
        """Largest placement of `path` inside (w, h), centered horizontally. For pictures not
        drawn at slot size (ArchR's PDFs, the workflow diagram, user images)."""
        iw, ih = Image.open(path).size
        sc = min(w / iw, h / ih)
        return self.pic(s, path, x + (w - iw * sc) / 2, y, iw * sc, ih * sc)

    @staticmethod
    def text(s, x, y, w, h, lines, pt=16, color=(0x22, 0x22, 0x22)):
        tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
        tf = tb.text_frame
        tf.word_wrap = True
        for i, line in enumerate(lines if isinstance(lines, list) else [lines]):
            para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
            r = para.add_run()
            r.text = line
            r.font.size = Pt(pt)
            r.font.color.rgb = RGBColor(*color)
        return tb

    @staticmethod
    def col_widths(df, total, pt, pad=0.25):
        """Column widths in inches: every VALUE fits on one line where the slide allows.

        Headers may wrap at word boundaries, so a header only needs its longest word; sizing
        to the full header let "median fragments per spot" take the width that the run IDs
        needed, and the run IDs wrapped. When the values cannot all fit, columns whose
        widest value is short (ids, counts, names) keep their width and the long free-text
        columns share the rest and wrap.
        """
        em = 0.58 * pt / 72
        val = [max([len(str(v)) for v in df[c]] or [0]) * em + pad for c in df.columns]
        head = [max(len(w) for w in str(c).split()) * em * 1.08 + pad for c in df.columns]
        need = [max(v, h) for v, h in zip(val, head)]
        if sum(need) <= total:
            extra = total - sum(need)
            return [n + extra * n / sum(need) for n in need]
        short_lim = 2.0
        fixed = sum(n for n in need if n <= short_lim)
        longs = [n for n in need if n > short_lim]
        rest = max(total - fixed, 1.2 * len(longs))
        return [n if n <= short_lim else rest * n / sum(longs) for n in need]

    def table(self, s, df, y, pt=14, row_h=None, col_w=None, max_bottom=BOTTOM):
        rows, cols = df.shape[0] + 1, df.shape[1]
        cw = col_w or self.col_widths(df, CONTENT_W, pt)
        cw = [c * CONTENT_W / sum(cw) for c in cw]
        if row_h is None:                           # fill the slot but never beyond ~3 lines
            row_h = min(0.80, max(0.36, (max_bottom - y) / rows))
        shp = s.shapes.add_table(rows, cols, Inches(X0), Inches(y), Inches(CONTENT_W), Inches(row_h * rows))
        t = shp.table
        for j, c in enumerate(cw):
            t.columns[j].width = Inches(c)
        for i in range(rows):
            t.rows[i].height = Inches(row_h)
        for j, name in enumerate(df.columns):
            cell = t.cell(0, j)
            cell.text = str(name)
            para = cell.text_frame.paragraphs[0]
            para.font.size = Pt(pt); para.font.bold = True
            para.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
            cell.fill.solid(); cell.fill.fore_color.rgb = D.HEADER_FILL
            cell.vertical_anchor = MSO_ANCHOR.MIDDLE
        for i, (_, r) in enumerate(df.iterrows(), start=1):
            for j, v in enumerate(r):
                cell = t.cell(i, j)
                cell.text = "" if v is None or (isinstance(v, float) and pd.isna(v)) else str(v)
                cell.text_frame.paragraphs[0].font.size = Pt(pt)
                cell.vertical_anchor = MSO_ANCHOR.MIDDLE
                cell.margin_top = cell.margin_bottom = Inches(0.03)
        return y + row_h * rows

    # ------------------------------------------------------------ slide types
    def s_samples(self, sp):
        path = sp.get("csv")
        if path:
            df = pd.read_csv(os.path.join(self.spec_dir, path))
        else:
            st = pd.read_csv(os.path.join(self.figdir, f"samples_{self.proj}.csv"))
            cols = {k: v for k, v in SAMPLE_COLUMNS.items() if k in st.columns}
            df = st[list(cols)].rename(columns=cols)
            if "Spots" in df:
                df["Spots"] = df["Spots"].map("{:,}".format)
        if sp.get("columns"):
            df = df[list(sp["columns"])].rename(columns=sp["columns"])
        t, sub = self.title_for("samples", sp)
        s, y = self.slide(t, sub, sp.get("notes"))
        pt = sp.get("pt", 16)
        bottom = self.table(s, df, y + 0.12, pt=pt, row_h=sp.get("row_h", 0.62), col_w=sp.get("col_w"))
        if sp.get("bullets"):
            self.text(s, X0, bottom + 0.35, CONTENT_W, BOTTOM - bottom - 0.4, sp["bullets"], pt=17)

    def s_table(self, sp):
        df = pd.read_csv(os.path.join(self.spec_dir, sp["csv"]))
        if sp.get("columns"):
            df = df[list(sp["columns"])].rename(columns=sp["columns"])
        s, y = self.slide(sp["title"], sp.get("subtitle"), sp.get("notes"))
        bottom = self.table(s, df, y + 0.12, pt=sp.get("pt", 14), row_h=sp.get("row_h"), col_w=sp.get("col_w"))
        if sp.get("bullets"):
            self.text(s, X0, bottom + 0.3, CONTENT_W, BOTTOM - bottom - 0.35, sp["bullets"], pt=16)

    def s_clustering(self, sp):
        a, b, c = self.need("cluster_spatial", "clustering"), self.fig("umap_cluster"), self.fig("umap_sample")
        if not a:
            return
        t, sub = self.title_for("clustering", sp)
        s, y = self.slide(t, sub, sp.get("notes"))
        w, h = SLOTS["cluster_spatial"]
        avail = BOTTOM - y
        k = min(1.0, avail / h)                    # a subtitle costs 0.4 in; shrink only then
        self.pic(s, a, X0, y, w * k, h * k)
        ux = X0 + w * k + GAP
        uw, uh = SLOTS["umap_cluster"]
        if b and c:
            self.pic(s, b, ux, y, uw * k, uh * k)
            self.pic(s, c, ux, y + (uh + GAP) * k, uw * k, uh * k)
        elif b:
            # Single-sample project: no UMAP by sample, so the cluster UMAP was drawn at
            # the full height of the right column (slot umap_cluster_tall) and goes in at
            # its own size -- stretching it here would distort it.
            tw, th = SLOTS["umap_cluster_tall"]
            self.pic(s, b, ux, y, tw * k, th * k)

    def s_celltype_table(self, sp):
        ct = self.spec.get("cell_types")
        if not ct and not sp.get("csv"):
            return
        if sp.get("csv"):
            return self.s_table({**sp, "title": sp.get("title", self.title_for("celltype_table", sp)[0]),
                                 "subtitle": self.title_for("celltype_table", sp)[1]})
        sizes = self.rep.get("cluster_sizes", {})
        main = self.rep.get("cluster_main_sample", {})
        ev = self.spec.get("cell_type_evidence", {})
        conf = self.spec.get("cell_type_confidence", {})
        top = self.rep.get("cluster_top_genes", {})
        one_sample = len(self.rep.get("samples", [])) < 2
        rows = [{"Cluster": c, "Putative cell type": ct.get(c, ""),
                 "Strongest genes": top.get(c, ""), "Spots": f"{sizes.get(c, 0):,}",
                 # One sample makes this column "D3145 (100%)" on every row.
                 "Main sample": "" if one_sample else main.get(c, ""),
                 "Supporting evidence": ev.get(c, ""),
                 "Confidence": conf.get(c, "")} for c in self.rep["clusters"]]
        df = pd.DataFrame(rows)
        df = df[[c for c in df.columns if df[c].astype(str).str.len().sum()]]
        t, sub = self.title_for("celltype_table", sp)
        s, y = self.slide(t, sub, sp.get("notes"))
        self.table(s, df, y + 0.12, pt=sp.get("pt", 14), col_w=sp.get("col_w"))

    def s_celltype_map(self, sp):
        a, b = self.need("celltype_spatial", "cell-type map"), self.fig("celltype_umap")
        if not a:
            return
        t, sub = self.title_for("celltype_map", sp)
        s, y = self.slide(t, sub, sp.get("notes"))
        w, h = SLOTS["celltype_spatial"]
        self.pic(s, a, X0, y, w, h)
        if b:
            self.pic(s, b, X0 + w + GAP, y, *SLOTS["celltype_umap"])

    def s_neighborhood(self, sp):
        a, b = self.need("neighborhood", "neighborhood"), self.fig("proportion")
        if not a:
            return
        t, sub = self.title_for("neighborhood", sp)
        s, y = self.slide(t, sub, sp.get("notes"))
        nh = BOTTOM - y
        w, h = SLOTS["neighborhood"]
        self.pic(s, a, X0, y, h=nh)
        if b:
            self.pic(s, b, X0 + w * nh / h + GAP, y, h=nh)

    def s_morphology(self, sp):
        subs = sp.get("subtitles", {})
        labels = self.rep["notes"].get("sample_labels", {})
        for run in self.rep["samples"]:
            p = self.need(f"morph_{run}", f"morphology {run}")
            if not p:
                continue
            t, sub = self.title_for("morphology", sp, label=labels.get(run, run))
            sub = subs.get(run, sub)
            has_he = self.rep["notes"].get(f"morph_{run}", {}).get("he")
            if sub is None:
                sub = "H&E is an adjacent section" if has_he else None
            s, y = self.slide(t, sub, sp.get("notes"))
            self.pic(s, p, X0, y, *SLOTS["morph"])

    def s_heatmaps(self, sp):
        hm = [self.fig(f"heatmap_p{i}") for i in (1, 2, 3)]
        hm = [h for h in hm if h]
        if not hm:
            self.left_out.append("marker heatmaps")
            return
        t, sub = self.title_for("heatmaps", sp)
        s, y = self.slide(t, sub, sp.get("notes"))
        cw = (CONTENT_W - GAP * (len(hm) - 1)) / len(hm)
        for i, h in enumerate(hm):
            self.fit_pic(s, h, X0 + i * (cw + GAP), y, cw, BOTTOM - y)

    def s_rows(self, sp):
        feats = sp["features"]
        slot = "row2" if len(feats) <= 2 else "row3"
        paths = [self.fig(f"row_{f}_{slot}") for f in feats]
        missing = [f for f, p in zip(feats, paths) if p is None]
        if missing:
            self.left_out.append(f"rows {missing}")
        paths = [p for p in paths if p]
        if not paths:
            return
        kinds = {self.rep["rows"].get(f, {}).get("kind", "ge") for f in feats}
        key = "rows_motif" if kinds == {"motif"} else "rows"
        t, sub = self.title_for(key, sp)
        s, y = self.slide(t, sub, sp.get("notes"))
        w, h = SLOTS[slot]
        gap = min(max((BOTTOM - y - h * len(paths)) / max(len(paths) - 1, 1), 0.0), 0.15)
        for k, p in enumerate(paths):
            self.pic(s, p, X0, y + k * (h + gap), w, h)

    def s_tracks(self, sp):
        genes = sp.get("genes") or self.rep.get("track_genes", [])
        cols = [self.fig(f"track_{g}") for g in genes[:3]]
        cols = [c for c in cols if c]
        if not cols:
            self.left_out.append("track browser")
            return
        t, sub = self.title_for("tracks", sp)
        s, y = self.slide(t, sub, sp.get("notes"))
        tw, th = SLOTS["track_col"]
        for k, c in enumerate(cols):
            self.pic(s, c, X0 + k * (tw + 0.1), y, tw, th)
        u = self.fig("track_umap")
        uw, uh = SLOTS["track_umap"]
        if u:
            self.pic(s, u, SLIDE_W - X0 - uw, y - 0.05, uw, uh)
        callouts = sp.get("callouts")
        if callouts is None:
            hl = {g: self.rep["notes"].get(f"tracks_{g}", {}).get("highlight", []) for g in genes[:3]}
            callouts = [f"{g}: {', '.join(c)}" for g, c in hl.items() if c]
        self.text(s, SLIDE_W - X0 - uw, y + uh + 0.05, uw, BOTTOM - y - uh - 0.1, callouts, pt=13)

    def s_trackset(self, sp):
        name = sp["name"]
        genes = self.spec.get("track_sets", {}).get(name, [])
        cols = [self.fig(f"trackset_{name}_{g}") for g in genes[:3]]
        cols = [c for c in cols if c]
        if not cols:
            self.left_out.append(f"track set {name}")
            return
        t, sub = self.title_for("trackset", sp, name=name)
        s, y = self.slide(t, sub, sp.get("notes"))
        cw, ch = SLOTS["track_col_wide"]
        step = (CONTENT_W - cw) / max(len(cols) - 1, 1) if len(cols) > 1 else 0
        for k, c in enumerate(cols):
            self.pic(s, c, X0 + k * step, y, cw, ch)

    def s_peaks(self, sp):
        a, b = self.need("peaks_cluster", "peaks"), self.fig("peaks_sample")
        if not a:
            return
        t, sub = self.title_for("peaks", sp)
        s, y = self.slide(t, sub, sp.get("notes"))
        w, h = SLOTS["peaks"]
        self.pic(s, a, X0, y, w, h)
        if b:
            self.pic(s, b, X0 + w + GAP, y, w, h)

    def s_qc(self, sp):
        qc = pd.read_csv(os.path.join(self.figdir, f"qc_table_{self.proj}.csv"))
        if "median fragments per spot" in qc:
            qc["median fragments per spot"] = qc["median fragments per spot"].map("{:,}".format)
        t, sub = self.title_for("qc", sp)
        notes = sp.get("notes", "The pipeline's obs column 'log10_nFrags' holds the natural log of nFrags "
                                "(archr_objs.R: log(metadata$nFrags)); the map is labeled ln(nFrags).")
        s, y = self.slide(t, sub, notes)
        self.table(s, qc, y + 0.05, pt=14, row_h=0.40)
        qw, qh = SLOTS["qc_map"]
        for k, stem in enumerate(("qc_TSS", "qc_frags")):
            p = self.fig(stem)
            if p:
                self.pic(s, p, X0 + k * (qw + GAP), BOTTOM - qh, qw, qh)

    def s_genome(self, sp):
        p = self.need("genome", "genome heatmap")
        if not p:
            return
        t, sub = self.title_for("genome", sp)
        s, y = self.slide(t, sub, sp.get("notes"))
        self.pic(s, p, X0, y, *SLOTS["genome"])

    def s_image(self, sp):
        s, y = self.slide(sp["title"], sp.get("subtitle"), sp.get("notes"))
        self.fit_pic(s, os.path.join(self.spec_dir, sp["path"]), X0, y, CONTENT_W, BOTTOM - y)

    def s_workflow(self, sp):
        p = self.fig("bioinformatics_workflow")
        if not p:
            return
        t, sub = self.title_for("workflow", sp)
        notes = sp.get("notes", "Not generated (Latch UI screenshots, supply from the workspace): "
                                "'Visualize and explore with AtlasXPlots' and 'Additional Data Outputs "
                                "available on LatchBio Workspace'.")
        s, y = self.slide(t, sub, notes)
        self.fit_pic(s, p, X0, y, CONTENT_W, BOTTOM - y)

    # ------------------------------------------------------------ plan
    def default_plan(self):
        rows = self.rep.get("rows", {})
        main = [f for f, r in rows.items() if not r.get("appendix")]
        app = [f for f, r in rows.items() if r.get("appendix")]
        ge = [f for f in main if rows[f]["kind"] == "ge"]
        mo = [f for f in main if rows[f]["kind"] != "ge"]
        plan = [{"type": "samples"}, {"type": "clustering"}]
        plan += [{"type": "neighborhood"}, {"type": "morphology"}, {"type": "heatmaps"}]
        plan += [{"type": "rows", "features": ge[i:i + 2]} for i in range(0, len(ge), 2)]
        if self.rep.get("track_genes"):
            plan.append({"type": "tracks"})
        plan += [{"type": "rows", "features": mo[i:i + 3]} for i in range(0, len(mo), 3)]
        plan += [{"type": "peaks"}, {"type": "section", "title": "Appendix"}, {"type": "qc"}]
        if self.spec.get("cell_types"):
            # Putative calls are supporting material, not the report's headline.
            plan += [{"type": "celltype_table"}, {"type": "celltype_map"}]
        plan += [{"type": "rows", "features": app[i:i + 2]} for i in range(0, len(app), 2)]
        if self.spec.get("genome_heatmap"):
            plan.append({"type": "genome"})
        for name in self.spec.get("track_sets", {}):
            plan.append({"type": "trackset", "name": name})
        plan.append({"type": "workflow"})
        return plan

    def build(self, out, draft=False):
        dk = self.spec.get("deck", {})
        s = self.prs.slides.add_slide(D._title_layout(self.prs))
        ph = next(iter(s.placeholders))
        ph.text_frame.text = dk.get("title", "Chart the unexplored depths of spatial biology")
        ph.text_frame.add_paragraph().text = (dk.get("subtitle", f"Data Report: {self.proj}")
                                              + (" (draft)" if draft else ""))
        for sp in dk.get("slides") or self.default_plan():
            kind = sp["type"]
            if kind == "section":
                self.section(sp["title"])
            else:
                getattr(self, f"s_{kind}")(sp)
        if draft:
            rows = list(self.rep.get("assumptions", []))
            if self.left_out:
                rows.append("Left out (failed, skipped or not applicable): " + ", ".join(self.left_out) + ".")
            rows.append("Not reviewed by a scientist. Check every slide before it goes to the customer, "
                        "and remove this slide.")
            sl, y = self.slide("Assumptions in this draft")
            self.table(sl, pd.DataFrame({"Assumption": rows}), y + 0.05, pt=14, row_h=0.5)
        self.prs.save(out)
        D.assert_deck_valid(out)
        return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("figure_dir")
    ap.add_argument("--spec")
    ap.add_argument("--out")
    ap.add_argument("--draft", action="store_true")
    ap.add_argument("--template", default=TEMPLATE)
    a = ap.parse_args()
    spec = json.load(open(a.spec)) if a.spec else {}
    spec_dir = os.path.dirname(os.path.abspath(a.spec)) if a.spec else "."
    d = Deck(a.figure_dir, spec, spec_dir, a.template)
    out = a.out or os.path.join(a.figure_dir, f"{d.proj}_Data_Report{'_draft' if a.draft else ''}.pptx")
    d.build(out, draft=a.draft)
    print(f"wrote {out}: {len(d.prs.slides)} slides" + (f"; left out: {d.left_out}" if d.left_out else ""))


if __name__ == "__main__":
    main()
