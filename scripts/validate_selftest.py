"""Self-test for atx_deck.validate_deck: every check must fire on a planted fault.

A validation pass nobody has seen fail is not evidence of anything. Each case below
builds one deliberately broken slide and asserts the matching code is reported, so a
future refactor that quietly stops checking something is caught here.

Run:  python validate_selftest.py <template.pptx>
"""
import sys

from pptx import Presentation
from pptx.util import Emu, Inches

import atx_deck as D

TEMPLATE = sys.argv[1] if len(sys.argv) > 1 else "atx_template.pptx"
IMG = sys.argv[2] if len(sys.argv) > 2 else None


def _fresh():
    prs = Presentation(TEMPLATE)
    for s in list(prs.slides._sldIdLst):
        prs.slides._sldIdLst.remove(s)
    return prs


def _slide(prs, title="Check", subtitle=None, chrome=True):
    lay = D._content_layout(prs)
    s = prs.slides.add_slide(lay)
    if chrome:
        D._inherit_chrome(s, lay)
    D._fill_placeholders(s, title, subtitle)
    return s


def codes(prs, path="_selftest.pptx"):
    prs.save(path)
    return {f["code"] for f in D.validate_deck(path)}


def expect(name, found, want):
    ok = want in found
    print(f"  {'ok  ' if ok else 'FAIL'} {name:22} -> {sorted(found) or 'clean'}")
    return ok


def main():
    results = []

    prs = _fresh(); s = _slide(prs, "Off without ext", "sub")
    ph = [p for p in s.placeholders if p.placeholder_format.idx == 13][0]
    ph.top = Emu(int(Inches(2.0)))                      # offset only, no extent
    ph._element.spPr.xfrm.remove(ph._element.spPr.xfrm.find(
        "{http://schemas.openxmlformats.org/drawingml/2006/main}ext"))
    results.append(expect("OFF_WITHOUT_EXT", codes(prs), "OFF_WITHOUT_EXT"))

    prs = _fresh(); s = _slide(prs, "Off slide")
    s.shapes.add_textbox(Emu(int(Inches(12.0))), Emu(int(Inches(1.0))),
                         Emu(int(Inches(3.0))), Emu(int(Inches(1.0)))).text_frame.text = "x"
    results.append(expect("OFF_SLIDE", codes(prs), "OFF_SLIDE"))

    prs = _fresh(); s = _slide(prs, "Crosses footer")
    s.shapes.add_textbox(Emu(int(Inches(1.0))), Emu(int(Inches(6.5))),
                         Emu(int(Inches(3.0))), Emu(int(Inches(0.8)))).text_frame.text = "x"
    results.append(expect("CROSSES_FOOTER", codes(prs), "CROSSES_FOOTER"))

    prs = _fresh(); s = _slide(prs, "Overlap")
    for dx in (0.0, 0.1):
        s.shapes.add_textbox(Emu(int(Inches(1.0 + dx))), Emu(int(Inches(2.0))),
                             Emu(int(Inches(3.0))), Emu(int(Inches(1.0)))
                             ).text_frame.text = "x"
    results.append(expect("SHAPE_OVERLAP", codes(prs), "SHAPE_OVERLAP"))

    prs = _fresh()
    lay = D._content_layout(prs)
    s = prs.slides.add_slide(lay)                        # no _inherit_chrome
    D._fill_placeholders(s, "No chrome", None)
    results.append(expect("MISSING_SLIDE_NUMBER", codes(prs), "MISSING_SLIDE_NUMBER"))

    prs = _fresh()
    lay = D._content_layout(prs)
    s = prs.slides.add_slide(lay)
    D._inherit_chrome(s, lay)                            # placeholders left empty
    results.append(expect("EMPTY_PLACEHOLDER", codes(prs), "EMPTY_PLACEHOLDER"))

    prs = _fresh(); s = _slide(prs, "Autofit", "sub")
    t_ph = [p for p in s.placeholders if p.placeholder_format.idx == 0][0]
    el = t_ph._element.find(
        ".//{http://schemas.openxmlformats.org/drawingml/2006/main}normAutofit")
    if el is not None:
        el.getparent().remove(el)                      # strip the shrink-to-fit flag
    results.append(expect("MISSING_AUTOFIT", codes(prs), "MISSING_AUTOFIT"))

    if IMG:
        prs = _fresh(); s = _slide(prs, "Aspect")
        s.shapes.add_picture(IMG, Emu(int(Inches(1.0))), Emu(int(Inches(2.0))),
                             width=Emu(int(Inches(4.0))), height=Emu(int(Inches(4.0))))
        results.append(expect("ASPECT_DISTORTED", codes(prs), "ASPECT_DISTORTED"))
    else:
        print("  skip ASPECT_DISTORTED        (no image path given)")

    print(f"\n{sum(results)}/{len(results)} checks fired as expected")
    return 0 if all(results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
