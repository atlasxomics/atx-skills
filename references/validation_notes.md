# Validation notes: regenerated figures vs a hand-built deck

Every rule in SKILL.md came from comparing generated output against a deck a bench
scientist had actually delivered. This records *what reproduces and what cannot*, so the
next person knows which differences are bugs and which are expected.

Deliberately no customer data here: no project or run identifiers, no sample conditions,
no measured values. The point is the reproducibility class of each panel, not one
project's numbers.

## Reproducibility classes

| Class | Meaning | What to do |
|---|---|---|
| **A** | The workflow already wrote the figure | Fetch it (`workflow_figure()`), do not redraw |
| **B** | Regenerable from the project tables, matches exactly | Generate and check against the delivered values |
| **C** | Drawn by Latch's proprietary viewer | Data reproduces; cosmetics will not match, and that is fine |
| **D** | Needs an input not stored in the project | Ask for it rather than inventing one |

## What was established panel by panel

- **Peaks per cluster / per sample** reproduces **exactly**, including the `C1(n = ...)`
  tick labels, the category order and the peak-annotation colour mapping. It comes
  straight from the ArchR peak-call summary, so any mismatch means the wrong table.
- **The appendix QC table** reproduces on all but one cell, and the single miss is a
  rounding transcription error *in the hand-built deck*, not in the generator — a second
  project showed the mirror-image miss, which is what proved the column choice correct.
- **Marker gene / peak / motif heatmaps** are class A: ArchR writes `heatmaps_all.pdf`
  and the deck uses its three pages. Do not attempt to redraw them.
- **Gene-activity panels** reproduce on data, with class-C cosmetics. The viewer's
  `Mean:` annotation reproduces here.
- **Spatial QC maps** reproduce the colourbar endpoints exactly but *not* the `Mean:`
  annotation — matching endpoints with a drifting mean points at a different denominator
  inside the viewer. Open, deferred.
- **Differential motifs** are class B once you stop trying to recompute them. The
  statistics come from an ArchR Compare workflow run; fetch its result table and use its
  own per-barcode group assignment. Re-thresholding to recreate the selection is how you
  get numbers that nearly agree.
- **Genome-browser tracks** reproduce when the view window matches. The window is a human
  IGV navigation and is recorded nowhere, so a differing shared maximum usually means a
  differing window, not a differing calculation.
- **Morphology overlays** need the registration transform. Where only the unregistered
  image exists, say so rather than approximating it.

## Two panels that cannot be generated at all

Screenshots of the Latch web UI (the Plots viewer and the workspace file listing). They
are class D by nature: ask for the screenshots, or drop the slides.

## The habit worth copying

Three transcription errors in hand-built decks surfaced only because generated numbers
were compared against them cell by cell. When a generated value disagrees with a
delivered one, the delivered one is not automatically right — check the source table
before "fixing" the generator.
