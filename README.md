# atx_skills — AtlasXomics report figure generation

A Claude Science skill that regenerates any figure or table from an AtlasXomics
spatial-epigenomics Data Report directly from a Latch project, and optionally assembles
them into a branded PowerPoint deck.

## Layout

```
SKILL.md                      the skill: role, the conversation, the four parts of a report
scripts/atx_project.py        Latch project -> table-first data access; sample table and
                              slide metadata; ranked marker genes
scripts/atx_plots_lib.py      GUI-free plotting functions lifted from atlasxomics/plots_epigenome
scripts/atx_figures.py        one generator per report figure
scripts/render_figures.py     renders every figure a report needs, in one pass on the box
scripts/atx_compare.py        differential analysis: fetch, check and plot compare clusters runs;
                              launch a new one after the scientist confirms
scripts/atx_deck.py           PowerPoint assembly (runs locally, not on the box)
scripts/atx_slots.py          slide slot geometry shared by the two dense-deck scripts
scripts/slide_figures.py      every deck figure drawn at its slot size (runs on the box)
scripts/dense_deck.py         the default deck: figures placed at 100%, driven by a JSON spec
scripts/build_report_deck.py  the standard deck for any project, from a render_figures output
                              folder (--draft adds the one-shot assumptions slide)
scripts/build_demo_deck.py    the Demo_001 deck, as a worked example of atx_deck
scripts/validate_selftest.py  checks that every deck validation rule fires on a planted fault
scripts/cleanup.py            end-of-session cleanup of the box's scratch and cache
references/figure_catalog.md  provenance of every figure, measured geometry, QC column rules
references/figures.md         per-figure detail: the four checks, marker genes, track browser,
                              differential motif slide, motif logos
references/differential_analysis.md   how compare clusters computes differential genes/peaks/motifs
references/pipeline.md        the ATX processing pipeline (epigenomics in depth; WT and co-profiling in brief)
references/deck.md            deck assembly: template, slide structures, validation, layout traps
references/validation_notes.md   regenerated vs delivered: what reproduces and what cannot
assets/atx_template.pptx      branded template: a delivered deck with its slides stripped
assets/bioinformatics_workflow.png   the workflow diagram slide
```

## Prerequisites

The only hard requirement is the **`ssh:analysis-box` compute target**. Large downloads and
all heavy computation run there: the box has a Latch token, and the objects are too large to
move (`combined_ge.h5ad` is ~18 GB). Nothing needs to be installed on a scientist's machine.

If a scientist's machine happens to have the `latch` CLI installed and signed in, the skill
also uses it for searching and reading small tables. It never asks a scientist to set that
up; without it, everything Latch-related runs on the box.

The box needs, and currently has:

- `/opt/analysis` venv with scanpy, squidpy, anndata, plotly, kaleido, snapatac2, python-pptx
- the 12 system libraries headless Chrome needs for static image export
- `~/.latch/token` and `~/.latch/workspace` set to 13502

Those live on the instance's root volume and will **not** survive a rebuild. Bake them
into the AMI if the instance is ever recreated.

## Provenance

`atx_plots_lib.py` is extracted verbatim from `atlasxomics/plots_epigenome`
`welcome/init.py` at commit `3d05a74f7683` (main as of 2026-09-01), with three documented
deviations listed in the module docstring. Sibling pins at extraction time were
`plots_wt c04d36d51a03` and `plots_copro eabed425ec64`. If Plots moves, re-extract and
re-run the smoke test rather than hand-patching.

## What is deliberately not here

- **Example figures.** The rendered PNGs run 1-2 MB each; they are artifacts of a run,
  not source. Add an `examples/` directory if you want them as documentation.
- **Private workflow repos.** `atlasXomics_pipeline_wf`, `rnaqc_wf`, `barcodeqc_wf`,
  `atx_glue`, `cram2frags`, `atx_convert` and `archr_genes` are not publicly readable;
  nothing here depends on them. `atx_glue` becomes relevant only for co-profiling panels.

`assets/atx_template.pptx` *is* committed despite being a binary, because branded output
is unreachable without it and it is only ~1 MB. Regenerate it by stripping the slides
from any delivered deck.
