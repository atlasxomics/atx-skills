# atx_skills — AtlasXomics report figure generation

A Claude Science skill that regenerates any figure or table from an AtlasXomics
spatial-epigenomics Data Report directly from a Latch project, and optionally assembles
them into a branded PowerPoint deck.

## Layout

```
SKILL.md                      the skill itself (instructions + triggering description)
scripts/atx_plots_lib.py      GUI-free plotting functions lifted from atlasxomics/plots_epigenome
scripts/atx_project.py        Latch project -> table-first data access
scripts/atx_figures.py        one generator per report figure
scripts/atx_deck.py           PowerPoint assembly (runs locally, not on the box)
references/figure_catalog.md  provenance of every figure, measured geometry, QC column rules
references/pipeline.md        the ATX processing pipeline (epigenomics in depth; WT and co-profiling in brief)
references/deck.md            deck assembly: template, slide structures, validation, layout traps
references/comparison_Demo_001.md   regenerated vs delivered, panel by panel
assets/atx_template.pptx      branded template: a delivered deck with its slides stripped
```

## Prerequisites

The only hard requirement is the **`ssh:analysis-box` compute target**. Every Latch
operation runs there: the token lives on that box, latch.bio is unreachable from a local
sandbox, and the objects are too large to move (`combined_ge.h5ad` is ~18 GB). Nothing
installs on a scientist's machine.

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
