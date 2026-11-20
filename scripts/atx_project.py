"""atx_project — resolve an AtlasXomics Latch project into the inputs the report
figures need, table-first.

Design rule (from Demo_001: combined_ge.h5ad is 19 GB against 64 GB RAM):
anything derivable from tables/ is read from tables/. The feature matrices are opened
backed, and only for gene/motif feature plots.
"""
from __future__ import annotations

import os
import shlex
import subprocess
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional

import numpy as np
import pandas as pd

LATCH = os.environ.get("LATCH_BIN", "/opt/analysis/bin/latch")

#: Morphology image folders, searched in order; the first that has <short>_fix.* wins.
IMAGE_ROOTS = ("latch://13502.account/spatials", "latch://atx-illumina.mount/Images")

#: Genome build -> species, for projects whose runs have no slide_metadata.json.
GENOME_SPECIES = {"hg38": "Homo sapiens", "hg19": "Homo sapiens",
                  "mm10": "Mus musculus", "mm39": "Mus musculus",
                  "rn6": "Rattus norvegicus", "rn7": "Rattus norvegicus",
                  "rn8": "Rattus norvegicus"}

#: slide_metadata.json fields worth showing, by their SLIMS display title. The condition
#: is usually in "Block experimental condition"; the slide-level field is often empty.
SLIDE_FIELDS = ("Sample Id", "Block name", "Species", "Organ", "Tissue type", "Disease",
                "Tissue slide experimental condition", "Block experimental condition",
                "Workflow", "CT Primary Antibody", "CT Secondary Antibody", "Source",
                "Comment")


def short_run_id(run_id):
    """RUN1_SEQ001 -> RUN1. Image folders drop the leading zeros (D02778 -> D2778)."""
    import re
    head = str(run_id).split("_")[0]
    m = re.match(r"^([A-Za-z]+)0*(\d+)$", head)
    return f"{m.group(1)}{int(m.group(2))}" if m else head


def parse_slide_metadata(doc: dict) -> dict:
    """Flatten a fastq2frags `run_metrics/slide_metadata.json` into {title: value}.

    The file holds one SLIMS record per object linked to the run (lysate, NGS libraries,
    sublibrary, tissue slide). The "Tissue slide" record is the fullest (it alone carries
    Disease and the CUT&Tag antibodies), so it wins; the others only fill gaps. A field
    that disagrees between records keeps the tissue-slide value."""
    def value(col):
        for k in ("displayValue", "joinedDisplayValue", "displayValues", "value"):
            v = col.get(k)
            if v not in (None, "", []):
                return "; ".join(map(str, v)) if isinstance(v, list) else v
        return None

    records = []
    for rec in doc.get("slides", []):
        flat = {}
        for col in rec.get("columns", []):
            title, v = col.get("title") or col.get("name"), value(col)
            if title and v is not None and title not in flat:
                flat[title] = v
        records.append(flat)
    records.sort(key=lambda r: r.get("Type") != "Tissue slide")
    out = {"Run Id": doc.get("run_id")}
    for r in records:
        for k, v in r.items():
            out.setdefault(k, v)
    slide = next((r for r in records if r.get("Type") == "Tissue slide"), None)
    out["Tissue slide Id"] = slide.get("Id") if slide else None
    return out

_LS_QUERY = """
query LDataChildren($argPath: String!) {
  ldataResolvePathData(argPath: $argPath) {
    finalLinkTarget {
      childLdataTreeEdges(filter: {child: {removed: {equalTo: false},
          pending: {equalTo: false}, copiedFrom: {isNull: true}}}) {
        nodes { child { name type } }
      }
    }
  }
}"""


def latch_ls(path: str, timeout: int = 300) -> List[str]:
    """Names in a Latch directory, directories with a trailing '/'.

    `latch ls` cuts names longer than 50 characters to 47 + '...' (Plots' compare_config
    files, morphology TIFs), so the names come from the Latch API when the latch Python
    package is importable, and from the CLI only as a fallback, which refuses to return
    a cut name rather than hand back one that matches nothing."""
    try:
        import gql
        import latch.ldata.path  # noqa: F401  (import first: avoids a latch_cli import cycle)
        from latch_cli.utils.path import normalize_path
        from latch_sdk_gql.execute import execute
    except ImportError:
        execute = None
    if execute is not None:
        data = execute(gql.gql(_LS_QUERY), {"argPath": normalize_path(path)})["ldataResolvePathData"]
        if data is None:
            raise FileNotFoundError(f"no such Latch directory: {path}")
        nodes = data["finalLinkTarget"]["childLdataTreeEdges"]["nodes"]
        return sorted(n["child"]["name"] + ("" if str(n["child"]["type"]).lower() == "obj" else "/")
                      for n in nodes)
    cmd = f"{shlex.quote(LATCH)} ls {shlex.quote(path)}"
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    if p.returncode != 0:
        if "no such" in (p.stdout + p.stderr).lower():
            raise FileNotFoundError(f"no such Latch directory: {path}")
        raise RuntimeError(f"latch ls failed for {path}\n{p.stderr[-800:]}")
    names = []
    for line in p.stdout.splitlines():
        line = line.rstrip()
        if not line or line.lstrip().startswith(("Size", "WARN", "Please")):
            continue
        name = line.split(None, 2)[-1] if line.lstrip().startswith("-") else line.split(None, 4)[-1]
        if name.rstrip("/").endswith("..."):
            raise RuntimeError(f"`latch ls` truncated a name in {path} ({name}); run where "
                               "the latch Python package is importable (analysis-box venv)")
        names.append(name)
    return names

ANALYSIS_DIRS = {"archr": "epi_analysis_archr", "snap": "epi_analysis_snap"}

# Small tables that are safe to pull wholesale (all well under 100 MB).
CELL_TABLES = ("obs.csv", "spatial.csv", "UMAPHarmony.csv")
SUMMARY_TABLES = ("medians.csv", "input_parameters.csv",
                  "Clusters_Peak-Call-Summary.csv", "Sample_Peak-Call-Summary.csv")


@dataclass
class Project:
    """A located analysis project plus its local cache directory."""
    name: str
    backend: str = "archr"
    cache: str = "/data/analysis_cache"
    _pulled: set = field(default_factory=set, repr=False)

    @property
    def latch_root(self) -> str:
        return f"latch:///{ANALYSIS_DIRS[self.backend]}/{self.name}"

    @property
    def local_root(self) -> str:
        return os.path.join(self.cache, self.name)

    # ---------- transfer ----------

    def _cp(self, remote: str, local: str, timeout: int = 1800) -> str:
        """latch cp one object; skip when already cached. Returns the local path."""
        if os.path.exists(local) and os.path.getsize(local) > 0:
            return local
        os.makedirs(os.path.dirname(local), exist_ok=True)
        cmd = f"{shlex.quote(LATCH)} cp {shlex.quote(remote)} {shlex.quote(local)}"
        p = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
        if not os.path.exists(local):
            raise RuntimeError(f"latch cp failed for {remote}\n{p.stdout[-800:]}\n{p.stderr[-800:]}")
        return local

    def table(self, name: str) -> str:
        """Local path to tables/<name>, fetching it on first use."""
        return self._cp(f"{self.latch_root}/tables/{name}",
                        os.path.join(self.local_root, "tables", name))

    def fetch_tables(self, names: Iterable[str]) -> Dict[str, str]:
        return {n: self.table(n) for n in names}

    # ---------- cell-level table-first loading ----------

    def cell_table(self) -> pd.DataFrame:
        """obs.csv joined to spatial.csv and UMAPHarmony.csv, ALIGNED ON THE BARCODE INDEX.

        The three files are not written in the same row order, so they must be joined on
        the index (`<run_id>#<barcode>-1`), never concatenated positionally.
        """
        obs = pd.read_csv(self.table("obs.csv"), index_col=0)
        spa = pd.read_csv(self.table("spatial.csv"), index_col=0)
        ump = pd.read_csv(self.table("UMAPHarmony.csv"), index_col=0)
        if not (len(obs) == len(spa) == len(ump)):
            raise ValueError(f"row counts differ: obs={len(obs)} spatial={len(spa)} umap={len(ump)}")
        spa = spa.reindex(obs.index)
        ump = ump.reindex(obs.index)
        missing = int(spa.iloc[:, 0].isna().sum() + ump.iloc[:, 0].isna().sum())
        if missing:
            raise ValueError(f"{missing} barcodes in obs.csv absent from spatial/UMAP tables")
        df = obs.copy()
        df[["spatial_1", "spatial_2"]] = spa.values
        df[["umap_1", "umap_2"]] = ump.values
        return df

    def cell_adata(self):
        """AnnData carrying obs + obsm only (no feature matrix).

        Sufficient for every cluster-, sample- and QC-driven figure: spatial cluster maps,
        proportions, neighborhood enrichment, spatial QC maps.
        """
        import anndata as ad
        df = self.cell_table()
        obsm = {
            "spatial": df[["spatial_1", "spatial_2"]].to_numpy(float),
            "X_umap": df[["umap_1", "umap_2"]].to_numpy(float),
        }
        obs = df.drop(columns=["spatial_1", "spatial_2", "umap_1", "umap_2"])
        for col in ("cluster", "sample", "condition"):
            if col in obs:
                obs[col] = obs[col].astype("category")
        a = ad.AnnData(obs=obs)
        for k, v in obsm.items():
            a.obsm[k] = v
        a.obsm["spatial_offset"] = obsm["spatial"].copy()
        return a

    def feature_adata(self, kind: str = "ge", backed: bool = True):
        """Open combined_ge.h5ad / combined_motifs.h5ad backed.

        Never opens the reduced *_sm_* objects: for ATAC their .X is float16, so any value
        read from them is imprecise.
        """
        fname = {"ge": "combined_ge.h5ad", "motifs": "combined_motifs.h5ad"}[kind]
        local = self._cp(f"{self.latch_root}/anndata/{fname}",
                         os.path.join(self.local_root, "anndata", fname))
        import anndata as ad
        return ad.read_h5ad(local, backed="r" if backed else None)

    # ---------- report tables ----------

    # ---------- samples and metadata ----------

    def slide_metadata(self, run_id: str) -> Optional[dict]:
        """`fastq2frags/<run_id>/run_metrics/slide_metadata.json`, flattened, or None.

        fastq2frags writes it only when SLIMS has records for the run, and only since
        2026-11-17, so older runs have none."""
        folder = f"latch:///fastq2frags/{run_id}/run_metrics"
        try:
            if "slide_metadata.json" not in latch_ls(folder):
                return None
        except FileNotFoundError:
            return None
        import json
        local = self._cp(f"{folder}/slide_metadata.json",
                         os.path.join(self.local_root, "metadata", f"{run_id}.json"))
        with open(local) as fh:
            return parse_slide_metadata(json.load(fh))

    def image_folder(self, run_id: str) -> Dict[str, Optional[str]]:
        """Where a run's morphology images are: the first IMAGE_ROOTS folder that exists,
        and its registered `_fix` image if there is one. H&E images are not tracked."""
        short = short_run_id(run_id)
        for root in IMAGE_ROOTS:
            try:
                names = latch_ls(f"{root}/{short}")
            except FileNotFoundError:
                continue
            fix = sorted(n for n in names if "_fix" in n and not n.endswith("/"))
            return {"image_folder": f"{root}/{short}",
                    "fix_image": f"{root}/{short}/{fix[0]}" if fix else None}
        return {"image_folder": None, "fix_image": None}

    def sample_table(self, images: bool = True) -> pd.DataFrame:
        """One row per run: cell count, the project's sample/condition labels, species,
        slide metadata, and image locations. The first thing to show a scientist.

        `metadata_source` says where the descriptive fields came from:
        "slide_metadata.json" (SLIMS, newer runs) or "project tables" (obs.csv and the
        genome in input_parameters.csv). Missing values are left empty for the scientist
        to fill, never guessed."""
        obs = pd.read_csv(self.table("obs.csv"), index_col=0)
        cond_cols = [c for c in obs.columns if c == "condition" or c.startswith("condition_")]
        genome = None
        try:
            ip = pd.read_csv(self.table("input_parameters.csv"))
            genome = str(ip["genome"].iloc[0]) if "genome" in ip else None
        except Exception:
            pass
        rows = []
        for run, sub in obs.groupby("sample", sort=True):
            row = {"run_id": run, "n_cells": len(sub)}
            if "SampleName" in sub:
                row["sample_name"] = "; ".join(map(str, sub["SampleName"].unique()))
            for c in cond_cols:
                row[c] = "; ".join(map(str, sub[c].unique()))
            row["genome"] = genome
            meta = self.slide_metadata(run)
            if meta:
                row["metadata_source"] = "slide_metadata.json"
                row["Tissue slide Id"] = meta.get("Tissue slide Id")
                for f in SLIDE_FIELDS:
                    row[f] = meta.get(f)
            else:
                row["metadata_source"] = "project tables"
                row["Species"] = GENOME_SPECIES.get(genome or "")
            if images:
                row.update(self.image_folder(run))
            rows.append(row)
        df = pd.DataFrame(rows)
        return df.dropna(axis=1, how="all")

    def ranked_genes(self, level: str = "cluster") -> pd.DataFrame:
        """`tables/ranked_genes_per_<level>.csv`: every gene, each group against all other
        cells (ArchR getMarkerFeatures on the GeneScoreMatrix). Columns: group_name, name,
        Log2FC, MeanDiff, FDR, Pval. For `condition`, the file may be named
        `_condition_1` or `_conditions1` (see Plots' resolve_heatmap_stats_table)."""
        names = {"cluster": ["ranked_genes_per_cluster.csv"],
                 "sample": ["ranked_genes_per_sample.csv"],
                 "condition": ["ranked_genes_per_condition.csv",
                               "ranked_genes_per_condition_1.csv",
                               "ranked_genes_per_conditions1.csv"]}[level]
        present = set(latch_ls(f"{self.latch_root}/tables"))
        for n in names:
            if n in present:
                return pd.read_csv(self.table(n))
        raise FileNotFoundError(f"no ranked_genes table for {level} (tried {', '.join(names)})")

    def top_ranked_genes(self, level: str = "cluster", n: int = 10, fdr: float = 0.05,
                         sort_by: str = "FDR", groups: Optional[List[str]] = None):
        """The top `n` genes higher in each group than in all other cells: FDR <= `fdr`
        and Log2FC > 0, ranked by `sort_by` (FDR ascending, or Log2FC / MeanDiff
        descending).

        FDR is the default because it gave the most interpretable lists on 'a test'
        (C7: Prkch, Lef1, Bcl11b, Tcf7, a T-cell cluster). Ranking by Log2FC puts genes
        with almost no gene score (Gm*, Mir*; MeanDiff 0.03-0.06) on top, and ranking by
        MeanDiff favours long genes (Csmd1, Dmd, Lrp1b) that collect gene score by size.
        Treat any list as candidates for the scientist, not as cell-type calls.

        Returns (table, rule). Quote the rule with any list you show: ArchR
        getMarkerFeatures on the GeneScoreMatrix, t-test, bias TSSEnrichment +
        log10(nFrags), each group against all other cells (archrproject/wf/archr_objs.R).
        That is a one-vs-rest marker ranking, not a pairwise comparison; for A vs B use
        compare clusters (atx_compare)."""
        if sort_by not in ("FDR", "Log2FC", "MeanDiff"):
            raise ValueError("sort_by must be FDR, Log2FC or MeanDiff")
        df = self.ranked_genes(level)
        df["group_name"] = df["group_name"].astype(str)
        if groups is not None:
            df = df[df["group_name"].isin([str(g) for g in groups])]
        sig = df[(df["FDR"] <= fdr) & (df["Log2FC"] > 0)]
        top = (sig.sort_values(["group_name", sort_by], ascending=[True, sort_by == "FDR"])
                  .groupby("group_name", sort=False).head(n)
                  [["group_name", "name", "Log2FC", "MeanDiff", "FDR"]]
                  .reset_index(drop=True))
        rule = (f"ranked_genes_per_{level}: ArchR getMarkerFeatures, GeneScoreMatrix, t-test, "
                f"bias TSSEnrichment + log10(nFrags), each {level} vs all other cells; "
                f"FDR <= {fdr}, Log2FC > 0, top {n} by {sort_by}")
        return top, rule

    def peak_call_summary(self, level: str = "cluster") -> pd.DataFrame:
        """Tidy the ArchR peak-call summary behind the 'Peaks per Cluster/Sample' charts.

        `Group` already carries the `C1(n = ...)` label and `Freq` is already in thousands.
        """
        name = {"cluster": "Clusters_Peak-Call-Summary.csv",
                "sample": "Sample_Peak-Call-Summary.csv"}[level]
        df = pd.read_csv(self.table(name))
        df = df.rename(columns={"Var1": "peak_anno", "Freq": "peaks_thousands"})
        df["is_union"] = df["Group"] == "UnionPeaks"
        df["n_cells"] = (df["Group"].str.extract(r"\(n = (\d+)\)")[0]
                         .astype("Int64"))
        df["label"] = df["Group"].str.replace(r"\(n = \d+\)", "", regex=True).str.strip()
        return df

    def qc_table(self, runs: Optional[List[str]] = None) -> pd.DataFrame:
        """The appendix 'Samples pass QCs for each assay' table.

        Joins the analysis project's medians.csv to each run's preprocessing summary.csv.
        Column provenance and formatting were verified against the delivered Demo_001
        and Demo_002 decks:

          TSS enrichment          medians.csv:median_TSS                 -> 2 dp
          median fragments/spot   medians.csv:median_fragments           -> nearest 100
          FRiP                    medians.csv:median_FRIP                -> 2 dp
          Number of reads         summary.csv:'Number of input read pairs'
                                    -> truncated to whole millions, or 1 dp in billions
          Fraction Duplicates     summary.csv:'Median Duplication Rate'  -> 2 dp
          Sample                  medians.csv:condition

        Caveat on duplicates: 'Median Duplication Rate' reproduces 5 of the 7 delivered rows
        exactly; two decks each showed one cell where the hand-typed value disagreed
        differ by one hundredth in opposite directions, which is consistent with manual
        transcription rather than a different source column.
        """
        med = pd.read_csv(self.table("medians.csv"))
        if runs is None:
            runs = med["run_id"].tolist()
        rows = []
        for rid in runs:
            local = self._cp(f"latch:///fastq2frags/{rid}/run_metrics/summary.csv",
                             os.path.join(self.local_root, "summaries", f"{rid}.csv"))
            s = pd.read_csv(local).iloc[0]
            m = med.loc[med["run_id"] == rid].iloc[0]
            rows.append({
                "Run ID": rid,
                "Sample": m["condition"],
                "TSS enrichment": round(float(m["median_TSS"]), 2),
                "median fragments per spot": int(round(float(m["median_fragments"]) / 100.0) * 100),
                "FRiP": round(float(m["median_FRIP"]), 2),
                "Number of reads": format_reads(float(s["Number of input read pairs"])),
                "Fraction Duplicates": round(float(s["Median Duplication Rate"]), 2),
            })
        return pd.DataFrame(rows)


def format_reads(n: float) -> str:
    """Match the decks: whole millions (truncated) below 1e9, else billions to 1 dp."""
    if n >= 1e9:
        return f"{n / 1e9:.1f} billion"
    return f"{int(n // 1e6)} million"
