"""atx_compare — read, check and plot AtlasXomics differential analysis (compare clusters).

Differential statistics on a Data Report come from the **compare clusters** Latch workflow
(repo `atlasxomics/combined_cluster_wf`, script `wf/compare_clusters.R`), which runs
`ArchR::getMarkerFeatures` on an ArchRProject between two cell groups, GroupA and GroupB.
It is normally launched from the Compare tab of the Epigenomics Plots app
(`plots_epigenome/compare/`). This module does not recompute any statistic: it fetches
the workflow's tables, checks them, and draws the same views Plots draws. When no run
exists it can prepare one and, only after the scientist confirms, launch it
(`describe_launch`, `check_launch`, `launch_compare`, `wait_for_run`).

See references/differential_analysis.md for the method and the traps. In short:

* GroupA rows and GroupB rows in `all_*.csv` are two SEPARATE tests (each group against
  the other, with its own bias-matched background), not mirror images. Pick one side
  with `one_direction()` and say which.
* `marker_motifs.csv` mislabels direction in its `Significance` column (workflow bug up
  to v0.10.8: it tests `mpvalue > 0` instead of `mmean > 0`). Use `fixed_marker_motifs()`.
* `all_peaks.csv` joins peak annotation on start/end only (up to v0.10.8), so peaks on
  different chromosomes with the same coordinates get crossed annotations.
  `flag_ambiguous_peaks()`. Both are fixed in v0.10.9; the guards are no-ops on fixed runs.

Latch access goes through `Project._cp` and the `latch` CLI named by `LATCH_BIN`. Compare
tables are small, so this module may run locally when the CLI works there; otherwise run
it on analysis-box.
"""
from __future__ import annotations

import json
import os
import re
import shlex
import subprocess
from dataclasses import dataclass, field
from datetime import datetime
from typing import Dict, Iterable, List, Optional

import pandas as pd

from atx_project import LATCH, Project, latch_ls

COMPARE_ROOT = "latch:///compare_outs"

#: feature -> (results subdirectory, all-features table, effect-size column)
FEATURES = {
    "gene": ("gene_results", "all_genes.csv", "Log2FC"),
    "motif": ("motif_results", "all_motifs.csv", "MeanDiff"),
    "peak": ("peak_results", "all_peaks.csv", "Log2FC"),
}

#: What each test is, from wf/compare_clusters.R. Quoted in reports so a number is never
#: shown without the test that produced it.
METHODS = {
    "gene": "ArchR getMarkerFeatures on GeneScoreMatrix, t-test, bias TSSEnrichment + log10(nFrags)",
    "peak": "ArchR getMarkerFeatures on PeakMatrix, Wilcoxon, bias TSSEnrichment + log10(nFrags)",
    "motif": ("ArchR getMarkerFeatures on MotifMatrix (chromVAR deviation z-scores), Wilcoxon, "
              "normBy none, bias TSSEnrichment + log10(nFrags); effect size is MeanDiff of z"),
    "motif_enrichment": ("ArchR peakAnnoEnrichment (hypergeometric) of motifs in differential "
                         "peaks, Pval <= 0.1 and Log2FC > 0 (up) or < 0 (down)"),
}

#: Plots' Compare tab defaults (plots_epigenome/compare/{gene,motif}_results.py).
PLOTS_FDR = 0.05
PLOTS_EFFECT = 0.01

#: plots_epigenome/compare/check_selection.py names the uploaded barcode file like this.
_CONFIG_RE = re.compile(
    r"^compare_config_(?P<slug>.+)_A(?P<n_a>\d+)_B(?P<n_b>\d+)_(?P<ts>\d{8}_\d{6})\.json$")


# ------------------------------------------------------------------ listing runs

def _latch_ls(path: str, timeout: int = 300) -> List[str]:
    """Full names in a Latch directory (directories keep their trailing '/')."""
    return latch_ls(path, timeout)


def list_compare_runs() -> List[str]:
    """Every run under latch:///compare_outs. List this BEFORE asking for a new run."""
    return sorted(n.rstrip("/") for n in _latch_ls(COMPARE_ROOT) if n.endswith("/"))


def parse_config_name(name: str) -> Optional[dict]:
    """Group sizes and time from a Plots `compare_config_*.json` file name, else None."""
    m = _CONFIG_RE.match(os.path.basename(name))
    if not m:
        return None
    return {"file": os.path.basename(name), "dataset": m["slug"],
            "n_groupA": int(m["n_a"]), "n_groupB": int(m["n_b"]),
            "created": datetime.strptime(m["ts"], "%Y%m%d_%H%M%S")}


def find_compare_configs(project: Project) -> pd.DataFrame:
    """The barcode files Plots wrote into the project root, one row per launch attempt.

    These are the only record of how many cells each group had *before* the workflow
    capped them, and of when a selection was made. They do not record the run name.
    """
    rows = [r for r in (parse_config_name(n) for n in _latch_ls(project.latch_root)) if r]
    return pd.DataFrame(rows).sort_values("created") if rows else pd.DataFrame(
        columns=["file", "dataset", "n_groupA", "n_groupB", "created"])


# ------------------------------------------------------------------ loading a run

@dataclass
class CompareRun:
    name: str
    tables: Dict[str, pd.DataFrame] = field(default_factory=dict)
    assignment: Optional[pd.DataFrame] = None
    missing: List[str] = field(default_factory=list)

    @property
    def group_sizes(self) -> Dict[str, int]:
        """Cells assigned to each group (before the max_cells cap inside the test)."""
        if self.assignment is None:
            return {}
        return self.assignment["UpdateClustName"].value_counts().astype(int).to_dict()


def _fetch(project: Project, run: str, rel: str) -> str:
    local = os.path.join(project.local_root, "compare", run, os.path.basename(rel))
    return project._cp(f"{COMPARE_ROOT}/{run}/{rel}", local)


def load_run(project: Project, run: str,
             features: Iterable[str] = ("gene", "motif", "peak"),
             extras: bool = True) -> CompareRun:
    """Fetch one run's tables. A feature whose table is absent is recorded in `.missing`
    rather than raised: the workflow skips peak and motif steps that fail and still
    succeeds, so a partial run is normal and must be reported, not hidden."""
    out = CompareRun(run)
    for feat in features:
        sub, fname, _ = FEATURES[feat]
        try:
            out.tables[feat] = pd.read_csv(_fetch(project, run, f"{sub}/{fname}"))
        except RuntimeError:
            out.missing.append(feat)
    if extras:
        for key, rel in [("marker_genes", "gene_results/marker_genes.csv"),
                         ("marker_motifs", "motif_results/marker_motifs.csv"),
                         ("motifs_up", "motif_results/upRegulated_motifs.csv"),
                         ("motifs_down", "motif_results/downRegulated_motifs.csv")]:
            try:
                out.tables[key] = pd.read_csv(_fetch(project, run, rel))
            except RuntimeError:
                out.missing.append(key)
    out.assignment = pd.read_csv(_fetch(project, run, "UpdateClustName_by_barcode.csv"))
    return out


def check_group_mapping(run: CompareRun, n_groupA: int, n_groupB: int) -> dict:
    """Confirm which group is which by cell count before naming them in a figure.

    GroupA/GroupB are positional. Compare the run's assignment counts with the counts the
    selection was supposed to have (from a compare_config file name, or the scientist).
    Raises when the counts match only after swapping, which is the mistake this catches.
    """
    got = run.group_sizes
    a, b = got.get("GroupA", 0), got.get("GroupB", 0)
    if (a, b) == (n_groupA, n_groupB):
        return {"ok": True, "GroupA": a, "GroupB": b}
    if (a, b) == (n_groupB, n_groupA):
        raise ValueError(f"{run.name}: groups are SWAPPED relative to the expected counts "
                         f"(GroupA={a}, GroupB={b}; expected A={n_groupA}, B={n_groupB})")
    return {"ok": False, "GroupA": a, "GroupB": b,
            "expected": {"GroupA": n_groupA, "GroupB": n_groupB}}


# ------------------------------------------------------------------ reading the tables

def one_direction(df: pd.DataFrame, group: str = "GroupA") -> pd.DataFrame:
    """Rows of the test for `group` against the other group.

    getMarkerFeatures tests each group against the rest separately, each with its own
    bias-matched and max_cells-capped background, so the GroupB rows are NOT the GroupA
    rows with the sign flipped (on a real run FDRs differ by up to ~1 between the two).
    Report from one side and name it. Plots and the delivered decks use GroupA.
    """
    if "group_name" not in df:
        return df
    sub = df[df["group_name"] == group]
    if sub.empty:
        raise ValueError(f"no rows for {group!r}; groups present: "
                         f"{sorted(df['group_name'].unique())}")
    return sub.copy()


def significant(df: pd.DataFrame, feature: str, group: str = "GroupA",
                fdr: float = PLOTS_FDR, effect: float = PLOTS_EFFECT) -> pd.DataFrame:
    """Features passing FDR <= fdr and |effect| >= effect, with a `higher_in` column.

    Defaults are the Plots Compare tab's (FDR 0.05, effect 0.01). The workflow's own
    marker_*.csv files use a different rule (raw p < 0.05, |effect| >= 0.4), so always
    state which threshold a count came from.
    """
    eff = FEATURES[feature][2]
    sub = one_direction(df, group)
    sub = sub[(sub["FDR"] <= fdr) & (sub[eff].abs() >= effect)].copy()
    other = "GroupB" if group == "GroupA" else "GroupA"
    sub["higher_in"] = sub[eff].gt(0).map({True: group, False: other})
    return sub.sort_values("FDR")


def summarize(run: CompareRun, group: str = "GroupA", fdr: float = PLOTS_FDR,
              effect: float = PLOTS_EFFECT, top_n: int = 10) -> dict:
    """Counts and top features per direction, with the test and thresholds attached."""
    out = {"run": run.name, "tested_side": group, "fdr": fdr, "effect": effect,
           "group_sizes": run.group_sizes, "missing": run.missing, "features": {}}
    for feat in ("gene", "motif", "peak"):
        if feat not in run.tables:
            continue
        sig = significant(run.tables[feat], feat, group, fdr, effect)
        eff = FEATURES[feat][2]
        label = "name" if feat != "peak" else "nearestGene"
        per = {}
        for side, part in sig.groupby("higher_in"):
            part = part.reindex(part[eff].abs().sort_values(ascending=False).index)
            per[side] = {"n": int(len(part)),
                         "top": part[label].astype(str).head(top_n).tolist()}
        out["features"][feat] = {"method": METHODS[feat], "effect_column": eff,
                                 "n_tested": int(len(one_direction(run.tables[feat], group))),
                                 "significant": per}
    return out


def fixed_marker_motifs(df: pd.DataFrame, p: float = 0.05, effect: float = 0.4) -> pd.DataFrame:
    """marker_motifs.csv with its Significance column recomputed correctly.

    compare_clusters.R labels direction with `mpvalue > 0`, which is always true, so every
    significant motif is called GroupA whatever the sign of MeanDiff. This reapplies the
    workflow's own rule (p < 0.05, |MeanDiff| >= 0.4) with the sign of `mmean`.
    """
    out = df.copy()
    sig = (out["mpvalue"] < p) & (out["mmean"].abs() >= effect)
    out["Significance"] = "Not significant"
    out.loc[sig & (out["mmean"] > 0), "Significance"] = "GroupA"
    out.loc[sig & (out["mmean"] < 0), "Significance"] = "GroupB"
    out["relabeled"] = out["Significance"].ne(df["Significance"])
    return out


def flag_ambiguous_peaks(df: pd.DataFrame) -> pd.DataFrame:
    """Mark peaks whose annotation may belong to a different chromosome.

    all_peaks.csv merges the peak annotation onto the statistics by start/end only, so two
    peaks on different chromosomes with the same coordinates are cross-joined: one peak
    appears twice, once with the other chromosome's nearestGene and distances. The
    statistics are fine; the annotation on those rows is not. Do not quote nearestGene for
    a row flagged here.
    """
    out = df.copy()
    key = ["group_name", "seqnames", "start", "end"]
    out["ambiguous_annotation"] = out.duplicated(key, keep=False)
    return out


# ------------------------------------------------------------------ figures

def compare_volcano(df: pd.DataFrame, feature: str, group: str = "GroupA",
                    fdr: float = PLOTS_FDR, effect: float = PLOTS_EFFECT,
                    width: int = 750, height: int = 640, top_n: int = 2):
    """The Plots Compare-tab volcano: FDR against Log2FC (genes) or MeanDiff (motifs)."""
    import atx_plots_lib as L
    eff = FEATURES[feature][2]
    sub = one_direction(df, group)
    other = "GroupB" if group == "GroupA" else "GroupA"
    fig = L.plot_volcano(sub, fdr, effect, group, other, pval_key="FDR", l2fc_key=eff,
                         names_key="name", plot_width=width, plot_height=height, top_n=top_n)
    sig = significant(df, feature, group, fdr, effect)
    return fig, {"feature": feature, "side": group, "fdr": fdr, "effect": effect,
                 "effect_column": eff, "n_significant": int(len(sig)),
                 "method": METHODS[feature]}


def compare_rank_plot(df: pd.DataFrame, feature: str, group: str = "GroupA",
                      y_col: Optional[str] = None, color_col: str = "FDR",
                      n_labels: int = 4, width: int = 760, height: int = 620):
    """The Plots Compare-tab rank plot. Genes rank by Log2FC, motifs by MeanDiff, colored
    by FDR with PuBu_r and marker size 6, as in plots_epigenome/compare/*_results.py."""
    import atx_plots_lib as L
    y_col = y_col or FEATURES[feature][2]
    sub = one_direction(df, group)
    fig = L.plot_ranked_feature_plotly(
        sub, y_col=y_col, x_col=None, n_labels=n_labels, label_col="name",
        color_col=color_col, colorscale="PuBu_r", marker_size=6, title="", y_label=y_col)
    fig.update_layout(width=width, height=height)
    return fig, {"feature": feature, "side": group, "y": y_col, "color": color_col,
                 "n": int(len(sub))}


# ------------------------------------------------------------------ preparing a new run

def write_groupings_file(obs: pd.DataFrame, column: str, group_a, group_b,
                         out_dir: str, dataset: str = "dataset") -> dict:
    """Write the barcode file the workflow's File mode reads, exactly as Plots does.

    `obs` is indexed by the `<run_id>#<barcode>-1` cell names the ArchRProject uses
    (obs.csv from the analysis project). Mirrors check_selection.py: a categorical column
    (fewer than 30 distinct values), two different values, both non-empty, no overlap.

    The file is written locally; `upload_groupings()` puts it in the project root where
    Plots would have put it, and `launch_compare()` starts the run once the scientist has
    confirmed it.
    """
    if column not in obs:
        raise KeyError(f"{column!r} not in obs")
    vals = obs[column]
    if not (isinstance(vals.dtype, pd.CategoricalDtype) or vals.nunique() < 30):
        raise ValueError(f"{column!r} looks continuous ({vals.nunique()} values); "
                         "threshold it into a categorical column first")
    if str(group_a) == str(group_b):
        raise ValueError("GroupA and GroupB must be different values")
    a = obs.index[vals.astype(str) == str(group_a)].astype(str).tolist()
    b = obs.index[vals.astype(str) == str(group_b)].astype(str).tolist()
    if not a or not b:
        raise ValueError(f"empty group: {group_a!r} has {len(a)} cells, {group_b!r} has {len(b)}")
    if set(a) & set(b):
        raise ValueError("the two groups overlap")
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", dataset).strip("._-") or "dataset"
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    name = f"compare_config_{slug}_A{len(a)}_B{len(b)}_{ts}.json"
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, name)
    with open(path, "w") as fh:
        json.dump({"groupA": a, "groupB": b}, fh)
    return {"path": path, "file": name, "column": column,
            "GroupA": str(group_a), "GroupB": str(group_b),
            "n_groupA": len(a), "n_groupB": len(b)}


def launch_parameters(project: Project, run_name: str, groupings_remote: str,
                      max_cells: int = 500, use_max_possible_cells: bool = False,
                      closest: bool = False) -> dict:
    """The compare clusters inputs for a File-mode run, as Plots builds them.

    `closest` defaults to False because that is what a Plots launch sends unless the
    "Closest cells" box is ticked; the workflow's own default is True. Results differ
    between the two, so pick one deliberately and record it with the run.

    Compare needs an ArchRProject, which only create ArchRProject (epi_analysis_archr)
    keeps; ATX_snap deletes its ArchR checkpoint after a successful run.
    """
    if project.backend != "archr":
        raise ValueError("compare clusters needs an ArchRProject; ATX_snap projects "
                         "(epi_analysis_snap) do not keep one")
    return {"project_name": run_name, "groupings": groupings_remote,
            "archrproject": f"{project.latch_root}/{project.name}_ArchRProject",
            "max_cells": max_cells, "use_max_possible_cells": use_max_possible_cells,
            "closest": closest}


# ------------------------------------------------------------------ launching (confirmed)
#
# Launching costs compute and writes to the scientist's workspace, so it happens only after
# the scientist has said yes to the exact summary `describe_launch()` produces. Pass
# `confirmed=True` only then. Runs from any interpreter that has the latch Python package
# (>= 2.62) and a signed-in token: analysis-box's /opt/analysis venv, or locally when the
# CLI works there.

COMPARE_WF = "wf.__init__.compare_workflow"
#: The version the Plots Compare tab launches (plots_epigenome/compare/launch_workflow.py).
#: It predates the v0.10.9 fixes, so its output still needs the guards above.
PLOTS_PINNED_VERSION = "0.10.6-c5165d-8acb83"


def registered_versions(timeout: int = 300) -> pd.DataFrame:
    """Every registered compare clusters version in this workspace, oldest first. `dev`
    marks `-wip-` builds (registered from uncommitted code); launch a release, not those.
    Pass the full `version` string (e.g. 0.10.6-c5165d-8acb83) to the launch helpers."""
    cmd = f"{shlex.quote(LATCH)} get-wf --name {shlex.quote(COMPARE_WF)}"
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    rows = []
    for line in p.stdout.splitlines():
        parts = [x.strip() for x in line.split("\t")]
        if len(parts) == 3 and parts[1] == COMPARE_WF and re.match(r"\d+\.\d+\.\d+", parts[2]):
            rows.append({"id": parts[0], "version": parts[2], "dev": "-wip-" in parts[2]})
    if not rows:
        raise RuntimeError(f"no versions of {COMPARE_WF} listed\n{p.stdout[-400:]}\n{p.stderr[-400:]}")
    df = pd.DataFrame(rows)
    df["semver"] = df["version"].map(_semver)
    return df.sort_values(["semver", "id"]).reset_index(drop=True)


def latest_release() -> str:
    """The newest registered non-dev version."""
    rel = registered_versions().query("~dev")
    return rel["version"].iloc[-1]


def upload_groupings(project: Project, info: dict, timeout: int = 600) -> str:
    """Upload a `write_groupings_file()` result to the project root, where Plots puts its
    compare_config files (so `find_compare_configs()` sees it). Returns the remote path."""
    remote = f"{project.latch_root}/{info['file']}"
    cmd = f"{shlex.quote(LATCH)} cp {shlex.quote(info['path'])} {shlex.quote(remote)}"
    p = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
    if p.returncode != 0 or info["file"] not in _latch_ls(project.latch_root):
        raise RuntimeError(f"upload of {info['file']} failed\n{p.stdout[-800:]}\n{p.stderr[-800:]}")
    return remote


def _semver(version: str) -> tuple:
    """(0, 10, 6) from '0.10.6-c5165d-8acb83'."""
    return tuple(int(x) for x in re.match(r"(\d+)\.(\d+)\.(\d+)", version).groups())


def describe_launch(params: dict, version: str, groups: Optional[dict] = None) -> str:
    """The summary to show the scientist before launching. `groups` is the
    `write_groupings_file()` result, for what each group is and its size."""
    lines = [f"Launch compare clusters {version} ({COMPARE_WF})",
             f"  output: {COMPARE_ROOT}/{params['project_name']}/",
             f"  ArchRProject: {params['archrproject']}",
             f"  groupings file: {params['groupings']}"]
    if groups:
        if groups["file"] != params["groupings"].rstrip("/").rsplit("/", 1)[-1]:
            raise ValueError(f"groups describe {groups['file']} but the launch uses "
                             f"{params['groupings']}")
        lines.append(f"  GroupA = {groups['column']} == {groups['GroupA']!r} ({groups['n_groupA']} cells)")
        lines.append(f"  GroupB = {groups['column']} == {groups['GroupB']!r} ({groups['n_groupB']} cells)")
    cap = ("min(nA, nB) cells per group" if params["use_max_possible_cells"]
           else f"at most {params['max_cells']} cells per group")
    lines.append(f"  cell cap: {cap}")
    lines.append(f"  closest: {params['closest']} "
                 f"({'workflow default' if params['closest'] else 'what Plots sends'})")
    if _semver(version) < (0, 10, 9):
        lines.append("  note: this version predates the v0.10.9 fixes (motif direction labels, "
                     "peak annotation join, multiple-condition matching)")
    lines.append("Uses Latch compute in the scientist's workspace.")
    return "\n".join(lines)


def _to_inputs(params: dict) -> dict:
    import latch.ldata.path  # noqa: F401  (import first: avoids a latch_cli import cycle)
    from latch.types.directory import LatchDir
    from latch.types.file import LatchFile
    from latch_cli.utils.path import normalize_path
    # Pin both to this workspace (latch://<id>.account/...) rather than relying on how
    # latch:/// resolves inside the execution.
    out = dict(params)
    out["groupings"] = LatchFile(normalize_path(params["groupings"]))
    out["archrproject"] = LatchDir(normalize_path(params["archrproject"]))
    return out


def check_launch(params: dict, version: str, groups: Optional[dict] = None) -> dict:
    """Everything `launch_compare()` checks, without launching: the workflow version is
    registered in this workspace, the inputs convert to its interface, the inputs exist and
    the output folder is free. Safe to call before asking the scientist.

    Pass `groups` (the `write_groupings_file()` result) before the upload: the groupings
    file is then checked locally instead of on Latch."""
    from latch.utils import current_workspace
    from latch_cli.services.launch.interface import get_workflow_interface
    from latch_cli.services.launch.type_converter import convert_inputs_to_literals
    import google.protobuf.json_format as gpjson
    from flyteidl.core import interface_pb2
    from flytekit.models.interface import VariableMap

    run = params["project_name"]
    if f"{run}/" in _latch_ls(COMPARE_ROOT):
        raise ValueError(f"{COMPARE_ROOT}/{run}/ already exists; pick another run name "
                         "(a second run into the same folder mixes outputs)")
    remote_keys = ["archrproject"]
    if groups is None:
        remote_keys.append("groupings")
    elif (groups["file"] != params["groupings"].rstrip("/").rsplit("/", 1)[-1]
          or not os.path.isfile(groups["path"])):
        raise ValueError(f"groupings {params['groupings']} does not match the local file "
                         f"{groups['path']}")
    for key in remote_keys:
        parent, name = params[key].rstrip("/").rsplit("/", 1)
        if not any(n.rstrip("/") == name for n in _latch_ls(parent)):
            raise FileNotFoundError(f"{key} not found: {params[key]}")
    ws = current_workspace()
    wf_id, interface, _ = get_workflow_interface(ws, COMPARE_WF, version)
    types = VariableMap.from_flyte_idl(
        gpjson.ParseDict(interface, interface_pb2.VariableMap())).variables
    convert_inputs_to_literals(params=_to_inputs(params), flyte_interface_types=types)
    return {"workspace": ws, "workflow_id": wf_id, "version": version}


def launch_compare(params: dict, version: str, *, confirmed: bool = False) -> dict:
    """Launch compare clusters. Only after the scientist has approved `describe_launch()`
    for these exact params and version; `confirmed` must say so.

    Returns a record to keep with the run (execution id, version, every input including
    `closest`, which the run's output folder does not record)."""
    if not confirmed:
        raise PermissionError("show describe_launch() to the scientist and launch only "
                              "after they say yes; then pass confirmed=True")
    check_launch(params, version)
    from latch_cli.services.launch.launch_v2 import launch
    ex = launch(wf_name=COMPARE_WF, params=_to_inputs(params), version=version)
    return {"execution_id": str(ex.id), "workflow": COMPARE_WF, "version": version,
            "launched_at": datetime.now().isoformat(timespec="seconds"),
            "output": f"{COMPARE_ROOT}/{params['project_name']}/", "params": dict(params)}


def wait_for_run(execution_id: str, poll_seconds: int = 60, timeout_hours: float = 12,
                 progress=print) -> str:
    """Poll a launched run until it finishes. Returns SUCCEEDED, FAILED or ABORTED; raises
    TimeoutError after `timeout_hours`. Run it in the background: compare clusters takes
    tens of minutes to hours. The workflow succeeds even when its optional peak or motif
    steps fail, so check `load_run(...).missing` afterwards."""
    import time
    from latch_cli.services.launch.launch_v2 import Execution
    ex = Execution(id=str(execution_id), python_outputs={})
    deadline = time.time() + timeout_hours * 3600
    last = None
    for _ in ex.poll():
        if ex.status != last:
            progress(f"{datetime.now():%H:%M:%S} execution {execution_id}: {ex.status}")
            last = ex.status
        if ex.status in {"SUCCEEDED", "FAILED", "ABORTED"}:
            return ex.status
        if time.time() > deadline:
            raise TimeoutError(f"execution {execution_id} still {ex.status} after {timeout_hours} h")
        time.sleep(poll_seconds)
