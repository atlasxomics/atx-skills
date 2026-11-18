"""End-of-session cleanup for the analysis box.

Two different things, with two different risk profiles:

* **Job scratch** (`/data/scratch/.claude-science/jobs/<id>/`) is pure residue -- every
  output was harvested into an artifact before the job was reported. It is always cleared,
  no question asked, except for directories younger than `--scratch-age` hours, which may
  belong to a job still running (this box is shared).

* **The heavy cache** (`anndata/`, `coverages/`) is 20 GB for a project like Demo_001 and
  re-downloads from Latch on demand. It is removed ONLY for a project named with
  `--finished`, i.e. only when the person doing the analysis has said they are done with
  that project. Everything else in the cache -- `tables/`, `spatial/`, `images/`,
  `figures/`, `compare/`, `summaries/` -- is ~150 MB, backs most of the figures, and stays.

    python cleanup.py                      # inventory + clear job scratch
    python cleanup.py --finished Demo_001 # the above, plus drop that project's heavy cache

THE BOX IS SHARED. Every destructive path here is pinned to one named project directory
and checked with realpath before anything is removed; there is deliberately no "all
projects" mode, because a blanket sweep of /data/analysis_cache would take out a
colleague's working set. Stopping the instance is NOT part of this script -- that is the
user's call, and the hourly charge is theirs to decide about.
"""
import argparse
import os
import shutil
import time

CACHE_ROOT = "/data/analysis_cache"
JOBS_ROOT = "/data/scratch/.claude-science/jobs"
HEAVY = ("anndata", "coverages")        # re-fetchable from Latch
KEEP = ("tables", "spatial", "images", "figures", "compare", "summaries")


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f}{unit}" if unit == "B" else f"{n:.1f}{unit}"
        n /= 1024.0


def dir_size(path):
    total = 0
    for root, _dirs, files in os.walk(path, onerror=lambda e: None):
        for f in files:
            try:
                total += os.lstat(os.path.join(root, f)).st_size
            except OSError:
                pass
    return total


def _safe_child(root, *parts):
    """Resolve root/parts and refuse anything that escapes root."""
    target = os.path.realpath(os.path.join(root, *parts))
    if not target.startswith(os.path.realpath(root) + os.sep):
        raise ValueError(f"refusing to touch {target!r}: outside {root!r}")
    return target


def inventory():
    print(f"== cache under {CACHE_ROOT} ==")
    if not os.path.isdir(CACHE_ROOT):
        print("  (absent)")
        return
    for proj in sorted(os.listdir(CACHE_ROOT)):
        pdir = os.path.join(CACHE_ROOT, proj)
        if not os.path.isdir(pdir):
            continue
        heavy = sum(dir_size(os.path.join(pdir, h))
                    for h in HEAVY if os.path.isdir(os.path.join(pdir, h)))
        total = dir_size(pdir)
        print(f"  {proj:20s} total {human(total):>9s}   heavy {human(heavy):>9s}"
              f"   keep {human(total - heavy):>9s}")


def clear_scratch(min_age_hours=2.0, dry_run=False):
    if not os.path.isdir(JOBS_ROOT):
        print(f"== job scratch: {JOBS_ROOT} absent ==")
        return 0
    cutoff = time.time() - min_age_hours * 3600
    freed = kept = removed = 0
    for name in sorted(os.listdir(JOBS_ROOT)):
        d = _safe_child(JOBS_ROOT, name)
        if not os.path.isdir(d):
            continue
        if os.path.getmtime(d) > cutoff:   # may still be running
            kept += 1
            continue
        size = dir_size(d)
        if not dry_run:
            shutil.rmtree(d, ignore_errors=True)
        freed += size
        removed += 1
    verb = "would free" if dry_run else "freed"
    print(f"== job scratch: {removed} dir(s) {verb} {human(freed)}; "
          f"{kept} kept (newer than {min_age_hours}h) ==")
    return freed


def drop_heavy(project, dry_run=False):
    pdir = _safe_child(CACHE_ROOT, project)
    if not os.path.isdir(pdir):
        raise SystemExit(f"no cache for project {project!r} under {CACHE_ROOT}")
    freed = 0
    for h in HEAVY:
        d = _safe_child(pdir, h)
        if not os.path.isdir(d):
            continue
        size = dir_size(d)
        if not dry_run:
            shutil.rmtree(d, ignore_errors=True)
        freed += size
        print(f"  {'would remove' if dry_run else 'removed'} {project}/{h}  {human(size)}")
    kept = [k for k in KEEP if os.path.isdir(os.path.join(pdir, k))]
    print(f"== {project}: {'would free' if dry_run else 'freed'} {human(freed)}; "
          f"kept {', '.join(kept) or 'nothing'} ==")
    return freed


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--finished", metavar="PROJECT",
                    help="project the user has said they are DONE with; its anndata/ and "
                         "coverages/ caches are removed")
    ap.add_argument("--scratch-age", type=float, default=2.0,
                    help="leave job scratch younger than this many hours (default 2)")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would be freed and remove nothing")
    a = ap.parse_args()

    inventory()
    freed = clear_scratch(a.scratch_age, dry_run=a.dry_run)
    if a.finished:
        freed += drop_heavy(a.finished, dry_run=a.dry_run)
    print(f"\nTOTAL {'would free' if a.dry_run else 'freed'}: {human(freed)}")
    if not a.finished:
        print("Heavy cache left in place (no project named with --finished).")
    print("The instance is still running; stop it separately if that is what you want.")


if __name__ == "__main__":
    main()
