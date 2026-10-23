# SPDX-License-Identifier: GPL-3.0-or-later
"""Bounded, resumable caches joining first-pass construction to overlap correction."""

import json
import os
import shutil
import time
from pathlib import Path

from sheet_overlaps import SheetIndex, checksum, correct_sheet, create_patch, digest


def save_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def trim_cache(root, limit, marker, pinned=()):
    entries = sorted(
        (p for p in root.iterdir() if p.is_dir() and (p / marker).is_file()),
        key=lambda p: p.stat().st_mtime_ns,
        reverse=True,
    )
    keep = set(pinned)
    keep.update(p.name for p in entries[:limit])
    for path in entries:
        if path.name not in keep:
            shutil.rmtree(path)


class OverlapProduction:
    def __init__(
        self,
        root,
        index,
        selected,
        fingerprint,
        build,
        *,
        width_m=2000,
        context_m=3600,
        reserve_gib=3,
        sheet_cache_limit=16,
        patch_cache_limit=64,
    ):
        self.root = Path(root)
        self.index = SheetIndex(index, selected)
        self.fingerprint = fingerprint
        self.build = build
        self.width_m, self.context_m = width_m, context_m
        self.reserve_gib = reserve_gib
        self.sheet_cache_limit, self.patch_cache_limit = (
            sheet_cache_limit,
            patch_cache_limit,
        )
        self.sheets = self.root / "initial-cache"
        self.results = self.root / "initial-results"
        self.patches = self.root / "patch-cache"
        self.active_sheets, self.active_patches = set(), set()
        for folder in (self.sheets, self.results, self.patches):
            folder.mkdir(exist_ok=True)

    def space_check(self):
        if shutil.disk_usage(self.root).free < self.reserve_gib * 1024**3:
            # Keep the current job's constraints until its COG has been checked.
            trim_cache(self.sheets, 0, "cache_record.json", self.active_sheets)
            trim_cache(self.patches, 0, "report.json", self.active_patches)
        if shutil.disk_usage(self.root).free < self.reserve_gib * 1024**3:
            raise OSError("Disk reserve reached during overlap preparation")

    def initial(self, sheet, *, require_diagnostics=False):
        folder = self.sheets / sheet
        record_path = folder / "cache_record.json"
        result_path = self.results / f"{sheet}.json"
        if not record_path.exists():
            if result_path.exists() and not require_diagnostics:
                previous = json.loads(result_path.read_text())
                if previous["fingerprint"] != self.fingerprint:
                    raise ValueError("Initial-result cache belongs to another run")
                if previous["outcome"]["status"] == "failed":
                    return None, previous["outcome"]
            self.space_check()
            if folder.exists():
                # A previous process may have stopped before its checkpoint.
                interrupted = (
                    self.root / "interrupted" / f"initial-{sheet}-{time.time_ns()}"
                )
                interrupted.parent.mkdir(exist_ok=True)
                folder.replace(interrupted)
            folder.mkdir()
            print(f"INITIAL {sheet}", flush=True)
            outcome = self.build(sheet, folder)
            files = {
                str(p.relative_to(folder)): checksum(p)
                for p in folder.rglob("*")
                if p.is_file()
            }
            record = {
                "fingerprint": self.fingerprint,
                "outcome": outcome,
                "files": files,
            }
            save_json(record_path, record)
            save_json(
                result_path, {"fingerprint": self.fingerprint, "outcome": outcome}
            )
        else:
            record = json.loads(record_path.read_text())
            if record["fingerprint"] != self.fingerprint:
                raise ValueError("Initial cache belongs to another run")
            for name, expected in record["files"].items():
                if checksum(folder / name) != expected:
                    raise ValueError(f"Cached first-pass file changed: {sheet}/{name}")
            outcome = record["outcome"]
        os.utime(folder, None)
        return folder, outcome

    def prepare(self, sheet, job):
        """Populate a fresh production job; source-conflict failures stay failures."""
        job = Path(job)
        self.active_sheets = set(self.index.neighbourhood(sheet))
        self.active_patches = set()
        initial, outcome = self.initial(sheet, require_diagnostics=True)
        shutil.copytree(initial, job, dirs_exist_ok=True)
        if outcome["status"] != "built":
            return outcome
        # Source geometry hashes stay the same; the validation path now points
        # at this job's independent working copy.
        report_path = job / "dem/report.json"
        report = json.loads(report_path.read_text())
        for source in report["source_files"]:
            path = Path(source["path"])
            source["path"] = str(job / path.relative_to(initial))
        save_json(report_path, report)
        folders, missing = {}, {}
        for neighbour in self.index.neighbourhood(sheet):
            directory, result = self.initial(neighbour)
            if result["status"] == "built":
                folders[neighbour] = directory / "dem"
            else:
                missing[neighbour] = result
        built, statuses = {}, {}

        def patch_for(key):
            if key in statuses:
                return built.get(key)
            specification = self.index.jobs[key]
            members = {s: folders[s] for s in specification["members"] if s in folders}
            absent = sorted(set(specification["members"]) - set(members))
            minimum = 2 if specification["kind"] == "edge" else 3
            if len(members) < minimum:
                statuses[key] = {
                    "status": "unavailable",
                    "missing_members": absent,
                    "reason": "Insufficient successfully built neighbours",
                }
                return None
            dependencies = []
            if specification["kind"] == "corner":
                for edge in self.index.incident_edges(specification):
                    if not set(self.index.jobs[edge]["members"]).issubset(members):
                        continue
                    found = patch_for(edge)
                    if found is None:
                        statuses[key] = {
                            "status": "blocked",
                            "reason": "Incident edge correction is unresolved",
                            "edge": edge,
                        }
                        return None
                    dependencies.append(found)
            signature = digest(
                {
                    "run": self.fingerprint,
                    "job": specification,
                    "inputs": {
                        s: {
                            n: checksum(p / n)
                            for n in ("dem.tif", "interpolation_inputs.npz")
                        }
                        for s, p in members.items()
                    },
                    "prior_edges": [checksum(p / "report.json") for p in dependencies],
                    "width_m": self.width_m,
                    "context_m": self.context_m,
                }
            )
            destination = self.patches / f"{key}-{signature[:16]}"
            try:
                if destination.exists():
                    saved = json.loads((destination / "report.json").read_text())
                    if saved["run_fingerprint"] != self.fingerprint:
                        raise RuntimeError("Cached patch belongs to another run")
                    for name in ("inputs", "patch"):
                        if (
                            checksum(destination / f"{name}.npz")
                            != saved[f"{name}_sha256"]
                        ):
                            raise RuntimeError("Cached overlap arrays changed")
                else:
                    self.space_check()
                    staging = self.patches / f"{key}-{signature[:16]}.partial"
                    if staging.exists():
                        shutil.rmtree(staging)
                    print(
                        f"OVERLAP {specification['kind']} {' '.join(members)}",
                        flush=True,
                    )
                    create_patch(
                        specification,
                        members,
                        self.index,
                        staging,
                        prior_patches=dependencies,
                        width_m=self.width_m,
                        context_m=self.context_m,
                        fingerprint=self.fingerprint,
                    )
                    staging.replace(destination)
                os.utime(destination, None)
                self.active_patches.add(destination.name)
                built[key] = destination
                statuses[key] = {
                    "status": "applied",
                    "missing_members": absent,
                    "report_sha256": checksum(destination / "report.json"),
                }
            except ValueError as exc:
                statuses[key] = {
                    "status": "blocked",
                    "reason": str(exc),
                    "missing_members": absent,
                }
            return built.get(key)

        for key in sorted(self.index.by_sheet[sheet]):
            patch_for(key)
        own_patches = {
            key: built[key] for key in self.index.by_sheet[sheet] if key in built
        }
        own_statuses = {key: statuses[key] for key in self.index.by_sheet[sheet]}
        summary = correct_sheet(
            sheet,
            job / "dem",
            self.index,
            own_patches,
            own_statuses,
            width_m=self.width_m,
        )
        summary["neighbour_failures"] = missing
        save_json(job / "overlap_neighbours.json", summary)
        return {**outcome, "overlap": summary}

    def trim(self):
        self.active_sheets, self.active_patches = set(), set()
        trim_cache(self.sheets, self.sheet_cache_limit, "cache_record.json")
        trim_cache(self.patches, self.patch_cache_limit, "report.json")
