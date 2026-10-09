#!/usr/bin/env python3
"""Copy unique finite-quarter runs into one provenance-preserving result tree."""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
from datetime import datetime, timezone
from pathlib import Path


CASE_NAME = "shiftx0.000_shifty0.000"
CATEGORIES = ("00_model", "01_results", "10_overview", "80_logs", "99_config")
TAG_PATTERN = re.compile(r"_(f\d+\.\d+_e\d+_\d{8}(?:_retry\d+)?)$")


def io_path(path: Path) -> Path:
    """Use the Windows extended path form for source trees beyond MAX_PATH."""
    absolute = str(path.resolve())
    return Path(f"\\\\?\\{absolute}") if os.name == "nt" else Path(absolute)


def tree_stats(path: Path) -> tuple[int, int]:
    files = [item for item in io_path(path).rglob("*") if item.is_file()]
    return len(files), sum(item.stat().st_size for item in files)


def source_tag(source: Path) -> str:
    match = TAG_PATTERN.search(source.name)
    if not match:
        raise ValueError(f"Cannot derive run tag from {source.name}")
    return match.group(1)


def destination(category: str, tag: str, symmetry: str | None = None) -> Path:
    prefix = Path(category)
    if category in {"10_overview", "99_config"}:
        prefix /= "source_runs"
    result = prefix / tag
    return result / symmetry if symmetry else result


def build_plan(sources: list[Path]) -> tuple[list[tuple[Path, Path]], list[dict[str, object]]]:
    plan: list[tuple[Path, Path]] = []
    details: list[dict[str, object]] = []
    overview_sources = 0

    for source in sources:
        tag = source_tag(source)
        case = source / CASE_NAME
        if not case.is_dir():
            raise FileNotFoundError(case)
        top_names = {item.name for item in source.iterdir()}
        symmetry_cases = sorted(item for item in case.iterdir() if item.is_dir())

        if top_names == {CASE_NAME}:
            if {item.name for item in symmetry_cases} != set(CATEGORIES):
                raise ValueError(f"Unexpected all-mode case layout: {case}")
            for category in CATEGORIES:
                child = case / category
                if tree_stats(child)[0]:
                    plan.append((child, destination(category, tag)))
            layout = "all_modes"
        elif top_names == {CASE_NAME, "overview", "four_mode_run_summary.json"}:
            overview_sources += 1
            plan.append((source / "overview", Path("10_overview")))
            for symmetry_dir in symmetry_cases:
                if not symmetry_dir.name.startswith("symmetry_"):
                    raise ValueError(f"Unexpected legacy symmetry directory: {symmetry_dir}")
                if {item.name for item in symmetry_dir.iterdir()} != set(CATEGORIES):
                    raise ValueError(f"Unexpected legacy case layout: {symmetry_dir}")
                for category in CATEGORIES:
                    child = symmetry_dir / category
                    if tree_stats(child)[0]:
                        plan.append((child, destination(category, tag, symmetry_dir.name)))
            plan.append(
                (
                    source / "four_mode_run_summary.json",
                    Path("99_config") / "source_runs" / tag / "four_mode_run_summary.json",
                )
            )
            layout = "legacy_four_symmetry"
        else:
            raise ValueError(f"Unexpected source layout: {source}")

        count, size = tree_stats(source)
        details.append(
            {"path": str(source), "tag": tag, "layout": layout, "file_count": count, "bytes": size}
        )

    if overview_sources != 1:
        raise ValueError(f"Expected exactly one legacy overview source, found {overview_sources}")
    if len({detail["tag"] for detail in details}) != len(details):
        raise ValueError("Run tags must be unique")
    plan.sort(key=lambda mapping: mapping[1] != Path("10_overview"))
    return plan, details


def copy_entry(source: Path, target: Path) -> None:
    if target.exists():
        raise FileExistsError(target)
    target.parent.mkdir(parents=True, exist_ok=True)
    if source.is_dir():
        shutil.copytree(io_path(source), io_path(target), copy_function=shutil.copy2)
        if tree_stats(source) != tree_stats(target):
            raise IOError(f"Copied tree verification failed: {source}")
    else:
        shutil.copy2(io_path(source), io_path(target))
        if source.stat().st_size != target.stat().st_size:
            raise IOError(f"Copied file verification failed: {source}")


def consolidate(target: Path, requested_sources: list[Path], dry_run: bool = False) -> dict[str, object]:
    target = target.resolve()
    resolved_inputs = [path.resolve(strict=True) for path in requested_sources]
    unique_sources = list(dict.fromkeys(resolved_inputs))
    duplicates = [str(path) for index, path in enumerate(resolved_inputs) if path in resolved_inputs[:index]]
    plan, source_details = build_plan(unique_sources)
    source_files = sum(int(detail["file_count"]) for detail in source_details)
    source_bytes = sum(int(detail["bytes"]) for detail in source_details)
    mapped_stats = [tree_stats(source) if source.is_dir() else (1, source.stat().st_size) for source, _ in plan]
    if (sum(item[0] for item in mapped_stats), sum(item[1] for item in mapped_stats)) != (
        source_files,
        source_bytes,
    ):
        raise ValueError("Copy plan does not cover every unique source file exactly once")

    summary: dict[str, object] = {
        "workflow": "finite_quarter_result_consolidation",
        "status": "dry_run" if dry_run else "complete",
        "target": str(target),
        "requested_sources": [str(path) for path in resolved_inputs],
        "unique_source_count": len(unique_sources),
        "duplicate_inputs_removed": duplicates,
        "deduplication_rule": "normalized source paths only; no frequency-based mode deletion",
        "source_file_count": source_files,
        "source_bytes": source_bytes,
        "expected_output_file_count": source_files + 1,
        "sources": source_details,
        "mappings": [
            {"source": str(source), "destination": str(relative), "file_count": stats[0], "bytes": stats[1]}
            for (source, relative), stats in zip(plan, mapped_stats, strict=True)
        ],
    }
    if dry_run:
        return summary
    if target.exists():
        raise FileExistsError(target)

    staging = target.parent / f".fq_integrate_staging_{os.getpid()}"
    if staging.exists():
        raise FileExistsError(staging)
    staging.mkdir()
    print(f"staging={staging}", flush=True)
    for index, (source, relative) in enumerate(plan, start=1):
        print(f"[{index}/{len(plan)}] {source} -> {relative}", flush=True)
        copy_entry(source, staging / relative)

    copied_files, copied_bytes = tree_stats(staging)
    if (copied_files, copied_bytes) != (source_files, source_bytes):
        raise IOError(
            f"Staging verification failed: {(copied_files, copied_bytes)} != {(source_files, source_bytes)}"
        )
    summary["created_at"] = datetime.now(timezone.utc).isoformat()
    manifest = staging / "99_config" / "integration_manifest.json"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if tree_stats(staging)[0] != source_files + 1:
        raise IOError("Manifest file-count verification failed")

    os.replace(staging, target)
    final_files, final_bytes = tree_stats(target)
    if final_files != source_files + 1 or final_bytes <= source_bytes:
        raise IOError("Final target verification failed")
    summary["output_file_count"] = final_files
    summary["output_bytes"] = final_bytes
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", type=Path)
    parser.add_argument("sources", nargs="+", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    print(json.dumps(consolidate(args.target, args.sources, args.dry_run), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
