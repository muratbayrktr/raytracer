#!/usr/bin/env python3
import csv
import json
import os
import re
from pathlib import Path
from typing import Dict, Any, List

RAW_ROOT = Path("../my_outputs_hw2/raw").resolve()
OUTPUT_ROOT = RAW_ROOT.parent
CSV_PATH = OUTPUT_ROOT / "benchmark_raw.csv"

SEQ_RE = re.compile(r"^(?P<base>.+)_(?P<frame>\d{3})_results$", re.IGNORECASE)
SINGLE_RE = re.compile(r"^(?P<base>.+)_results$", re.IGNORECASE)


def find_result_files(root: Path) -> List[Path]:
    """Find all *_results.json files under root."""
    result_files = []
    for dirpath, _, filenames in os.walk(root):
        for name in filenames:
            if not name.endswith(".json"):
                continue
            if not name.endswith("_results.json"):
                continue
            result_files.append(Path(dirpath) / name)
    return result_files


def group_key_for_result(path: Path) -> str:
    """
    Produce a logical group key for a results file.

    Examples (relative to RAW_ROOT):
      raven/davids_000_results.json -> "raven/davids"
      raven/davids_001_results.json -> "raven/davids"
      grass/grass_desert_results.json -> "grass/grass_desert"
      dragon_metal_results.json -> "dragon_metal"
    """
    rel = path.relative_to(RAW_ROOT)
    rel_dir = rel.parent.as_posix()  # e.g. "raven", "raven/dragon", "."
    stem = path.stem  # filename without .json, e.g. "davids_000_results"

    m_seq = SEQ_RE.match(stem)
    if m_seq:
        base = m_seq.group("base")
    else:
        m_single = SINGLE_RE.match(stem)
        base = m_single.group("base") if m_single else stem

    if rel_dir == ".":
        return base
    else:
        return f"{rel_dir}/{base}"


def parse_result_json(path: Path) -> Dict[str, Any]:
    """Parse the results JSON and normalize fields."""
    with path.open("r") as f:
        data = json.load(f)

    # Normalize timing fields (default to 0.0 if missing)
    def get_float(key: str) -> float:
        val = data.get(key)
        try:
            return float(val)
        except (TypeError, ValueError):
            return 0.0

    # Normalize bool-ish fields
    def get_bool_or_none(key: str):
        if key not in data:
            return None
        val = data[key]
        if isinstance(val, bool):
            return val
        if isinstance(val, str):
            if val.lower() in ("true", "1", "yes"):
                return True
            if val.lower() in ("false", "0", "no"):
                return False
        return None

    result = {
        "sceneName": data.get("sceneName"),
        "preprocessingTimeMs": get_float("preprocessingTimeMs"),
        "renderTimeMs": get_float("renderTimeMs"),
        "totalTimeMs": get_float("totalTimeMs"),
        "enableBackFaceCulling": get_bool_or_none("enableBackFaceCulling"),
        "isMultiThreaded": get_bool_or_none("isMultiThreaded"),
        "useBVH": get_bool_or_none("useBVH"),
        "_path": path,
        "_raw": data,
    }
    return result


def aggregate_groups(files: List[Path]) -> List[Dict[str, Any]]:
    """Group result files and average sequences."""
    groups: Dict[str, List[Dict[str, Any]]] = {}

    for path in files:
        key = group_key_for_result(path)
        groups.setdefault(key, []).append(parse_result_json(path))

    summaries = []

    for key, entries in groups.items():
        n = len(entries)
        # Average timings
        avg_pre = sum(e["preprocessingTimeMs"] for e in entries) / n
        avg_render = sum(e["renderTimeMs"] for e in entries) / n
        avg_total = sum(e["totalTimeMs"] for e in entries) / n

        scene_names = [e["sceneName"] for e in entries if e["sceneName"]]
        sorted_scene_names = sorted(scene_names)
        scene_name = sorted_scene_names[0] if sorted_scene_names else ""
        if scene_name.endswith("_000.png"):
            scene_name_parts = scene_name.split(".png")[0].split("_")
            scene_name_parts = scene_name_parts[:-1]
            scene_name = "_".join(scene_name_parts) + ".mp4"


        def summarize_bool(field: str) -> str:
            vals = [e[field] for e in entries if e[field] is not None]
            if not vals:
                return ""
            if all(v is True for v in vals):
                return "true"
            if all(v is False for v in vals):
                return "false"
            return "mixed"

        summary = {
            "group": key,
            "sceneName": scene_name,
            "num_results": n,
            "avg_preprocessingTimeMs": avg_pre,
            "avg_renderTimeMs": avg_render,
            "avg_totalTimeMs": avg_total,
            "enableBackFaceCulling": summarize_bool("enableBackFaceCulling"),
            "isMultiThreaded": summarize_bool("isMultiThreaded"),
            "useBVH": summarize_bool("useBVH"),
        }
        summaries.append(summary)


    summaries.sort(key=lambda r: r["avg_totalTimeMs"], reverse=True)
    return summaries


def write_csv(summaries: List[Dict[str, Any]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "group",
        "sceneName",
        "num_results",
        "avg_preprocessingTimeMs",
        "avg_renderTimeMs",
        "avg_totalTimeMs",
        "enableBackFaceCulling",
        "isMultiThreaded",
        "useBVH",
    ]

    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in summaries:
            writer.writerow(row)


def print_table(summaries: List[Dict[str, Any]]) -> None:
    """Print a simple ASCII table to stdout."""
    headers = [
        "group",
        "sceneName",
        "num",
        "pre(ms)",
        "render(ms)",
        "total(ms)",
        "cull",
        "mt",
        "bvh",
    ]

    rows = []
    for r in summaries:
        rows.append([
            str(r["group"]),
            str(r["sceneName"]),
            str(r["num_results"]),
            f"{r['avg_preprocessingTimeMs']:.1f}",
            f"{r['avg_renderTimeMs']:.1f}",
            f"{r['avg_totalTimeMs']:.1f}",
            r["enableBackFaceCulling"],
            r["isMultiThreaded"],
            r["useBVH"],
        ])

    col_widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            col_widths[i] = max(col_widths[i], len(cell))

    def fmt_row(row_cells):
        return " | ".join(cell.ljust(col_widths[i]) for i, cell in enumerate(row_cells))

    print(fmt_row(headers))
    print("-+-".join("-" * w for w in col_widths))

    for row in rows:
        print(fmt_row(row))


def main() -> None:
    print(f"[INFO] Benchmarking results under: {RAW_ROOT}")
    result_files = find_result_files(RAW_ROOT)
    print(f"[INFO] Found {len(result_files)} *_results.json files.")

    if not result_files:
        print("[WARN] No result files found. Nothing to benchmark.")
        return

    summaries = aggregate_groups(result_files)
    write_csv(summaries, CSV_PATH)
    print(f"[INFO] Wrote CSV: {CSV_PATH}")

    print("\n[INFO] Benchmark table:")
    print_table(summaries)


if __name__ == "__main__":
    main()