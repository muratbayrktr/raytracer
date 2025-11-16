#!/usr/bin/env python3
import json
import os
import re
import shutil
import subprocess
from pathlib import Path
from typing import Optional, Tuple, List


# --- Configuration ---
INPUT_ROOT = Path("../inputs_hw2").resolve()
OUTPUT_ROOT = Path("../my_outputs_hw2").resolve()
RAW_ROOT = OUTPUT_ROOT / "raw"
PRETTY_ROOT = OUTPUT_ROOT / "pretty"
RAYTRACER_BIN = Path("./raytracer").resolve()


def find_scene_files(root: Path) -> List[Path]:
    """Recursively find all .json files under root that are likely scene files."""
    scene_files = []
    for dirpath, _, filenames in os.walk(root):
        for name in filenames:
            if not name.lower().endswith(".json"):
                continue
            # Skip output metadata jsons if any are inside inputs by mistake
            if name.endswith("_results.json"):
                continue
            scene_files.append(Path(dirpath) / name)
    return scene_files


def parse_image_name(scene_path: Path) -> Optional[str]:
    """Parse the scene JSON and return the ImageName (e.g. ellipsoids.png)."""
    try:
        with scene_path.open("r") as f:
            data = json.load(f)
    except Exception as e:
        print(f"[WARN] Failed to parse JSON {scene_path}: {e}")
        return None

    scene = data.get("Scene") or {}
    cameras = scene.get("Cameras") or {}
    cam = cameras.get("Camera")

    if cam is None:
        print(f"[WARN] No Camera in Scene for {scene_path}")
        return None

    # Camera can be a dict or list of dicts
    if isinstance(cam, list):
        cam = cam[0] if cam else None

    if not isinstance(cam, dict):
        print(f"[WARN] Unexpected Camera format in {scene_path}")
        return None

    img_name = cam.get("ImageName")
    if not img_name:
        print(f"[WARN] No ImageName in Camera for {scene_path}")
        return None

    return str(img_name).strip()


def expected_output_names(image_name: str) -> Tuple[str, list]:
    """
    From a PNG image name (e.g. 'ellipsoids.png') derive:
    - ppm file name: 'ellipsoids.ppm'
    - candidate result json file names.
    """
    base = Path(image_name).name  # strip any path inside ImageName
    stem, _ext = os.path.splitext(base)  # 'ellipsoids', '.png'

    ppm_name = f"{stem}.ppm"

    # We don't know exactly how the raytracer names the results JSON,
    # so we try a few reasonable possibilities:
    candidates = [
        f"{stem}_results.json",
        f"{stem}.ppm_results.json",
        f"{stem}.png_results.json",
    ]

    return ppm_name, candidates


def find_file_anywhere(root: Path, filename: str) -> Optional[Path]:
    """Search for a file with given filename anywhere under root."""
    for path in root.rglob(filename):
        if path.is_file():
            return path
    return None


def find_results_anywhere(root: Path, candidates: list) -> Optional[Path]:
    """Search for any of the candidate result filenames anywhere under root."""
    for cand in candidates:
        found = find_file_anywhere(root, cand)
        if found:
            return found
    return None


def ensure_rendered_and_placed(scene_path: Path) -> None:
    """
    Ensure the scene is rendered and its outputs are placed in my_outputs_hw2/raw/
    mirroring the input folder structure.
    """
    rel_dir = scene_path.relative_to(INPUT_ROOT).parent
    raw_dir = RAW_ROOT / rel_dir
    raw_dir.mkdir(parents=True, exist_ok=True)

    img_name = parse_image_name(scene_path)
    if img_name is None:
        return

    ppm_name, result_candidates = expected_output_names(img_name)

    # Check if ppm already exists anywhere under my_outputs_hw2
    ppm_anywhere = find_file_anywhere(OUTPUT_ROOT, ppm_name)

    if ppm_anywhere is None:
        # Need to render
        print(f"[INFO] Rendering scene: {scene_path}")
        try:
            subprocess.run(
                [str(RAYTRACER_BIN), str(scene_path)],
                check=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"[ERROR] raytracer failed for {scene_path}: {e}")
            return

        # After rendering, look for the ppm again
        ppm_anywhere = find_file_anywhere(OUTPUT_ROOT, ppm_name)
        if ppm_anywhere is None:
            print(f"[WARN] Could not find expected output ppm '{ppm_name}' after rendering {scene_path}")
            return
    else:
        print(f"[INFO] Skipping render for {scene_path}, found existing {ppm_name}")

    # Find results JSON (optional but useful)
    results_anywhere = find_results_anywhere(OUTPUT_ROOT, result_candidates)
    if results_anywhere is None:
        print(f"[WARN] No results JSON found for image {img_name} (candidates: {result_candidates})")

    # Copy them into raw/ mirrored directory
    target_ppm = raw_dir / ppm_name
    if not target_ppm.exists():
        print(f"[INFO] Copying {ppm_anywhere} -> {target_ppm}")
        shutil.copy2(ppm_anywhere, target_ppm)
    else:
        print(f"[INFO] Raw ppm already in place: {target_ppm}")

    if results_anywhere is not None:
        target_results = raw_dir / results_anywhere.name
        if not target_results.exists():
            print(f"[INFO] Copying {results_anywhere} -> {target_results}")
            shutil.copy2(results_anywhere, target_results)
        else:
            print(f"[INFO] Raw results already in place: {target_results}")


def ppm_to_png(ppm_path: Path, png_path: Path) -> None:
    """Convert a .ppm file to .png using ffmpeg."""
    if png_path.exists():
        return

    png_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"[INFO] Converting {ppm_path} -> {png_path}")
    try:
        subprocess.run(
            ["ffmpeg", "-y", "-i", str(ppm_path), str(png_path)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] ffmpeg failed while converting {ppm_path}: {e}")


def build_pretty_outputs() -> None:
    """
    From my_outputs_hw2/raw/, generate PNGs and videos into my_outputs_hw2/pretty/
    with the same folder structure.

    - Every *.ppm -> *.png in the mirrored pretty/ directory.
    - If we detect numbered sequences like name_000.ppm, name_001.ppm, ...
      we also build name.mp4 from the corresponding PNG frames.
    """
    for dirpath, _dirnames, filenames in os.walk(RAW_ROOT):
        raw_dir = Path(dirpath)
        rel = raw_dir.relative_to(RAW_ROOT)
        pretty_dir = PRETTY_ROOT / rel
        pretty_dir.mkdir(parents=True, exist_ok=True)

        ppm_files = [f for f in filenames if f.lower().endswith(".ppm")]
        if not ppm_files:
            continue

        # First, convert all PPM to PNG
        numbered_groups = {}  # prefix -> list of (frame_index, png_path)
        frame_regex = re.compile(r"^(?P<prefix>.+)_(?P<frame>\d{3})\.ppm$", re.IGNORECASE)

        for fname in ppm_files:
            ppm_path = raw_dir / fname
            stem = ppm_path.stem  # e.g., 'davids_000' or 'ellipsoids'
            png_path = pretty_dir / f"{stem}.png"

            ppm_to_png(ppm_path, png_path)

            m = frame_regex.match(fname)
            if m:
                prefix = m.group("prefix")  # 'davids'
                frame_index = int(m.group("frame"))
                numbered_groups.setdefault(prefix, []).append((frame_index, png_path))

        # Then, create videos for sequences with multiple frames
        for prefix, frames in numbered_groups.items():
            if len(frames) <= 1:
                # Treat single frame as just a PNG; no video needed
                continue

            frames.sort(key=lambda t: t[0])  # sort by frame index

            video_path = pretty_dir / f"{prefix}.mp4"
            if video_path.exists():
                print(f"[INFO] Video already exists: {video_path}")
                continue

            print(f"[INFO] Creating video {video_path} from frames {prefix}_###.png")
            try:
                # We rely on ffmpeg's pattern matching; assume frames named prefix_000.png, etc.
                subprocess.run(
                    [
                        "ffmpeg",
                        "-y",
                        "-framerate",
                        "30",
                        "-i",
                        f"{prefix}_%03d.png",
                        "-c:v",
                        "libx264",
                        "-pix_fmt",
                        "yuv420p",
                        str(video_path),
                    ],
                    check=True,
                    cwd=str(pretty_dir),
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )
            except subprocess.CalledProcessError as e:
                print(f"[ERROR] ffmpeg failed while making video {video_path}: {e}")


def main() -> None:
    print(f"[INFO] Input root:  {INPUT_ROOT}")
    print(f"[INFO] Output root: {OUTPUT_ROOT}")
    print(f"[INFO] Raytracer:   {RAYTRACER_BIN}")

    RAW_ROOT.mkdir(parents=True, exist_ok=True)
    PRETTY_ROOT.mkdir(parents=True, exist_ok=True)

    # 1. Traverse inputs and ensure all scenes are rendered + placed into raw/
    scene_files = find_scene_files(INPUT_ROOT)
    print(f"[INFO] Found {len(scene_files)} JSON scene candidate(s).")

    for scene in sorted(scene_files):
        ensure_rendered_and_placed(scene)

    # 2. Build pretty/ outputs: PNGs + videos
    build_pretty_outputs()
    print("[INFO] Done.")


if __name__ == "__main__":
    main()