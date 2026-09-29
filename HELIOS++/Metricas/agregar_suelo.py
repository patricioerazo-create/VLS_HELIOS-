#!/usr/bin/env python3
"""Add synthetic ground grid for 3DFin CSF.

Sets:
  las     -> Metricas/Nubes Terreno  -> Metricas/3DFin/las/
  helios  -> Metricas/Nubes HELIOS    -> Metricas/3DFin/helios/
            Ground Z from Nubes Reales (same id) or DEM; not cloud z_min when
            the trunk base is missing (helios z_min above reference).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import laspy
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from alinear_las_dem import base_xy_zmin
from dem_grid import DEFAULT_DEM, DemGrid

_METRICAS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = _METRICAS_DIR.parent
OUTPUT_ROOT = _METRICAS_DIR / "3DFin"

SET_CONFIG: dict[str, dict] = {
    "las": {
        "input": _METRICAS_DIR / "Nubes Terreno",
        "output": OUTPUT_ROOT / "las",
        "ground_z": "cloud",
    },
    "helios": {
        "input": _METRICAS_DIR / "Nubes HELIOS",
        "output": OUTPUT_ROOT / "helios",
        "ground_z": "reference",
    },
}

RADIO_PARCHE_MIN_M = 4.0
RADIO_PARCHE_FACTOR = 1.2
PASO_REJILLA_M = 0.05
GROSOR_BASE_M = 0.15
Z_OFFSET_M = -0.02
RADIO_INTERIOR_M = 0.5
FUSTE_GAP_M = 0.25
CLASS_SYNTH_GROUND = 20
REFERENCE_LAS = _METRICAS_DIR / "Nubes Reales"



def reference_z_min(tree_id: str, bx: float, by: float, dem: DemGrid | None) -> float:
    for ext in (".las", ".laz"):
        ref_path = REFERENCE_LAS / f"{tree_id}{ext}"
        if ref_path.is_file():
            ref = laspy.read(ref_path)
            return float(np.min(ref.z))
    if dem is None:
        raise FileNotFoundError(
            f"No reference LAS for {tree_id} and DEM not loaded"
        )
    return dem.z_at(bx, by)


def ground_z_for_cloud(
    coords: np.ndarray,
    *,
    mode: str,
    tree_id: str,
    dem: DemGrid | None,
) -> tuple[float, float, float, str]:
    x, y, z = coords[:, 0], coords[:, 1], coords[:, 2]
    bx, by, z_cloud_min = base_xy_zmin(x, y, z)

    if mode == "cloud":
        z_ground = z_cloud_min + Z_OFFSET_M
        return bx, by, z_ground, "cloud_zmin"

    z_ref = reference_z_min(tree_id, bx, by, dem)
    gap = z_cloud_min - z_ref
    if gap > FUSTE_GAP_M:
        note = f"reference (fuste gap {gap:+.2f} m)"
        z_ground = z_ref + Z_OFFSET_M
    else:
        note = "reference"
        z_ground = z_ref + Z_OFFSET_M
    return bx, by, z_ground, note


def patch_radius(coords: np.ndarray, x_centro: float, y_centro: float) -> float:
    r_xy = np.hypot(coords[:, 0] - x_centro, coords[:, 1] - y_centro)
    return max(RADIO_PARCHE_MIN_M, float(r_xy.max()) * RADIO_PARCHE_FACTOR)


def generate_ground(x_centro: float, y_centro: float, z_suelo: float, radio: float) -> np.ndarray:
    half = radio + PASO_REJILLA_M
    xs = np.arange(x_centro - half, x_centro + half + PASO_REJILLA_M, PASO_REJILLA_M)
    ys = np.arange(y_centro - half, y_centro + half + PASO_REJILLA_M, PASO_REJILLA_M)
    xx, yy = np.meshgrid(xs, ys)
    rr = np.hypot(xx - x_centro, yy - y_centro)
    mask = (rr >= RADIO_INTERIOR_M) & (rr <= radio)
    n = int(mask.sum())
    return np.column_stack(
        [xx[mask].ravel(), yy[mask].ravel(), np.full(n, z_suelo, dtype=np.float64)]
    )


def save_las_from_source(
    src_las: laspy.LasData,
    ground_coords: np.ndarray,
    out_path: Path,
) -> None:
    n_tree = len(src_las.points)
    n_ground = len(ground_coords)
    n_total = n_tree + n_ground

    out_header = laspy.LasHeader(
        point_format=src_las.header.point_format.id,
        version=src_las.header.version,
    )
    out_header.offsets = src_las.header.offsets.copy()
    out_header.scales = src_las.header.scales.copy()
    for vlr in src_las.header.vlrs:
        out_header.vlrs.append(vlr)

    las = laspy.LasData(out_header)
    las.x = np.concatenate([np.asarray(src_las.x), ground_coords[:, 0]])
    las.y = np.concatenate([np.asarray(src_las.y), ground_coords[:, 1]])
    las.z = np.concatenate([np.asarray(src_las.z), ground_coords[:, 2]])

    if "classification" in las.point_format.dimension_names:
        cls = np.zeros(n_total, dtype=np.uint8)
        if "classification" in src_las.point_format.dimension_names:
            cls[:n_tree] = np.asarray(src_las.classification)
        cls[n_tree:] = CLASS_SYNTH_GROUND
        las.classification = cls

    out_path.parent.mkdir(parents=True, exist_ok=True)
    las.write(str(out_path))


def process_file(path: Path, cfg: dict, dem: DemGrid | None) -> None:
    src_las = laspy.read(path)
    coords = np.column_stack([src_las.x, src_las.y, src_las.z]).astype(np.float64)
    tree_id = path.stem

    x_centro, y_centro, z_suelo, note = ground_z_for_cloud(
        coords, mode=cfg["ground_z"], tree_id=tree_id, dem=dem
    )
    radio = patch_radius(coords, x_centro, y_centro)
    ground = generate_ground(x_centro, y_centro, z_suelo, radio)

    out_path = cfg["output"] / f"{tree_id}_con_suelo.las"
    save_las_from_source(src_las, ground, out_path)

    rel = out_path.relative_to(_METRICAS_DIR)
    print(
        f"OK {path.name} -> {rel} | {len(coords):,}+{len(ground):,} pts | "
        f"z_suelo={z_suelo:.3f} ({note}) | R≈{radio:.1f} m"
    )


def collect_inputs(input_dir: Path, only: str | None) -> list[Path]:
    files = sorted(input_dir.glob("*.las")) + sorted(input_dir.glob("*.laz"))
    files = sorted({p.stem: p for p in files}.values(), key=lambda p: int(p.stem))
    if only:
        files = [p for p in files if p.stem == only]
    return files


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--set", choices=sorted(SET_CONFIG), default="helios")
    ap.add_argument("--id", help="Single tree id, e.g. 3001")
    ap.add_argument("--input", type=Path, help="Input file")
    ap.add_argument("--input-dir", type=Path, help="Input directory override")
    ap.add_argument("--all", action="store_true", help="Process all trees in input dir")
    ap.add_argument("--dem", type=Path, default=DEFAULT_DEM)
    args = ap.parse_args()

    cfg = SET_CONFIG[args.set].copy()
    if args.input_dir:
        cfg["input"] = args.input_dir

    dem: DemGrid | None = None
    if cfg["ground_z"] == "reference":
        dem = DemGrid.from_tif(args.dem)

    if args.input:
        process_file(args.input.resolve(), cfg, dem)
        return 0

    if args.id:
        path = None
        for ext in (".las", ".laz"):
            p = Path(cfg["input"]) / f"{args.id}{ext}"
            if p.is_file():
                path = p
                break
        if path is None:
            print(f"Not found: {args.id} in {cfg['input']}", file=sys.stderr)
            return 1
        process_file(path, cfg, dem)
        return 0

    if args.all or (not args.id and not args.input):
        files = collect_inputs(Path(cfg["input"]), None)
        if not files:
            print(f"No LAS in {cfg['input']}", file=sys.stderr)
            return 1
        for path in files:
            process_file(path, cfg, dem)
        print(f"\nDone: {len(files)} -> {cfg['output']}")
        return 0

    ap.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
