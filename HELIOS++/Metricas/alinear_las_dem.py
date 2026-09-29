#!/usr/bin/env python3
"""Shift tree LAS/LAZ Z so z_min matches DEM 2m at tree base XY (bilinear)."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import laspy
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dem_grid import DEFAULT_DEM, DemGrid

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = PROJECT / "Metricas/Nubes Terreno"
DEFAULT_OUTPUT = PROJECT / "Metricas/Nubes Reales"
DEFAULT_OUT_CSV = PROJECT / "Metricas/las_dem_alignment.csv"
BASE_TOL_M = 0.05


def find_input(input_dir: Path, tid: str) -> Path | None:
    for ext in (".las", ".laz"):
        p = input_dir / f"{tid}{ext}"
        if p.is_file():
            return p
    return None


def base_xy_zmin(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, float, float]:
    zmin = float(z.min())
    mask = z <= zmin + BASE_TOL_M
    if not np.any(mask):
        mask = z == z.min()
    return float(x[mask].mean()), float(y[mask].mean()), zmin


def shift_las(in_path: Path, out_path: Path, delta_z: float) -> tuple[int, float, float]:
    las = laspy.read(in_path)
    n = len(las.points)
    z_before = float(las.z.min()) if n else float("nan")
    if abs(delta_z) >= 1e-9 and n:
        las.z = las.z + delta_z
    z_after = float(las.z.min()) if n else float("nan")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    las.write(out_path)
    return n, z_before, z_after


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--dem", type=Path, default=DEFAULT_DEM)
    ap.add_argument("--out-csv", type=Path, default=DEFAULT_OUT_CSV)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    dem = DemGrid.from_tif(args.dem)
    rows: list[dict[str, float | str]] = []
    processed = 0
    skipped = 0
    gaps: list[float] = []

    for tid in range(3001, 3080):
        sid = str(tid)
        in_path = find_input(args.input, sid)
        if in_path is None:
            print(f"skip missing LAS {sid}", file=sys.stderr)
            skipped += 1
            continue

        las = laspy.read(in_path)
        if len(las.points) == 0:
            skipped += 1
            continue
        x = np.asarray(las.x, dtype=np.float64)
        y = np.asarray(las.y, dtype=np.float64)
        z = np.asarray(las.z, dtype=np.float64)
        bx, by, las_zmin = base_xy_zmin(x, y, z)
        dem_z = dem.z_at(bx, by)
        delta = dem_z - las_zmin

        rows.append(
            {
                "id": sid,
                "x_sample": bx,
                "y_sample": by,
                "dem_z_m": dem_z,
                "las_zmin_before_m": las_zmin,
                "delta_z_m": delta,
            }
        )

        out_path = args.output / f"{sid}.las"
        if args.dry_run:
            print(f"{sid}: base=({bx:.2f},{by:.2f}) dem={dem_z:.3f} zmin={las_zmin:.3f} d={delta:+.3f}")
        else:
            n, z0, z1 = shift_las(in_path, out_path, delta)
            gap = z1 - dem_z
            gaps.append(gap)
            print(
                f"{sid}: {n:,} pts | z_min {z0:.3f}->{z1:.3f} | "
                f"dem={dem_z:.3f} gap={gap:+.3f} | d={delta:+.3f}"
            )
        processed += 1

    if rows:
        args.out_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.out_csv.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)
        d = np.array([float(r["delta_z_m"]) for r in rows])
        print(
            f"\nDone: {processed} trees, {skipped} skipped. "
            f"delta_z mean={d.mean():+.3f} std={d.std():.3f} "
            f"range=[{d.min():+.3f}, {d.max():+.3f}]"
        )
        if gaps:
            g = np.array(gaps)
            print(f"Post-shift gap (z_min - dem@base): mean={g.mean():+.4f} max|.|={np.abs(g).max():.4f} m")
        print(f"CSV: {args.out_csv}")
    return 0 if skipped == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
