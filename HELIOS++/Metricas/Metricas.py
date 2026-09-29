#!/usr/bin/env python3
"""Crown metrics from tree point clouds (LAS, LAZ, XYZ).

Center XY from arboles_plot_metrics_xy.csv.
Outputs: crown radii (40/130/220/310 deg), projected crown area (10 cm grid).
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import sys
from pathlib import Path

import numpy as np

CROWN_AZIMUTHS_DEG = (40.0, 130.0, 220.0, 310.0)
CROWN_SECTOR_HALF_DEG = 25.0

_METRICAS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = _METRICAS_DIR.parent
DEFAULT_CENTERS_CSV = _METRICAS_DIR / "arboles_plot_metrics_xy.csv"

SET_CONFIG: dict[str, dict] = {
    "terreno": {
        "label": "terreno",
        "default_input": _METRICAS_DIR / "Nubes Terreno",
        "extensions": (".las", ".laz", ".xyz", ".txt"),
        "output_dir": _METRICAS_DIR / "terreno",
    },
    "helios": {
        "label": "HELIOS",
        "default_input": _METRICAS_DIR / "Nubes HELIOS",
        "extensions": (".las", ".laz", ".xyz", ".txt"),
        "output_dir": _METRICAS_DIR / "helios",
    },
}

CSV_FIELDS = [
    "id_arbol",
    "altura_m",
    "superficie_copa_m2",
    "r_40_m",
    "r_130_m",
    "r_220_m",
    "r_310_m",
]

PRESERVE_FIELDS = ("altura_m",)


def load_centers_csv(path: Path) -> dict[str, dict]:
    """Return {tree_id: {x, y, nota} | missing coords}."""
    out: dict[str, dict] = {}
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            tid = str(row.get("id", "")).strip()
            if not tid:
                continue
            xs = str(row.get("x", "")).strip()
            ys = str(row.get("y", "")).strip()
            nota = str(row.get("nota", "")).strip()
            if xs.upper() == "X" or ys.upper() == "X":
                out[tid] = {"x": None, "y": None, "nota": nota or "sin coordenadas"}
                continue
            try:
                out[tid] = {"x": float(xs), "y": float(ys), "nota": nota}
            except ValueError:
                out[tid] = {"x": None, "y": None, "nota": nota or "coordenadas invalidas"}
    return out


def load_xyz_text(path: Path) -> np.ndarray:
    xs, ys, zs = [], [], []
    opener = gzip.open if path.suffix.lower() == ".gz" else path.open
    mode = "rt" if path.suffix.lower() == ".gz" else "r"
    with opener(path, mode, encoding="utf-8", errors="ignore") as f:  # type: ignore[arg-type]
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("//"):
                continue
            parts = line.replace(",", " ").split()
            if len(parts) < 3:
                continue
            try:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
                zs.append(float(parts[2]))
            except ValueError:
                continue
    if not xs:
        raise ValueError(f"No XYZ points in {path}")
    return np.column_stack(
        [
            np.asarray(xs, dtype=np.float64),
            np.asarray(ys, dtype=np.float64),
            np.asarray(zs, dtype=np.float64),
        ]
    )


def load_las(path: Path) -> np.ndarray:
    try:
        import laspy
    except ImportError as exc:
        raise ImportError(
            "Reading LAS/LAZ requires laspy. Install with: pip install laspy[lazrs]"
        ) from exc
    with laspy.open(path) as reader:
        las = reader.read()
    return np.column_stack(
        [
            np.asarray(las.x, dtype=np.float64),
            np.asarray(las.y, dtype=np.float64),
            np.asarray(las.z, dtype=np.float64),
        ]
    )


def load_points(path: Path) -> np.ndarray:
    suffix = path.suffix.lower()
    if suffix in {".xyz", ".txt"}:
        return load_xyz_text(path)
    if suffix in {".las", ".laz"}:
        return load_las(path)
    raise ValueError(f"Unsupported format: {path.suffix}")


def tree_id_from_path(path: Path) -> str:
    return path.stem


def azimuth_unit(az_deg: float) -> tuple[float, float]:
    a = math.radians(az_deg)
    return math.sin(a), math.cos(a)


def _sector_mask(ang: np.ndarray, az_deg: float, half_sec_rad: float) -> np.ndarray:
    a0 = math.radians(az_deg)
    d = (ang - a0 + math.pi) % (2.0 * math.pi) - math.pi
    return np.abs(d) <= half_sec_rad


def measure_crown_radii(
    verts: np.ndarray,
    *,
    cx: float,
    cy: float,
    crown_azs_deg: tuple[float, ...] = CROWN_AZIMUTHS_DEG,
    sector_half_deg: float = CROWN_SECTOR_HALF_DEG,
) -> dict:
    x, y = verts[:, 0], verts[:, 1]
    ang = np.arctan2(x - cx, y - cy)
    half_sec = math.radians(sector_half_deg)
    radii: dict[float, float] = {}

    for az in crown_azs_deg:
        dx, dy = azimuth_unit(az)
        sel = _sector_mask(ang, az, half_sec)
        if not sel.any():
            radii[az] = float("nan")
            continue
        t_crown = (x[sel] - cx) * dx + (y[sel] - cy) * dy
        radii[az] = float(t_crown.max())

    return {
        "centro_x": cx,
        "centro_y": cy,
        "r_40_m": radii[40.0],
        "r_130_m": radii[130.0],
        "r_220_m": radii[220.0],
        "r_310_m": radii[310.0],
        "sector_aperture_deg": 2.0 * sector_half_deg,
    }


def find_inputs(inputs: list[Path], extensions: tuple[str, ...]) -> list[Path]:
    paths: list[Path] = []
    ext_set = {e.lower() for e in extensions}
    for p in inputs:
        if p.is_file() and p.suffix.lower() in ext_set:
            paths.append(p)
        elif p.is_dir():
            for ext in extensions:
                paths.extend(sorted(p.glob(f"*{ext}")))
                paths.extend(sorted(p.glob(f"**/*{ext}")))
    seen: set[Path] = set()
    out: list[Path] = []
    for p in paths:
        rp = p.resolve()
        if rp not in seen:
            seen.add(rp)
            out.append(p)
    return out


def metrics_row(path: Path, metrics: dict) -> dict[str, str | float]:
    def _fmt_radius(v: float) -> str | float:
        if isinstance(v, float) and math.isnan(v):
            return ""
        return round(float(v), 3)

    return {
        "id_arbol": tree_id_from_path(path),
        "altura_m": metrics.get("altura_m", ""),
        "superficie_copa_m2": round(metrics["superficie_copa_m2"], 2),
        "r_40_m": _fmt_radius(metrics.get("r_40_m", "")),
        "r_130_m": _fmt_radius(metrics.get("r_130_m", "")),
        "r_220_m": _fmt_radius(metrics.get("r_220_m", "")),
        "r_310_m": _fmt_radius(metrics.get("r_310_m", "")),
    }


def error_row(path: Path, message: str) -> dict[str, str | float]:
    return {
        "id_arbol": tree_id_from_path(path),
        "altura_m": "",
        "superficie_copa_m2": "",
        "r_40_m": "",
        "r_130_m": "",
        "r_220_m": "",
        "r_310_m": "",
    }


def load_existing_preserved(csv_path: Path) -> dict[str, dict[str, str]]:
    if not csv_path.is_file():
        return {}
    with csv_path.open(newline="", encoding="utf-8") as f:
        return {
            r["id_arbol"]: {k: r.get(k, "") for k in PRESERVE_FIELDS}
            for r in csv.DictReader(f)
            if r.get("id_arbol")
        }


def resolve_set_paths(set_name: str, args: argparse.Namespace) -> tuple[Path, Path]:
    cfg = SET_CONFIG[set_name]
    input_dir = args.input or cfg["default_input"]
    output_dir = args.output_dir or cfg["output_dir"]
    csv_path = args.output or (output_dir / "metricas.csv")
    return input_dir, csv_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Crown radii per tree (center from CSV)")
    parser.add_argument(
        "--set",
        choices=sorted(SET_CONFIG),
        default="terreno",
        help="Measurement set (defines default input/output paths)",
    )
    parser.add_argument("inputs", nargs="*", type=Path, help="Input file(s) or directories")
    parser.add_argument("--input", type=Path, default=None, help="Input directory override")
    parser.add_argument("-o", "--output", type=Path, default=None, help="Output CSV path")
    parser.add_argument("--output-dir", type=Path, default=None, help="Output directory for CSV")
    parser.add_argument(
        "--centers-csv",
        type=Path,
        default=DEFAULT_CENTERS_CSV,
        help=f"Tree center XY CSV (default: {DEFAULT_CENTERS_CSV.name})",
    )
    parser.add_argument("--only", type=str, default=None, help="Process only this tree id")
    args = parser.parse_args(argv)

    if not args.centers_csv.is_file():
        print(f"Centers CSV not found: {args.centers_csv}", file=sys.stderr)
        return 1

    centers = load_centers_csv(args.centers_csv)
    cfg = SET_CONFIG[args.set]
    default_input, csv_path = resolve_set_paths(args.set, args)
    input_paths = args.inputs if args.inputs else [default_input]
    existing_heights = load_existing_preserved(csv_path)

    files = find_inputs(input_paths, cfg["extensions"])
    if args.only:
        files = [p for p in files if p.name == args.only or p.stem == args.only]
    if not files:
        print(f"No input files found for set '{args.set}' in {input_paths}", file=sys.stderr)
        return 1

    if str(_METRICAS_DIR) not in sys.path:
        sys.path.insert(0, str(_METRICAS_DIR))
    from superficie_copa import measure_crown_surface

    rows: list[dict[str, str | float]] = []
    ok_count = 0
    for path in files:
        tree_id = tree_id_from_path(path)
        print(f"Processing [{args.set}] {path} ...", flush=True)
        center = centers.get(tree_id)
        has_center = center is not None and center.get("x") is not None and center.get("y") is not None

        try:
            verts = load_points(path)
            xy = verts[:, :2]
            surface = measure_crown_surface(xy)

            if has_center:
                cx, cy = center["x"], center["y"]
                crown = measure_crown_radii(verts, cx=cx, cy=cy)
            else:
                crown = {
                    "r_40_m": 0.0,
                    "r_130_m": 0.0,
                    "r_220_m": 0.0,
                    "r_310_m": 0.0,
                }

            m = {
                **crown,
                "superficie_copa_m2": surface["area_grid_010_m2"],
                "altura_m": "",
            }
            row = metrics_row(path, m)
            tid = row["id_arbol"]
            if tid in existing_heights:
                row.update(existing_heights[tid])
            rows.append(row)
            ok_count += 1
            if has_center:
                print(
                    f"  centro=({cx:.2f}, {cy:.2f})  "
                    f"r40={row['r_40_m']} r130={row['r_130_m']} "
                    f"r220={row['r_220_m']} r310={row['r_310_m']}  "
                    f"superficie={row['superficie_copa_m2']}",
                    flush=True,
                )
            else:
                print(
                    f"  sin centro CSV  superficie={row['superficie_copa_m2']} radios=0",
                    flush=True,
                )
        except Exception as exc:  # noqa: BLE001
            msg = str(exc)
            print(f"  ERROR: {msg}", file=sys.stderr, flush=True)
            rows.append(error_row(path, msg))

    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CSV_FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {csv_path} ({ok_count}/{len(rows)} ok)", flush=True)
    return 0 if ok_count else 2


if __name__ == "__main__":
    raise SystemExit(main())
