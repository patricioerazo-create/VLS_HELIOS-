#!/usr/bin/env python3
"""Tree height as simple Z range: z_max - z_min (all points in cloud).

Updates altura_m in Metricas/<set>/metricas.csv.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

_METRICAS_DIR = Path(__file__).resolve().parent

SET_CONFIG: dict[str, dict] = {
    "terreno": {
        "label": "terreno",
        "default_input": _METRICAS_DIR / "Nubes Terreno",
        "output_dir": _METRICAS_DIR / "terreno",
    },
    "helios": {
        "label": "HELIOS",
        "default_input": _METRICAS_DIR / "Nubes HELIOS",
        "output_dir": _METRICAS_DIR / "helios",
    },
}

EXTENSIONS = (".las", ".laz", ".xyz", ".txt")
METRICAS_FIELDS = [
    "id_arbol",
    "altura_m",
    "superficie_copa_m2",
    "r_40_m",
    "r_130_m",
    "r_220_m",
    "r_310_m",
]


def load_las_z(path: Path) -> np.ndarray:
    import laspy

    with laspy.open(path) as reader:
        las = reader.read()
    return np.asarray(las.z, dtype=np.float64)


def load_xyz_z(path: Path) -> np.ndarray:
    zs: list[float] = []
    with path.open(encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.replace(",", " ").split()
            if len(parts) < 3:
                continue
            try:
                zs.append(float(parts[2]))
            except ValueError:
                continue
    if not zs:
        raise ValueError(f"No Z values in {path}")
    return np.asarray(zs, dtype=np.float64)


def load_z(path: Path) -> np.ndarray:
    suffix = path.suffix.lower()
    if suffix in {".las", ".laz"}:
        return load_las_z(path)
    if suffix in {".xyz", ".txt"}:
        return load_xyz_z(path)
    raise ValueError(f"Unsupported format: {path.suffix}")


def tree_id_from_path(path: Path) -> str:
    stem = path.stem
    if stem.endswith("_con_suelo"):
        return stem[: -len("_con_suelo")]
    return stem


def find_inputs(input_dir: Path) -> list[Path]:
    paths: list[Path] = []
    for ext in EXTENSIONS:
        paths.extend(input_dir.glob(f"*{ext}"))
    seen: set[str] = set()
    out: list[Path] = []
    for p in sorted(
        paths,
        key=lambda x: int(tree_id_from_path(x))
        if tree_id_from_path(x).isdigit()
        else tree_id_from_path(x),
    ):
        tid = tree_id_from_path(p)
        if tid in seen:
            continue
        seen.add(tid)
        out.append(p)
    return out


def measure_height(z: np.ndarray) -> float:
    return float(z.max() - z.min())


def load_metricas_rows(csv_path: Path) -> dict[str, dict[str, str]]:
    if not csv_path.is_file():
        return {}
    with csv_path.open(newline="", encoding="utf-8") as f:
        return {r["id_arbol"]: r for r in csv.DictReader(f) if r.get("id_arbol")}


def write_metricas(csv_path: Path, rows: list[dict[str, str | float | int]]) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=METRICAS_FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--set", choices=sorted(SET_CONFIG), required=True)
    ap.add_argument("--input", type=Path, help="Input directory override")
    ap.add_argument("-o", "--output", type=Path, help="Output metricas.csv path")
    ap.add_argument("--only", help="Process single tree id")
    args = ap.parse_args(argv)

    cfg = SET_CONFIG[args.set]
    input_dir = args.input or cfg["default_input"]
    output_csv = args.output or (cfg["output_dir"] / "metricas.csv")

    if not input_dir.is_dir():
        print(f"Input directory not found: {input_dir}", file=sys.stderr)
        return 1

    files = find_inputs(input_dir)
    if args.only:
        files = [p for p in files if tree_id_from_path(p) == args.only]
    if not files:
        print(f"No files in {input_dir}", file=sys.stderr)
        return 1

    existing = load_metricas_rows(output_csv)
    updated: dict[str, dict[str, str | float | int]] = {}
    ok = 0

    for path in files:
        tid = tree_id_from_path(path)
        try:
            altura = round(measure_height(load_z(path)), 3)
            base = dict(existing.get(tid, {}))
            base.update({"id_arbol": tid, "altura_m": altura})
            for field in METRICAS_FIELDS:
                base.setdefault(field, "")
            updated[tid] = base
            ok += 1
            print(f"[{args.set}] {tid}: altura={altura:.3f} m", flush=True)
        except Exception as exc:
            print(f"[{args.set}] {tid}: ERROR {exc}", file=sys.stderr, flush=True)
            base = dict(existing.get(tid, {}))
            base.update({"id_arbol": tid, "altura_m": ""})
            for field in METRICAS_FIELDS:
                base.setdefault(field, "")
            updated[tid] = base

    for tid, row in existing.items():
        if tid not in updated:
            updated[tid] = row

    out_rows = [
        updated[tid]
        for tid in sorted(updated.keys(), key=lambda x: int(x) if x.isdigit() else x)
    ]
    write_metricas(output_csv, out_rows)
    print(f"\nUpdated {output_csv} ({ok}/{len(files)} heights ok, {len(out_rows)} rows total)")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
