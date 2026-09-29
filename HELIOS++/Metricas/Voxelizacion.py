#!/usr/bin/env python3
"""Voxelización de nubes LAS/LAZ → archivos .vox para HELIOS++ (detailedvoxels).

Resumen
-------
Convierte cada nube de puntos de un árbol en una rejilla 3D de voxels ocupados
y escribe el formato ``VOXEL SPACE`` que carga el filtro ``detailedvoxels`` de
HELIOS++ (escena Ruil_MLS: ``data/sceneparts/vox/30[0-9][0-9].vox``).

Entrada por defecto: ``Metricas/Nubes Reales/{id}.las`` (Z absoluta, alineada DEM).
Salida por defecto: ``data/sceneparts/vox/{id}.vox``.

Algoritmo (por árbol)
---------------------
1. Leer todos los puntos (X, Y, Z) en float64.
2. Origen de la rejilla: ``min(x,y,z) - vsize/2`` (convención tipo Open3D).
3. Índice de celda: ``trunc((punto - origen) / vsize)`` → enteros (i, j, k).
4. Celdas únicas: si varios puntos caen en la misma celda, cuenta una sola.
5. Escribir .vox: cabecera + una fila por celda ocupada con ``PadBVTotal=1``
   (voxel totalmente opaco / transmisión 1 en HELIOS).

Parámetros CLI
--------------
--vsize   Tamaño del cubo (m). Producción Ruil_MLS: 0.035 (3.5 cm).
--pad     PadBVTotal por voxel (default 1.0).
--input   Carpeta con LAS/LAZ.
--output  Carpeta destino .vox.
--only    IDs específicos, ej. ``--only 3001 3043``.

Librerías
---------
- ``laspy``   Lectura LAS/LAZ (coordenadas x, y, z).
- ``numpy``   Arrays, min/max, índices de voxel, ``unique``.
- ``pathlib``, ``argparse``, ``sys``, ``time`` (stdlib).

Uso
---
    python3 Metricas/Voxelizacion.py
    python3 Metricas/Voxelizacion.py --vsize 0.035 --only 3001
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import laspy
import numpy as np

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_INPUT = PROJECT / "Metricas/Nubes Reales"
DEFAULT_OUTPUT = PROJECT / "data/sceneparts/vox"

HEADER_COLS = (
    "i j k PadBVTotal angleMean bsEntering bsIntercepted bsPotential "
    "ground_distance lMeanTotal lgTotal nbEchos nbSampling transmittance "
    "attenuation attenuationBiasCorrection"
)


def read_xyz(las_path: Path) -> np.ndarray:
    """Load point coordinates as N×3 float64 array (x, y, z)."""
    pc = laspy.read(las_path)
    return np.column_stack(
        (
            np.asarray(pc.x, dtype=np.float64),
            np.asarray(pc.y, dtype=np.float64),
            np.asarray(pc.z, dtype=np.float64),
        )
    )


def voxelize_to_vox(
    xyz: np.ndarray, vsize: float, out_path: Path, pad: float = 1.0
) -> dict:
    """Discretize points into voxels and write HELIOS .vox file."""
    las_min = xyz.min(axis=0)
    origin = las_min - vsize / 2.0
    idx = np.trunc((xyz - origin) / vsize).astype(np.int64)
    unique_idx = np.unique(idx, axis=0)

    min_corner = origin
    max_corner = origin + (unique_idx.max(axis=0) + 1) * vsize
    split = unique_idx.max(axis=0) + 1

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fout:
        fout.write("VOXEL SPACE\n")
        fout.write(
            f"#min_corner: {min_corner[0]:.4f} {min_corner[1]:.4f} {min_corner[2]:.4f}\n"
        )
        fout.write(
            f"#max_corner: {max_corner[0]:.4f} {max_corner[1]:.4f} {max_corner[2]:.4f}\n"
        )
        fout.write(f"#split: {int(split[0])} {int(split[1])} {int(split[2])}\n")
        fout.write(f"#res: {vsize:.4f}\n")
        fout.write(f"{HEADER_COLS}\n")
        for i, j, k in unique_idx:
            fout.write(
                f"{int(i)} {int(j)} {int(k)} {pad:g} 0 0 0 0 0 0 0 0 0 0 0 0\n"
            )

    return {
        "voxel_count": len(unique_idx),
        "points": len(xyz),
        "min_corner": min_corner.tolist(),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--vsize", type=float, default=0.035, help="Voxel edge length (m)")
    ap.add_argument("--pad", type=float, default=1.0, help="PadBVTotal per voxel")
    ap.add_argument("--only", nargs="*", help="Tree ids e.g. 3043 3001")
    args = ap.parse_args(argv)

    files = sorted(args.input.glob("*.las")) + sorted(args.input.glob("*.laz"))
    files = sorted({p.stem: p for p in files}.values(), key=lambda p: int(p.stem))
    if args.only:
        wanted = {s.replace(".las", "").replace(".laz", "") for s in args.only}
        files = [p for p in files if p.stem in wanted]
    if not files:
        print(f"No LAS/LAZ in {args.input}", file=sys.stderr)
        return 1

    t0 = time.time()
    print(
        f"Voxelizing {len(files)} files | vsize={args.vsize} m pad={args.pad} "
        f"-> {args.output}\n"
    )

    for n, las_path in enumerate(files, 1):
        out_path = args.output / f"{las_path.stem}.vox"
        info = voxelize_to_vox(read_xyz(las_path), args.vsize, out_path, args.pad)
        mb = out_path.stat().st_size / (1024 * 1024)
        print(
            f"[{n:02d}/{len(files)}] {las_path.stem} | "
            f"{info['points']:,} pts | {info['voxel_count']:,} voxels | "
            f"{mb:.2f} MB",
            flush=True,
        )

    print(f"\nDone in {time.time() - t0:.1f}s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
