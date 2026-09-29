"""GeoTIFF DEM grid for bilinear elevation sampling."""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np

PROJECT = Path(__file__).resolve().parent.parent
DEFAULT_DEM = PROJECT / "data/sceneparts/tif/DEM_Ruiles_2m.tif"


class DemGrid:
    """GeoTIFF DEM as a regular grid for bilinear sampling."""

    def __init__(self, array: np.ndarray, gt: tuple[float, float, float, float, float, float]):
        self.array = array
        self.gt = gt
        self.ny, self.nx = array.shape

    @classmethod
    def from_tif(cls, dem_path: Path) -> "DemGrid":
        info = subprocess.run(
            ["bash", "-lc", f"module load gdal/3.11.0-zen4-k; gdalinfo '{dem_path}'"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
        size_line = next(l for l in info.splitlines() if l.startswith("Size is "))
        nx, ny = map(int, size_line.replace("Size is ", "").split(", "))
        origin_line = next(l for l in info.splitlines() if l.startswith("Origin = "))
        ox, oy = map(float, origin_line.split("(")[1].split(")")[0].split(",")[:2])
        pix_line = next(l for l in info.splitlines() if l.startswith("Pixel Size = "))
        px, py = map(float, pix_line.split("(")[1].split(")")[0].split(",")[:2])
        gt = (ox, px, 0.0, oy, 0.0, py)

        xyz_path = dem_path.with_suffix(".tmp.xyz")
        subprocess.run(
            [
                "bash",
                "-lc",
                f"module load gdal/3.11.0-zen4-k; gdal_translate -of XYZ '{dem_path}' '{xyz_path}'",
            ],
            check=True,
        )
        pts = np.loadtxt(xyz_path)
        xyz_path.unlink(missing_ok=True)
        z = pts[:, 2].reshape(ny, nx)
        return cls(z, gt)

    def z_at(self, x: float, y: float) -> float:
        gt = self.gt
        px = (x - gt[0]) / gt[1]
        py = (y - gt[3]) / gt[5]
        if px < 0 or py < 0 or px > self.nx - 1 or py > self.ny - 1:
            raise ValueError(f"({x:.3f}, {y:.3f}) outside DEM extent")
        x0, y0 = int(np.floor(px)), int(np.floor(py))
        x1, y1 = min(x0 + 1, self.nx - 1), min(y0 + 1, self.ny - 1)
        tx, ty = px - x0, py - y0
        z00, z10 = self.array[y0, x0], self.array[y0, x1]
        z01, z11 = self.array[y1, x0], self.array[y1, x1]
        return float(
            (1 - tx) * (1 - ty) * z00
            + tx * (1 - ty) * z10
            + (1 - tx) * ty * z01
            + tx * ty * z11
        )
