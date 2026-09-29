"""Crown projected area (m²) from full point cloud, zenithal view.

Methods: convex hull, alpha shape (Delaunay triangles), and occupied grid.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy.spatial import ConvexHull, Delaunay, QhullError, cKDTree

GRID_RES_M = 0.10
ALPHA_SAMPLE_MAX = 15000
DRAW_POINTS_MAX = 180000


def _load_fonts() -> tuple:
    try:
        return (
            ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf", 15),
            ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans.ttf", 11),
            ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSans.ttf", 10),
        )
    except OSError:
        d = ImageFont.load_default()
        return d, d, d


def area_grid(xy: np.ndarray, res: float = GRID_RES_M) -> tuple[float, np.ndarray, float, float]:
    x0, y0 = float(xy[:, 0].min()), float(xy[:, 1].min())
    ix = np.floor((xy[:, 0] - x0) / res).astype(np.int64)
    iy = np.floor((xy[:, 1] - y0) / res).astype(np.int64)
    cells = np.unique(np.column_stack([ix, iy]), axis=0)
    return len(cells) * res * res, cells, x0, y0


def _alpha_edges_and_area(points: np.ndarray, alpha: float) -> tuple[float, set[tuple[int, int]]]:
    tri = Delaunay(points)
    edges: set[tuple[int, int]] = set()
    tri_area = 0.0

    def add_edge(i: int, j: int) -> None:
        e = (min(i, j), max(i, j))
        rev = (e[1], e[0])
        if e in edges:
            edges.remove(e)
        elif rev in edges:
            edges.remove(rev)
        else:
            edges.add(e)

    for simplex in tri.simplices:
        pa, pb, pc = points[simplex[0]], points[simplex[1]], points[simplex[2]]
        a = float(np.linalg.norm(pa - pb))
        b = float(np.linalg.norm(pb - pc))
        c = float(np.linalg.norm(pc - pa))
        s = (a + b + c) / 2.0
        area_sq = max(s * (s - a) * (s - b) * (s - c), 0.0)
        circum_r = (a * b * c) / (4.0 * area_sq) if area_sq > 0 else float("inf")
        if circum_r < alpha:
            tri_area += float(area_sq**0.5)
            add_edge(int(simplex[0]), int(simplex[1]))
            add_edge(int(simplex[1]), int(simplex[2]))
            add_edge(int(simplex[2]), int(simplex[0]))
    return tri_area, edges


def _edges_to_polygon(points: np.ndarray, edges: set[tuple[int, int]]) -> np.ndarray:
    if not edges:
        return np.empty((0, 2))
    adj: dict[int, list[int]] = {}
    for i, j in edges:
        adj.setdefault(i, []).append(j)
        adj.setdefault(j, []).append(i)
    start = min(adj)
    loop = [start]
    prev, cur = -1, start
    for _ in range(len(edges) + 2):
        nbs = [n for n in adj[cur] if n != prev]
        if not nbs:
            break
        nxt = nbs[0]
        if nxt == start and len(loop) > 2:
            break
        loop.append(nxt)
        prev, cur = cur, nxt
    return points[np.array(loop)]


def _polygon_area(poly: np.ndarray) -> float:
    if len(poly) < 3:
        return float("nan")
    x, y = poly[:, 0], poly[:, 1]
    return float(0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))))


def _stratified_sample(xy: np.ndarray, max_pts: int, cell: float) -> np.ndarray:
    ix = np.floor((xy[:, 0] - xy[:, 0].min()) / cell).astype(np.int64)
    iy = np.floor((xy[:, 1] - xy[:, 1].min()) / cell).astype(np.int64)
    keys = ix * 1_000_003 + iy
    _, idx = np.unique(keys, return_index=True)
    if len(idx) <= max_pts:
        return xy[idx]
    pick = idx[np.linspace(0, len(idx) - 1, max_pts, dtype=int)]
    return xy[pick]


def _auto_alpha_radius(sample: np.ndarray, grid_res: float) -> float:
    nn_d = cKDTree(sample).query(sample, k=min(2, len(sample)))[0]
    spacing = float(np.median(nn_d if nn_d.ndim == 1 else nn_d[:, -1]))
    return float(np.clip(max(grid_res * 45.0, spacing * 25.0), 3.5, 15.0))


def _point_in_triangle(p: np.ndarray, a: np.ndarray, b: np.ndarray, c: np.ndarray) -> bool:
    v0, v1, v2 = c - a, b - a, p - a
    dot00 = float(np.dot(v0, v0))
    dot01 = float(np.dot(v0, v1))
    dot02 = float(np.dot(v0, v2))
    dot11 = float(np.dot(v1, v1))
    dot12 = float(np.dot(v1, v2))
    denom = dot00 * dot11 - dot01 * dot01
    if abs(denom) < 1e-12:
        return False
    inv = 1.0 / denom
    u = (dot11 * dot02 - dot01 * dot12) * inv
    v = (dot00 * dot12 - dot01 * dot02) * inv
    return u >= -1e-9 and v >= -1e-9 and (u + v) <= 1.0 + 1e-9


def _alpha_occupancy_cells(sample: np.ndarray, alpha: float, res: float) -> tuple[np.ndarray, float, float]:
    tri = Delaunay(sample)
    x0 = float(sample[:, 0].min())
    y0 = float(sample[:, 1].min())
    occupied: set[tuple[int, int]] = set()
    for simplex in tri.simplices:
        pa, pb, pc = sample[simplex[0]], sample[simplex[1]], sample[simplex[2]]
        a = float(np.linalg.norm(pa - pb))
        b = float(np.linalg.norm(pb - pc))
        c = float(np.linalg.norm(pc - pa))
        s = (a + b + c) / 2.0
        area_sq = max(s * (s - a) * (s - b) * (s - c), 0.0)
        circum_r = (a * b * c) / (4.0 * area_sq) if area_sq > 0 else float("inf")
        if circum_r >= alpha:
            continue
        xmin = float(min(pa[0], pb[0], pc[0]))
        xmax = float(max(pa[0], pb[0], pc[0]))
        ymin = float(min(pa[1], pb[1], pc[1]))
        ymax = float(max(pa[1], pb[1], pc[1]))
        ix0 = int(np.floor((xmin - x0) / res))
        ix1 = int(np.floor((xmax - x0) / res))
        iy0 = int(np.floor((ymin - y0) / res))
        iy1 = int(np.floor((ymax - y0) / res))
        for ix in range(ix0, ix1 + 1):
            for iy in range(iy0, iy1 + 1):
                cx = x0 + (ix + 0.5) * res
                cy = y0 + (iy + 0.5) * res
                if _point_in_triangle(np.array([cx, cy]), pa, pb, pc):
                    occupied.add((ix, iy))
    cells = np.array(sorted(occupied), dtype=np.int64) if occupied else np.empty((0, 2), dtype=np.int64)
    return cells, x0, y0


def area_convex_hull(xy: np.ndarray) -> tuple[float, np.ndarray]:
    if len(xy) < 3:
        return float("nan"), np.empty((0, 2))
    try:
        hull = ConvexHull(xy)
    except QhullError:
        return float("nan"), np.empty((0, 2))
    poly = xy[hull.vertices]
    area = _polygon_area(poly)
    if math.isnan(area) or area <= 0:
        area = float(hull.volume)
    return area, poly


def area_alpha_shape(xy: np.ndarray) -> tuple[float, float, np.ndarray, float, float]:
    sample = _stratified_sample(xy, ALPHA_SAMPLE_MAX, GRID_RES_M)
    alpha = _auto_alpha_radius(sample, GRID_RES_M)
    tri_area, edges = _alpha_edges_and_area(sample, alpha)
    poly = _edges_to_polygon(sample, edges)
    area = tri_area
    if area <= 0 and len(poly) >= 3:
        area = _polygon_area(poly)
    alpha_cells, ax0, ay0 = _alpha_occupancy_cells(sample, alpha, GRID_RES_M)
    return area, alpha, alpha_cells, ax0, ay0


def measure_crown_surface(xy: np.ndarray) -> dict:
    a_grid, grid_cells, grid_x0, grid_y0 = area_grid(xy, GRID_RES_M)
    a_alpha, alpha_param, alpha_cells, alpha_x0, alpha_y0 = area_alpha_shape(xy)
    a_convex, convex_poly = area_convex_hull(xy)
    return {
        "area_convex_m2": a_convex,
        "area_alpha_m2": a_alpha,
        "area_grid_010_m2": a_grid,
        "alpha_param_m": alpha_param,
        "grid_res_m": GRID_RES_M,
        "_grid_cells": grid_cells,
        "_grid_origin": (grid_x0, grid_y0),
        "_alpha_cells": alpha_cells,
        "_alpha_origin": (alpha_x0, alpha_y0),
        "_convex_hull_poly": convex_poly,
    }


def _draw_surface_figure(
    xy: np.ndarray,
    *,
    tree_id: str,
    method: str,
    area_m2: float,
    extra: str,
    out_path: Path,
    center: tuple[float, float] | None,
    fill_rgba: tuple[int, int, int, int],
    line_rgb: tuple[int, int, int],
    grid_cells: np.ndarray,
    grid_origin: tuple[float, float],
    grid_res: float,
) -> None:
    cx_ref = center[0] if center else float(xy[:, 0].mean())
    cy_ref = center[1] if center else float(xy[:, 1].mean())
    r = float(np.hypot(xy[:, 0] - cx_ref, xy[:, 1] - cy_ref).max()) * 1.08
    xmin, span = cx_ref - r, 2 * r

    plot, band, side = 920, 70, 20
    w, h = plot + 2 * side, plot + 2 * band
    px0, py0 = side, band

    def px(xv: float) -> float:
        return px0 + (xv - xmin) / span * plot

    def py(yv: float) -> float:
        return py0 + plot - (yv - (cy_ref - r)) / span * plot

    font_b, font_s, _ = _load_fonts()
    img = Image.new("RGBA", (w, h), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.text((w / 2, 14), f"{tree_id} — superficie copa ({method})", fill=(33, 37, 41, 255), font=font_b, anchor="mt")
    draw.text(
        (w / 2, 36),
        f"nube completa · vista cenital · {extra} · área = {area_m2:.2f} m²",
        fill=(108, 117, 125, 255),
        font=font_s,
        anchor="mt",
    )
    draw.rectangle((px0, py0, px0 + plot - 1, py0 + plot - 1), fill=(250, 250, 252, 255), outline=(210, 210, 215, 255))

    idx = np.linspace(0, len(xy) - 1, min(len(xy), DRAW_POINTS_MAX), dtype=int)
    for xv, yv in xy[idx]:
        draw.point((int(px(xv)), int(py(yv))), fill=(70, 115, 80, 200))

    x0, y0 = grid_origin
    overlay = Image.new("RGBA", (w, h), (0, 0, 0, 0))
    od = ImageDraw.Draw(overlay, "RGBA")
    for cx, cy in grid_cells:
        x1 = x0 + float(cx) * grid_res
        y1 = y0 + float(cy) * grid_res
        x2, y2 = x1 + grid_res, y1 + grid_res
        od.rectangle((px(x1), py(y2), px(x2), py(y1)), fill=fill_rgba, outline=(*line_rgb, 180))
    img = Image.alpha_composite(img, overlay)
    draw = ImageDraw.Draw(img, "RGBA")

    if center:
        uc, vc = px(center[0]), py(center[1])
        draw.ellipse((uc - 6, vc - 6, uc + 6, vc + 6), fill=(244, 162, 97, 255), outline=(33, 37, 41, 255), width=2)

    ref = r * 0.09
    draw.line((px(cx_ref), py(cy_ref), px(cx_ref), py(cy_ref + ref)), fill=(29, 53, 87, 255), width=2)
    draw.line((px(cx_ref), py(cy_ref), px(cx_ref + ref), py(cy_ref)), fill=(29, 53, 87, 255), width=2)
    draw.text((px(cx_ref), py(cy_ref + ref) - 8), "N", fill=(29, 53, 87, 255), font=font_s, anchor="mm")
    draw.text((px(cx_ref + ref) + 8, py(cy_ref)), "E", fill=(29, 53, 87, 255), font=font_s, anchor="lm")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    img.convert("RGB").save(out_path)


def _draw_convex_figure(
    xy: np.ndarray,
    *,
    tree_id: str,
    area_m2: float,
    hull_poly: np.ndarray,
    out_path: Path,
    center: tuple[float, float] | None,
) -> None:
    cx_ref = center[0] if center else float(xy[:, 0].mean())
    cy_ref = center[1] if center else float(xy[:, 1].mean())
    r = float(np.hypot(xy[:, 0] - cx_ref, xy[:, 1] - cy_ref).max()) * 1.08
    xmin, span = cx_ref - r, 2 * r

    plot, band, side = 920, 70, 20
    w, h = plot + 2 * side, plot + 2 * band
    px0, py0 = side, band

    def px(xv: float) -> float:
        return px0 + (xv - xmin) / span * plot

    def py(yv: float) -> float:
        return py0 + plot - (yv - (cy_ref - r)) / span * plot

    font_b, font_s, _ = _load_fonts()
    img = Image.new("RGBA", (w, h), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img, "RGBA")
    draw.text(
        (w / 2, 14),
        f"{tree_id} — superficie copa (convex hull)",
        fill=(33, 37, 41, 255),
        font=font_b,
        anchor="mt",
    )
    draw.text(
        (w / 2, 36),
        f"nube completa · vista cenital · envolvente convexa · área = {area_m2:.2f} m²",
        fill=(108, 117, 125, 255),
        font=font_s,
        anchor="mt",
    )
    draw.rectangle((px0, py0, px0 + plot - 1, py0 + plot - 1), fill=(250, 250, 252, 255), outline=(210, 210, 215, 255))

    idx = np.linspace(0, len(xy) - 1, min(len(xy), DRAW_POINTS_MAX), dtype=int)
    for xv, yv in xy[idx]:
        draw.point((int(px(xv)), int(py(yv))), fill=(70, 115, 80, 200))

    if len(hull_poly) >= 3:
        overlay = Image.new("RGBA", (w, h), (0, 0, 0, 0))
        od = ImageDraw.Draw(overlay, "RGBA")
        pts = [(px(float(xv)), py(float(yv))) for xv, yv in hull_poly]
        od.polygon(pts, fill=(230, 111, 81, 45), outline=(230, 111, 81, 220))
        for i in range(len(pts)):
            od.line([pts[i], pts[(i + 1) % len(pts)]], fill=(230, 111, 81, 255), width=2)
        img = Image.alpha_composite(img, overlay)
        draw = ImageDraw.Draw(img, "RGBA")

    if center:
        uc, vc = px(center[0]), py(center[1])
        draw.ellipse((uc - 6, vc - 6, uc + 6, vc + 6), fill=(244, 162, 97, 255), outline=(33, 37, 41, 255), width=2)

    ref = r * 0.09
    draw.line((px(cx_ref), py(cy_ref), px(cx_ref), py(cy_ref + ref)), fill=(29, 53, 87, 255), width=2)
    draw.line((px(cx_ref), py(cy_ref), px(cx_ref + ref), py(cy_ref)), fill=(29, 53, 87, 255), width=2)
    draw.text((px(cx_ref), py(cy_ref + ref) - 8), "N", fill=(29, 53, 87, 255), font=font_s, anchor="mm")
    draw.text((px(cx_ref + ref) + 8, py(cy_ref)), "E", fill=(29, 53, 87, 255), font=font_s, anchor="lm")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    img.convert("RGB").save(out_path)


def render_surface_figures(
    xy: np.ndarray,
    surface: dict,
    *,
    tree_id: str,
    center: tuple[float, float] | None,
    figuras_dir: Path,
) -> None:
    _draw_surface_figure(
        xy,
        tree_id=tree_id,
        method="alpha shape",
        area_m2=surface["area_alpha_m2"],
        extra=f"alpha ≈ {surface['alpha_param_m']:.2f} m · muestra estratificada ≤{ALPHA_SAMPLE_MAX}",
        out_path=figuras_dir / f"{tree_id}_superficie_alpha.png",
        center=center,
        fill_rgba=(69, 123, 157, 55),
        line_rgb=(69, 123, 157),
        grid_cells=surface["_alpha_cells"],
        grid_origin=surface["_alpha_origin"],
        grid_res=GRID_RES_M,
    )
    _draw_surface_figure(
        xy,
        tree_id=tree_id,
        method=f"rejilla {GRID_RES_M:.2f} m",
        area_m2=surface["area_grid_010_m2"],
        extra=f"celdas ocupadas × {GRID_RES_M}² · todos los puntos",
        out_path=figuras_dir / f"{tree_id}_superficie_grid.png",
        center=center,
        fill_rgba=(42, 157, 143, 55),
        line_rgb=(42, 157, 143),
        grid_cells=surface["_grid_cells"],
        grid_origin=surface["_grid_origin"],
        grid_res=GRID_RES_M,
    )
    _draw_convex_figure(
        xy,
        tree_id=tree_id,
        area_m2=surface["area_convex_m2"],
        hull_poly=surface["_convex_hull_poly"],
        out_path=figuras_dir / f"{tree_id}_superficie_convex.png",
        center=center,
    )
