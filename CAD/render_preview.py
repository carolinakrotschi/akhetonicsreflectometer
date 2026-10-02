"""Schattierte Vorschau der Fiber-Box (PNG) zur Kontrolle (Software-Render, kein OpenGL)."""
import numpy as np, cadquery as cq, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import fiberbox as fb

GREY, PINK, WHITE, ORANGE = (0.80, 0.80, 0.82), (0.78, 0.13, 0.42), (0.97, 0.97, 0.97), (0.95, 0.6, 0.25)

def tris(shape, tol):
    v, t = shape.tessellate(tol, 0.4)
    v = np.array([(p.x, p.y, p.z) for p in v]); return v[np.array(t)]

bx = fb.make_box().val(); inl = [w.val() for _, w in fb.make_inlays()]
h = fb.make_holder().val(); lid = fb.make_lid().val(); tray = fb.load_tray().val()

def scene(lid_on, inside, out_dir):
    parts = [(tris(bx, 0.3), GREY, 1.0)]
    # Inlays minimal nach aussen versetzt, damit sie vor der Wand liegen
    parts += [(tris(s.translate(cq.Vector(0, 0.08 * out_dir if s.BoundingBox().ymin < 1 else -0.08 * out_dir, 0)), 0.05), PINK, 1.0) for s in inl]
    if lid_on: parts.append((tris(lid, 0.3), GREY, 1.0))
    if inside:
        for cy in fb.HOLDER_CYS:
            for i in range(fb.N_HOLDERS):
                parts.append((tris(h.moved(cq.Location((fb.tray_cx, cy, fb.FLOOR + i * fb.HOLDER_PITCH))), 0.3), PINK, 1.0))
    kp = fb.make_keeper().val()
    for loc in fb.keeper_locations():
        parts.append((tris(kp.moved(loc), 0.2), PINK, 1.0))
        for i in range(fb.N_TRAYS):
            parts.append((tris(tray.moved(cq.Location((fb.tray_cx, fb.tray_y_front + fb.TRAY_ZMAX, fb.FLOOR + i * fb.TRAY_PITCH))), 0.6), WHITE, 1.0))
    return parts

def shot(name, az, el, parts):
    T = np.concatenate([p[0] for p in parts])
    base = np.concatenate([np.tile(p[1], (len(p[0]), 1)) for p in parts])
    n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]); n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    L = np.array([-0.5, -0.8, 0.9]); L /= np.linalg.norm(L)
    col = np.clip(base * (0.5 + 0.5 * np.abs(n @ L))[:, None], 0, 1)
    fig = plt.figure(figsize=(12, 9)); ax = fig.add_subplot(projection="3d")
    ax.add_collection3d(Poly3DCollection(T, facecolors=col, linewidths=0))
    ax.set_xlim(0, 235); ax.set_ylim(-50, 185); ax.set_zlim(-55, 180); ax.set_box_aspect((1, 1, 1))
    ax.view_init(el, az); ax.set_axis_off()
    plt.savefig(f"preview_{name}.png", dpi=100, bbox_inches="tight"); plt.close(fig); print(name)

shot("front", -65, 18, scene(True, False, -1))
shot("rear", 115, 18, scene(True, False, +1))
shot("inside", -60, 50, scene(False, True, -1))
