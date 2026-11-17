"""
view_3d_mesh.py - draws the 3D half-pipe mesh as an SVG picture (no Fluent, no plotting library needed).

    python 02_geometry_mesh\\3d\\view_3d_mesh.py      -> 02_geometry_mesh\\3d\\pipe3d_B3D_view.svg (open in any browser)

Panels: (A) 3D view of the half pipe (length compressed), (B) the cross-section mesh,
(C) zoom on the wall: inflation layers and the 1 mm copper wall, (D) side view of the axial cell spacing.
"""
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_3d_mesh as m3  # noqa: E402

P = m3.LEVELS["B3D"]
xy, quads, zone, info = m3.cross_section(P)
x_nodes = m3.m2.axial_nodes(P["dx_u"], P["dx_h"])
R_I, R_O, L, XS = m3.R_I, m3.R_O, m3.L_TOTAL, m3.X_SPLIT

FLUID, SOLID, LINE = "#cfe3f7", "#e3a97a", "#33506e"
SECTOR = {1: "#d1495b", 2: "#edae49", 3: "#00798c", 4: "#6a4c93"}
UNHEATED = "#b8bec6"
W, H = 1100, 1240
out = []


def poly(pts, fill, stroke="none", sw=0.0, op=1.0):
    d = " ".join(f"{a:.2f},{b:.2f}" for a, b in pts)
    out.append(f'<polygon points="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}" fill-opacity="{op}"/>')


def text(x, y, s, size=14, anchor="start", weight="normal", fill="#1d2733"):
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    out.append(f'<text x="{x:.1f}" y="{y:.1f}" font-size="{size}" text-anchor="{anchor}" '
               f'font-weight="{weight}" fill="{fill}">{s}</text>')


def line(x0, y0, x1, y1, stroke=LINE, sw=0.6):
    out.append(f'<line x1="{x0:.2f}" y1="{y0:.2f}" x2="{x1:.2f}" y2="{y1:.2f}" stroke="{stroke}" stroke-width="{sw}"/>')


# ---------------------------------------------------------------- (A) 3D view
def panel_3d(ox, oy, w, h):
    text(ox, oy - 12, "A  3D half pipe (length compressed 40x; real pipe is 4 m long, 25 mm across)", 15, weight="bold")
    stretch = 1.0 / 40.0                                   # visual compression of x
    az, el = math.radians(-35), math.radians(25)

    def proj(x, y, z):
        X = x * stretch
        # rotate about the vertical (y) axis, then tilt
        xr = X * math.cos(az) + z * math.sin(az)
        zr = -X * math.sin(az) + z * math.cos(az)
        yr = y * math.cos(el) - zr * math.sin(el)
        depth = y * math.sin(el) + zr * math.cos(el)
        return xr, yr, depth

    ncirc = info["ncirc"]
    th = np.linspace(0, math.pi, ncirc + 1)
    polys = []
    xs = list(np.linspace(0, XS, 9)) + list(np.linspace(XS, L, 9))[1:]
    for a in range(ncirc):
        sec = a // P["nsec"] + 1
        for b in range(len(xs) - 1):
            x0, x1 = xs[b], xs[b + 1]
            col = UNHEATED if x1 <= XS + 1e-9 else SECTOR[sec]
            pts3 = [(x0, R_O * math.cos(th[a]), R_O * math.sin(th[a])), (x1, R_O * math.cos(th[a]), R_O * math.sin(th[a])),
                    (x1, R_O * math.cos(th[a + 1]), R_O * math.sin(th[a + 1])),
                    (x0, R_O * math.cos(th[a + 1]), R_O * math.sin(th[a + 1]))]
            pr = [proj(*p) for p in pts3]
            polys.append((np.mean([p[2] for p in pr]), [(p[0], p[1]) for p in pr], col, 0.25))
    # inlet end face (x = 0): fluid + solid quads
    for q, zn in zip(quads, zone):
        pr = [proj(0.0, xy[n][0], xy[n][1]) for n in q]
        polys.append((min(p[2] for p in pr) - 1.0, [(p[0], p[1]) for p in pr], FLUID if zn == 0 else SOLID, 0.15))
    # scale to the panel
    allp = np.array([p for _, ps, _, _ in polys for p in ps])
    mn, mx = allp.min(0), allp.max(0)
    sc = min((w - 20) / (mx[0] - mn[0]), (h - 20) / (mx[1] - mn[1]))
    tr = lambda p: (ox + 10 + (p[0] - mn[0]) * sc, oy + h - 10 - (p[1] - mn[1]) * sc)
    for depth, ps, col, sw in sorted(polys, key=lambda t: -t[0]):
        poly([tr(p) for p in ps], col, LINE, sw)
    # labels
    a_in = tr(proj(0, R_O * 1.25, 0)[:2])
    text(a_in[0] - 10, a_in[1] - 6, "inlet (x = 0)", 13, "end")
    a_mid = tr(proj(XS, R_O * 1.1, R_O * 0.2)[:2])
    text(a_mid[0], a_mid[1] - 34, "x = 2 m: heating starts", 13, "middle")
    a_out = tr(proj(L, R_O * 1.2, 0)[:2])
    text(a_out[0] + 6, a_out[1] - 4, "outlet (x = 4 m)", 13)
    text(ox + w / 2, oy + h + 16, "flat underside = symmetry plane z = 0 (other half is the mirror image)", 13, "middle")
    lx, ly = ox + w - 230, oy + 8
    for i, (lab, col) in enumerate([("unheated wall, x < 2 m", UNHEATED)] +
                                   [(f"heated sector q{k} (45 deg)", SECTOR[k]) for k in range(1, 5)] +
                                   [("fluid (water)", FLUID), ("copper wall, 1 mm", SOLID)]):
        poly([(lx, ly + 20 * i), (lx + 14, ly + 20 * i), (lx + 14, ly + 20 * i + 14), (lx, ly + 20 * i + 14)], col, LINE, 0.5)
        text(lx + 20, ly + 20 * i + 12, lab, 13)


# ---------------------------------------------------------------- (B) cross-section
def panel_section(ox, oy, size):
    text(ox, oy - 12, f"B  Cross-section mesh: {len(quads)} cells (O-grid, 6 cells per 45 deg sector)", 15, weight="bold")
    sc = (size - 40) / (2 * R_O)
    cx, cy = ox + size / 2, oy + size * 0.78
    tr = lambda p: (cx + p[0] * sc, cy - p[1] * sc)
    for q, zn in zip(quads, zone):
        poly([tr(xy[n]) for n in q], FLUID if zn == 0 else SOLID, LINE, 0.35)
    for k in range(1, 5):
        t = math.radians(45 * (k - 1) + 22.5)
        p = tr((1.12 * R_O * math.cos(t), 1.12 * R_O * math.sin(t)))
        text(p[0], p[1] + 5, f"q{k}", 15, "middle", "bold", SECTOR[k])
        t0 = math.radians(45 * k)
        a, b = tr((R_O * math.cos(t0), R_O * math.sin(t0))), tr((1.06 * R_O * math.cos(t0), 1.06 * R_O * math.sin(t0)))
        line(*a, *b, "#1d2733", 1.4)
    text(cx, cy + 26, "symmetry plane z = 0", 13, "middle")
    text(cx, cy + 44, "(y axis to the right, z axis up; heated sectors q1..q4 on the outer copper surface)", 12, "middle", fill="#5a6573")
    return tr, sc


# ---------------------------------------------------------------- (C) wall zoom
def panel_zoom(ox, oy, w, h):
    text(ox, oy - 12, "C  Zoom on the wall at the top of the pipe (2.5 mm tall window)", 15, weight="bold")
    r0, r1 = R_O - 2.5e-3, R_O + 0.1e-3
    t_half = math.radians(9.5)
    sc = h / (r1 - r0)
    cx = ox + w / 2
    tr = lambda p: (cx + p[0] * sc, oy + h - (math.hypot(*p) * 0 + p[1] - r0) * sc)
    out.append(f'<clipPath id="zc"><rect x="{ox}" y="{oy}" width="{w}" height="{h}"/></clipPath><g clip-path="url(#zc)">')
    for q, zn in zip(quads, zone):
        pts = [xy[n] for n in q]
        if max(math.hypot(*p) for p in pts) < r0 or max(abs(math.atan2(p[0], p[1])) for p in pts) > t_half:
            continue
        poly([tr(p) for p in pts], FLUID if zn == 0 else SOLID, LINE, 0.4)
    out.append("</g>")
    out.append(f'<rect x="{ox}" y="{oy}" width="{w}" height="{h}" fill="none" stroke="#1d2733" stroke-width="1"/>')
    yi = tr((0, R_I))[1]
    yo = tr((0, R_O))[1]
    text(ox + w + 8, yo + 4, "outer surface (heat flux q)", 12)
    text(ox + w + 8, (yo + yi) / 2 + 4, "copper wall 1 mm, 4 cells", 12)
    text(ox + w + 8, yi + 4, "inner wall: first cell 5 um", 12)
    text(ox + w + 8, yi + 20, "22 layers growing x1.2, then 8", 12)
    text(ox + w + 8, yi + 36, "transition cells into the core", 12)
    text(ox + w + 8, yi + 52, "(y+ about 0.07: no wall functions)", 12, fill="#5a6573")


# ---------------------------------------------------------------- (D) axial spacing
def panel_axial(ox, oy, w, h):
    nax = len(x_nodes) - 1
    text(ox, oy - 12, f"D  Side view, cells along the pipe: {nax} (20 mm unheated, 5 mm heated)", 15, weight="bold")
    sc = w / L
    poly([(ox, oy), (ox + XS * sc, oy), (ox + XS * sc, oy + h), (ox, oy + h)], "#eef1f4", "#1d2733", 0.8)
    poly([(ox + XS * sc, oy), (ox + w, oy), (ox + w, oy + h), (ox + XS * sc, oy + h)], "#fbe3d6", "#1d2733", 0.8)
    for xv in x_nodes:
        line(ox + xv * sc, oy, ox + xv * sc, oy + h, LINE, 0.35)
    for xv, lab in ((0, "0"), (1, "1 m"), (2, "2 m"), (3, "3 m"), (4, "4 m")):
        text(ox + xv * sc, oy + h + 16, lab, 12, "middle")
    poly([(ox + 3.0 * sc, oy - 4), (ox + 3.85 * sc, oy - 4), (ox + 3.85 * sc, oy), (ox + 3.0 * sc, oy)], "#2f855a")
    text(ox + 3.425 * sc, oy + h + 32, "Nu evaluation window 3.0 to 3.85 m", 12, "middle", fill="#2f855a")
    text(ox + 1.0 * sc, oy + h + 32, "unheated: flow develops", 12, "middle", fill="#5a6573")


out.append(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" '
           'font-family="Segoe UI, Arial, sans-serif">')
out.append(f'<rect width="{W}" height="{H}" fill="#ffffff"/>')
text(30, 36, f"3D half-pipe mesh B3D: {len(quads) * (len(x_nodes) - 1):,} cells (student limit 512,000)", 20, weight="bold")
panel_3d(30, 80, 1040, 390)
panel_section(30, 540, 520)
panel_zoom(600, 560, 220, 360)
panel_axial(40, 1120, 1020, 50)
out.append("</svg>")
dst = HERE / "pipe3d_B3D_view.svg"
dst.write_text("\n".join(out), encoding="utf-8")
print("written", dst)
