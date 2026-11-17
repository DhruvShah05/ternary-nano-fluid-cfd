"""
make_2d_mesh.py  -  Phase 2: scripted 2D axisymmetric mesh for the heated pipe.

Writes a Fluent ASCII .msh file (structured quads, conformal fluid/solid interface).
Axis of the pipe = x (y is the radial direction), so use Fluent's 2D "Axisymmetric".

Geometry (spec Sec. 3):  L = 4 m, heated for x in [2, 4],
  fluid  0      <= r <= 0.0115
  solid  0.0115 <= r <= 0.0125

Zones written (names are exact):
  cell zones : fluid, solid
  boundaries : inlet (velocity-inlet), outlet (pressure-outlet), axis,
               outer_unheated (wall), outer_heated (wall),
               solid_end_in (wall), solid_end_out (wall)
  interface  : interface_fs (interior between fluid and solid cells)
               -> must become a coupled wall pair in Fluent (checked in Phase 2 script)

Mesh levels follow spec Section 4 (M1/M2/M3).  Usage:
    python make_2d_mesh.py M2 pipe2d_M2.msh
Wall-function mesh (standard k-epsilon, standard wall functions):
    python make_2d_mesh.py M2 pipe2d_M2WF.msh --first-layer 2.2e-3 --growth 1.15 --layers 2
"""
import sys
import math
import numpy as np

R_FLUID = 0.0115
R_OUTER = 0.0125
L_TOTAL = 4.0
X_SPLIT = 2.0

LEVELS = {
    #       growth, n_layers, dx_unheated, dx_heated, n_solid
    "M1": dict(growth=1.20, n_layers=30, dx_u=8e-3,   dx_h=5e-3, n_solid=4),
    "M2": dict(growth=1.15, n_layers=35, dx_u=5e-3,   dx_h=3e-3, n_solid=5),
    "M3": dict(growth=1.12, n_layers=42, dx_u=3.5e-3, dx_h=2e-3, n_solid=7),
    # 2D twin of the planned 3D Gate B half-pipe mesh (512k-cell student limit): M1 wall layers and heated
    # spacing, but 20 mm cells in the unheated part (nothing is measured there)
    "B3D": dict(growth=1.20, n_layers=30, dx_u=20e-3, dx_h=5e-3, n_solid=4),
}
FIRST_LAYER = 5e-6
AXIAL_TRANSITION_RATIO = 1.1   # spec: <= 1.1 across x = 2


# ----------------------------------------------------------------------------
def radial_fluid(first, growth, n_layers, R):
    """Radial node positions axis->wall: uniform core, then geometric inflation at the wall."""
    infl = [first * growth ** k for k in range(n_layers)]        # wall -> inward
    t_infl = sum(infl)
    rem = R - t_infl
    if rem <= 0:
        raise ValueError("Inflation layers thicker than the pipe radius")
    cap = infl[-1]
    n_core = max(1, round(rem / cap))
    core = [rem / n_core] * n_core
    sizes_wall_to_axis = infl + core
    r = [R]
    for s in sizes_wall_to_axis:
        r.append(r[-1] - s)
    r = np.array(r[::-1])
    r[0] = 0.0
    return r


def axial_nodes(dx_u, dx_h):
    """x positions: uniform heated part, smooth (ratio 1.1) transition on the unheated side of x=2."""
    n_h = round((L_TOTAL - X_SPLIT) / dx_h)
    xh = [X_SPLIT + (L_TOTAL - X_SPLIT) * i / n_h for i in range(n_h + 1)]
    dxh = (L_TOTAL - X_SPLIT) / n_h
    # transition cells going upstream from x = 2
    trans, s = [], dxh
    while s * AXIAL_TRANSITION_RATIO < dx_u:
        s *= AXIAL_TRANSITION_RATIO
        trans.append(s)
    x_trans_start = X_SPLIT - sum(trans)
    n_u = max(1, round(x_trans_start / dx_u))
    xu = [x_trans_start * i / n_u for i in range(n_u + 1)]
    xt = [x_trans_start]
    for s in reversed(trans):
        xt.append(xt[-1] + s)
    xt[-1] = X_SPLIT
    x = np.array(xu[:-1] + xt[:-1] + xh)
    return x


def build_mesh(level="M2", first=None, growth=None, n_layers=None):
    """first / growth / n_layers override the level's wall-layer settings (used for wall-function meshes)."""
    p = dict(LEVELS[level])
    if growth is not None:
        p["growth"] = growth
    if n_layers is not None:
        p["n_layers"] = n_layers
    x = axial_nodes(p["dx_u"], p["dx_h"])
    rf = radial_fluid(FIRST_LAYER if first is None else first, p["growth"], p["n_layers"], R_FLUID)
    rs = np.linspace(R_FLUID, R_OUTER, p["n_solid"] + 1)[1:]
    y = np.concatenate([rf, rs])
    Nx, Nf, Ns = len(x) - 1, len(rf) - 1, p["n_solid"]
    return x, y, Nx, Nf, Ns


# ----------------------------------------------------------------------------
def write_msh(path, level="M2", first=None, growth=None, n_layers=None):
    x, y, Nx, Nf, Ns = build_mesh(level, first, growth, n_layers)
    Ny = Nf + Ns
    node = lambda i, j: 1 + j * (Nx + 1) + i          # 1-based node id
    nfc = Nx * Nf
    ncells = Nx * Ny

    def cell(i, j):                                    # 1-based cell id; fluid first, then solid
        return 1 + j * Nx + i                          # rows j=0..Nf-1 fluid, then solid (contiguous)

    # ---- faces, grouped by zone -------------------------------------------
    Z = {k: [] for k in ["int_f", "int_s", "interface", "inlet", "outlet", "axis",
                         "outer_unheated", "outer_heated", "solid_end_in", "solid_end_out"]}
    # internal build uses c0 on the right of (a->b); the file writes the nodes reversed (b->a), so Fluent sees c0 on the LEFT
    # axial-direction faces (edges along +x at row j)
    for j in range(Ny + 1):
        for i in range(Nx):
            a, b = node(i, j), node(i + 1, j)
            if j == 0:
                Z["axis"].append((b, a, cell(i, 0), 0))
            elif j == Ny:
                zone = "outer_unheated" if x[i + 1] <= X_SPLIT + 1e-12 else "outer_heated"
                Z[zone].append((a, b, cell(i, Ny - 1), 0))
            else:
                zone = "int_f" if j < Nf else ("interface" if j == Nf else "int_s")
                Z[zone].append((a, b, cell(i, j - 1), cell(i, j)))
    # radial-direction faces (edges along +y at column i)
    for i in range(Nx + 1):
        for j in range(Ny):
            a, b = node(i, j), node(i, j + 1)
            if i == 0:
                Z["inlet" if j < Nf else "solid_end_in"].append((a, b, cell(0, j), 0))
            elif i == Nx:
                Z["outlet" if j < Nf else "solid_end_out"].append((b, a, cell(Nx - 1, j), 0))
            else:
                Z["int_f" if j < Nf else "int_s"].append((a, b, cell(i, j), cell(i - 1, j)))

    # zone table: key, zone id, bc code, type name, display name
    ztab = [("int_f", 4, 2, "interior", "interior_fluid"),
            ("int_s", 5, 2, "interior", "interior_solid"),
            ("interface", 6, 2, "interior", "interface_fs"),
            ("inlet", 7, 10, "velocity-inlet", "inlet"),
            ("outlet", 8, 5, "pressure-outlet", "outlet"),
            ("axis", 9, 37, "axis", "axis"),
            ("outer_unheated", 10, 3, "wall", "outer_unheated"),
            ("outer_heated", 11, 3, "wall", "outer_heated"),
            ("solid_end_in", 12, 3, "wall", "solid_end_in"),
            ("solid_end_out", 13, 3, "wall", "solid_end_out")]
    nfaces = sum(len(v) for v in Z.values())
    h = lambda v: format(v, "x")

    with open(path, "w") as f:
        f.write('(0 "2D axisymmetric heated pipe, written by make_2d_mesh.py")\n')
        f.write("(2 2)\n")
        nn = (Nx + 1) * (Ny + 1)
        f.write(f"(10 (0 1 {h(nn)} 0 2))\n")
        f.write(f"(10 (1 1 {h(nn)} 1 2)(\n")
        for j in range(Ny + 1):
            for i in range(Nx + 1):
                f.write(f"{x[i]:.10e} {y[j]:.10e}\n")
        f.write("))\n")
        f.write(f"(12 (0 1 {h(ncells)} 0 0))\n")
        f.write(f"(12 (2 1 {h(nfc)} 1 3))\n")
        f.write(f"(12 (3 {h(nfc + 1)} {h(ncells)} 1 3))\n")
        f.write(f"(13 (0 1 {h(nfaces)} 0 0))\n")
        start = 1
        for key, zid, bc, _, _ in ztab:
            faces = Z[key]
            if not faces:
                continue
            f.write(f"(13 ({h(zid)} {h(start)} {h(start + len(faces) - 1)} {h(bc)} 2)(\n")
            for n0, n1, c0, c1 in faces:
                f.write(f"{h(n1)} {h(n0)} {h(c0)} {h(c1)}\n")  # reversed order: run 1 showed ALL faces left-handed with the other order
            f.write("))\n")
            start += len(faces)
        f.write("(45 (2 fluid fluid)())\n")
        f.write("(45 (3 solid solid)())\n")
        for key, zid, _, tname, dname in ztab:
            if Z[key]:
                f.write(f"(45 ({zid} {tname} {dname})())\n")
    return dict(x=x, y=y, Nx=Nx, Nf=Nf, Ns=Ns, ncells=ncells, nodes=nn, nfaces=nfaces,
                zones={k: len(v) for k, v in Z.items()}, faces=Z)


# ----------------------------------------------------------------------------
def validate(m):
    """Independent geometric checks on the generated mesh (no Fluent needed)."""
    x, y, Nx, Nf, Ns = m["x"], m["y"], m["Nx"], m["Nf"], m["Ns"]
    Ny = Nf + Ns
    node_xy = lambda n: ((n - 1) % (Nx + 1), (n - 1) // (Nx + 1))
    cen = lambda c: (0.5 * (x[(c - 1) % Nx] + x[(c - 1) % Nx + 1]),
                     0.5 * (y[(c - 1) // Nx] + y[(c - 1) // Nx + 1]))
    bad = 0
    cell_faces = np.zeros(m["ncells"] + 1, int)
    for key, faces in m["faces"].items():
        for n0, n1, c0, c1 in faces:
            (i0, j0), (i1, j1) = node_xy(n1), node_xy(n0)   # as written to file
            dx, dy = x[i1] - x[i0], y[j1] - y[j0]
            for c, side in ((c0, +1), (c1, -1)):       # c0 on LEFT (cross>0), c1 on RIGHT (cross<0)
                if c == 0:
                    continue
                cx, cy = cen(c)
                cross = dx * (cy - y[j0]) - dy * (cx - x[i0])
                if cross * side <= 0:
                    bad += 1
                cell_faces[c] += 1
    ok_counts = bool(np.all(cell_faces[1:] == 4))
    # areas and axisymmetric volume
    dxs, dys = np.diff(x), np.diff(y)
    area = np.outer(dys, dxs).sum()
    vol = (np.outer(0.5 * (y[1:] + y[:-1]) * dys, dxs)).sum() * 2 * math.pi
    exact = math.pi * R_OUTER ** 2 * L_TOTAL
    ar = np.outer(dys, 1.0 / dxs)
    aspect = np.maximum(ar, 1.0 / ar)
    return dict(orientation_errors=bad, every_cell_has_4_faces=ok_counts,
                volume_error_pct=100 * (vol - exact) / exact,
                max_aspect=float(aspect.max()), first_layer_um=1e6 * (y[Nf] - y[Nf - 1]),
                x_split_present=bool(np.any(np.isclose(x, X_SPLIT))),
                max_axial_ratio=float(np.max(np.maximum(dxs[1:] / dxs[:-1], dxs[:-1] / dxs[1:]))),
                max_radial_ratio_fluid=float(np.max(np.maximum(dys[1:Nf] / dys[:Nf - 1], dys[:Nf - 1] / dys[1:Nf]))) if Nf > 1 else 1.0,
                min_dx_mm=1e3 * dxs.min(), max_dx_mm=1e3 * dxs.max())


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Write a 2D axisymmetric pipe mesh (Fluent .msh).")
    ap.add_argument("level", nargs="?", default="M2", choices=sorted(LEVELS))
    ap.add_argument("out", nargs="?", default=None)
    ap.add_argument("--first-layer", type=float, default=None,
                    help="first fluid cell height at the wall in metres (default 5e-6). "
                         "For STANDARD WALL FUNCTIONS use about 2.2e-3 (y+ roughly 30 to 70 at Re 10k to 25k)")
    ap.add_argument("--growth", type=float, default=None, help="wall-layer growth ratio override")
    ap.add_argument("--layers", type=int, default=None, help="number of wall layers override")
    a = ap.parse_args()
    out = a.out or f"pipe2d_{a.level}.msh"
    m = write_msh(out, a.level, a.first_layer, a.growth, a.layers)
    v = validate(m)
    print(f"Level {a.level}: {m['ncells']} cells ({m['Nx']} axial x {m['Nf']} fluid + {m['Ns']} solid radial), "
          f"{m['nodes']} nodes, {m['nfaces']} faces -> {out}")
    print("Zone face counts:", m["zones"])
    for k, val in v.items():
        print(f"  {k}: {val}")
    print("  radial fluid cell heights (mm, wall -> axis):",
          [round(float(1e3 * d), 3) for d in np.diff(m["y"][:m["Nf"] + 1])[::-1]])
