"""
make_3d_mesh.py  -  Phase 3: scripted 3D HALF-PIPE mesh (structured hex O-grid) for Gate B and the sweep.

Writes a Fluent ASCII .msh file directly (no SpaceClaim / ANSYS Meshing needed), like make_2d_mesh.py.

Geometry (spec Sec. 3): axis +X, L = 4 m, D_i = 23 mm, D_o = 25 mm, heated for x in [2, 4].
HALF MODEL (student licence: 512k cells): z >= 0 only, symmetry plane z = 0.
  Angle theta is measured in the Y-Z plane from +Y toward +Z (spec). With theta0 = 22.5 deg
  (DEVIATION from the spec default 0, allowed because gravity is off and only the circular order of
  the sectors matters, handoff decision 2) sector k spans [(k-1)*45, k*45] deg, so the half model
  contains sectors 1..4 and the symmetry plane z = 0 maps sector k onto sector 9-k
  (q8 = q1, q7 = q2, q6 = q3, q5 = q4). Valid for uniform flux and for the placeholder pattern.

Cross-section (O-grid, all quads):
  core   : rectangle y in [-a, a], z in [0, a], NC x NC/2 cells, nodes on the rays of the ring
  ring   : NCIRC cells over 180 deg (NSEC cells per 45 deg sector) x (N_INFL + N_TRANS) fluid layers
           (first layer 5 um, growth 1.2, then a geometric transition into the core) + N_SOLID solid layers
Axial: same node distribution as the 2D level B3D (20 mm unheated, 5 mm heated, ratio <= 1.1 at x = 2).

Zones written (exact names):
  cell zones : fluid, solid
  boundaries : inlet, outlet, solid_end_in, solid_end_out, outer_unheated, outer_q1..outer_q4,
               symmetry_fluid, symmetry_solid
  interface  : interface_fs (interior between fluid and solid; becomes a coupled wall pair after the
               Fluent slit command, exactly as in 2D)
A sidecar JSON (same name, .json) records the exact faceted areas so the Fluent driver can set the
inlet velocity for the intended mass flow.

Usage (PowerShell, project root, venv active):
    python 02_geometry_mesh\\3d\\make_3d_mesh.py               # level B3D -> pipe3d_B3D.msh
    python 02_geometry_mesh\\3d\\make_3d_mesh.py --tiny        # 10 axial cells, for a quick Fluent read test
"""
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "2d"))
import make_2d_mesh as m2  # noqa: E402  (axial node distribution and geometry constants)

R_I, R_O = m2.R_FLUID, m2.R_OUTER
L_TOTAL, X_SPLIT = m2.L_TOTAL, m2.X_SPLIT

LEVELS = {
    # nsec: cells per 45 deg sector; a: half-width of the core rectangle (m)
    "B3D": dict(nsec=6, a=5.0e-3, first=5e-6, growth=1.20, n_infl=22, n_trans=8, n_solid=4,
                dx_u=20e-3, dx_h=5e-3),
    # Gate D mesh study (B3D is the FINE mesh). Every spacing is scaled by r: first cell x r, growth 1.2^r,
    # layer counts / r, cells per sector / r, axial dx x r. Heated dx divides the 5 cm Nu bands exactly
    # (5.0 / 6.25 / 8.33 mm = 10 / 8 / 6 cells per band). Medium r = 1.25, coarse r = 1.667 (vs fine).
    "D3M": dict(nsec=5, a=5.0e-3, first=6.25e-6, growth=1.20 ** 1.25, n_infl=18, n_trans=6, n_solid=3,
                dx_u=25e-3, dx_h=6.25e-3),
    "D3C": dict(nsec=4, a=5.0e-3, first=8.333e-6, growth=1.20 ** (5 / 3), n_infl=13, n_trans=5, n_solid=2,
                dx_u=33.33e-3, dx_h=2.0 / 240),
}
N_SECTORS_HALF = 4

# Fluent zone table: key -> (zone id, bc type code, type name)
BC = {"interior": 2, "wall": 3, "pressure-outlet": 5, "symmetry": 7, "velocity-inlet": 10}


# ----------------------------------------------------------------------------- cross-section
def solve_ratio(first, n, span):
    """Growth ratio r with first*(1 + r + ... + r^(n-1)) = span."""
    if abs(first * n - span) < 1e-15:
        return 1.0
    lo, hi = (1.0 + 1e-12, 3.0) if span > first * n else (0.2, 1.0 - 1e-12)
    for _ in range(200):
        r = 0.5 * (lo + hi)
        s = first * (r ** n - 1.0) / (r - 1.0)
        if (s > span) == (span > first * n):
            hi = r
        else:
            lo = r
    return 0.5 * (lo + hi)


def ray_levels(rho, p):
    """Radii along one ray from the core boundary (rho) to R_O: fluid transition, inflation, solid."""
    infl = [p["first"] * p["growth"] ** k for k in range(p["n_infl"])]          # wall -> inward
    t_infl = sum(infl)
    span = (R_I - t_infl) - rho
    if span <= 0:
        raise ValueError("core too large for the inflation layers")
    first_t = infl[-1] * p["growth"]
    r = solve_ratio(first_t, p["n_trans"], span)
    trans = [first_t * r ** k for k in range(p["n_trans"])]                      # inflation edge -> core
    sizes_core_to_wall = trans[::-1] + infl[::-1]
    radii = [rho]
    for s_ in sizes_core_to_wall:
        radii.append(radii[-1] + s_)
    radii[-1] = R_I
    radii += list(np.linspace(R_I, R_O, p["n_solid"] + 1)[1:])
    return np.array(radii), r


def cross_section(p):
    """Returns node coords (N2, 2) in (y, z), quads (n, 4) CCW in the (y, z) plane, quad zone (0 fluid, 1 solid),
    and bookkeeping."""
    nsec, a = p["nsec"], p["a"]
    ncirc = N_SECTORS_HALF * nsec                       # cells around 180 deg
    nc_y, nc_z = 2 * nsec, nsec                         # core cells (y, z)
    dth = math.pi / ncirc
    theta = np.arange(ncirc + 1) * dth                  # ring rays, 0..pi
    # core node coordinates on the rays (tan / cot distribution, corners at 45 and 135 deg)
    yt = np.array([a / math.tan(theta[nsec + (nc_y - i)]) if abs(math.tan(theta[nsec + (nc_y - i)])) > 1e-14 else 0.0
                   for i in range(nc_y + 1)])          # top row, i = 0 (y=-a) .. nc_y (y=+a)
    yt[0], yt[-1], yt[nc_y // 2] = -a, a, 0.0
    zr = np.array([a * math.tan(theta[j]) for j in range(nc_z + 1)])
    zr[-1] = a
    nodes, idx = [], {}

    def add(key, yz):
        idx[key] = len(nodes)
        nodes.append(yz)

    for j in range(nc_z + 1):
        for i in range(nc_y + 1):
            add(("c", i, j), (yt[i], zr[j]))

    # perimeter index p_ (ray) -> core node key
    def core_key(pr):
        if pr <= nsec:                        # right side, theta 0..45: y = +a, z = a tan(theta)
            return ("c", nc_y, pr)
        if pr <= 3 * nsec:                    # top, theta 45..135
            return ("c", nc_y - (pr - nsec), nc_z)
        return ("c", 0, ncirc - pr)           # left side, theta 135..180

    # ring nodes: level 0 = core boundary node (shared)
    levels, ratios = [], []
    for pr in range(ncirc + 1):
        y0, z0 = nodes[idx[core_key(pr)]]
        rho = math.hypot(y0, z0)
        rad, rr = ray_levels(rho, p)
        ratios.append(rr)
        levels.append(len(rad) - 1)
        c, s_ = math.cos(theta[pr]), math.sin(theta[pr])
        for lv in range(1, len(rad)):
            add(("r", pr, lv), (rad[lv] * c, rad[lv] * s_))
    nlev = levels[0]
    n_fluid_layers = p["n_trans"] + p["n_infl"]

    def rkey(pr, lv):
        return core_key(pr) if lv == 0 else ("r", pr, lv)

    quads, zone = [], []
    for j in range(nc_z):                                  # core quads, CCW in (y, z)
        for i in range(nc_y):
            quads.append([idx[("c", i, j)], idx[("c", i + 1, j)], idx[("c", i + 1, j + 1)], idx[("c", i, j + 1)]])
            zone.append(0)
    for pr in range(ncirc):                                # ring quads
        for lv in range(nlev):
            q = [idx[rkey(pr, lv)], idx[rkey(pr, lv + 1)], idx[rkey(pr + 1, lv + 1)], idx[rkey(pr + 1, lv)]]
            quads.append(q)
            zone.append(0 if lv < n_fluid_layers else 1)
    xy = np.array(nodes)
    quads = np.array(quads)
    # make every quad CCW in (y, z)
    for k, q in enumerate(quads):
        P = xy[q]
        area2 = np.sum(P[:, 0] * np.roll(P[:, 1], -1) - np.roll(P[:, 0], -1) * P[:, 1])
        if area2 < 0:
            quads[k] = q[::-1]
    return xy, quads, np.array(zone), dict(ncirc=ncirc, n_fluid_layers=n_fluid_layers,
                                           trans_ratio_min=min(ratios), trans_ratio_max=max(ratios))


def classify_edges(xy, quads, zone):
    """2D edges: interior (fluid/solid/interface) or boundary (outer sector k / symmetry fluid|solid)."""
    emap = {}
    for c, q in enumerate(quads):
        for k in range(4):
            a, b = q[k], q[(k + 1) % 4]
            emap.setdefault((min(a, b), max(a, b)), []).append((c, a, b))
    edges = []          # (n_a, n_b, c_left, c_right_or_-1, kind)
    for (lo, hi), lst in emap.items():
        if len(lst) == 2:
            (c0, a, b), (c1, _, _) = lst
            z0, z1 = zone[c0], zone[c1]
            kind = "interface" if z0 != z1 else ("int_f" if z0 == 0 else "int_s")
            edges.append((a, b, c0, c1, kind))
        else:
            c0, a, b = lst[0]
            pa, pb = xy[a], xy[b]
            ra, rb = np.hypot(*pa), np.hypot(*pb)
            if abs(pa[1]) < 1e-12 and abs(pb[1]) < 1e-12:
                kind = "sym_f" if zone[c0] == 0 else "sym_s"
            elif abs(ra - R_O) < 1e-9 and abs(rb - R_O) < 1e-9:
                mid = 0.5 * (pa + pb)
                th = math.degrees(math.atan2(mid[1], mid[0]))
                kind = f"outer_{int(th // 45) + 1}"
            else:
                raise RuntimeError(f"unclassified boundary edge {pa} {pb}")
            edges.append((a, b, c0, -1, kind))
    return edges


# ----------------------------------------------------------------------------- 3D assembly
def build(level="B3D", tiny=False):
    p = LEVELS[level]
    xy, quads, zone, info = cross_section(p)
    x = m2.axial_nodes(p["dx_u"], p["dx_h"])
    if tiny:
        x = np.linspace(0.0, L_TOTAL, 11)
    nx = len(x) - 1
    n2 = len(xy)
    fl = np.where(zone == 0)[0]
    so = np.where(zone == 1)[0]
    nF, nS = len(fl), len(so)
    # 3D cell id (1-based): fluid cells first (all slabs), then solid
    cid2 = np.empty(len(quads), int)
    cid2[fl] = np.arange(nF)
    cid2[so] = np.arange(nS)

    def cell(q2, k):
        q2 = np.asarray(q2)
        base = np.where(zone[q2] == 0, k * nF + cid2[q2], nF * nx + k * nS + cid2[q2])
        return base + 1

    nid = lambda n, k: k * n2 + np.asarray(n) + 1            # 1-based node id
    edges = classify_edges(xy, quads, zone)
    Z = {}

    def put(key, arr):
        Z.setdefault(key, []).append(arr)

    ks = np.arange(nx)
    # lateral faces: 2D edge x axial interval
    for a, b, c0, c1, kind in edges:
        f = np.stack([nid(a, ks), nid(b, ks), nid(b, ks + 1), nid(a, ks + 1),
                      cell(np.full(nx, c0), ks), cell(np.full(nx, c1), ks) if c1 >= 0 else np.zeros(nx, int)], 1)
        if kind.startswith("outer_"):
            xm = 0.5 * (x[:-1] + x[1:])
            heated = xm > X_SPLIT
            if np.any(~heated):
                put("outer_unheated", f[~heated])
            if np.any(heated):
                put(f"outer_q{kind.split('_')[1]}", f[heated])
        else:
            put(kind, f)
    # axial faces: quad x station
    qn = quads
    for k in range(nx + 1):
        f = np.zeros((len(quads), 6), int)
        f[:, :4] = nid(qn, k)
        if k == 0:
            f[:, 4] = cell(np.arange(len(quads)), 0)
        elif k == nx:
            f[:, 4] = cell(np.arange(len(quads)), nx - 1)
        else:
            f[:, 4] = cell(np.arange(len(quads)), k - 1)
            f[:, 5] = cell(np.arange(len(quads)), k)
        if k in (0, nx):
            end = "in" if k == 0 else "out"
            put("inlet" if end == "in" else "outlet", f[zone == 0])
            put(f"solid_end_{end}", f[zone == 1])
        else:
            put("int_f", f[zone == 0])
            put("int_s", f[zone == 1])
    Z = {k: np.concatenate(v) for k, v in Z.items()}
    # node coordinates (x, y, z)
    X = np.repeat(x, n2)
    YZ = np.tile(xy, (nx + 1, 1))
    coords = np.column_stack([X, YZ])
    # cell centroids (for orientation)
    ncell = (nF + nS) * nx
    cen = np.zeros((ncell + 1, 3))
    qc = xy[quads].mean(axis=1)
    xc = 0.5 * (x[:-1] + x[1:])
    for k in range(nx):
        ids = cell(np.arange(len(quads)), k)
        cen[ids, 0] = xc[k]
        cen[ids, 1:] = qc
    return dict(x=x, xy=xy, quads=quads, zone=zone, nx=nx, n2=n2, nF=nF, nS=nS, ncell=ncell,
                Z=Z, coords=coords, cen=cen, info=info, level=level)


def orient(m, normal_toward="c0"):
    """Reorder face nodes so the right-hand normal (n1-n0)x(n2-n0) points toward c1 (or toward c0).
    For boundary faces (c1 = 0) 'toward c1' means pointing OUT of the domain."""
    flips = 0
    for key, f in m["Z"].items():
        P = m["coords"][f[:, :4] - 1]
        nrm = np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0])
        fc = P.mean(axis=1)
        d = fc - m["cen"][f[:, 4]]                       # c0 -> face: points away from c0
        s = np.einsum("ij,ij->i", nrm, d)
        want_pos = normal_toward == "c1"
        bad = (s < 0) if want_pos else (s > 0)
        f[bad, :4] = f[bad][:, ::-1][:, 2:6]             # reverse node order (n3 n2 n1 n0)
        flips += int(bad.sum())
    return flips


ZTAB = [("int_f", 4, "interior", "interior_fluid"), ("int_s", 5, "interior", "interior_solid"),
        ("interface", 6, "interior", "interface_fs"), ("inlet", 7, "velocity-inlet", "inlet"),
        ("outlet", 8, "pressure-outlet", "outlet"), ("solid_end_in", 9, "wall", "solid_end_in"),
        ("solid_end_out", 10, "wall", "solid_end_out"), ("outer_unheated", 11, "wall", "outer_unheated"),
        ("outer_q1", 12, "wall", "outer_q1"), ("outer_q2", 13, "wall", "outer_q2"),
        ("outer_q3", 14, "wall", "outer_q3"), ("outer_q4", 15, "wall", "outer_q4"),
        ("sym_f", 16, "symmetry", "symmetry_fluid"), ("sym_s", 17, "symmetry", "symmetry_solid")]


def face_area(m, key):
    f = m["Z"][key]
    P = m["coords"][f[:, :4] - 1]
    return float(0.5 * np.linalg.norm(np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 1]), axis=1).sum())


def write(m, path):
    h = lambda v: format(int(v), "x")
    nn = len(m["coords"])
    nfc = m["nF"] * m["nx"]
    nfaces = sum(len(v) for v in m["Z"].values())
    with open(path, "w", newline="\n") as fh:
        fh.write('(0 "3D half-pipe O-grid, written by make_3d_mesh.py")\n(2 3)\n')
        fh.write(f"(10 (0 1 {h(nn)} 0 3))\n(10 (1 1 {h(nn)} 1 3)(\n")
        np.savetxt(fh, m["coords"], fmt="%.10e")
        fh.write("))\n")
        fh.write(f"(12 (0 1 {h(m['ncell'])} 0 0))\n")
        fh.write(f"(12 (2 1 {h(nfc)} 1 4))\n(12 (3 {h(nfc + 1)} {h(m['ncell'])} 1 4))\n")
        fh.write(f"(13 (0 1 {h(nfaces)} 0 0))\n")
        start = 1
        for key, zid, tname, _ in ZTAB:
            f = m["Z"].get(key)
            if f is None or not len(f):
                continue
            fh.write(f"(13 ({h(zid)} {h(start)} {h(start + len(f) - 1)} {h(BC[tname])} 4)(\n")
            np.savetxt(fh, f, fmt="%x")
            fh.write("))\n")
            start += len(f)
        fh.write("(45 (2 fluid fluid)())\n(45 (3 solid solid)())\n")
        for key, zid, tname, dname in ZTAB:
            if key in m["Z"]:
                fh.write(f"(45 ({zid} {tname} {dname})())\n")


def validate(m, normal_toward="c0"):
    """Independent checks (no Fluent): every hex has 6 faces, closed cells, exact areas/volume, orientation."""
    sgn = 1.0 if normal_toward == "c1" else -1.0     # +1: right-hand normal points away from c0
    cnt = np.zeros(m["ncell"] + 1, int)
    net = np.zeros((m["ncell"] + 1, 3))                 # sum of outward area vectors per cell -> 0 if closed
    vol = np.zeros(m["ncell"] + 1)
    for key, f in m["Z"].items():
        P = m["coords"][f[:, :4] - 1]
        A = sgn * 0.5 * np.cross(P[:, 2] - P[:, 0], P[:, 3] - P[:, 1])  # outward area vector of c0
        fc = P.mean(axis=1)
        np.add.at(cnt, f[:, 4], 1)
        np.add.at(net, f[:, 4], A)                       # normal points away from c0 (after orient)
        np.add.at(vol, f[:, 4], np.einsum("ij,ij->i", fc, A) / 3.0)
        inner = f[:, 5] > 0
        np.add.at(cnt, f[inner, 5], 1)
        np.add.at(net, f[inner, 5], -A[inner])
        np.add.at(vol, f[inner, 5], -np.einsum("ij,ij->i", fc[inner], A[inner]) / 3.0)
    cnt, net, vol = cnt[1:], net[1:], vol[1:]
    exact_vol = 0.5 * math.pi * R_O ** 2 * L_TOTAL
    ncirc = m["info"]["ncirc"]
    poly = 0.5 * ncirc * math.sin(math.pi / ncirc) * 2 / math.pi    # faceted / exact area ratio for a half circle
    sec_exact = R_O * math.pi / 4 * (L_TOTAL - X_SPLIT)
    return dict(
        cells=m["ncell"], fluid_cells=m["nF"] * m["nx"], solid_cells=m["nS"] * m["nx"], nodes=len(m["coords"]),
        cross_section_cells=len(m["quads"]), axial_cells=m["nx"],
        every_cell_has_6_faces=bool(np.all(cnt == 6)),
        cells_closed=bool(np.abs(net).max() < 1e-12),
        min_cell_volume=float(vol.min()), negative_volumes=int((vol <= 0).sum()),
        volume_error_pct=100 * (vol.sum() / (exact_vol * poly) - 1),
        faceting_area_deficit_pct=100 * (1 - poly),
        inlet_area=face_area(m, "inlet"), inlet_area_exact=0.5 * math.pi * R_I ** 2,
        sector_area={k: face_area(m, f"outer_q{k}") for k in range(1, 5)}, sector_area_exact=sec_exact,
        interface_heated_area_exact=0.5 * math.pi * 2 * R_I * (L_TOTAL - X_SPLIT),
        transition_ratio_range=(m["info"]["trans_ratio_min"], m["info"]["trans_ratio_max"]),
        zones={k: len(v) for k, v in m["Z"].items()},
        under_512k=m["ncell"] <= 512_000)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("level", nargs="?", default="B3D", choices=sorted(LEVELS))
    ap.add_argument("--tiny", action="store_true", help="10 axial cells only (Fluent read / orientation test)")
    ap.add_argument("--normal-toward", default="c0", choices=["c1", "c0"],
                    help="face-normal convention. VERIFIED with Fluent 2026 R1 mesh check on the tiny mesh: Fluent needs the "
                         "right-hand normal (n1-n0)x(n2-n0) to point TOWARD c0 (the other choice gave 100 %% left-handed faces)")
    a = ap.parse_args()
    m = build(a.level, a.tiny)
    flips = orient(m, a.normal_toward)
    name = f"pipe3d_{a.level}{'_tiny' if a.tiny else ''}"
    out = HERE / f"{name}.msh"
    v = validate(m, a.normal_toward)
    for k, val in v.items():
        print(f"  {k}: {val}")
    if not (v["every_cell_has_6_faces"] and v["cells_closed"] and v["negative_volumes"] == 0):
        print("VALIDATION FAILED - mesh not written")
        return 1
    write(m, out)
    side = dict(level=a.level, tiny=a.tiny, normal_toward=a.normal_toward, half_model=True,
                theta0_deg=22.5, params=LEVELS[a.level], **{k: v[k] for k in
                ("cells", "inlet_area", "inlet_area_exact", "sector_area", "sector_area_exact",
                 "interface_heated_area_exact", "faceting_area_deficit_pct", "axial_cells", "cross_section_cells")})
    (HERE / f"{name}.json").write_text(json.dumps(side, indent=2))
    print(f"written {out} ({out.stat().st_size / 1e6:.1f} MB), {flips} faces re-ordered for orientation")
    return 0


if __name__ == "__main__":
    sys.exit(main())
