"""
run_gateD.py - Gate D mesh study in 3D: 3 meshes x 4 Reynolds numbers, then Nu table + GCI.

    .venv\\Scripts\\python.exe 03_fluent_setup\\run_gateD.py              # run everything not done yet, then analyse
    .venv\\Scripts\\python.exe 03_fluent_setup\\run_gateD.py --analyse    # only rebuild the summary from finished runs

Meshes (02_geometry_mesh/3d, make_3d_mesh.py): D3C coarse 107,712 cells, D3M medium 239,540, B3D fine 449,328.
Each spacing is scaled by the same factor (first cell, growth, layers, cells per sector, axial dx).
Runs go coarse -> medium -> fine (fast feedback first). A run whose result JSON exists is skipped, so the
script can be stopped and restarted. Each run is the normal 3D driver (fluent_gateB.py, uniform flux).

Output: 05_results/gateD_summary.csv and gateD_summary.txt
GCI method: Celik et al. (2008), J. Fluids Eng. 130, 078001 (Richardson extrapolation, Fs = 1.25).
Spec Gate D pass: Nu change between the two finest meshes below 1 %.
"""
import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import fluent_gateA as ga  # noqa: E402

RES = ROOT / "05_results"
MESH_DIR = ROOT / "02_geometry_mesh" / "3d"
MESHES = ["D3C", "D3M", "B3D"]          # coarse, medium, fine
RE_LIST = [12000, 14000, 16000, 25000]
PREFIX = "gateD"
VOLUME = 0.5 * math.pi * 0.0125 ** 2 * 4.0   # half pipe incl. wall (only ratios of h matter)


def result_file(re, mesh):
    return RES / f"{PREFIX}_Re{re}_{mesh}.json"


def run_all():
    py = sys.executable
    for mesh in MESHES:
        for re in RE_LIST:
            out = result_file(re, mesh)
            if out.exists():
                print(f"skip  Re {re} {mesh}: {out.name} exists", flush=True)
                continue
            print(f"\n=== Re {re}  mesh {mesh}   started {time.strftime('%H:%M')}", flush=True)
            t0 = time.time()
            rc = subprocess.call([py, str(HERE / "fluent_gateB.py"), "--re", str(re), "--mesh", mesh,
                                  "--prefix", PREFIX])
            print(f"=== Re {re}  mesh {mesh}   exit code {rc}   {(time.time() - t0) / 60:.1f} min", flush=True)


def gci(phi, n_cells):
    """phi, n_cells ordered FINE, MEDIUM, COARSE. Returns dict (Celik 2008)."""
    h = [(VOLUME / n) ** (1 / 3) for n in n_cells]
    r21, r32 = h[1] / h[0], h[2] / h[1]
    e21, e32 = phi[1] - phi[0], phi[2] - phi[1]
    out = dict(r21=r21, r32=r32, change_fine_medium_pct=100 * abs(e21) / abs(phi[0]),
               change_medium_coarse_pct=100 * abs(e32) / abs(phi[1]))
    if abs(e21) < 1e-12 or abs(e32) < 1e-12:
        out.update(note="differences too small for Richardson extrapolation")
        return out
    s = math.copysign(1.0, e32 / e21)
    p = 2.0
    for _ in range(200):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p_new = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
        if abs(p_new - p) < 1e-10:
            break
        p = 0.5 * (p + p_new)
    phi_ext = (r21 ** p * phi[0] - phi[1]) / (r21 ** p - 1)
    out.update(monotonic=e21 * e32 > 0, order_p=p, Nu_extrapolated=phi_ext,
               fine_error_vs_extrap_pct=100 * abs((phi_ext - phi[0]) / phi_ext),
               GCI_fine_pct=100 * 1.25 * abs(e21 / phi[0]) / (r21 ** p - 1),
               GCI_medium_pct=100 * 1.25 * abs(e32 / phi[1]) / (r32 ** p - 1))
    return out


def analyse():
    cells = {m: json.loads((MESH_DIR / f"pipe3d_{m}.json").read_text())["cells"] for m in MESHES}
    rows, lines = [], []
    pr = ga.WATER["mu"] * ga.WATER["cp"] / ga.WATER["k"]
    lines.append("Gate D: 3D mesh study, water, uniform flux, Nu = fully developed window (x = 3.0-3.85 m)")
    lines.append("meshes: " + ", ".join(f"{m} {cells[m]:,} cells" for m in MESHES))
    lines.append("")
    for re in RE_LIST:
        nu, extra = {}, {}
        for m in MESHES:
            f = result_file(re, m)
            if f.exists():
                d = json.loads(f.read_text())
                nu[m] = d["Nu_window"]
                extra[m] = dict(conv=d.get("converged"), its=d.get("iterations"), eb=d.get("energy_balance_pct"),
                                yp=d.get("yplus_max"))
        gn = ga.gnielinski(re, pr)
        db = 0.023 * re ** 0.8 * pr ** 0.4
        lines.append(f"Re {re}:  Gnielinski {gn:.2f}   Dittus-Boelter {db:.2f}")
        for m in MESHES:
            if m in nu:
                e = extra[m]
                lines.append(f"   {m:4s} {cells[m]:>8,} cells  Nu {nu[m]:9.4f}  vs Gn {100 * (nu[m] / gn - 1):+6.2f} %  "
                             f"converged={e['conv']} its={e['its']} energy bal {e['eb']:.2e} %  y+max {e['yp']:.3f}")
            else:
                lines.append(f"   {m:4s} not run yet")
        row = dict(Re=re, Nu_gnielinski=gn, Nu_dittus_boelter=db, **{f"Nu_{m}": nu.get(m) for m in MESHES})
        if len(nu) == 3:
            g = gci([nu["B3D"], nu["D3M"], nu["D3C"]], [cells["B3D"], cells["D3M"], cells["D3C"]])
            row.update(g)
            row["gateD_pass"] = g["change_fine_medium_pct"] < 1.0
            lines.append(f"   change medium->fine {g['change_fine_medium_pct']:.3f} %, coarse->medium "
                         f"{g['change_medium_coarse_pct']:.3f} %   Gate D (<1 %): {'PASS' if row['gateD_pass'] else 'FAIL'}")
            if "order_p" in g:
                lines.append(f"   GCI: order p {g['order_p']:.2f} ({'monotonic' if g['monotonic'] else 'OSCILLATORY'}), "
                             f"Nu extrapolated {g['Nu_extrapolated']:.3f}, GCI fine {g['GCI_fine_pct']:.3f} %, "
                             f"GCI medium {g['GCI_medium_pct']:.3f} %")
            else:
                lines.append(f"   GCI: {g['note']}")
        lines.append("")
        rows.append(row)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k != "Re", k))
    with open(RES / f"{PREFIX}_summary.csv", "w") as fh:
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join("" if r.get(k) is None else str(r[k]) for k in keys) + "\n")
    (RES / f"{PREFIX}_summary.txt").write_text("\n".join(lines), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--analyse", action="store_true", help="only build the summary from finished runs")
    a = ap.parse_args()
    if not a.analyse:
        run_all()
    analyse()
