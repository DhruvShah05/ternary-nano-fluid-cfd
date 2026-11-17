"""
probe_wall_temperature.py - READ-ONLY check of which Fluent field gives the true inner-wall temperature.
It does NOT solve. It reads saved Gate A case+data files and, for two bands in the evaluation window,
prints the bulk temperature and several candidate wall-temperature fields, with the Nusselt number each
one would imply:   Nu = q_in * D_i / (k * (T_wall - T_bulk)).

Why: with standard wall functions the first cell is about 2 mm thick, so the 'temperature' field on the
clipped wall band may be a near-wall cell value, not the wall face temperature (Nu came out at about 1200).

Run (PowerShell, project root, venv active; close any Fluent GUI first):
    python 03_fluent_setup\\probe_wall_temperature.py --tag Re10000_M2WF_ke_std_swf --tag Re10000_M2
Output: 05_results\\probe_wall_temperature.txt  (upload it)
"""
import argparse
import re
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import fluent_gateA as ga  # noqa: E402

CANDIDATES = ["temperature", "wall-temperature", "wall-adjacent-temperature", "reference-temperature-at-y+"]
BANDS = (20, 30)
K_WATER = ga.WATER["k"]


def implied_nu(q, tw, tb):
    d = tw - tb
    return q * ga.D_I / (K_WATER * d) if abs(d) > 1e-9 else float("nan")


def analyse(post, bands, log):
    """post needs .mw(surface, field), .aw(surface, field), .fluid_wall, .shadow_wall.
    Returns a list of result rows (also logged)."""
    rows = []
    for i in bands:
        tb = post.mw(f"pl_{i:02d}", "temperature")
        q = abs(post.aw(f"bd_{i:02d}", "heat-flux"))
        log(f"\n  band {i}: bulk T (plane, mass-weighted) = {tb:.4f} K   q_in = {q:.2f} W/m2")
        for side, surf in (("fluid side", f"bd_{i:02d}"), ("solid side", f"bs_{i:02d}")):
            for fld in CANDIDATES:
                try:
                    tw = post.aw(surf, fld)
                except Exception as e:  # noqa: BLE001
                    log(f"    {side:10s} {fld:30s} FAILED: {str(e)[:110]}")
                    continue
                nu = implied_nu(q, tw, tb)
                log(f"    {side:10s} {fld:30s} T = {tw:.4f} K   T-Tb = {tw - tb:+.4f} K   implied Nu = {nu:.2f}")
                rows.append(dict(band=i, side=side, field=fld, T=tw, dT=tw - tb, Nu=nu))
    return rows


def read_case(s, case, log):
    try:
        s.file.read_case_data(file_name=str(case))
        log("  read via file.read_case_data")
        return
    except Exception as e:  # noqa: BLE001
        log(f"  read_case_data failed ({str(e)[:120]}); trying read_case + read_data")
    s.file.read_case(file_name=str(case))
    s.file.read_data(file_name=str(case).replace(".cas.h5", ".dat.h5"))
    log("  read via read_case + read_data")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", action="append", required=True, help="e.g. Re10000_M2WF_ke_std_swf (repeatable)")
    a = ap.parse_args()
    out = ROOT / "05_results" / "probe_wall_temperature.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    lines = []

    def log(m=""):
        m = str(m)
        print(m, flush=True)
        lines.append(m)
        out.write_text("\n".join(lines), encoding="utf-8")

    import ansys.fluent.core as pyfluent
    solver = pyfluent.launch_fluent(mode="solver", dimension=2, precision="double", processor_count=ga.PROCS)
    try:
        s = solver.settings
        for tag in a.tag:
            log(f"\n=========== {tag} ===========")
            case = ROOT / "04_cases" / "gateA" / f"gateA_{tag}.cas.h5"
            if not case.exists():
                log(f"  case file not found: {case}")
                continue
            m = re.match(r"Re(\d+)_", tag)
            re_ = float(m.group(1)) if m else 10000.0
            nums = ga.case_numbers(re_)
            read_case(s, case, log)
            post = ga.Post(solver, Path(tempfile.gettempdir()) / "_probe_tmp.txt", log)
            post.detect_fluid_wall()
            post.shadow_wall = "interface_fs-shadow" if post.fluid_wall == "interface_fs" else "interface_fs"
            post.build_surfaces()
            for i in BANDS:     # same x-range clip on the solid-side wall
                o = post._make(post.surf.iso_clip, f"bs_{i:02d}")
                o.field = "x-coordinate"
                o.surfaces = [post.shadow_wall]
                o.range.minimum = ga.X_HEAT + ga.BAND * i
                o.range.maximum = ga.X_HEAT + ga.BAND * (i + 1)
            log(f"  Gnielinski Nu = {nums['nu_gnielinski']:.2f}; expected wall-to-bulk dT at that Nu = "
                f"{nums['q_in_expected'] * ga.D_I / (K_WATER * nums['nu_gnielinski']):.2f} K")
            analyse(post, BANDS, log)
        log("\nDONE. Upload 05_results\\probe_wall_temperature.txt")
    finally:
        try:
            solver.exit()
        except Exception:  # noqa: BLE001
            pass


if __name__ == "__main__":
    main()
