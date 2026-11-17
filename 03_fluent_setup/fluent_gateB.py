"""
Gate B: 3D HALF-PIPE, all heated sectors at the same flux, water, compared with the 2D axisymmetric result.
Pass (spec Sec. 8): window Nu within 2 % of the 2D result. Fluent 2026 R1 / PyFluent 0.43.0, Student licence.

Usage (PowerShell, project root, venv active; close any Fluent GUI first):
    python 02_geometry_mesh\\3d\\make_3d_mesh.py                       # writes pipe3d_B3D.msh (once)
    python 03_fluent_setup\\fluent_gateB.py --mesh B3D_tiny --test     # Fluent read + plumbing check, minutes
    python 03_fluent_setup\\fluent_gateB.py --re 10000                 # the Gate B run (hours)

Everything physical is identical to fluent_gateA.py (same water, flux per m2, turbulence model
realizable k-epsilon + enhanced wall treatment, second-order discretisation, relaxation 0.5,
pseudo-time 0.1, flow-first start, convergence tests, band / plane / window definitions). Only differences:
  * 3D half model (z >= 0) with symmetry plane z = 0 (zones symmetry_fluid, symmetry_solid).
  * Heated wall is split into outer_q1..outer_q4 (theta0 = 22.5 deg, see make_3d_mesh.py).
  * Inlet velocity is set from the FACETED inlet area in the mesh sidecar JSON so the mass flow equals the
    2D case exactly (the polygon cross-section is 0.29 % smaller than the circle).
  * Fluent reports HALF the mass flow and heat of the full pipe; the bulk temperature rise is unchanged.
"""
import argparse
import json
import math
import shutil
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import fluent_gateA as ga  # noqa: E402

MESH_DIR = ROOT / "02_geometry_mesh" / "3d"
SECTORS = (1, 2, 3, 4)
GATE_B_TOL = 0.02
REF_2D = {"M2": "gateA_Re{re}_M2_ke_rlz_ewt.json", "B3D": "gateA_Re{re}_B3D_ke_rlz_ewt.json"}


def case_numbers_3d(re, side, fluxes=None):
    """Full-pipe numbers from Gate A, converted to the half model with the mesh's faceted areas."""
    n = ga.case_numbers(re)
    p = ga.WATER
    mdot_half = n["mdot"] / 2.0
    v = mdot_half / (p["rho"] * side["inlet_area"])
    q = {k: (n["q_out"] if fluxes is None else fluxes[k - 1]) for k in SECTORS}
    q_total = sum(q[k] * side["sector_area"][str(k)] for k in SECTORS)
    return dict(n, v=v, v_circle=n["v"], mdot=mdot_half, q_sector=q, q_total=q_total,
                dTb_expected=q_total / (mdot_half * p["cp"]),
                q_in_expected=sum(q.values()) / len(q) * ga.D_O / ga.D_I)


def make_tag(re_, mesh, test=False):
    return f"Re{int(re_)}_{mesh}" + ("_TEST" if test else "")


# ---------------------------------------------------------------------------- Fluent side
def prepare_mesh(solver, mesh, log, transcript=None):
    s = solver.settings
    s.file.read_case(file_name=str(mesh))
    solver.tui.mesh.modify_zones.slit_interior_between_diff_materials()
    s.setup.models.energy.enabled = True
    log("Mesh read (3D), interface slit, energy on.")
    try:
        solver.tui.mesh.check()
        log("  mesh check run (see the transcript for its report)")
    except Exception as e:  # noqa: BLE001
        log(f"  mesh check call failed: {str(e)[:150]}")
    if transcript is not None and Path(transcript).exists():
        txt = Path(transcript).read_text(errors="replace")
        bad = [ln.strip() for ln in txt.splitlines()
               if "left-handed faces" in ln or "non-positive volume" in ln or "wrong node order" in ln]
        if bad or "Mesh check failed" in txt:
            raise RuntimeError("Fluent mesh check FAILED (" + "; ".join(bad[:3]) + "). Regenerate the mesh; do not solve.")
        log("  Fluent mesh check: no left-handed faces, no non-positive volumes")


def apply_setup(solver, nums, log, relax, scale1, turb):
    s = solver.settings
    bc = s.setup.boundary_conditions
    mats = s.setup.materials
    st = ga.Steps(log)
    p = ga.WATER
    st.run(f"viscous model {turb}", lambda: ga.configure_turbulence(s.setup.models.viscous, turb, None, None, log))

    def fluid():
        if "nf" not in list(mats.fluid.keys()):
            mats.fluid.create("nf")
        m = mats.fluid["nf"]
        for prop, val in (("density", p["rho"]), ("specific_heat", p["cp"]),
                          ("thermal_conductivity", p["k"]), ("viscosity", p["mu"])):
            obj = getattr(m, prop)
            obj.option = "value"
            obj.value = val
    st.run("material nf (water, constant)", fluid)

    def copper():
        if "copper" not in list(mats.solid.keys()):
            mats.database.copy_by_name(type="solid", name="copper")
    st.run("material copper from database", copper)
    st.run("fluid zone -> nf",
           lambda: setattr(s.setup.cell_zone_conditions.fluid["fluid"].general, "material", "nf"))
    st.run("solid zone -> copper",
           lambda: setattr(s.setup.cell_zone_conditions.solid["solid"].general, "material", "copper"))

    def inlet():
        vi = bc.velocity_inlet["inlet"]
        ga.setv(vi.momentum, "velocity_magnitude", nums["v"])
        ga.setv(vi.thermal, "temperature", ga.T_IN)
        vi.turbulence.turbulence_specification = "Intensity and Hydraulic Diameter"
        vi.turbulence.turbulent_intensity = ga.TURB_INTENSITY
        ga.setv(vi.turbulence, "hydraulic_diameter", ga.D_I)
    st.run(f"velocity inlet {nums['v']:.6f} m/s", inlet)

    def outlet():
        po = bc.pressure_outlet["outlet"]
        ga.setv(po.momentum, "gauge_pressure", 0.0)
        po.turbulence.turbulence_specification = "Intensity and Hydraulic Diameter"
        po.turbulence.backflow_turbulent_intensity = ga.TURB_INTENSITY
        ga.setv(po.turbulence, "backflow_hydraulic_diameter", ga.D_I)
        ga.setv(po.thermal, "backflow_total_temperature", ga.T_IN)
    st.run("pressure outlet", outlet)

    def heated():
        for k in SECTORS:
            th = bc.wall[f"outer_q{k}"].thermal
            th.thermal_condition = "Heat Flux"
            ga.setv(th, "heat_flux", nums["q_sector"][k])
    st.run("outer_q1..q4 heat flux " + ", ".join(f"{nums['q_sector'][k]:.1f}" for k in SECTORS), heated)

    def adiabatic():
        for n in ("outer_unheated", "solid_end_in", "solid_end_out"):
            th = bc.wall[n].thermal
            th.thermal_condition = "Heat Flux"
            ga.setv(th, "heat_flux", 0.0)
    st.run("adiabatic walls", adiabatic)

    def symmetry():
        names = list(bc.symmetry.keys())
        if sorted(names) != ["symmetry_fluid", "symmetry_solid"]:
            raise RuntimeError(f"symmetry zones are {names}")
    st.run("symmetry plane zones present", symmetry)

    def relax_fn():
        pv = s.solution.controls.p_v_controls
        pv.explicit_momentum_under_relaxation = relax
        pv.explicit_pressure_under_relaxation = relax
    st.run(f"explicit relaxation {relax}", relax_fn)

    def residuals():
        eq = s.solution.monitor.residual.equations
        second = "omega" if ga.TURB_PRESETS[turb]["model"] == "k-omega" else "epsilon"
        for n, t in (("continuity", 1e-6), ("x-velocity", 1e-6), ("y-velocity", 1e-6), ("z-velocity", 1e-6),
                     ("k", 1e-6), (second, 1e-6), ("energy", 1e-9)):
            eq[n].absolute_criteria = t
            eq[n].check_convergence = False
    st.run("residual targets, auto-stop OFF", residuals)

    def discretisation():
        ds = s.solution.methods.spatial_discretization.discretization_scheme
        for eq in list(ds.get_state()):
            ds[eq] = "second-order" if eq == "pressure" else "second-order-upwind"
        got = ds.get_state()
        bad = {k: v for k, v in got.items() if v != ("second-order" if k == "pressure" else "second-order-upwind")}
        if bad:
            raise RuntimeError(f"discretisation not second order: {bad}")
        log(f"    verified discretisation: {got}")
    st.run("second-order discretisation (all equations)", discretisation)
    st.run(f"pseudo-time scale factor {scale1}",
           lambda: setattr(s.solution.run_calculation.pseudo_time_settings.time_step_method,
                           "time_step_size_scale_factor", scale1))
    return st


class Post3D(ga.Post):
    def energy_balance(self):
        t_in = self.mw("inlet", "temperature")
        t_out = self.mw("outlet", "temperature")
        mdot = abs(self.flux("mass_flow", "inlet"))
        q_wall = sum(abs(self.flux("heat_transfer", f"outer_q{k}")) for k in SECTORS)
        q_out = mdot * ga.WATER["cp"] * (t_out - t_in)
        return dict(T_in_meas=t_in, T_out=t_out, mdot_fluent=mdot, Q_wall_fluent=q_wall,
                    Q_enthalpy=q_out, eb=(q_wall - q_out) / q_wall if q_wall else float("nan"))

    def wall_t_spread(self, band=30):
        """Circumferential inner-wall temperature spread (max - min) on the band starting at x = 3.5 m."""
        si = self.rep.surface_integrals
        fmax, fmin = getattr(si, "facet_max", None), getattr(si, "facet_min", None)
        if fmax is None or fmin is None:
            return None
        name = f"bd_{band:02d}"
        return (self._scalar(fmax, surface_names=[name], report_of=self.tw_field)
                - self._scalar(fmin, surface_names=[name], report_of=self.tw_field))


def finalise(post, nums, total, converged, tag, results_dir, log, re_, extra, prefix="gateB"):
    is_gate = prefix == "gateB"
    prof = [post.band_values(i) for i in range(ga.N_BANDS)]
    win = [r for r in prof if ga.WIN_FIRST <= r["i"] <= ga.WIN_LAST]
    nu_win = sum(r["Nu"] for r in win) / len(win)
    spread = (max(r["Nu"] for r in win) - min(r["Nu"] for r in win)) / nu_win
    bal = post.energy_balance()
    dp = post.aw("pl_f12", "pressure") - post.aw("pl_f20", "pressure")
    f_cfd = 2.0 * ga.D_I * dp / (ga.WATER["rho"] * nums["v"] ** 2 * 0.8)
    qin_dev = max(abs(r["q"] - nums["q_in_expected"]) / nums["q_in_expected"] for r in win)
    yp = post.y_plus_stats()
    try:
        wts = post.wall_t_spread()
    except Exception as e:  # noqa: BLE001
        log(f"  wall-T spread failed: {str(e)[:150]}")
        wts = None
    refs = {}
    for key, pat in REF_2D.items():
        fp = results_dir / pat.format(re=int(re_))
        if fp.exists():
            nu2 = json.loads(fp.read_text())["Nu_window"]
            refs[key] = dict(file=fp.name, Nu_2D=nu2, diff_pct=100 * (nu_win / nu2 - 1))
    main_ref = refs.get("M2") if is_gate else None
    sector = {}
    for k in SECTORS:   # outer-surface temperature of each heated sector (area average, x = 2..4 m)
        try:
            sector[f"T_outer_q{k}_avg"] = post.aw(f"outer_q{k}", "temperature")
        except Exception as e:  # noqa: BLE001
            log(f"  sector q{k} temperature failed: {str(e)[:120]}")
    res = dict(tag=tag, Re=re_, converged=converged, iterations=total,
               Nu_window=nu_win, Nu_heated_avg=sum(r["Nu"] for r in prof) / len(prof),
               Nu_gnielinski=nums["nu_gnielinski"], gnielinski_error_pct=100 * (nu_win / nums["nu_gnielinski"] - 1),
               reference_2D=refs, gateB_diff_pct=main_ref["diff_pct"] if main_ref else None,
               gateB_pass=(abs(main_ref["diff_pct"]) <= 100 * GATE_B_TOL) if main_ref else None,
               window_spread_pct=100 * spread, window_flat_flag=spread <= 0.02,
               energy_balance_pct=100 * bal["eb"], fluent_mdot_half=bal["mdot_fluent"],
               fluent_mdot_over_expected=bal["mdot_fluent"] / nums["mdot"],
               T_out=bal["T_out"], dTb_expected=nums["dTb_expected"], dTb_cfd=bal["T_out"] - bal["T_in_meas"],
               q_in_max_dev_from_expected_pct=100 * qin_dev, f_cfd=f_cfd, f_correlation=nums["f"],
               delta_p_1p2_to_2p0=dp, wall_T_spread_x3p5=wts, inlet_velocity=nums["v"],
               q_outer_sector_Wm2={f"q{k}": nums["q_sector"][k] for k in SECTORS},
               Q_total_half_W=nums["q_total"], **sector, **yp, **extra)
    (results_dir / f"{prefix}_{tag}.json").write_text(json.dumps(res, indent=2))
    with open(results_dir / f"{prefix}_{tag}_profile.csv", "w") as fh:
        fh.write("band,x_centre,Tb,Tw,q_in,Nu\n")
        for r in prof:
            fh.write(f"{r['i']},{r['x']:.4f},{r['Tb']:.6f},{r['Tw']:.6f},{r['q']:.4f},{r['Nu']:.4f}\n")
    log("\n==== RESULT ====")
    for k, v in res.items():
        log(f"  {k}: {v}")
    if main_ref:
        log(f"\nGate B ({tag}): {'PASS' if res['gateB_pass'] else 'FAIL'}  (Nu_3D {nu_win:.3f} vs Nu_2D(M2) "
            f"{main_ref['Nu_2D']:.3f}, {main_ref['diff_pct']:+.2f} %, limit 2 %)")
    elif is_gate:
        log("\nGate B: no 2D reference file found; run fluent_gateA.py --re <Re> --mesh M2 first.")
    else:
        log(f"\nNon-uniform flux case ({tag}): Nu (circumferential average, window) {nu_win:.3f}; "
            f"uniform-flux Gate B gave 84.286 for comparison. Not a validation gate.")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--re", type=float, default=10000.0)
    ap.add_argument("--mesh", default="B3D", help="mesh level: B3D or B3D_tiny (file 02_geometry_mesh/3d/pipe3d_<mesh>.msh)")
    ap.add_argument("--turb", default=ga.DEFAULT_TURB, choices=sorted(ga.TURB_PRESETS))
    ap.add_argument("--test", action="store_true", help="short plumbing check only")
    ap.add_argument("--relax", type=float, default=ga.DEFAULT_RELAX)
    ap.add_argument("--scale1", type=float, default=ga.DEFAULT_SCALE1)
    ap.add_argument("--suns", default=None,
                    help="per-sector outer flux in suns for q1,q2,q3,q4, e.g. 1,1,20,20 (mirrored: q8,q7,q6,q5). "
                         "Omit for the uniform Gate B flux.")
    ap.add_argument("--sun", type=float, default=1000.0, help="W/m2 per sun (default 1000)")
    a = ap.parse_args()

    fluxes = None
    if a.suns:
        suns = [float(t) for t in a.suns.split(",")]
        if len(suns) != 4 or min(suns) < 0:
            print("--suns needs four non-negative numbers, e.g. 1,1,20,20")
            return 1
        fluxes = [s_ * a.sun for s_ in suns]
    tag = make_tag(a.re, a.mesh, a.test)
    prefix = "gateB"
    if fluxes:
        tag += "_sun" + "-".join(f"{s_:g}" for s_ in suns)
        prefix = "sun3d"
    results_dir = ROOT / "05_results"
    results_dir.mkdir(parents=True, exist_ok=True)
    logfile = results_dir / f"{prefix}_{tag}_log.txt"
    lines = []

    def log(msg=""):
        print(msg, flush=True)
        lines.append(str(msg))
        logfile.write_text("\n".join(lines), encoding="utf-8")

    mesh = MESH_DIR / f"pipe3d_{a.mesh}.msh"
    side_f = mesh.with_suffix(".json")
    if not mesh.exists() or not side_f.exists():
        log(f"Mesh or sidecar not found: {mesh}. Run: python 02_geometry_mesh\\3d\\make_3d_mesh.py"
            + (" --tiny" if "tiny" in a.mesh else ""))
        return 1
    if a.re < ga.RE_MIN:
        log(f"STOP: Re = {a.re:g} is below {ga.RE_MIN:g}.")
        return 1
    if ga.TURB_PRESETS[a.turb]["wall_function"]:
        log("STOP: wall-function models need a coarse-wall mesh; this 3D mesh has a 5 um first layer.")
        return 1
    side = json.loads(side_f.read_text())
    if side["cells"] > 512_000:
        log(f"STOP: mesh has {side['cells']} cells, above the 512k student limit.")
        return 1
    nums = case_numbers_3d(a.re, side, fluxes)
    log(f"{prefix} {tag}  mesh {mesh.name} ({side['cells']} cells, half model)   turbulence {a.turb}   "
        f"relax={a.relax} scale1={a.scale1}")
    log("Expected (half model): " + ", ".join(f"{k}={v:.5g}" for k, v in nums.items() if isinstance(v, float)))
    dp_darcy = nums["f"] * (4.0 / ga.D_I) * ga.WATER["rho"] * nums["v"] ** 2 / 2.0
    log(f"Expected inlet static pressure about {dp_darcy:.0f} Pa (Darcy)")

    import ansys.fluent.core as pyfluent
    t0 = time.time()
    solver = pyfluent.launch_fluent(mode="solver", dimension=3, precision="double", processor_count=ga.PROCS)
    tmp_t = Path(tempfile.gettempdir()) / f"fluent_transcript_gateB_{tag}.txt"
    try:
        tmp_t.unlink(missing_ok=True)   # Fluent appends to an existing transcript; a stale one fooled the mesh-check guard
        solver.tui.file.start_transcript(f'"{tmp_t}"' if " " in str(tmp_t) else str(tmp_t))
        prepare_mesh(solver, mesh, log, tmp_t)
        log("Applying setup:")
        steps = apply_setup(solver, nums, log, a.relax, a.scale1, a.turb)
        if steps.fails:
            log(f"\nSETUP FAILED for: {steps.fails}\nNo iterations were run.")
            return 2
        s = solver.settings
        log("Initialising (flow-first: temperature equation OFF):")
        ga.set_temperature_equation(s, solver, False, log)
        ga.initialise(solver, log)
        post = Post3D(solver, results_dir / "_tmp_report.txt", log)
        iterate = lambda n: s.solution.run_calculation.iterate(iter_count=int(n))  # noqa: E731
        log("Warm-up 50 iterations, then post-processing setup:")
        iterate(50)
        post.detect_fluid_wall()
        post.build_surfaces()
        total = 50
        log("Flow stage (temperature equation OFF):")
        if a.test:
            iterate(50)
            total += 50
            log(f"  test flow probe: {post.flow_probe()}")
            fl_its = 50
        else:
            fl_its, status, hist = ga.run_flow_stage(iterate, post.flow_probe, log, p_limit=max(1e4, 20.0 * dp_darcy))
            total += fl_its
            if status != "settled":
                log(f"\nFLOW STAGE ended as '{status}' after {fl_its} its; stopping (no ladder in 3D yet). "
                    f"p_in history {[round(h, 1) for h in hist[-6:]]}")
                return 4
            log(f"  flow settled after {fl_its} its")
        log("Energy stage (temperature equation ON):")
        ga.set_temperature_equation(s, solver, True, log)
        if a.test:
            iterate(50)
            snap = post.snapshot()
            log(f"  snapshot: Nu_win={snap['Nu']:.3f} energy balance={snap['eb']*100:.3f}% T_out={snap['T'][-1]:.4f}")
            log(f"  balance detail: {snap['bal']}")
            log(f"  sample band: {snap['rows'][0]}")
            log(f"  y+: {post.y_plus_stats()}")
            log(f"  wall-T spread at x=3.5: {post.wall_t_spread()}")
            log(f"TEST DONE in {(time.time()-t0)/60:.1f} min.")
            return 0
        e_total, converged, _ = ga.run_iterations(iterate, lambda: None, post.snapshot, None, log,
                                                  warm=0, probe=post.probe)
        total += e_total
        log(f"Finished: {total} iterations (flow {fl_its} + energy {e_total}), converged={converged}, "
            f"wall time {(time.time()-t0)/60:.1f} min")
        finalise(post, nums, total, converged, tag, results_dir, log, a.re,
                 extra=dict(mesh=a.mesh, turb=a.turb, cells=side["cells"], flow_iterations=fl_its,
                            recipe_relax=a.relax, recipe_scale1=a.scale1), prefix=prefix)
        case_dir = ROOT / "04_cases" / prefix
        case_dir.mkdir(parents=True, exist_ok=True)
        s.file.write_case_data(file_name=str(case_dir / f"{prefix}_{tag}.cas.h5"))
        log(f"Saved case and data to {case_dir}")
        return 0
    except Exception as e:  # noqa: BLE001
        log(f"\nABORTED: {type(e).__name__}: {e}")
        return 3
    finally:
        try:
            solver.tui.file.stop_transcript()
        except Exception:  # noqa: BLE001
            pass
        try:
            solver.exit()
        except Exception:  # noqa: BLE001
            pass
        if tmp_t.exists():
            shutil.copy(tmp_t, results_dir / f"{prefix}_{tag}_fluent_transcript.txt")


if __name__ == "__main__":
    sys.exit(main())
