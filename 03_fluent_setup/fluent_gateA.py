"""
Gate A: 2D axisymmetric, uniform flux, water, compared with Gnielinski.
Fluent 2026 R1 / PyFluent 0.43.0, Student licence (max 4 processes).

Usage (PowerShell, from the project root, venv active):
    python 03_fluent_setup\\fluent_gateA.py --test              # plumbing check, ~2 min
    python 03_fluent_setup\\fluent_gateA.py --re 10000 --mesh M2
    python 03_fluent_setup\\fluent_gateA.py --re 25000 --mesh M2

Turbulence model: DEFAULT is realizable k-epsilon with enhanced wall treatment (ke_rlz_ewt).
  DEVIATION from spec 6.3 (k-omega SST): SST failed Gate A (+16.7 % at Re 10,000, +10.9 % at Re 25,000),
  realizable k-epsilon EWT passed at Re 25,000 (+6.9 %). Pending the supervisor's agreement.
  Use --turb sst to reproduce the spec model.
    --turb ke_rlz_ewt (default) | sst | ke_std_swf | ke_rlz_swf | ke_std_ewt
        ke_*_swf = k-epsilon with STANDARD WALL FUNCTIONS -> needs the coarse-wall mesh:  --mesh M2WF
        ke_*_ewt = k-epsilon with ENHANCED WALL TREATMENT -> works on the resolved meshes (M1/M2/M3)
    --prt X [--prt-wall Y]   DIAGNOSTIC ONLY: change the constant turbulent Prandtl number (spec: 0.85)

ALWAYS run --test first. Several Fluent calls below are UNVERIFIED (see the
list in the chat reply). Setup failures are collected and reported together
before any long iteration is started.

Definitions I had to choose where the spec is silent (all reported in output):
  * Band i (i = 0..39) is the inner-wall strip x in [2+0.05 i, 2.05+0.05 i].
    Its bulk temperature is taken from the plane at the band centre.
  * The 3.0 to 3.85 m window = bands i = 20..36 (17 bands).
  * Spec 6.8(a) wall-temperature stations x = 3.0, 3.5, 3.85 are monitored with
    the bands i = 20, 30, 36 (starting at 3.0, 3.5 and ending at 3.85).
  * FLOW-FIRST (deviation from spec 6.9, constant-property fluids only): the
    temperature equation is switched off while the flow converges, then
    switched on. With constant properties the two are exactly decoupled.
    This is NOT valid for temperature-dependent fluids (oils, salts).
  * DEVIATION: explicit relaxation 0.5 and pseudo-time scale 0.1 are used instead of the spec's
    0.75 and 1.0, because the spec values gave reverse flow at the outlet and divergence.
    The scale factor is NOT raised afterwards (0.1 -> 0.5 made the flow diverge); use --scale2 only to experiment.
"""
import argparse
import json
import math
import re
import shutil
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
PROCS = 4  # Student licence limit (a request for 6 killed Fluent at launch)

# ------------------------------------------------------------ fixed spec values
D_I, D_O = 0.023, 0.025
R_I = D_I / 2.0
X_HEAT = 2.0
T_IN = 300.0
WATER = dict(rho=996.5, cp=4181.0, k=0.613, mu=8.55e-4)
DT_TARGET = 5.0        # mean wall-to-bulk dT used to scale the flux (handoff decision 3)
TURB_INTENSITY = 0.05
BAND = 0.05
N_BANDS = 40
WIN_FIRST, WIN_LAST = 20, 36            # bands inside 3.0 to 3.85 m
MONITOR_BANDS = (20, 30, 36)
BLOCK = 100
FIRST_STAGE = 2000          # (kept for reference) spec: first 2000 iterations
ENERGY_FIRST = 500          # energy-stage iterations before convergence checking starts (T_out was steady by 500)
FLOW_MIN, FLOW_MAX, FLOW_TOL = 500, 6000, 1e-3
FLOW_MASS_TOL = 5e-3          # |m_out/m_in - 1| needed to call the flow settled
LADDER = [(0.3, 0.03), (0.2, 0.01)]   # fallbacks after the default recipe: (explicit relaxation, pseudo-time scale)
DEFAULT_RELAX, DEFAULT_SCALE1 = 0.5, 0.1   # the spec's 0.75 / 1.0 DIVERGE (reverse flow at the outlet) in this setup
MAX_ITER = 8000
TOL_T, TOL_NU, TOL_EB = 1e-4, 5e-4, 2e-3   # spec 6.8: 0.01 %, 0.05 %, 0.2 %


# Turbulence presets (setting paths verified from introspect_dump_ke.txt, Fluent 2026 R1):
#   viscous.model 'k-omega' | 'k-epsilon'; viscous.k_epsilon_model 'standard'|'realizable'|'rng';
#   viscous.near_wall_treatment.wall_treatment 'standard-wall-fn'|'enhanced-wall-treatment'|...
TURB_PRESETS = {
    "sst":        dict(model="k-omega",   variant=None,         wall=None,                      wall_function=False),
    "ke_std_swf": dict(model="k-epsilon", variant="standard",   wall="standard-wall-fn",        wall_function=True),
    "ke_rlz_swf": dict(model="k-epsilon", variant="realizable", wall="standard-wall-fn",        wall_function=True),
    "ke_std_ewt": dict(model="k-epsilon", variant="standard",   wall="enhanced-wall-treatment", wall_function=False),
    "ke_rlz_ewt": dict(model="k-epsilon", variant="realizable", wall="enhanced-wall-treatment", wall_function=False),
}
# Wall temperature field. VERIFIED by probe_wall_temperature.py: on a clipped wall band the 'temperature' field equals
# 'wall-adjacent-temperature' (the near-wall CELL value), NOT the wall face temperature. That is harmless on a resolved
# wall (5 um cell: about 1 % effect on Nu) but wrong by a factor of about 16 with wall functions (2 mm cell).
TW_FIELD = "wall-temperature"
WF_MIN_YPLUS = 11.225     # Fluent switches to the log-law at y* = 11.225; 30 to 300 is the usual recommendation
DEFAULT_TURB = "ke_rlz_ewt"   # realizable k-epsilon + enhanced wall treatment (SST failed Gate A; see docstring)
RE_MIN = 3000.0               # Gnielinski lower limit; stops typos such as --re 1000
WF_MESH = "M2WF"          # coarse-wall mesh for standard wall functions (y+ about 30 to 70)


def make_tag(re_, mesh, turb=None, prt=None, prt_wall=None, test=False, tw_field=TW_FIELD):
    """Output tag. SST keeps the original bare tag (Re10000_M2) so earlier files are not renamed; every other
    model, including the default ke_rlz_ewt, gets a suffix (Re10000_M2_ke_rlz_ewt)."""
    turb = DEFAULT_TURB if turb is None else turb
    tag = f"Re{int(re_)}_{mesh}"
    if turb != "sst":
        tag += f"_{turb}"
    if prt is not None:
        tag += f"_prt{prt:g}"
    if prt_wall is not None:
        tag += f"_prtw{prt_wall:g}"
    if tw_field != TW_FIELD:
        tag += "_legacyTw"
    return tag + ("_TEST" if test else "")


def configure_turbulence(v, turb, prt, prt_wall, log):
    """Set the viscous model and read the state back; raise if anything did not take."""
    cfg = TURB_PRESETS[turb]
    if cfg["model"] == "k-omega":
        stt = v.get_state()
        if stt.get("model") != "k-omega":
            v.model = "k-omega"
        if stt.get("k_omega_model") != "sst":
            v.k_omega_model = "sst"
    else:
        v.model = "k-epsilon"
        v.k_epsilon_model = cfg["variant"]
        v.near_wall_treatment.wall_treatment = cfg["wall"]
    if prt is not None:
        v.options.turbulent_prandtl_numbers.energy_prandtl_number.value = prt
    if prt_wall is not None:
        v.options.turbulent_prandtl_numbers.wall_prandtl_number.value = prt_wall
    got = v.get_state()
    problems = []
    if got.get("model") != cfg["model"]:
        problems.append(f"model is {got.get('model')!r}, wanted {cfg['model']!r}")
    if cfg["model"] == "k-omega":
        if got.get("k_omega_model") != "sst":
            problems.append(f"k_omega_model is {got.get('k_omega_model')!r}")
    else:
        if got.get("k_epsilon_model") != cfg["variant"]:
            problems.append(f"k_epsilon_model is {got.get('k_epsilon_model')!r}, wanted {cfg['variant']!r}")
        wt = (got.get("near_wall_treatment") or {}).get("wall_treatment")
        if wt != cfg["wall"]:
            problems.append(f"wall_treatment is {wt!r}, wanted {cfg['wall']!r}")
    pt = ((got.get("options") or {}).get("turbulent_prandtl_numbers")) or {}
    e_val = (pt.get("energy_prandtl_number") or {}).get("value")
    w_val = (pt.get("wall_prandtl_number") or {}).get("value")
    if prt is not None and e_val != prt:
        problems.append(f"energy Pr_t is {e_val!r}, wanted {prt!r}")
    if prt_wall is not None and w_val != prt_wall:
        problems.append(f"wall Pr_t is {w_val!r}, wanted {prt_wall!r}")
    if problems:
        raise RuntimeError("; ".join(problems))
    wt = (got.get("near_wall_treatment") or {}).get("wall_treatment")
    log(f"    verified: model={got.get('model')} variant={got.get('k_epsilon_model') or got.get('k_omega_model')} "
        f"wall_treatment={wt} Pr_t(energy)={e_val} Pr_t(wall)={w_val}")


# ------------------------------------------------------------------- pure maths
def friction_factor(re):
    return (1.84 * math.log10(re) - 1.64) ** -2


def gnielinski(re, pr):
    f = friction_factor(re)
    return (f / 8.0) * (re - 1000.0) * pr / (1.0 + 12.7 * math.sqrt(f / 8.0) * (pr ** (2.0 / 3.0) - 1.0))


def case_numbers(re):
    p = WATER
    pr = p["mu"] * p["cp"] / p["k"]
    v = re * p["mu"] / (p["rho"] * D_I)
    area = math.pi * R_I ** 2
    mdot = p["rho"] * v * area
    nu_g = gnielinski(re, pr)
    h = nu_g * p["k"] / D_I
    q_in = DT_TARGET * h                 # flux on the inner wall
    q_out = q_in * D_I / D_O             # flux to apply on the outer wall
    q_total = q_out * math.pi * D_O * 2.0
    return dict(re=re, pr=pr, v=v, mdot=mdot, f=friction_factor(re), nu_gnielinski=nu_g,
                h_expected=h, q_in_expected=q_in, q_out=q_out, q_total=q_total,
                dTb_expected=q_total / (mdot * p["cp"]))


_FLOAT = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def to_number(x):
    """Pull one float out of whatever a Fluent report command returned."""
    if x is None or isinstance(x, bool):
        return None
    if isinstance(x, (int, float)):
        return float(x)
    if isinstance(x, dict):
        for v in reversed(list(x.values())):
            n = to_number(v)
            if n is not None:
                return n
        return None
    if isinstance(x, (list, tuple)):
        for v in reversed(x):
            n = to_number(v)
            if n is not None:
                return n
        return None
    if isinstance(x, str):
        m = _FLOAT.findall(x)
        return float(m[-1]) if m else None
    return None


def parse_report_file(path):
    """Last numeric token of the last parseable line of a Fluent report file."""
    for line in reversed(Path(path).read_text(errors="ignore").splitlines()):
        toks = line.split()
        if not toks:
            continue
        try:
            return float(toks[-1])
        except ValueError:
            continue
    return None


def check_converged(hist, n_blocks=3):
    """hist: one dict per 100-iteration block with keys T (list), Nu, eb.
    Spec 6.8: temperatures < 0.01 %, Nu < 0.05 % over 300 its, energy balance < 0.2 %."""
    if len(hist) <= n_blocks:
        return False, {}
    a, b = hist[-1], hist[-1 - n_blocks]
    dT = max(abs(x - y) / abs(y) for x, y in zip(a["T"], b["T"]))
    dNu = abs(a["Nu"] - b["Nu"]) / abs(b["Nu"])
    ok = dT < TOL_T and dNu < TOL_NU and abs(a["eb"]) < TOL_EB
    return ok, dict(dT_rel=dT, dNu_rel=dNu, eb=a["eb"])


def run_iterations(iterate, after_warm, snap, set_scale, log,
                   first=ENERGY_FIRST, max_iter=MAX_ITER, block=BLOCK, warm=50, probe=None, scale2=None):
    """Stage 1: scale factor 1 up to `first` its (snapshots only in the last 3 blocks).
    Stage 2: scale factor 5, blocks of 100 until converged or max_iter."""
    total = 0
    if warm:
        iterate(warm)
        total += warm
    after_warm()
    hist = []
    while total < first:
        n = min(block, first - total)
        iterate(n)
        total += n
        if probe:
            probe(total)
        else:
            log(f"  iterations: {total}")
    if scale2 is not None:      # raising the pseudo-time factor destabilised the flow in the first Gate A run
        set_scale(scale2)
    hist = []
    converged = False
    while total < max_iter:
        iterate(block)
        total += block
        s = snap()
        hist.append(s)
        ok, info = check_converged(hist)
        log(f"  it {total}: Nu_win={s['Nu']:.4f} Tout={s['T'][-1]:.5f} eb={s['eb']*100:.4f}% {info}")
        if ok:
            converged = True
            break
    return total, converged, hist


# ----------------------------------------------------------------- Fluent side
def setv(parent, name, val):
    """Set a scalar that may be a plain number or an {'option','value'} group."""
    child = getattr(parent, name)
    st = None
    try:
        st = child.get_state()
    except Exception:  # noqa: BLE001
        pass
    if isinstance(st, dict) and "value" in st:
        child.value = val
    else:
        setattr(parent, name, val)


def run_flow_stage(iterate, flow_probe, log, block=BLOCK, min_its=FLOW_MIN, max_its=FLOW_MAX,
                   tol=FLOW_TOL, mass_tol=FLOW_MASS_TOL, p_limit=1e4):
    """Iterate with the temperature equation off. Returns (its, status, hist);
    status is 'settled', 'diverged' or 'slow' (stable but not settled by max_its).
    Settled = inlet static pressure steady (< tol over 300 its) AND outlet/inlet mass flow within mass_tol.
    Diverged = non-finite or |p_in| > 1e5, or after 500 its |p_in| > p_limit (default 1e4; the driver passes 20 x the Darcy estimate), or after 500 its
    the ratio -m_out/m_in is off by more than 50 %.
    SIGN CONVENTION (verified from Fluent's transcript): Fluent counts flow INTO the domain as positive, so a
    healthy outlet reads NEGATIVE (m_out = -m_in) and a POSITIVE m_out means reverse flow at the outlet."""
    hist, its = [], 0
    while its < max_its:
        iterate(block)
        its += block
        r = flow_probe()
        hist.append(r["p_in"])
        ratio = -r["m_out"] / r["m_in"] if r["m_in"] else float("nan")   # healthy = +1
        log(f"  flow it {its}: p_in={r['p_in']:.3f} Pa  m_in={r['m_in']:.5f}  m_out={r['m_out']:.5f}"
            f"  (out/in={ratio:.3f})  T_out={r['T_out']:.4f}")
        if not math.isfinite(r["p_in"]) or abs(r["p_in"]) > 1e5:
            return its, "diverged", hist
        if its >= 500 and (abs(r["p_in"]) > p_limit or not math.isfinite(ratio) or abs(ratio - 1.0) > 0.5):
            return its, "diverged", hist
        if its >= min_its and len(hist) > 3 and abs(hist[-4]) > 0:
            if abs(hist[-1] - hist[-4]) / abs(hist[-4]) < tol and abs(ratio - 1.0) < mass_tol:
                return its, "settled", hist
    return its, "slow", hist


def set_stability(s, relax, scale):
    pv = s.solution.controls.p_v_controls
    pv.explicit_momentum_under_relaxation = relax
    pv.explicit_pressure_under_relaxation = relax
    s.solution.run_calculation.pseudo_time_settings.time_step_method.time_step_size_scale_factor = scale


def set_temperature_equation(s, solver, on, log):
    """Switch the energy (temperature) equation solve on/off. Several routes; first that verifies wins."""
    eq = s.solution.controls.equations
    attempts = (
        ("equations['temperature'] = ...", lambda: eq.__setitem__("temperature", on)),
        ("equations.set_state", lambda: eq.set_state({"temperature": on})),
        ("equations.temperature = ...", lambda: setattr(eq, "temperature", on)),
        ("TUI solve/set/equations/temperature",
         lambda: solver.tui.solve.set.equations.temperature("yes" if on else "no")),
    )
    for name, fn in attempts:
        try:
            fn()
            stt = eq.get_state()
            if isinstance(stt, dict) and stt.get("temperature") == on:
                log(f"  temperature equation -> {'ON' if on else 'OFF'} (via {name})")
                return
            log(f"  attempt '{name}' ran but state is {stt}")
        except Exception as e:  # noqa: BLE001
            log(f"  attempt '{name}' failed: {str(e)[:150]}")
    raise RuntimeError("could not switch the temperature equation; see attempts above")


class Steps:
    def __init__(self, log):
        self.log = log
        self.fails = []

    def run(self, name, fn):
        try:
            fn()
            self.log(f"  OK    {name}")
        except Exception as e:  # noqa: BLE001
            self.log(f"  FAIL  {name}: {type(e).__name__}: {str(e)[:300]}")
            self.fails.append(name)


def prepare_mesh(solver, mesh, log):
    s = solver.settings
    s.file.read_case(file_name=str(mesh))
    s.setup.general.solver.two_dim_space = "axisymmetric"
    solver.tui.mesh.modify_zones.slit_interior_between_diff_materials()
    s.setup.models.energy.enabled = True
    log("Mesh read, axisymmetric, interface slit, energy on.")


def apply_setup(solver, nums, log, relax=0.75, scale1=1.0, turb=None, prt=None, prt_wall=None):
    turb = DEFAULT_TURB if turb is None else turb
    s = solver.settings
    bc = s.setup.boundary_conditions
    mats = s.setup.materials
    st = Steps(log)
    p = WATER

    def viscous():
        configure_turbulence(s.setup.models.viscous, turb, prt, prt_wall, log)
    st.run("viscous model " + ("k-omega SST" if turb == "sst" else turb), viscous)

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
        setv(vi.momentum, "velocity_magnitude", nums["v"])
        setv(vi.thermal, "temperature", T_IN)
        vi.turbulence.turbulence_specification = "Intensity and Hydraulic Diameter"
        vi.turbulence.turbulent_intensity = TURB_INTENSITY
        setv(vi.turbulence, "hydraulic_diameter", D_I)
    st.run("velocity inlet", inlet)

    def outlet():
        po = bc.pressure_outlet["outlet"]
        setv(po.momentum, "gauge_pressure", 0.0)
        po.turbulence.turbulence_specification = "Intensity and Hydraulic Diameter"
        po.turbulence.backflow_turbulent_intensity = TURB_INTENSITY
        setv(po.turbulence, "backflow_hydraulic_diameter", D_I)
        setv(po.thermal, "backflow_total_temperature", T_IN)
    st.run("pressure outlet", outlet)

    def heated():
        th = bc.wall["outer_heated"].thermal
        th.thermal_condition = "Heat Flux"
        setv(th, "heat_flux", nums["q_out"])
    st.run("outer_heated heat flux", heated)

    def adiabatic():
        for n in ("outer_unheated", "solid_end_in", "solid_end_out"):
            th = bc.wall[n].thermal
            th.thermal_condition = "Heat Flux"
            setv(th, "heat_flux", 0.0)
    st.run("adiabatic walls", adiabatic)

    relax_val = relax

    def relax_fn():
        pv = s.solution.controls.p_v_controls
        pv.explicit_momentum_under_relaxation = relax_val
        pv.explicit_pressure_under_relaxation = relax_val
    st.run(f"explicit relaxation {relax_val}", relax_fn)

    def residuals():
        eq = s.solution.monitor.residual.equations
        second = "omega" if TURB_PRESETS[turb]["model"] == "k-omega" else "epsilon"
        for n, t in (("continuity", 1e-6), ("x-velocity", 1e-6), ("y-velocity", 1e-6),
                     ("k", 1e-6), (second, 1e-6), ("energy", 1e-9)):
            eq[n].absolute_criteria = t
            eq[n].check_convergence = False   # we control stopping, not Fluent
    st.run("residual targets, auto-stop OFF", residuals)

    def discretisation():
        # Spec 6.7: second order for pressure, second-order upwind for everything else. Fluent's DEFAULT for k and
        # epsilon is FIRST-order upwind (seen in case_settings_Re10000_M2_ke_rlz_ewt.txt), so set every equation
        # explicitly and read it back.
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


def initialise(solver, log):
    s = solver.settings
    s.solution.initialization.hybrid_initialize()
    patch = s.solution.initialization.patch.calculate_patch
    tries = ({"cell_zones": ["fluid", "solid"], "variable": "temperature", "value": T_IN},
             {"domain": "mixture", "cell_zones": ["fluid", "solid"], "variable": "temperature",
              "value": T_IN})
    for kw in tries:
        try:
            patch(**kw)
            log(f"  patch T={T_IN} OK with {sorted(kw)}")
            return
        except Exception as e:  # noqa: BLE001
            log(f"  patch attempt {sorted(kw)} failed: {str(e)[:200]}")
    log("  WARNING: temperature patch failed; continuing with hybrid-init temperatures.")


class Post:
    def __init__(self, solver, tmpfile, log, tw_field=TW_FIELD):
        self.tw_field = tw_field
        self.s = solver.settings
        self.surf = self.s.results.surfaces
        self.rep = self.s.results.report
        self.tmp = str(tmpfile)
        self.log = log
        self.fluid_wall = None
        self.first_return_logged = False
        self.file_mode = False

    # -- generic report call with file fallback
    def _scalar(self, cmd, **kw):
        v = None
        if not self.file_mode:
            r = cmd(**kw)
            if not self.first_return_logged:
                self.log(f"  [first report return: type={type(r).__name__}, value={str(r)[:120]!r}]")
                self.first_return_logged = True
            v = to_number(r)
            if v is None:
                self.file_mode = True   # Fluent 2026 R1 returns None: read via report file from now on
        if v is None:
            cmd(**kw, write_to_file=True, file_name=self.tmp, append_data=False)
            v = parse_report_file(self.tmp)
        if v is None or not math.isfinite(v):
            raise RuntimeError(f"could not read a number from report {kw}")
        return v

    def aw(self, surface, field):
        return self._scalar(self.rep.surface_integrals.area_weighted_avg,
                            surface_names=[surface], report_of=field)

    def mw(self, surface, field):
        return self._scalar(self.rep.surface_integrals.mass_weighted_avg,
                            surface_names=[surface], report_of=field)

    def flux(self, which, zone):
        cmd = getattr(self.rep.fluxes, which)
        last = None
        for kw in ({"zones": [zone]}, {"locations": [zone]}):
            try:
                return self._scalar(cmd, **kw)
            except Exception as e:  # noqa: BLE001
                last = e
        raise last

    # -- surfaces
    def _make(self, cont, name):
        if name in list(cont.keys()):
            return cont[name]
        cont.create(name)
        return cont[name]

    def make_plane(self, name, x):
        o = self._make(self.surf.iso_surface, name)
        o.field = "x-coordinate"
        o.iso_values = [x]
        o.zones = ["fluid"]

    def make_band(self, name, x0, x1):
        o = self._make(self.surf.iso_clip, name)
        o.field = "x-coordinate"
        o.surfaces = [self.fluid_wall]
        try:
            o.range.minimum = x0
            o.range.maximum = x1
        except Exception:  # noqa: BLE001
            o.range = {"minimum": x0, "maximum": x1}

    def detect_fluid_wall(self):
        yp = {n: self.aw(n, "y-plus") for n in ("interface_fs", "interface_fs-shadow")}
        self.log(f"  area-avg y+ on interface walls: {yp}")
        best = max(yp, key=lambda k: yp[k])
        if yp[best] <= 0:
            raise RuntimeError("y+ is zero on both interface walls; cannot tell the fluid side")
        self.fluid_wall = best
        self.log(f"  fluid-side wall = {best}")

    def build_surfaces(self):
        for i in range(N_BANDS):
            xc = X_HEAT + BAND * i + BAND / 2
            self.make_plane(f"pl_{i:02d}", xc)
            self.make_band(f"bd_{i:02d}", X_HEAT + BAND * i, X_HEAT + BAND * (i + 1))
        self.make_plane("pl_f12", 1.2)
        self.make_plane("pl_f20", 2.0)
        self.log(f"  created {N_BANDS} planes, {N_BANDS} bands, 2 friction planes")

    # -- quantities
    def band_values(self, i):
        tb = self.mw(f"pl_{i:02d}", "temperature")
        tw = self.aw(f"bd_{i:02d}", self.tw_field)
        q = abs(self.aw(f"bd_{i:02d}", "heat-flux"))
        nu = (q / (tw - tb)) * D_I / WATER["k"] if abs(tw - tb) > 1e-9 else float("nan")
        return dict(i=i, x=X_HEAT + BAND * i + BAND / 2, Tb=tb, Tw=tw, q=q, Nu=nu)

    def probe(self, total):
        """Cheap stage-1 progress line and divergence guard."""
        t_out = self.mw("outlet", "temperature")
        self.log(f"  iterations: {total}   T_out = {t_out:.4f} K")
        if total >= 500 and not (T_IN - 2.0 < t_out < T_IN + 30.0):
            raise RuntimeError(f"T_out = {t_out:.2f} K at {total} its is unphysical -> diverged or badly stuck")

    def flow_probe(self):
        return dict(p_in=self.aw("inlet", "pressure"), m_in=self.flux("mass_flow", "inlet"),
                    m_out=self.flux("mass_flow", "outlet"), T_out=self.mw("outlet", "temperature"))

    def energy_balance(self):
        t_in = self.mw("inlet", "temperature")
        t_out = self.mw("outlet", "temperature")
        mdot = abs(self.flux("mass_flow", "inlet"))
        q_wall = abs(self.flux("heat_transfer", "outer_heated"))
        q_out = mdot * WATER["cp"] * (t_out - t_in)
        return dict(T_in_meas=t_in, T_out=t_out, mdot_fluent=mdot, Q_wall_fluent=q_wall,
                    Q_enthalpy=q_out, eb=(q_wall - q_out) / q_wall if q_wall else float("nan"))

    def snapshot(self):
        """Cheap per-block monitor: window bands only."""
        rows = [self.band_values(i) for i in range(WIN_FIRST, WIN_LAST + 1)]
        eb = self.energy_balance()
        tw = {r["i"]: r["Tw"] for r in rows}
        snap = dict(T=[tw[i] for i in MONITOR_BANDS] + [eb["T_out"]],
                    Nu=sum(r["Nu"] for r in rows) / len(rows), eb=eb["eb"], rows=rows, bal=eb)
        for v in snap["T"] + [snap["Nu"], snap["eb"]]:
            if not math.isfinite(v):
                raise RuntimeError("non-finite value in monitor -> run has diverged")
        if not (250.0 < snap["T"][-1] < 400.0):
            raise RuntimeError(f"outlet temperature {snap['T'][-1]} K out of range -> diverged")
        return snap

    def y_plus_stats(self):
        out = {"yplus_avg": self.aw(self.fluid_wall, "y-plus"), "yplus_max": None, "yplus_max_window": None}
        fm = getattr(self.rep.surface_integrals, "facet_max", None)
        if fm is not None:
            try:
                out["yplus_max"] = self._scalar(fm, surface_names=[self.fluid_wall], report_of="y-plus")
            except Exception as e:  # noqa: BLE001
                self.log(f"  facet_max for y+ failed: {str(e)[:150]}")
            try:
                out["yplus_max_window"] = max(
                    self._scalar(fm, surface_names=[f"bd_{i:02d}"], report_of="y-plus")
                    for i in range(WIN_FIRST, WIN_LAST + 1))
            except Exception as e:  # noqa: BLE001
                self.log(f"  facet_max for y+ on window bands failed: {str(e)[:150]}")
        else:
            self.log("  surface_integrals.facet_max not found: y+ max not available.")
        return out


# ----------------------------------------------------------------------- driver
def finalise(post, nums, total, converged, tag, results_dir, log, extra=None):
    prof = [post.band_values(i) for i in range(N_BANDS)]
    win = [r for r in prof if WIN_FIRST <= r["i"] <= WIN_LAST]
    nu_win = sum(r["Nu"] for r in win) / len(win)
    nu_heated = sum(r["Nu"] for r in prof) / len(prof)
    spread = (max(r["Nu"] for r in win) - min(r["Nu"] for r in win)) / nu_win
    bal = post.energy_balance()
    dp = post.aw("pl_f12", "pressure") - post.aw("pl_f20", "pressure")
    f_cfd = 2.0 * D_I * dp / (WATER["rho"] * nums["v"] ** 2 * 0.8)
    qin_dev = max(abs(r["q"] - nums["q_in_expected"]) / nums["q_in_expected"] for r in win)
    yp = post.y_plus_stats()
    err = nu_win / nums["nu_gnielinski"] - 1.0
    res = dict(tag=tag, Re=nums["re"], converged=converged, iterations=total,
               Nu_window=nu_win, Nu_heated_avg=nu_heated, Nu_gnielinski=nums["nu_gnielinski"],
               gateA_error_pct=100 * err, gateA_pass=abs(err) <= 0.10,
               window_spread_pct=100 * spread, window_flat_flag=spread <= 0.02,
               energy_balance_pct=100 * bal["eb"], fluent_mdot_over_expected=bal["mdot_fluent"] / nums["mdot"],
               T_out=bal["T_out"], dTb_expected=nums["dTb_expected"], dTb_cfd=bal["T_out"] - bal["T_in_meas"],
               q_in_max_dev_from_expected_pct=100 * qin_dev, f_cfd=f_cfd, f_correlation=nums["f"],
               delta_p_1p2_to_2p0=dp, **yp, **(extra or {}))
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / f"gateA_{tag}.json").write_text(json.dumps(res, indent=2))
    with open(results_dir / f"gateA_{tag}_profile.csv", "w") as fh:
        fh.write("band,x_centre,Tb,Tw,q_in,Nu\n")
        for r in prof:
            fh.write(f"{r['i']},{r['x']:.4f},{r['Tb']:.6f},{r['Tw']:.6f},{r['q']:.4f},{r['Nu']:.4f}\n")
    log("\n==== RESULT ====")
    for k, v in res.items():
        log(f"  {k}: {v}")
    log(f"\nGate A ({tag}): {'PASS' if res['gateA_pass'] else 'FAIL'}  "
        f"(Nu_CFD {nu_win:.2f} vs Gnielinski {nums['nu_gnielinski']:.2f}, {100*err:+.2f} %)")
    ratio = nu_win / nums["nu_gnielinski"]
    if not (0.5 < ratio < 2.0):
        log(f"WARNING: Nu is {ratio:.2f} x Gnielinski. That is not a turbulence-model difference; suspect the "
            f"post-processing (wall temperature field '{post.tw_field}') before trusting this number.")
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--re", type=float, default=10000.0)
    ap.add_argument("--mesh", default="M2", choices=["M1", "M2", "M3", "B3D", WF_MESH])
    ap.add_argument("--tw-field", default=TW_FIELD, choices=[TW_FIELD, "temperature"],
                    help="field used for the inner-wall temperature. Default 'wall-temperature' is correct for every "
                         "model; 'temperature' reproduces the OLD (wrong for wall functions) results")
    ap.add_argument("--turb", default=DEFAULT_TURB, choices=sorted(TURB_PRESETS),
                    help=f"turbulence model preset (default {DEFAULT_TURB}; 'sst' = the spec's model). "
                         "ke_*_swf needs --mesh M2WF")
    ap.add_argument("--prt", type=float, default=None,
                    help="DIAGNOSTIC: constant energy turbulent Prandtl number (spec: 0.85)")
    ap.add_argument("--prt-wall", type=float, default=None, help="DIAGNOSTIC: constant wall turbulent Prandtl number")
    ap.add_argument("--test", action="store_true", help="short plumbing check only")
    ap.add_argument("--no-ladder", action="store_true", help="use only --relax/--scale1, no automatic retries")
    ap.add_argument("--no-transcript", action="store_true", help="do not start a Fluent transcript")
    ap.add_argument("--relax", type=float, default=DEFAULT_RELAX,
                    help="explicit momentum/pressure relaxation (spec says 0.75, which diverges here)")
    ap.add_argument("--scale1", type=float, default=DEFAULT_SCALE1,
                    help="pseudo-time scale factor, first stages (spec says 1, which diverges here)")
    ap.add_argument("--scale2", type=float, default=None,
                    help="pseudo-time factor to switch to after the first energy block. Default: do not change it "
                         "(switching 0.1 -> 0.5 made the first Gate A run diverge)")
    a = ap.parse_args()

    tag = make_tag(a.re, a.mesh, a.turb, a.prt, a.prt_wall, a.test, a.tw_field)
    results_dir = ROOT / "05_results"
    results_dir.mkdir(parents=True, exist_ok=True)
    logfile = results_dir / f"gateA_{tag}_log.txt"
    lines = []

    def log(msg=""):
        print(msg, flush=True)
        lines.append(str(msg))
        logfile.write_text("\n".join(lines), encoding="utf-8")

    mesh = ROOT / "02_geometry_mesh" / "2d" / f"pipe2d_{a.mesh}.msh"
    if not mesh.exists():
        log(f"Mesh not found: {mesh}")
        return 1
    if a.re < RE_MIN:
        log(f"STOP: Re = {a.re:g} is below {RE_MIN:g}. Gate A needs turbulent flow (Gnielinski is valid from "
            f"Re about 3000; the spec asks for 10000 and 25000). Did you mean --re {a.re * 10:g}?")
        return 1
    cfg = TURB_PRESETS[a.turb]
    if cfg["wall_function"] and a.mesh != WF_MESH:
        log(f"STOP: {a.turb} uses wall functions, which need y+ of about 30 or more, but mesh {a.mesh} has a 5 um "
            f"first layer (y+ < 1). Use --mesh {WF_MESH} (python make_2d_mesh.py M2 pipe2d_{WF_MESH}.msh "
            f"--first-layer 2.2e-3 --growth 1.15 --layers 2).")
        return 1
    if a.mesh == WF_MESH and not cfg["wall_function"]:
        log(f"WARNING: mesh {WF_MESH} is a coarse-wall mesh (y+ about 30 to 70) but {a.turb} expects a resolved wall; "
            "results are not meaningful for SST / enhanced wall treatment validation.")
    nums = case_numbers(a.re)
    log(f"Turbulence: {a.turb} {cfg}" + (f"  Pr_t energy={a.prt}" if a.prt is not None else "")
        + (f"  Pr_t wall={a.prt_wall}" if a.prt_wall is not None else ""))
    log(f"Gate A {tag}   mesh {mesh.name}   relax={a.relax} scale1={a.scale1} scale2={a.scale2 if a.scale2 is not None else '(unchanged)'}")
    log("Expected: " + ", ".join(f"{k}={v:.5g}" for k, v in nums.items()))
    dp_darcy = nums["f"] * (4.0 / D_I) * WATER["rho"] * nums["v"] ** 2 / 2.0
    log(f"Expected inlet static pressure (gauge, outlet = 0): about {dp_darcy:.0f} Pa (Darcy; developing flow adds a bit)")

    import ansys.fluent.core as pyfluent  # lazy so the maths can be unit-tested anywhere
    transcript = {}
    t0 = time.time()
    solver = pyfluent.launch_fluent(mode="solver", dimension=2, precision="double", processor_count=PROCS)
    try:
        # TUI string arguments are passed RAW (split at spaces), and the project path contains
        # spaces, so write the transcript to a no-space temp path and copy it afterwards.
        if not a.no_transcript:
            try:
                tmp_t = Path(tempfile.gettempdir()) / f"fluent_transcript_{tag}.txt"
                arg = f'"{tmp_t}"' if " " in str(tmp_t) else str(tmp_t)
                solver.tui.file.start_transcript(arg)
                transcript["tmp"] = tmp_t
                log(f"Fluent transcript -> {tmp_t} (copied to 05_results at the end)")
            except Exception as e:  # noqa: BLE001
                log(f"(transcript not started: {str(e)[:150]}). If later steps fail with "
                    "'An error or interrupt occurred', rerun with --no-transcript.")

        prepare_mesh(solver, mesh, log)
        log("Applying setup:")
        steps = apply_setup(solver, nums, log, relax=a.relax, scale1=a.scale1,
                            turb=a.turb, prt=a.prt, prt_wall=a.prt_wall)
        if steps.fails:
            log(f"\nSETUP FAILED for: {steps.fails}\nNo iterations were run. Send me this log.")
            return 2
        s = solver.settings
        log("Initialising (flow-first: temperature equation OFF):")
        set_temperature_equation(s, solver, False, log)
        initialise(solver, log)
        post = Post(solver, results_dir / "_tmp_report.txt", log, tw_field=a.tw_field)
        iterate = lambda n: s.solution.run_calculation.iterate(iter_count=int(n))  # noqa: E731

        def set_scale(x):
            s.solution.run_calculation.pseudo_time_settings.time_step_method.time_step_size_scale_factor = x
            log(f"  pseudo-time scale factor -> {x}; convergence history restarted")

        log("Warm-up 50 iterations, then post-processing setup:")
        iterate(50)
        post.detect_fluid_wall()
        post.build_surfaces()
        total = 50

        log("Flow stage (temperature equation OFF):")
        used = (a.relax, a.scale1)
        if a.test:
            iterate(100)
            total += 100
            log(f"  test flow probe: {post.flow_probe()}")
            fl_its = 100
        else:
            recipes = [(a.relax, a.scale1)]
            if not a.no_ladder:
                recipes += [r for r in LADDER if r[0] <= a.relax and r[1] < a.scale1]
            used, fl_its = None, 0
            for k, (rx, sc) in enumerate(recipes):
                log(f"Flow attempt {k+1}/{len(recipes)}: explicit relaxation {rx}, pseudo-time scale {sc}")
                if k > 0:
                    set_stability(s, rx, sc)
                    initialise(solver, log)   # fresh start; temperature equation is still OFF
                its, status, hist = run_flow_stage(iterate, post.flow_probe, log, p_limit=max(1e4, 20.0 * dp_darcy))
                total += its
                fl_its = its
                if status == "settled":
                    used = (rx, sc)
                    log(f"  flow settled after {its} its with relax={rx}, scale={sc}")
                    p_last = hist[-1]
                    if not (0.5 * dp_darcy < p_last < 2.0 * dp_darcy):
                        log(f"  WARNING: settled p_in = {p_last:.1f} Pa is far from the Darcy estimate "
                            f"{dp_darcy:.0f} Pa; check the setup before trusting Nu.")
                    break
                log(f"  attempt {k+1} ended as '{status}' after {its} its "
                    f"(p_in history: {[round(h, 1) for h in hist[-6:]]})")
            if used is None:
                log("\nALL FLOW ATTEMPTS FAILED. The flow solver does not reach a steady solution even with "
                    "conservative settings, so suspect the mesh or a boundary condition. "
                    "Send this log and the transcript file.")
                return 4

        log("Energy stage (temperature equation ON):")
        set_temperature_equation(s, solver, True, log)

        if a.test:
            iterate(100)
            snap = post.snapshot()
            log(f"  snapshot: Nu_win={snap['Nu']:.3f}  energy balance={snap['eb']*100:.3f}%  T_out={snap['T'][-1]:.4f}")
            log(f"  balance detail: {snap['bal']}")
            log(f"  sample band: {snap['rows'][0]}")
            log(f"  y+: {post.y_plus_stats()}")
            log("TEST DONE.")
            return 0

        e_total, converged, _ = run_iterations(iterate, lambda: None, post.snapshot, set_scale, log,
                                               warm=0, probe=post.probe,
                                               scale2=None if a.scale2 is None else min(a.scale2, 5.0 * used[1]))
        total += e_total
        log(f"Finished: {total} iterations (flow {fl_its} + energy {e_total}), converged={converged}, "
            f"wall time {(time.time()-t0)/60:.1f} min")
        res = finalise(post, nums, total, converged, tag, results_dir, log,
                 extra=dict(recipe_relax=used[0], recipe_scale1=used[1], flow_iterations=fl_its,
                          turb=a.turb, mesh=a.mesh, prt_energy=a.prt, prt_wall=a.prt_wall,
                          tw_field=a.tw_field))
        if cfg["wall_function"]:
            yp_avg = res.get("yplus_avg")
            if yp_avg is not None and not (WF_MIN_YPLUS <= yp_avg <= 300.0):
                log(f"WARNING: average wall y+ = {yp_avg:.1f} is outside {WF_MIN_YPLUS} to 300 (log-law range "
                    "for standard wall functions); this result is NOT valid. Tell me this number.")
            else:
                note = "  (below the usual 30 but inside the log-law range)" if yp_avg < 30.0 else ""
                log(f"Wall-function check: average wall y+ = {yp_avg:.1f}, max = {res.get('yplus_max')}{note}")
        case_dir = ROOT / "04_cases" / "gateA"
        case_dir.mkdir(parents=True, exist_ok=True)
        s.file.write_case_data(file_name=str(case_dir / f"gateA_{tag}.cas.h5"))
        log(f"Saved case and data to {case_dir}")
        return 0
    except Exception as e:  # noqa: BLE001
        log(f"\nABORTED: {type(e).__name__}: {e}")
        return 3
    finally:
        if transcript.get("tmp"):
            try:
                solver.tui.file.stop_transcript()
            except Exception:  # noqa: BLE001
                pass
        try:
            solver.exit()
        except Exception:  # noqa: BLE001
            pass
        if transcript.get("tmp") and Path(transcript["tmp"]).exists():
            try:
                shutil.copy(transcript["tmp"], results_dir / f"gateA_{tag}_fluent_transcript.txt")
            except Exception:  # noqa: BLE001
                pass


if __name__ == "__main__":
    sys.exit(main())