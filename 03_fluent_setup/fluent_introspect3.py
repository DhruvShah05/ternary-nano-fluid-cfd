"""
Round 3 introspection (Fluent 2026 R1, PyFluent 0.43.0).

Purpose: resolve the names and values that are still unknown before writing
the Gate A setup script. Nothing here is a real case setup; it only probes.

Run from the project root (PowerShell):
    .\.venv\Scripts\Activate.ps1
    python 03_fluent_setup\fluent_introspect3.py
Optional: pass a specific mesh file as the first argument.

Output: 03_fluent_setup\introspect_dump3.txt  (rewritten after every line,
so a crash still leaves a usable log). Upload that file.
"""
import sys
import pprint
from pathlib import Path

import ansys.fluent.core as pyfluent

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OUT = HERE / "introspect_dump3.txt"
_lines = []


def log(msg=""):
    print(msg)
    _lines.append(str(msg))
    OUT.write_text("\n".join(_lines), encoding="utf-8")


def section(title, fn):
    log(f"\n##### {title}")
    try:
        fn()
    except Exception as e:  # noqa: BLE001
        log(f"FAILED: {type(e).__name__}: {str(e)[:400]}")


def get(obj, name):
    """Return attribute; call it if it is a method."""
    v = getattr(obj, name, None)
    if callable(v):
        try:
            v = v()
        except Exception as e:  # noqa: BLE001
            v = f"<call failed: {e}>"
    return v


def describe(obj, state=True):
    log(f"child_names: {get(obj, 'child_names')}")
    log(f"command_names: {get(obj, 'command_names')}")
    if state:
        try:
            log(pprint.pformat(obj.get_state(), width=140))
        except Exception as e:  # noqa: BLE001
            log(f"get_state failed: {e}")


def allowed(label, obj):
    log(f"{label} allowed_values: {get(obj, 'allowed_values')}")


def find_mesh():
    if len(sys.argv) > 1:
        return Path(sys.argv[1])
    cands = sorted((ROOT / "02_geometry_mesh" / "2d").rglob("*.msh"),
                   key=lambda p: p.stat().st_mtime)
    if not cands:
        raise SystemExit("No .msh found under 02_geometry_mesh\\2d. "
                         "Pass the mesh path as the first argument.")
    return cands[-1]


# ----------------------------------------------------------------- launch
mesh = find_mesh()
log(f"Mesh: {mesh}")

try:
    solver = pyfluent.launch_fluent(mode="solver", dimension=2,
                                    precision="double", processor_count=4)
except Exception as e:  # noqa: BLE001
    log(f"launch with dimension=2 failed ({e}); trying Dimension.TWO")
    solver = pyfluent.launch_fluent(mode="solver",
                                    dimension=pyfluent.Dimension.TWO,
                                    precision="double", processor_count=4)

s = solver.settings

# Sequence verified by smoke test v5: read -> axisymmetric -> slit.
s.file.read_case(file_name=str(mesh))
s.setup.general.solver.two_dim_space = "axisymmetric"
solver.tui.mesh.modify_zones.slit_interior_between_diff_materials()
s.setup.models.energy.enabled = True
bc = s.setup.boundary_conditions
log("Setup sequence done (read, axisymmetric, slit, energy on).")

# ------------------------------------------------- 1. boundary conditions
section("velocity_inlet turbulence_specification",
        lambda: allowed("turbulence_specification",
                        bc.velocity_inlet["inlet"].turbulence.turbulence_specification))
section("velocity_inlet momentum",
        lambda: describe(bc.velocity_inlet["inlet"].momentum))
section("pressure_outlet turbulence",
        lambda: describe(bc.pressure_outlet["outlet"].turbulence))
section("pressure_outlet momentum",
        lambda: describe(bc.pressure_outlet["outlet"].momentum))
section("wall outer_heated thermal_condition",
        lambda: allowed("thermal_condition",
                        bc.wall["outer_heated"].thermal.thermal_condition))
section("wall solid_end_out thermal",
        lambda: log(pprint.pformat(bc.wall["solid_end_out"].thermal.get_state())))
section("wall zone names",
        lambda: log(f"walls: {list(bc.wall.keys())}"))

# ---------------------------------------------------------- 2. materials
mats = s.setup.materials


def _mat_containers():
    log("fluid container:")
    describe(mats.fluid, state=False)
    log("solid container:")
    describe(mats.solid, state=False)
    log(f"existing fluids: {list(mats.fluid.keys())}")
    log(f"existing solids: {list(mats.solid.keys())}")


section("materials containers", _mat_containers)
section("materials.database copy_by_name args",
        lambda: log(f"argument_names: {get(mats.database.copy_by_name, 'argument_names')}"))


def _copy_copper():
    mats.database.copy_by_name(type="solid", name="copper")
    log(f"solids after copy: {list(mats.solid.keys())}")
    log(pprint.pformat(mats.solid["copper"].get_state(), width=140))


section("try: copy copper from database", _copy_copper)


def _create_nf():
    mats.fluid.create("nf")
    log(f"fluids after create: {list(mats.fluid.keys())}")
    log(pprint.pformat(mats.fluid["nf"].get_state(), width=140))
    allowed("density.option", mats.fluid["nf"].density.option)


section("try: create fluid 'nf'", _create_nf)

# ------------------------------------------------------------ 3. zones
section("cell_zone_conditions.solid['solid']",
        lambda: describe(s.setup.cell_zone_conditions.solid["solid"]))
section("cell zone names",
        lambda: log(f"fluid zones: {list(s.setup.cell_zone_conditions.fluid.keys())}; "
                    f"solid zones: {list(s.setup.cell_zone_conditions.solid.keys())}"))

# ---------------------------------------------------------- 4. solution
section("p_v_controls (flow_courant_number etc.)",
        lambda: log(pprint.pformat(s.solution.controls.p_v_controls.get_state())))
section("solution.methods.axisymmetric",
        lambda: describe(s.solution.methods.axisymmetric))
section("initialization.patch.calculate_patch args",
        lambda: log(f"argument_names: "
                    f"{get(s.solution.initialization.patch.calculate_patch, 'argument_names')}"))
section("file.export",
        lambda: describe(s.file.export, state=False))
section("residual criterion_type allowed",
        lambda: allowed("criterion_type",
                        s.solution.monitor.residual.options.criterion_type))

# ------------------------------------------------ 5. post-processing
section("hybrid_initialize (needed so field data exists)",
        lambda: s.solution.initialization.hybrid_initialize())

section("results.surfaces",
        lambda: describe(s.results.surfaces, state=False))


def _try_surface(kind):
    cont = getattr(s.results.surfaces, kind, None)
    if cont is None:
        log(f"{kind}: NOT PRESENT in results.surfaces")
        return
    name = f"tmp_{kind}"
    try:
        cont.create(name)
    except Exception as e:  # noqa: BLE001
        log(f"{kind}.create failed ({e}); trying item assignment")
        cont[name] = {}
    obj = cont[name]
    describe(obj)


for _kind in ("iso_surface", "iso_clip", "line_surface", "point_surface", "plane_surface"):
    section(f"surface type: {_kind}", lambda k=_kind: _try_surface(k))


def _report_fields():
    cmd = s.results.report.surface_integrals.area_weighted_avg
    ro = getattr(cmd, "report_of", None)
    vals = get(ro, "allowed_values") if ro is not None else None
    log(f"report_of is None: {ro is None}")
    if isinstance(vals, (list, tuple)):
        keys = ("y-plus", "yplus", "y-star", "ystar", "heat-flux", "temperature",
                "x-coordinate", "density", "pressure", "wall-shear")
        hits = [v for v in vals if any(k in str(v).lower() for k in keys)]
        log(f"total fields: {len(vals)}")
        log(f"matching fields: {hits}")
    else:
        log(f"allowed_values: {vals}")
    sn = getattr(cmd, "surface_names", None)
    log(f"surface_names allowed: {get(sn, 'allowed_values') if sn is not None else None}")


section("report_of / surface_names allowed values", _report_fields)

log("\nDONE. Upload introspect_dump3.txt")
try:
    solver.exit()
except Exception:  # noqa: BLE001
    pass