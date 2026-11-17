"""
fluent_introspect.py - dump the exact setting names/values Fluent 2026 R1 exposes through PyFluent,
so the Gate A setup script can be written with real names instead of guesses.
Changes nothing in the case. Run:
    python fluent_introspect.py pipe2d_M2.msh 4
Output: introspect_dump.txt  (UPLOAD THIS FILE, do not paste it).
"""
import sys
import pprint
import traceback
from fluent_common import *

arg = sys.argv[1] if len(sys.argv) > 1 else MESH_NAME_DEFAULT
nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
mesh = find_mesh(arg, __file__)
preflight(mesh)
solver = launch(nproc, gui=False)
load_2d_mesh(solver, mesh)
s = solver.settings
out = open("introspect_dump.txt", "w", encoding="utf-8")


def dump(title, getter, state=True):
    out.write(f"\n##### {title}\n")
    try:
        obj = getter()
        for attr in ("child_names", "command_names"):
            try:
                out.write(f"{attr}: {list(getattr(obj, attr))}\n")
            except Exception:
                pass
        if state:
            st = obj.get_state() if hasattr(obj, "get_state") else obj
            out.write(pprint.pformat(st, width=150, compact=True) + "\n")
    except Exception as e:
        out.write(f"FAILED: {type(e).__name__}: {str(e)[:300]}\n")
    out.flush()
    print("dumped:", title)


dump("setup.general", lambda: s.setup.general)
dump("setup.models", lambda: s.setup.models)
dump("setup.models.viscous", lambda: s.setup.models.viscous)
dump("setup.models.energy", lambda: s.setup.models.energy)
dump("setup.materials (names only)", lambda: s.setup.materials, state=False)
dump("materials.fluid['air']", lambda: s.setup.materials.fluid["air"])
dump("materials.solid['aluminum']", lambda: s.setup.materials.solid["aluminum"])
dump("setup.cell_zone_conditions (names only)", lambda: s.setup.cell_zone_conditions, state=False)
dump("cell_zone_conditions.fluid['fluid']", lambda: s.setup.cell_zone_conditions.fluid["fluid"])
dump("cell_zone_conditions.solid['solid']", lambda: s.setup.cell_zone_conditions.solid["solid"])
dump("setup.boundary_conditions (names only)", lambda: s.setup.boundary_conditions, state=False)
dump("bc.velocity_inlet['inlet']", lambda: s.setup.boundary_conditions.velocity_inlet["inlet"])
dump("bc.pressure_outlet['outlet']", lambda: s.setup.boundary_conditions.pressure_outlet["outlet"])
dump("bc.wall['outer_heated']", lambda: s.setup.boundary_conditions.wall["outer_heated"])
dump("bc.wall['outer_unheated']", lambda: s.setup.boundary_conditions.wall["outer_unheated"])
dump("bc.wall['interface_fs']", lambda: s.setup.boundary_conditions.wall["interface_fs"])
dump("bc.wall['interface_fs-shadow']", lambda: s.setup.boundary_conditions.wall["interface_fs-shadow"])
dump("bc.axis['axis']", lambda: s.setup.boundary_conditions.axis["axis"])
dump("solution (names only)", lambda: s.solution, state=False)
dump("solution.methods", lambda: s.solution.methods)
dump("solution.controls", lambda: s.solution.controls)
dump("solution.initialization (names only)", lambda: s.solution.initialization, state=False)
dump("solution.run_calculation (names only)", lambda: s.solution.run_calculation, state=False)
dump("solution.monitor (names only)", lambda: s.solution.monitor, state=False)
dump("results (names only)", lambda: s.results, state=False)
dump("results.surfaces (names only)", lambda: s.results.surfaces, state=False)
dump("results.report (names only)", lambda: s.results.report, state=False)
dump("results.report.surface_integrals (names only)", lambda: s.results.report.surface_integrals, state=False)
dump("results.report.fluxes (names only)", lambda: s.results.report.fluxes, state=False)
dump("file (names only)", lambda: s.file, state=False)
dump("tui boundary_conditions.wall children",
     lambda: [n for n in dir(solver.tui.define.boundary_conditions) if not n.startswith("_")], state=True)

out.close()
print("\nDone. Upload introspect_dump.txt")
solver.exit()
