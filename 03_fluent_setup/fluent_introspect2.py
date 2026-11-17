"""
fluent_introspect2.py - round 2: turn the energy equation ON, then dump the thermal settings,
solver/run controls, monitors and post-processing commands. Changes nothing else.
    python fluent_introspect2.py pipe2d_M2.msh 4
Output: introspect_dump2.txt  (UPLOAD this file).
"""
import sys
import pprint
from fluent_common import *

arg = sys.argv[1] if len(sys.argv) > 1 else MESH_NAME_DEFAULT
nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
mesh = find_mesh(arg, __file__)
preflight(mesh)
solver = launch(nproc, gui=False)
load_2d_mesh(solver, mesh)
s = solver.settings
s.setup.models.energy.enabled = True
out = open("introspect_dump2.txt", "w", encoding="utf-8")


def dump(title, getter, state=True, args=False):
    out.write(f"\n##### {title}\n")
    try:
        obj = getter()
        for attr in ("child_names", "command_names") + (("argument_names",) if args else ()):
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


dump("setup.models.energy", lambda: s.setup.models.energy)
dump("materials.fluid['air'] (energy on)", lambda: s.setup.materials.fluid["air"])
dump("materials.solid['aluminum'] (energy on)", lambda: s.setup.materials.solid["aluminum"])
dump("materials.fluid.air.density", lambda: s.setup.materials.fluid["air"].density)
dump("materials.fluid.air.viscosity", lambda: s.setup.materials.fluid["air"].viscosity)
dump("materials.database command names", lambda: s.setup.materials.database, state=False)
dump("cell_zone_conditions.fluid['fluid']", lambda: s.setup.cell_zone_conditions.fluid["fluid"])
dump("bc.velocity_inlet['inlet'].thermal", lambda: s.setup.boundary_conditions.velocity_inlet["inlet"].thermal)
dump("bc.velocity_inlet['inlet'].turbulence", lambda: s.setup.boundary_conditions.velocity_inlet["inlet"].turbulence)
dump("bc.pressure_outlet['outlet'].thermal", lambda: s.setup.boundary_conditions.pressure_outlet["outlet"].thermal)
dump("bc.wall['outer_heated'].thermal", lambda: s.setup.boundary_conditions.wall["outer_heated"].thermal)
dump("bc.wall['outer_unheated'].thermal", lambda: s.setup.boundary_conditions.wall["outer_unheated"].thermal)
dump("bc.wall['solid_end_in'].thermal", lambda: s.setup.boundary_conditions.wall["solid_end_in"].thermal)
dump("bc.wall['interface_fs'].thermal", lambda: s.setup.boundary_conditions.wall["interface_fs"].thermal)
dump("bc.wall['interface_fs-shadow'].thermal", lambda: s.setup.boundary_conditions.wall["interface_fs-shadow"].thermal)
dump("bc.wall['interface_fs'] full", lambda: s.setup.boundary_conditions.wall["interface_fs"])
dump("models.viscous (energy on)", lambda: s.setup.models.viscous)
dump("solution.methods (energy on)", lambda: s.solution.methods)
dump("solution.controls.equations", lambda: s.solution.controls.equations)
dump("solution.controls.under_relaxation", lambda: s.solution.controls.under_relaxation)
dump("solution.controls.pseudo_time_method_local_time_step", lambda: s.solution.controls.pseudo_time_method_local_time_step)
dump("solution.controls.pseudo_time_explicit_relaxation_factor", lambda: s.solution.controls.pseudo_time_explicit_relaxation_factor)
dump("solution.controls.p_v_controls", lambda: s.solution.controls.p_v_controls)
dump("solution.controls.courant_number", lambda: s.solution.controls.courant_number)
dump("solution.run_calculation.pseudo_time_settings", lambda: s.solution.run_calculation.pseudo_time_settings)
dump("solution.run_calculation.parameters", lambda: s.solution.run_calculation.parameters)
dump("solution.run_calculation.iterate args", lambda: s.solution.run_calculation.iterate, state=False, args=True)
dump("solution.run_calculation.calculate args", lambda: s.solution.run_calculation.calculate, state=False, args=True)
dump("solution.initialization", lambda: s.solution.initialization)
dump("solution.initialization.patch (names)", lambda: s.solution.initialization.patch, state=False, args=True)
dump("solution.initialization.hybrid_initialize args", lambda: s.solution.initialization.hybrid_initialize, state=False, args=True)
dump("solution.monitor.residual", lambda: s.solution.monitor.residual)
dump("solution.monitor.report_files (names)", lambda: s.solution.monitor.report_files, state=False)
dump("solution.monitor.convergence_conditions", lambda: s.solution.monitor.convergence_conditions)
dump("solution.report_definitions (names)", lambda: s.solution.report_definitions, state=False)
dump("results.report.surface_integrals.area_weighted_avg", lambda: s.results.report.surface_integrals.area_weighted_avg, state=False, args=True)
dump("results.report.surface_integrals.mass_weighted_avg", lambda: s.results.report.surface_integrals.mass_weighted_avg, state=False, args=True)
dump("results.report.surface_integrals.integral", lambda: s.results.report.surface_integrals.integral, state=False, args=True)
dump("results.report.fluxes.heat_transfer", lambda: s.results.report.fluxes.heat_transfer, state=False, args=True)
dump("results.report.fluxes.mass_flow", lambda: s.results.report.fluxes.mass_flow, state=False, args=True)
dump("file.write_case_data args", lambda: s.file.write_case_data, state=False, args=True)

# temporary surfaces, to learn how they are described
for kind in ("plane_surface", "iso_surface", "iso_clip", "zone_surface"):
    try:
        container = getattr(s.results.surfaces, kind)
        container.create(name="tmp_" + kind)
        dump(f"surfaces.{kind}['tmp_{kind}']", lambda: container["tmp_" + kind])
    except Exception as e:
        out.write(f"\n##### create {kind}\nFAILED: {type(e).__name__}: {str(e)[:300]}\n")
        out.flush()

out.close()
print("\nDone. Upload introspect_dump2.txt")
solver.exit()
