"""
fluent_common.py - shared helpers for all Fluent scripts in this project.
Verified working on Fluent 2026 R1 / PyFluent 0.43.0 (smoke test run 5):
  read mesh via settings.file.read_case, axisymmetric via setup.general.solver.two_dim_space,
  interface via tui.mesh.modify_zones.slit_interior_between_diff_materials().
"""
import os
import re
import sys

MESH_NAME_DEFAULT = "pipe2d_M2.msh"


def find_mesh(arg, script_file):
    here = os.path.dirname(os.path.abspath(script_file))
    cands = [arg, os.path.join(os.getcwd(), arg), os.path.join(here, arg),
             os.path.join(here, "..", "02_geometry_mesh", "2d", os.path.basename(arg))]
    for c in cands:
        if os.path.exists(c):
            return os.path.abspath(c)
    print("ERROR: mesh file not found. Looked at:")
    for c in cands:
        print("   ", os.path.abspath(c))
    sys.exit(1)


def preflight(mesh, say=print):
    """Reject damaged .msh files (e.g. pasted terminal text) before launching Fluent."""
    ok_hdr = re.compile(r"^\((0|2|10|12|13|45)\s")
    ok_data = re.compile(r"^[0-9a-fA-F\.\-\+eE\s]+$")
    bad = []
    with open(mesh, "r", errors="replace") as fh:
        for ln, line in enumerate(fh, 1):
            t = line.strip()
            if not t or ok_hdr.match(t) or ok_data.match(t) or t in (")", "))"):
                continue
            bad.append((ln, t[:80]))
    if bad:
        say(f"ERROR: {len(bad)} invalid line(s) in the mesh file; first: {bad[:5]}")
        sys.exit(1)
    say("Pre-flight OK, mesh size (bytes):", os.path.getsize(mesh))


def launch(nproc=4, gui=False):
    import ansys.fluent.core as pyfluent
    return pyfluent.launch_fluent(mode="solver", dimension=2, precision="double",
                                  processor_count=nproc, ui_mode="gui" if gui else "no_gui")


def load_2d_mesh(solver, mesh):
    """Read mesh, set axisymmetric, turn interface_fs into a coupled wall pair."""
    solver.settings.file.read_case(file_name=mesh)
    solver.settings.setup.general.solver.two_dim_space = "axisymmetric"
    solver.tui.mesh.modify_zones.slit_interior_between_diff_materials()
