"""
fluent_smoke_test.py  -  Phase 2, step 1 (v2): launch Fluent 2D, read our mesh, report zones.

Run from anywhere:
    python fluent_smoke_test.py pipe2d_M2.msh 4
Arguments: mesh file name or path (default pipe2d_M2.msh), processes (default 4).
If the mesh is not at the given path, the script also looks in the current folder and
..\\02_geometry_mesh\\2d relative to this script, then stops with a clear message.
Output is also written to smoke_test_log.txt (in the folder you run from).
"""
import sys
import os
import traceback

arg = sys.argv[1] if len(sys.argv) > 1 else "pipe2d_M2.msh"
nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 4
here = os.path.dirname(os.path.abspath(__file__))
candidates = [arg,
              os.path.join(os.getcwd(), arg),
              os.path.join(here, arg),
              os.path.join(here, "..", "02_geometry_mesh", "2d", os.path.basename(arg))]
mesh = next((os.path.abspath(c) for c in candidates if os.path.exists(c)), None)

log = open("smoke_test_log.txt", "w")


def say(*a):
    msg = " ".join(str(x) for x in a)
    print(msg)
    log.write(msg + "\n")
    log.flush()


def step(name, fn):
    say(f"\n=== {name} ===")
    try:
        r = fn()
        if r is not None:
            say(r)
        say(f"[OK] {name}")
        return True
    except Exception:
        say(f"[FAILED] {name}")
        say(traceback.format_exc(limit=2))
        return False


if mesh is None:
    say("ERROR: mesh file not found. Looked at:")
    for c in candidates:
        say("   ", os.path.abspath(c))
    say("Put pipe2d_M2.msh next to this script or pass its full path. Nothing launched.")
    sys.exit(1)


# ---- Pre-flight: make sure the .msh file contains only valid Fluent mesh lines -------------
import re
_ok_hdr = re.compile(r"^\((0|2|10|12|13|45)\s")
_ok_data = re.compile(r"^[0-9a-fA-F\.\-\+eE\s]+$")
bad = []
with open(mesh, "r", errors="replace") as fh:
    for ln, line in enumerate(fh, 1):
        t = line.strip()
        if not t or _ok_hdr.match(t) or _ok_data.match(t) or t in (")", "))"):
            continue
        bad.append((ln, t[:80]))
say("Mesh file size (bytes):", os.path.getsize(mesh))
if bad:
    say(f"ERROR: {len(bad)} line(s) in the mesh file are not valid mesh data. First ones:")
    for ln, t in bad[:8]:
        say(f"   line {ln}: {t}")
    say("The file is damaged (probably text pasted into it). Re-download pipe2d_M2.msh and replace it.")
    say("Expected size of the correct pipe2d_M2.msh: 4130717 bytes.")
    sys.exit(1)
say("Pre-flight OK: mesh file contains only valid mesh lines.")

import ansys.fluent.core as pyfluent

say("PyFluent version:", getattr(pyfluent, "__version__", "unknown"))
say("Mesh file:", mesh)

solver = pyfluent.launch_fluent(mode="solver", dimension=2, precision="double",
                                processor_count=nproc, ui_mode="gui")
say("Launched:", solver.get_fluent_version())

# Paths with spaces break the text interface unless quoted, so try the settings API first.
def read_mesh():
    errors = []
    for label, fn in [
        ("settings.file.read_case", lambda: solver.settings.file.read_case(file_name=mesh)),
        ("settings.file.read_mesh", lambda: solver.settings.file.read_mesh(file_name=mesh)),
        ("tui read_case (quoted)", lambda: solver.tui.file.read_case('"' + mesh.replace("\\", "/") + '"')),
    ]:
        try:
            fn()
            return f"read OK via {label}"
        except Exception as e:
            errors.append(f"{label}: {str(e)[:200]}")
    raise RuntimeError("all read methods failed:\n" + "\n".join(errors))


step("Read mesh", read_mesh)
step("Mesh check", lambda: solver.tui.mesh.check())

def verbose_check():
    out = []
    for label, fn in [
        ("tui.mesh.check_verbosity(2)", lambda: solver.tui.mesh.check_verbosity(2)),
        ("settings.mesh.check_verbosity=2", lambda: setattr(solver.settings.mesh, "check_verbosity", 2)),
    ]:
        try:
            fn()
            out.append("set verbosity via " + label)
            break
        except Exception as e:
            out.append(f"{label} failed: {str(e)[:120]}")
    return "\n".join(out)

step("Raise mesh-check verbosity", verbose_check)
step("Mesh check (verbose)", lambda: solver.tui.mesh.check())
step("Report poor elements", lambda: solver.tui.mesh.repair_improve.report_poor_elements())
step("Children of settings.mesh", lambda: sorted(n for n in dir(solver.settings.mesh) if not n.startswith("_")))
step("Mesh size info", lambda: solver.tui.mesh.size_info())
step("List boundary zones", lambda: solver.tui.define.boundary_conditions.list_zones())
step("Set axisymmetric",
     lambda: setattr(solver.settings.setup.general.solver, "two_dim_space", "axisymmetric"))
step("Read back 2D space setting", lambda: solver.settings.setup.general.solver.two_dim_space())


step("Try slit_interior_between_diff_materials (makes fluid/solid interface a coupled wall?)",
     lambda: solver.tui.mesh.modify_zones.slit_interior_between_diff_materials())
step("List boundary zones AFTER slit", lambda: solver.tui.define.boundary_conditions.list_zones())
step("Mesh check AFTER slit", lambda: solver.tui.mesh.check())

def list_children(obj):
    return sorted(n for n in dir(obj) if not n.startswith("_"))

step("Commands under tui.define.boundary_conditions.modify_zones",
     lambda: list_children(solver.tui.define.boundary_conditions.modify_zones))
step("Commands under tui.mesh.modify_zones",
     lambda: list_children(solver.tui.mesh.modify_zones))
step("Mesh Interfaces in settings", lambda: list_children(solver.settings.setup.mesh_interfaces))

say("\nLook in the Fluent window: Outline View > Setup > Boundary Conditions.")
say("Press Enter here to close Fluent and finish.")
try:
    input()
except EOFError:
    pass
solver.exit()
