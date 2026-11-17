"""
view_gateB_in_fluent.py - opens the finished Gate B 3D result in the Fluent window, with ready-made pictures.

    .venv\\Scripts\\python.exe 03_fluent_setup\\view_gateB_in_fluent.py

Fluent opens (GUI, 4 cores), reads 04_cases\\gateB\\gateB_Re10000_B3D.cas.h5 + .dat.h5 and creates these
contour objects (Outline tree > Results > Graphics > Contours; double-click one to show it):
    temp_section_x3p5  temperature on a cross-section at x = 3.5 m (water + copper wall)
    vel_section_x3p5   velocity magnitude on the same cross-section
    temp_section_x1    temperature at x = 1.0 m (unheated part: uniform 300 K)
    temp_outer_wall    temperature on the outer copper surface (heated sectors + unheated part)
    temp_symmetry      temperature on the symmetry plane z = 0 (long side view of the pipe)
    vel_symmetry       velocity magnitude on the symmetry plane

The case/data files are only READ, never overwritten. Fluent stays open when this script ends;
close it yourself (File > Exit) when you are done.
"""
from pathlib import Path

import ansys.fluent.core as pyfluent

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "04_cases" / "gateB" / "gateB_Re10000_B3D.cas.h5"


def step(name, fn):
    try:
        fn()
        print(f"  ok    {name}")
    except Exception as e:  # a missing picture must not stop the others
        print(f"  SKIP  {name}: {type(e).__name__}: {str(e)[:160]}")


def main():
    if not CASE.exists():
        raise SystemExit(f"case not found: {CASE}")
    print("starting Fluent (GUI) ...")
    solver = pyfluent.launch_fluent(mode="solver", dimension=3, precision="double", processor_count=4,
                                    ui_mode="gui", cleanup_on_exit=False)
    s = solver.settings
    print("reading case + data ...")
    s.file.read_case_data(file_name=str(CASE))

    surfs = s.results.surfaces
    gfx = s.results.graphics

    def plane(name, x):
        def fn():
            surfs.plane_surface.create(name)
            p = surfs.plane_surface[name]
            p.method = "yz-plane"
            p.x = x
        return fn

    step("plane x = 3.5 m", plane("x_3p5", 3.5))
    step("plane x = 1.0 m", plane("x_1p0", 1.0))

    def contour(name, field, surfaces, rng=None):
        def fn():
            gfx.contour.create(name)
            c = gfx.contour[name]
            c.field = field
            c.surfaces_list = surfaces
            c.node_values = True
            if rng:
                c.range_options.auto_range = False
                c.range_options.minimum, c.range_options.maximum = rng
        return fn

    outer = ["outer_unheated", "outer_q1", "outer_q2", "outer_q3", "outer_q4"]
    step("temp_section_x1", contour("temp_section_x1", "temperature", ["x_1p0"]))
    step("temp_outer_wall", contour("temp_outer_wall", "temperature", outer))
    step("temp_symmetry", contour("temp_symmetry", "temperature", ["symmetry_fluid", "symmetry_solid"]))
    step("vel_symmetry", contour("vel_symmetry", "velocity-magnitude", ["symmetry_fluid"]))
    step("vel_section_x3p5", contour("vel_section_x3p5", "velocity-magnitude", ["x_3p5"]))
    step("temp_section_x3p5", contour("temp_section_x3p5", "temperature", ["x_3p5"]))

    # show the most telling picture first: the heated cross-section, looking down the pipe axis
    step("display temp_section_x3p5", lambda: gfx.contour["temp_section_x3p5"].display())
    step("view along the pipe axis", lambda: gfx.views.restore_view(view_name="right"))
    step("fit to window", lambda: gfx.views.auto_scale())

    print("\nFluent is open. In the Outline tree go to Results > Graphics > Contours and double-click a contour.")
    print("Mouse in the graphics window: left-drag = rotate, middle-drag = pan, wheel / right-drag = zoom.")
    print("The case files were only read. Close Fluent with File > Exit when done (no need to save).")


if __name__ == "__main__":
    main()
