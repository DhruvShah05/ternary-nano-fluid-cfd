"""
fluent_introspect_ke.py - READ-ONLY discovery of the k-epsilon / wall-treatment / turbulent-Prandtl settings.
It does NOT solve anything. It reads the M2 mesh, then prints which settings exist and which values are allowed,
first with the current k-omega SST model (turbulent Prandtl options), then after switching to k-epsilon
(standard and realizable). Output: 03_fluent_setup\\introspect_dump_ke.txt  (upload that file).

Run (PowerShell, project root, venv active; close any Fluent GUI first):
    python 03_fluent_setup\\fluent_introspect_ke.py
"""
import pprint
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fluent_common as fc  # noqa: E402

OUT = HERE / "introspect_dump_ke.txt"
LINES = []


def P(*a):
    t = " ".join(str(x) for x in a)
    print(t, flush=True)
    LINES.append(t)
    OUT.write_text("\n".join(LINES), encoding="utf-8")


def safe(label, fn):
    try:
        return fn()
    except Exception as e:  # noqa: BLE001
        P(f"  [{label}] FAILED: {type(e).__name__}: {str(e)[:220]}")
        return None


def allowed(obj):
    f = getattr(obj, "allowed_values", None)
    if f is None:
        return None
    try:
        return f() if callable(f) else f
    except Exception:  # noqa: BLE001
        return None


def state_str(obj, limit=1800):
    try:
        st = obj.get_state()
    except Exception as e:  # noqa: BLE001
        return f"<get_state failed: {str(e)[:100]}>"
    t = pprint.pformat(st, width=140)
    return t if len(t) <= limit else t[:limit] + " ...[truncated]"


def dump(obj, name):
    """One level: child names, command names, state, and allowed values of every child that has them."""
    P(f"\n##### {name}")
    av = allowed(obj)
    if av is not None:
        P("allowed_values:", av)
    kids = safe("child_names", lambda: list(obj.child_names)) or []
    P("child_names:", kids)
    cmds = safe("command_names", lambda: list(obj.command_names))
    if cmds:
        P("command_names:", cmds)
    P("state:", state_str(obj))
    for k in kids:
        child = safe(f"getattr {k}", lambda k=k: getattr(obj, k))
        if child is None:
            continue
        a = allowed(child)
        if a is not None:
            P(f"  .{k} allowed_values:", a)


def main():
    mesh = fc.find_mesh("pipe2d_M2.msh", __file__)
    fc.preflight(mesh, P)
    solver = fc.launch(nproc=4)
    try:
        fc.load_2d_mesh(solver, mesh)
        s = solver.settings
        s.setup.models.energy.enabled = True
        P("Mesh read, axisymmetric, slit, energy on. Currently k-omega SST.")
        v = s.setup.models.viscous

        P("\n================ PART A: current model (k-omega SST) ================")
        dump(v.model, "viscous.model")
        safe("A1", lambda: dump(v.options, "viscous.options"))
        safe("A2", lambda: dump(v.options.turbulent_prandtl_numbers, "viscous.options.turbulent_prandtl_numbers"))
        safe("A3", lambda: dump(v.options.turbulent_prandtl_numbers.energy_prandtl_number,
                                "...energy_prandtl_number"))
        safe("A4", lambda: dump(v.options.turbulent_prandtl_numbers.wall_prandtl_number,
                                "...wall_prandtl_number"))
        safe("A5", lambda: dump(v.turbulence_expert, "viscous.turbulence_expert"))
        safe("A6", lambda: dump(v.near_wall_treatment, "viscous.near_wall_treatment (k-omega)"))

        P("\n================ PART B: switch to k-epsilon ================")
        safe("set model k-epsilon", lambda: setattr(v, "model", "k-epsilon"))
        P("model now:", safe("model", lambda: v.model.get_state() if hasattr(v.model, 'get_state') else v.model))
        P("viscous state:", state_str(v, 3000))
        dump(v.k_epsilon_model, "viscous.k_epsilon_model")
        for variant in ("standard", "realizable"):
            P(f"\n---------------- k_epsilon_model = {variant} ----------------")
            safe(f"set {variant}", lambda variant=variant: setattr(v, "k_epsilon_model", variant))
            P("k_epsilon_model now:", safe("get", lambda: v.k_epsilon_model.get_state()))
            safe("B1", lambda: dump(v.near_wall_treatment, f"viscous.near_wall_treatment ({variant})"))
            safe("B2", lambda: dump(v.k_epsilon, f"viscous.k_epsilon ({variant})"))
            safe("B3", lambda: dump(v.options, f"viscous.options ({variant})"))
            safe("B4", lambda: dump(v.options.turbulent_prandtl_numbers, f"...turbulent_prandtl_numbers ({variant})"))

        P("\n================ PART C: other things the gate script touches, now with k-epsilon ================")
        safe("C1", lambda: P("residual equations:", list(s.solution.monitor.residual.equations.keys())))
        safe("C2", lambda: P("controls.equations state:", s.solution.controls.equations.get_state()))
        safe("C3", lambda: dump(s.solution.controls.pseudo_time_explicit_relaxation_factor,
                                "pseudo_time_explicit_relaxation_factor"))
        bc = s.setup.boundary_conditions
        safe("C4", lambda: dump(bc.velocity_inlet["inlet"].turbulence, "velocity_inlet inlet turbulence"))
        safe("C5", lambda: dump(bc.pressure_outlet["outlet"].turbulence, "pressure_outlet outlet turbulence"))
        safe("C6", lambda: dump(bc.wall["interface_fs"].momentum, "wall interface_fs momentum"))
        P("\nDONE. Upload 03_fluent_setup\\introspect_dump_ke.txt")
    finally:
        try:
            solver.exit()
        except Exception:  # noqa: BLE001
            pass


if __name__ == "__main__":
    main()
