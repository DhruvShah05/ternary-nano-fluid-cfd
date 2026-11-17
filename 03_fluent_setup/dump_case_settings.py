"""
dump_case_settings.py - READ-ONLY audit of a saved Gate A case. It does NOT iterate or save anything.
Reads 04_cases/gateA/gateA_<tag>.cas.h5 (+ .dat.h5), prints every setting that affects the answer
(general, models, materials, cell zones, boundary conditions, methods, controls) and runs Fluent's mesh check.

Run (PowerShell, project root, venv active; close any Fluent GUI first):
    python 03_fluent_setup\\dump_case_settings.py --tag Re10000_M2_ke_rlz_ewt
Output: 05_results\\case_settings_<tag>.txt
"""
import argparse
import pprint
import shutil
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="Re10000_M2_ke_rlz_ewt")
    a = ap.parse_args()
    cas = ROOT / "04_cases" / "gateA" / f"gateA_{a.tag}.cas.h5"
    out = ROOT / "05_results" / f"case_settings_{a.tag}.txt"
    lines = []

    def put(title, obj):
        lines.append(f"\n==================== {title}")
        lines.append(obj if isinstance(obj, str) else pprint.pformat(obj, width=120, sort_dicts=False))
        out.write_text("\n".join(lines), encoding="utf-8")

    import ansys.fluent.core as pyfluent
    solver = pyfluent.launch_fluent(mode="solver", dimension=2, precision="double", processor_count=4)
    tmp_t = Path(tempfile.gettempdir()) / f"fluent_dump_{a.tag}.txt"
    try:
        solver.tui.file.start_transcript(f'"{tmp_t}"' if " " in str(tmp_t) else str(tmp_t))
        s = solver.settings
        s.file.read_case_data(file_name=str(cas))
        put("case", str(cas))
        sections = {
            "general": lambda: s.setup.general.get_state(),
            "models.viscous": lambda: s.setup.models.viscous.get_state(),
            "models.energy": lambda: s.setup.models.energy.get_state(),
            "materials": lambda: s.setup.materials.get_state(),
            "cell_zone_conditions": lambda: s.setup.cell_zone_conditions.get_state(),
            "boundary_conditions": lambda: s.setup.boundary_conditions.get_state(),
            "solution.methods": lambda: s.solution.methods.get_state(),
            "solution.controls": lambda: s.solution.controls.get_state(),
            "pseudo_time_settings": lambda: s.solution.run_calculation.pseudo_time_settings.get_state(),
        }
        for name, fn in sections.items():
            try:
                put(name, fn())
            except Exception as e:  # noqa: BLE001
                put(name, f"(could not read: {type(e).__name__}: {str(e)[:300]})")
        for name, fn in (("mesh check", lambda: solver.tui.mesh.check()),
                         ("mesh quality", lambda: solver.tui.mesh.quality())):
            try:
                fn()
                put(name, "(printed to the transcript, appended below)")
            except Exception as e:  # noqa: BLE001
                put(name, f"(failed: {str(e)[:300]})")
        try:
            solver.tui.file.stop_transcript()
        except Exception:  # noqa: BLE001
            pass
        if tmp_t.exists():
            put("Fluent transcript", tmp_t.read_text(encoding="utf-8", errors="replace"))
    finally:
        solver.exit()
    print(f"written {out}")


if __name__ == "__main__":
    main()
