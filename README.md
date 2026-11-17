# Ternary nanofluid pipe CFD (ANSYS Fluent)

Scripted ANSYS Fluent 2026 R1 (Student) study of turbulent heat transfer in a heated copper pipe,
building towards Nusselt-number data for ternary nanofluids. Everything is driven from Python (PyFluent).

Pipe: 23 mm inner / 25 mm outer diameter copper, 4 m long. The first 2 m are unheated (flow develops);
the last 2 m are heated through 8 circumferential sectors (q1..q8). Model: realizable k-epsilon with
enhanced wall treatment, conjugate heat transfer through the wall, second-order discretisation.

## Folders

| Folder | Contents |
|---|---|
| `00_docs/teacher_case` | settings extracted from the supervisor's reference Fluent case |
| `02_geometry_mesh` | mesh generators: `2d/make_2d_mesh.py` (axisymmetric), `3d/make_3d_mesh.py` (half pipe, 449k cells) |
| `03_fluent_setup` | Fluent drivers: `fluent_gateA.py` (2D), `fluent_gateB.py` (3D, incl. `--suns` per-sector flux), viewers |
| `05_results` | result JSON, Nu profiles (CSV), run logs and Fluent transcripts |
| `nf_properties.py` | nanofluid property calculator (+ `test_nf_properties.py`) |

Meshes (`*.msh`) and Fluent case/data files (`*.cas.h5`, `*.dat.h5`) are not in the repository:
regenerate them with the scripts below.

## Run (PowerShell, from this folder)

```
python -m venv .venv
.venv\Scripts\pip install -r requirements.txt
.venv\Scripts\python 02_geometry_mesh\2d\make_2d_mesh.py
.venv\Scripts\python 02_geometry_mesh\3d\make_3d_mesh.py
.venv\Scripts\python 03_fluent_setup\fluent_gateA.py --re 10000 --mesh M2      # 2D, minutes
.venv\Scripts\python 03_fluent_setup\fluent_gateB.py --re 10000                 # 3D, ~50 min
.venv\Scripts\python 03_fluent_setup\view_gateB_in_fluent.py                    # open the 3D result in Fluent
```

## Validation status (water, Re 10,000)

| Gate | Test | Result |
|---|---|---|
| A | 2D vs Gnielinski (within 10 %) | Nu 84.07, +15.6 % at Re 10k (model bias, same in supervisor's case); +6.9 % at Re 25k |
| B | 3D half pipe vs 2D (within 2 %) | Nu 84.29, +0.26 %: pass |
| C | Energy balance below 0.2 % | 0.0001 %: pass |
| D-H | 3D mesh study, nanofluid check, repeatability, flux, pipe material | to do |
