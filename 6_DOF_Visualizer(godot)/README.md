# Hopper 6-DOF Telemetry Visualizer

A Godot 4 visualizer for a rocket simulation or live telemetry stream. It accepts MATLAB/Simulink NED data and provides 3D flight, playback, plots, cameras, and export controls.

## Project layout

The active scene is [Scenes/main.tscn](Scenes/main.tscn), which attaches [Scripts/main.gd](Scripts/main.gd). The main script is a composition root; it only creates modules and connects their signals.

| Script | Responsibility |
| --- | --- |
| `main.gd` | Connects modules and handles source/load/export commands. |
| `telemetry_database.gd` | Loads, normalizes, archives, exports, and retains telemetry rows. |
| `telemetry_schema.gd` | The authoritative CSV column table. |
| `attitude_transform.gd` | Converts MATLAB quaternions to Godot mesh orientation. |
| `playback_controller.gd` | Playback clock, speed, stepping, and scrubbing. |
| `flight_world.gd` | Ground, lighting, rocket, traces, landing aids, and cameras. |
| `rocket_visual.gd` | CG/CP markers, propellant levels, and vector arrows. |
| `visualizer_dashboard.gd` | Responsive draggable UI, graphs, HUD, and file picker. |
| `telemetry_graph.gd` | Reusable graph card with hover inspection and pop-out support. |
| `demo_telemetry.gd` | Starter flight shown before a CSV is loaded. |

No active script uses `class_name`; every dependency is explicitly preloaded by its consumer. Modules therefore stay private to this visualizer instead of becoming project-wide globals.

## MATLAB / Simulink telemetry contract

The canonical input is the 18-column output produced by `simToCSV.m`:

| Index | Field | Meaning |
| ---: | --- | --- |
| 0 | `time` | Simulation time, seconds |
| 1–3 | `north`, `east`, `up` | Model position, m; the exported `z` column is up-positive altitude |
| 4 | `thrust` | Force magnitude, N |
| 5–8 | `q0`, `q1`, `q2`, `q3` | MATLAB quaternion, scalar first |
| 9–10 | `fuel_mass`, `ox_mass` | Propellant mass, kg |
| 11 | `control` | Normalized TVC/RCS control proxy |
| 12–14 | `cg_x`, `cg_y`, `cg_z` | Body-frame CG, m |
| 15–17 | `vn`, `ve`, `vd` | NED velocity, m/s |

Header rows are optional. The reader also accepts the old 17-column Euler format and normalizes it internally, but new data should use the quaternion format.

### Quaternion and coordinate conversion

`attitude_transform.gd` assumes the common MATLAB/Aerospace convention:

- `[A, B, C, D]` means scalar-first `[q0, q1, q2, q3]`.
- The quaternion rotates body coordinates into NED coordinates.
- MATLAB body axes are forward/right/down (`+X/+Y/+Z`), while the supplied mesh has its nose along local `+Y`.
- The visualizer uses native Godot axes: +X east/right, +Y up, and +Z south/back. Exported position maps as `(north, east, up) → (east, up, -north)`. NED velocity maps as `(vn, ve, vd) → (ve, -vd, -vn)`.

Godot receives `NED_TO_GODOT × quaternion_basis × BODY_FROM_MODEL`. Position and velocity use the same NED conversion, so attitude, traces, and HUD values agree. The ground is Godot's normal XZ plane and +Y is vertical. Zero or negative altitude rests the engine on the ground through the mesh-base offset in `flight_world.gd`.

If the source block instead outputs NED-to-body quaternions, invert the quaternion basis in `attitude_transform.gd`. Confirm the source convention with a known horizontal/north-facing attitude and a vertical launch state before interpreting flight data.

## UI

- Load a CSV with **Load CSV**, or start newline-delimited TCP CSV intake using **TCP :8080**.
- Playback supports pause, stepping, scrubbing, and 0.25×–2× speed.
- Hold right mouse and drag to orbit. Fly-by and engine landing cameras remain available.
- Graph cards can be collapsed, popped out into a resizable OS window, and inspected with a time/value hover readout.
- Drag a card by its dotted title bar. **Reset layout** restores the responsive default location.
- The UI uses proportional card layouts and wrapping control rows to fit different window sizes without its default cards overlapping.

## Verification

The active scripts were parsed by Godot 4.6 and the scene was started headlessly after this refactor. The supplied `Test_Data/godotData.csv` was also read through the active parser: it contains 27,462 valid unit quaternions, and its initial attitude maps the mesh nose to +Y (up) in native Godot coordinates. There were no project script parse or runtime errors.
