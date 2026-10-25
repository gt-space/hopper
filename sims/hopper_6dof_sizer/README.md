# Hopper 6DOF Vehicle Sizer

MATLAB/Simulink sims used for preliminary vehicle sizing of the Hopper (YJSP's Propulsive Lander), modeling wind, slosh, resonant behavior, & mass/inertia variation.

## Quick start (MATLAB R2026a)

1. Open the project: double-click `HopperSim.prj`, or run `openProject('sims/hopper_6dof_sizer')`.
   This puts every source folder on the path and sends Simulink caches to `work/`. No `addpath` or `cd` is needed.
2. Run `sim_setup`, then simulate `hopper_6dof_NED_v2`.
3. Monte Carlo: `mc_runner(10)` (serial) or `mc_parallelization` (parallel, batched). Every run calls CoolProp through Python (`pip install CoolProp`).

Save models in **R2026a** only. A model saved in a newer release can't be opened by everyone else.

## Layout

| Folder | Contents |
|---|---|
| `models/` | Top-level models (`hopper_6dof_NED_v2.slx` is the master). `dev/` holds work in progress (LQI); `standalone/` holds component models. |
| `subsystems/` | Subsystem files referenced by the models (FDIR, actuators, propulsion sensors/valves, navigation filters). |
| `src/` | MATLAB code: `setup`, `sizing`, `dynamics`, `gnc`, `propulsion`, `slosh`, `fdir`, `monte_carlo`, `analysis`, `util`. |
| `inputs/` | `mission_inputs.m` and `mc_params.json`. |
| `data/` | Lookup tables the sims read (`mdot_lookup.xlsx`, `wind_vectors2.mat`, `cg_I_LUT.mat`, `baffle_lut.*`, `engine_contour.xlsx`). |
| `results/` | Tracked Monte Carlo and sweep results (`.mat`). |
| `tests/`, `validation/` | Automated tests and validation studies. |
| `fsw_bridge/` | C codegen and Rust SIL bridge to flight software. |
| `tools/refactor/` | Scripts that build or convert models reproducibly (for example `build_v2_fluid_integration`). |
| `docs/` | Model notes (`SLOSH_MODEL.md`). |
| `archive/` | Superseded models and scripts, kept for reference and **not** on the path. |
| `outputs/`, `work/` | Generated figures, videos, CSVs and build caches. Gitignored. |

If you move or add files with `git mv` or Explorer instead of the Project panel, run `syncProjectFiles` so the project's file list matches the folders.

Scripts that write files use `hopperRoot('data', ...)`, `hopperRoot('results', ...)` or `hopperRoot('outputs', ...)`, so output lands in the right folder whatever the current folder is.

## Diffing and merging models

`.slx` files are binary to git. `.gitattributes` routes them to MATLAB's merge driver (`merge=mlAutoMerge`). Set it up once per machine, using your own MATLAB install path:

```bash
git config --global merge.mlAutoMerge.driver '"C:/Program Files/MATLAB/R2026a/bin/win64/mlAutoMerge.bat" %O %A %B %A'
git config --global difftool.mlDiff.cmd '"C:/Program Files/MATLAB/R2026a/bin/win64/mlDiff.bat" $LOCAL $REMOTE'
git config --global mergetool.mlMerge.cmd '"C:/Program Files/MATLAB/R2026a/bin/win64/mlMerge.bat" $BASE $LOCAL $REMOTE $MERGED'
```

- **Diffs:** `git difftool -t mlDiff -- <model>.slx` opens Simulink's Comparison Tool.
- **Conflicts:** `mlAutoMerge` resolves what it can. For the rest, `git mergetool -t mlMerge` opens the three-way model merge.

Avoid having two people edit the same `.slx`. Shared pieces are moving into separate subsystem files so changes rarely overlap.
