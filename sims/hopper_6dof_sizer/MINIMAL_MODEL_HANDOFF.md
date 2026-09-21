# 6DOF environments: original, toggleable, minimal

Three versions of the hopper 6DOF sim, for trading fidelity against speed.

| file | what it is | active blocks | time per step |
|---|---|---|---|
| `hopper_6dof_NED_v2.slx` | **original / master**, full physics. Unchanged. | 1011 | baseline |
| `hopper_6dof_NED_v2_simplified.slx` | **toggleable**: Slosh, Wind Forces and Gusts are Variant Subsystems switched by `enable_slosh` / `enable_wind` | 479 with both off | ~85% with both off |
| `hopper_6dof_NED_v2_minimal.slx` | **minimal**: slosh, wind, gusts, the fluid network and the ground block deleted; mass, inertia and CG frozen at wet-mass values | 265 | ~78% |

Time per step was measured in two separate sessions (78% / 79% for minimal); run-to-run machine noise is about ±10%.

## Using them

```matlab
sim_setup_cached                    % or sim_setup; builds the workspace
enable_slosh.Value = 0;             % toggleable only; both default to 1 (full physics)
sim('hopper_6dof_NED_v2_simplified')
```

Turn a feature off for **one run only**, without touching the workspace:

```matlab
in = Simulink.SimulationInput('hopper_6dof_NED_v2_simplified');
out = sim(in.setVariable('enable_slosh', 0));
```

`compare_sims` runs all four configurations (original, toggleable on, toggleable off, minimal), reports time per simulated step and the max difference per logged signal, and plots overlays.

### Safety nets built into the models
- **Toggleable model `InitFcn`:** if `enable_slosh` / `enable_wind` are not defined, they default to **1**. The model compiles under any setup script and never silently drops physics. This is what makes the three variant subsystems safe to paste into the master: no shared setup script needs to change.
- **Warnings:** the `Disabled` choice of Slosh and of Wind Forces each print a warning when active, so a reduced-physics run is never silent.
- **Minimal model `InitFcn`:** computes its frozen values (`I_fixed = MoI_init`, `dIdt_fixed = zeros(3)`, `cg_fixed`, `engine_cg_fixed = [IN.mount_cg IN.off_axis_x IN.off_axis_y]`, `ox_mass_fixed`, `fu_mass_fixed`) from the current workspace on every run, so it works under `sim_setup`, `sim_setup_cached` or `mc_sim_setup` and follows Monte Carlo mass and inertia dispersions.

## Monte Carlo

Both MC scripts accept any of the three models via the shared helper `mc_sim_input.m`:

```matlab
mc_runner(20, [], [], 'minimal')                  % [] keeps a default
mc_runner(20, [], [], 'toggleable', [0 1])        % flags = [enable_slosh enable_wind]
```

In `mc_parallelization.m`, set `mc_model` and `mc_flags` at the top. The default (`'original'`) behaves exactly as before. Results for other models go to `mc_results_parallel_<model>.csv` / `mc_results_<model>.mat`, so runs don't overwrite each other.

**Measured per-worker cost** (one MATLAB process running one full master simulation): **2.43 GB peak**, **~63 s** per scenario warm, **~113 s** on a worker's first (compiling) run. 8 workers need ~20 GB.

**Requirements** not met on every machine: Parallel Computing Toolbox installed (not just licensed), and the CoolProp Python package in the Python MATLAB uses (`pyenv`), which `prop_system` calls through `mc_sim_setup`.

## Limitations
- **Minimal model:** CG is a fixed number (`cg_fixed`, the CG/MOI script's t = 0 output), so it does **not** respond to `cg_factor` dispersions, and it ignores slosh and wind dispersions entirely. Use the toggleable model with features off if CG dispersion matters. Because mass never depletes, it needs more thrust late in flight than the original; z tracks the original within 0.34 m of 52 m.
- The original's gust generators are unseeded, so two runs of the *original* differ (x ~0.38, thrust ~1300). Treat differences of that size as noise; roll / yaw differences near 180 or 360 are angle wrapping.

## How the models were built and verified
- **Toggleable:** Slosh, Wind Forces and the gust generators (grouped into `Gusts`) each wrapped in a Variant Subsystem with `Enabled` (original contents) and `Disabled` (zeros) choices. With both flags at 1 it matches the original to within the original's own run-to-run noise.
- **Minimal:** built from the toggleable model with both flags off, then the disabled features, `Fluid Network2` and the commented-out `Ground ` block deleted, force / moment rerouted from `Control Inputs` straight into the 6DOF block (verified bit-identical to toggleable-off). Then the CG/MOI script and the propellant-depletion chain were replaced by the frozen constants; every signal matches the pre-change model exactly at t = 0.
- `cg_fixed` (1.2566) is the CG/MOI script's own t = 0 output. It differs from the sizer's `cg_init` (1.1399) by 0.117 m; the model's value was kept so behaviour at t = 0 is unchanged.

## Notes for editing these models
- A `delete_block` + "remove dangling lines" pass must run **before** adding new lines: a branched line's segments report no source, so the cleanup can delete live branches (this once silently severed `Mass -> Divide` and produced NaN at t = 0).
- Constant blocks default to `SampleTime = inf`, which makes a downstream To Workspace log only 1–3 points; the constants feeding logs use `-1`.
- A To Workspace block keeps only the **last 1000 points** by default.
- Block names with special characters: `CG // MOI Script v3` (double the slash in paths), `Ground ` (trailing space).
