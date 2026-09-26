# FSW bridge: Rust flight software in the loop with the 6DOF sim

The plant ("environment": actuators, propulsion mass, CG/MOI, slosh, wind,
6DOF dynamics) is the `Environment` subsystem of
`../hopper_6dof_NED_v2_fswBridge.slx`. It is turned into C with Embedded
Coder, and a Rust program steps it in lockstep with luna's flight-software
controller.

```
  Rust (hopper_sil)                        C (generated from Simulink)
  luna LqrController  --u_cmd[4]-->  env_update()   advance 1 ms
         ^                           env_output()   state at time t
         +--------- x_true[13] -----------+
```

- `u_cmd = [thrust, tvc_pitch, tvc_yaw, rcs]`, applied on the next step
  (`Environment/u_delay`), as on the vehicle.
- `x_true = [pos NED (3), vel (3), body rates (3), quaternion q0..q3 (4)]`,
  the true state. Sensor models are not in yet.

## Workflow

In MATLAB, from `sims/hopper_6dof_sizer` (fresh session each time; the build
refuses if `hopper_env` is already loaded):

```matlab
addpath fsw_bridge
build_env_c     % Environment -> C in build/hopper_env_ert_rtw (~2 min)
export_gains    % LQR gain tables -> build/gains/*.csv
verify_env_c    % optional: C vs Simulink replay check (gusts off)
```

Then the Rust loop (needs Rust with the `x86_64-pc-windows-gnu` toolchain and
the luna repo cloned next to hopper):

```bash
cd fsw_bridge/sil
HOPPER_CC=<MinGW>/bin/gcc.exe cargo run --release -- --matlab-k2-sign
```

`HOPPER_CC` is the MinGW gcc MATLAB uses (`mex -setup C` shows it). Other
build settings: `HOPPER_ENV_CODE`, `MATLAB_ROOT`, `LUNA_DIR` (see `build.rs`).
The run writes `build/sil_run.csv`: `t, x_true[13], u_cmd[4], thrust, z,
ox_mass, fuel_mass` per step.

To compare a Rust run with Simulink step for step, build against the
gust-free code (`HOPPER_ENV_CODE=../build/verify_nogust/hopper_env_ert_rtw`)
and run `compare_sil('build/sil_run.csv')` in MATLAB.

## Things to know

- **K2 sign.** luna `flight2/src/control/lqr.rs` computes
  `u = unom - K1*dx + K2`; the Simulink controller uses `- K2`.
  `--matlab-k2-sign` negates the table so the Rust controller matches
  Simulink. Without it the flight separates from Simulink at t = 0.18 s and
  lands ~0.6 s later.
- **Gusts.** The generated code draws different random numbers than
  Simulink's gust blocks, so runs with gusts on do not match Simulink step
  for step. The gust on/off generator's seed is `Auto`, not controllable yet.
- **One environment per process.** The generated code keeps its state in
  globals.
- **Parameters are inlined.** Changing a workspace value (mass, wind, ...)
  needs `build_env_c` again.
- **Landing is very sensitive.** With identical logic, 1e-12 differences
  grow to metres in the last ~1.5 s before touchdown, in Simulink too.
