# FSW bridge: Rust flight software in the loop with the 6DOF sim

The plant ("environment": actuators, propulsion mass, CG/MOI, slosh, wind,
6DOF dynamics, sensor models) is the `Environment` subsystem of
`../hopper_6dof_NED_v2_fswBridge.slx`. It is turned into C with Embedded
Coder, and a Rust program steps it in lockstep with luna's flight-software
controller.

```
  Rust (hopper_sil)                        C (generated from Simulink)
  luna LqrController  --u_cmd[4]-->  env_update()   advance 1 ms
         ^                           env_output()   state + sensors at time t
         +--- x_true[13], sensors --------+
```

- `u_cmd = [thrust, tvc_pitch, tvc_yaw, rcs]`, applied on the next step
  (`Environment/u_delay`), as on the vehicle.
- `x_true = [pos NED (3), vel (3), body rates (3), quaternion q0..q3 (4)]`,
  the true state. The controller still flies on it until an estimator runs
  on the sensor data.
- Sensors (`Environment/Sensors`), in luna's units:

  | Sensor | Rate | Output |
  |---|---|---|
  | IMU (ADIS16500) | 1 kHz | accel m/s^2, gyro deg/s |
  | Magnetometer (LIS2MDL) | 100 Hz | field, Gauss |
  | Barometer (MS5611) | 100 Hz | pressure Pa, temperature degC |
  | GPS (ZED-F9P) | 5 Hz, 0.1 s latency | lat/lon deg, alt m, NED velocity, fix flag, satellites |
  | LiDAR (4x TF-Luna) | 100 Hz | range per beam m, -1 when invalid |

  The models follow the Notion pages (Dynamics & GNC > Navigation >
  Measurement Models). Parameters are in `sensor_params.m`; values marked
  PLACEHOLDER are datasheet figures or guesses to replace with test data.
  LiDAR is a simple flat-ground stand-in for the ray-traced model on Notion.

## Workflow

In MATLAB, from `sims/hopper_6dof_sizer` (fresh session each time; the build
refuses if `hopper_env` is already loaded):

```matlab
addpath fsw_bridge
add_sensors     % only after editing sensor_blocks/*.m: rebuilds Environment/Sensors
check_sensors   % optional: measurement errors vs truth
build_env_c     % Environment -> C in build/hopper_env_ert_rtw (~3 min)
export_gains    % LQR gain tables -> build/gains/*.csv
verify_env_c    % optional: C vs Simulink replay check (gusts off)
```

The sensor block code lives in `sensor_blocks/*.m`; edit it there and run
`add_sensors`, not the blocks in the model.

Then the Rust loop (needs Rust with the `x86_64-pc-windows-gnu` toolchain and
the luna repo cloned next to hopper):

```bash
cd fsw_bridge/sil
HOPPER_CC=<MinGW>/bin/gcc.exe cargo run --release
```

`HOPPER_CC` is the MinGW gcc MATLAB uses (`mex -setup C` shows it). Other
build settings: `HOPPER_ENV_CODE`, `MATLAB_ROOT`, `LUNA_DIR` (see `build.rs`).
The run writes `build/sil_run.csv` with a header row: time, true state,
command, truth logs, then every sensor reading, one row per 1 ms step.

To compare a Rust run with Simulink step for step, build against the
gust-free code (`HOPPER_ENV_CODE=../build/verify_nogust/hopper_env_ert_rtw`)
and run `compare_sil('build/sil_run.csv')` in MATLAB.

## Things to know

- **K2 sign.** luna `flight2/src/control/lqr.rs` computes
  `u = unom - K1*dx + K2`; the Simulink controller uses `- K2`. The team
  uses `- K2`, so hopper_sil negates the table by default. `--luna-k2-sign`
  flies lqr.rs as written (separates from Simulink at t = 0.18 s, lands
  ~0.6 s later).
- **Noise is repeatable.** Every sensor draws from its own xorshift stream
  seeded by `SENS.seed` (`sensor_params(seed)`), so a seed reproduces a run
  exactly, in Simulink and in the C code alike.
- **Gusts.** The gust blocks use Simulink's own random blocks, which draw
  different numbers in generated code, and the gust on/off generator's seed
  is `Auto`. Runs with gusts on do not match Simulink step for step.
- **One environment per process.** The generated code keeps its state in
  globals.
- **Parameters are inlined.** Changing a workspace value (mass, wind, sensor
  seed, ...) needs `build_env_c` again.
- **Landing is very sensitive.** With identical logic, 1e-12 differences
  grow to metres in the last ~1.5 s before touchdown, in Simulink too.
