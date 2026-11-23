# FSW bridge: Rust flight software in the loop with the 6DOF sim

## The harness route (CONTROLLER=fsw)

The plant now comes straight from the variant harness,
`models/hopper_6dof.slx`, with `CONTROLLER=fsw` (docs/VARIANTS.md,
*Flight software in the loop*): the whole vehicle minus the controller,
which the flight software replaces through the root ports `u_fsw`
(thrust N, TVC pitch rad, TVC yaw rad, RCS torque N m) and `x_fsw` (position
NED, velocity NED, body rates, quaternion). No model is copied.

It also puts out the harness's shared sensor model (`src/sensors`), so the
flight software can fly its own navigation filter: `sens_meas` (26),
`sens_valid` (9, 1/0) and `sens_stamp` (9, s), laid out as
`src/sensors/sensor_layout.m` and taken after the fault gates. A sensor's
stamp changes when it takes a new sample. `plant.json` records the sensor
parameters baked in (`sensor_params`, after the vehicle config and any
`'SensorParams'` override), the layout, `cg_height` (the full-tank CG
above the base on the pad, m: sensor positions in `sensor_params` are from
the base, while the GPS and LiDAR outputs are referred to the CG), the
channel names under `outputs`,
`fault_injection` and the fault table (`'Faults'`, used only with
`FAULT_INJECTION=on` in `'Variants'`). Parameters are inlined: a different
GPS sigma or fault means a rebuild, e.g.
`build_plant_c('SensorParams', struct('gps', struct('sigma', [1.5; 1.61; 4.43])))`.

IMU placement: by default the accelerometer reads the specific force at
the CG (`imu_from_cg` = 0 in `plant.json`). Setting `imu.mount` (m from
the base, body axes, like the other sensor mounts) puts it there, and the
accelerometer then also reads alpha × r + ω × (ω × r), with r from the CG
and alpha the true rate's backward difference over the IMU step:
`build_plant_c('SensorParams', struct('imu', struct('mount', [0.45; 0.03; 0])))`.
`plant.json`'s `imu_from_cg` is that mount minus the pad CG (`cg_height`
along body x).

Ground: `build_plant_c('Ground', 'strut')` builds the plant with leg
springs (`GROUND=strut`) and no touchdown stop, so the vehicle lands,
comes to rest and stays flying in the plant until the flight software
stops; the IMU sees the impact and the rest. `plant.json` says `ground`
and `touchdown_stop`. Without it (the default) the plant stops at the first
descending contact, but only after the vehicle has been 0.1 m up: one that
never climbs that far falls through the pad and does not stop. With no
ground the vehicle also sinks about 0.1 m while thrust builds, so fly the
strut plant for liftoff. On the struts the plant's own engine cutoff
(`IN.mission.engine_cutoff_alt`, 5 mm on the CTV) still cuts thrust and
does not latch, and `u_fsw` thrust below the throttle floor (7 N on the
CTV) is raised to the floor. The CTV's leg stiffness is the sizer's
(2.19 MN/m per leg, critically damped per leg) and its feet are an
ASSUMPTION (4 at 0.15 m radius): use it to exercise touchdown logic, not
for landing loads. Only the feet touch the ground: a vehicle that tips
over falls through the pad. `verify_plant_c('Ground', 'strut')` checks the
C against Simulink through contact and 2 s past it (its open-loop replay
crashes and tips, which is fine for an equality check). `build_plant_c`
writes the strut build to its own folder if you pass one, e.g.
`build_plant_c('Ground', 'strut', 'Folder', 'fsw_bridge/build/plant_strut')`.

| Script | What it does |
|---|---|
| `build_plant_c` | Embedded Coder C for the harness with `CONTROLLER=fsw` (default `VEHICLE=ctv`, ARM64), packed with every header into `build/plant/` with `sources.txt`, `defines.txt` and `plant.json` |
| `export_plant_gains` | The scheduled LQR tables the vehicle's controller flies (for the CTV, `CONTROLLER=lqr_ctv_tuning`: `K1flat_ctv`/`K2grid_ctv`) as CSV in `build/gains/<vehicle>/`, plus `gains.json` with the controller, the K1 layout (row-major), the law, a SHA-256 per CSV and the git commit |
| `export_controller_vectors` | Golden `(t, x) -> u` cases for the flight software's LQR, next to the gains: `controller_vectors.csv` and `.json`. `u` comes from the model's own `Gain_Selector` code run on the exported tables. The cases cover the nominal flight every 10 ms, grid points, dispersed states with half the quaternions negated, and times past both ends. |
| `verify_plant_c` | Replays the LQR flight's commands at 100 Hz through Simulink (`CONTROLLER=fsw`) and through the generated C (`plant_replay.c`, gcc in WSL), and compares every state and every sensor port. Gusts off. |

The flight software's side, the Rust wrapper and flight2's `--sim-plant`
bench mode, live in luna (`sitl/hopper_plant`, flight2 feature
`sim-plant`); its README has the bench procedure.

Versions: `plant.json` and `gains.json` both carry `git_commit`. The
flight software should refuse a gain set whose commit differs from the
plant's, and should check that `plant.json`'s `controller` is `fsw`
(by name: the index is 5 in `variantAxes`, but it was 4 on the unmerged
dev-fsw-harness branch, which `lqr_ctv_tuning` also took).

Timing: each 10 ms period the flight software reads `x_fsw` at the period's
start, and its command reaches the vehicle one plant step (1 ms) later
(`fsw_latency`), then holds for the rest of the period. Simulink and the
generated code do the same, so they match bit for bit.

## The original route (hopper_6dof_NED_v2_fswBridge.slx)

Kept until the harness route replaces it. The plant ("environment": actuators, propulsion mass, CG/MOI, slosh, wind,
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
