# FDIR and fault injection in the 6DOF sim

How the sim models the Hopper's fault detection, isolation and recovery, and how it injects the faults that exercise it.

**Sources:** the Hopper FDIR CoDR (GNC section, slides 19–28) and the Hopper GNC PDR (flight modes, slides 10–11; FDIR, slides 77–86). Where the two differ, the PDR wins.

Both features are variant axes of `models/hopper_6dof.slx` (see [VARIANTS.md](VARIANTS.md)):
- **`FAULT_INJECTION`:** `off` / `on`.
- **`FDIR`:** `off` / `on`.

With both `off`, the harness flies exactly as v2 does.

## Frame

The PDR's per-frame order runs at 100 Hz, and the sim runs it in the same order:

1. **Sensor ingest:** `Sensors`, then the fault-injection gates.
2. **Nav propagate and update:** the truth state until NAVIGATION is wired in.
3. **G&C feasibility.**
4. **FDIR executive.**
5. **Mode manager.**
6. **Actuator output:** the actuator gates, then the plant.

A mode change takes effect on the next frame.

```
truth ─► Sensors ─► [sensor gates] ─► (nav) ─► Guidance ─► Control ─► [actuator gates] ─► Actuators ─► plant
                         │                         ▲            ▲              │
                         ▼                         │            │              ▼
                    FDIR monitors ◄── estimator port (NAVIGATION) ◄── actuator feedback
                         │
                         ▼
          persistence ► isolation ► latched fault ► response table ► Mode manager
                                                                       │ active mode, response level,
                                                                       │ masks, gain set, abort reference,
                                                                       ▼ safing commands
```

## Mode manager

The flight-mode state machine follows PDR slides 10–11.

| Mode | Phase mask | Commanded behavior | Nominal exit |
|---|---|---|---|
| SAFE | PAD | actuation inhibited | initialization requested |
| INITIALIZE | PAD | nav origin, biases, references | initialization succeeds |
| READY | PAD | evaluate readiness, including the IBIT checks | arm command, arm conditions met |
| ARMED | PAD | mission configuration locked | start command |
| ENGINE START | STARTUP | engine startup; held on the pad | stable, flight-capable engine |
| PAD CHECKOUT | STARTUP | final launch checks | launch authorized, conditions met |
| ASCENT | ASCENT | track the ascent profile | enter the hover region |
| HOVER | HOVER | hold the hover reference, time the hover | hover complete |
| DESCENT | DESCENT | track the descent profile | terminal-descent entry |
| TERMINAL DESCENT | DESCENT | final approach, tighter limits | thrust cut at `IN.mission.engine_cutoff_alt` |
| TOUCHDOWN CONFIRM | LANDED | confirm contact, coordinate shutdown | touchdown and shutdown confirmed |
| LANDED | LANDED | safe post-flight state | ground reset |
| ABORT LAND (R2) | by altitude | abort-to-land reference | confirmed touchdown |
| NAV-DEGRADED ABORT (R3) | by altitude | shadow estimator, then land | confirmed touchdown |
| SAFED (R4) | — | throttle closed, gimbal centered, RCS off | ground reset |

**Ground commands in the sim:** the sim authorizes init, arm, start and launch automatically, and the pre-launch modes pass through in a configurable pad time. A test can instead issue them from a command schedule, including the manual abort (R2) and terminate (R4), which are always accepted.

**Transitions:** they use the PDR's conditions (hover capture limits, terminal-descent entry) and run on the sim's state, not on the reference timing.

## Fault injection

The CoDR places injection in "gates that restrict / modify signals between blocks". Each gate reads one fault table, `faults`, built by `faultScenario`. A row gives `target`, `mode`, `t_on`, `t_off`, `magnitude`, plus a seed for random modes. A Monte Carlo can disperse it, and it applies to the CTV as well as the hopper.

| Gate | Targets | Modes |
|---|---|---|
| Sensor | IMU accel, IMU gyro, lidar beams 1–4, GPS, baro, mag | loss (NaN / invalid flag), bias, drift (ramp), saturation, stuck (hold), noise burst, jump, stale timestamp / age, poor fix (GPS quality) |
| Gimbal | pitch, yaw | stuck, slow (rate cut), hardover, bias |
| RCS fan | fans 1–4 | dead, reduced thrust |
| Throttle | main engine | stuck thrust (T/W below 0.9, at 1, above 1.1), reduced authority |
| Estimator | reserved for NAVIGATION | divergence, NaN |

## Sensors

The harness's shared `Sensors` block (`src/sensors/sensor_frame.m`: the sensor model, then the sensor fault gate) runs at the IMU rate, 500 Hz on the hopper. It is active when FDIR is on or NAVIGATION isn't `none`. FDIR reads it through rate transitions to its 100 Hz frame, and the navigation filter reads it at 500 Hz, so FDIR, the filter and fault injection all see the same measurements.

Channels (`sensor_layout`):
- IMU specific force and rates
- barometric altitude
- four lidar beam ranges
- GPS position and velocity of antenna 1, with fix quality, satellite count and age
- the differential-GPS baseline, antenna 1 to 2 in NED
- magnetometer

**Geometry:** lidar mounts and beams and the GPS antennas come from the vehicle config (`VC.sensors`), placed from the true CG. The same geometry goes to the filter (`lidar_geometry`, `nav_measurements`).

**Validity and timing:** every channel carries a valid flag and a timestamp, which the raw checks need. Each sensor samples at its own rate and holds between samples. FDIR's IMU timing check expects one sample per frame, since it sees the 500 Hz IMU at 100 Hz.

**Noise:** seeded noise sources, so runs repeat exactly. The sigmas are `sim_sensor_params`; the navigation adapter hands the same sigmas to the filter as each measurement's covariance.

## FDIR executive

The executive runs the CoDR and PDR ladder: monitor → persistence → isolation → latched fault flag → response table → mode manager.

**Fault table:** `inputs/fdir_fault_catalog.csv`, readable as text and diffable. Each row has these fields:
- `id` and `fta_row`
- `monitor`, the function that evaluates it
- `filter`: `N_OF_M`, `CONSECUTIVE` or `CUSUM`, with N, M and the CUSUM limit
- `clear_policy`: `NEVER_IN_FLIGHT`, `AUTO_CLEAR` or `GROUND_ONLY`
- `phase_mask`: any of PAD, STARTUP, ASCENT, HOVER, DESCENT, LANDED
- `response`: R0–R4, or NO_ARM for IBIT checks
- `target`: the task the response acts on
- `isolated_action`, for example hold last sample
- `escalate_to`, for example R1 then R2

The GNC rows below are filled in, with thresholds marked provisional. The PDR says numeric thresholds come from Monte Carlo false-alarm against missed-detection trades.

**Fault state:** INACTIVE → PENDING → LATCHED, with an M-bit history, first-latch frame, fail count and an inhibited flag. Once latched, a fault clears only according to its clear policy.

**Arbitration:** the active response is the maximum over all latched faults, and it only escalates. R0 and R1 stack without changing mode.

**R4 permission:** always allowed on the pad. In flight, it's allowed only below the safe height `h_cut` (≈1–2 m, set by the legs) or when the vehicle is uncontrollable.

## GNC monitors

**Sensors** (CoDR 20–22, PDR 80):

| Sensor | Checks | Response |
|---|---|---|
| IMU | saturation above 0.95×full scale; stuck (windowed variance); rate of change (jerk bound); timestamp dt within [0.5, 1.5]× nominal | isolated: hold last sample; persistent: gyro R4, accel R3 |
| Lidar, per beam | valid flag and signal strength; range bounds with the pad blind zone masked; rate bound | invalid sample; persistent R0 |
| GPS | 3D fix and more than 6 satellites; age below ~200 ms; jump test against the velocity estimate; pad multipath mask | don't fuse; persistent jump R0 |
| Barometer | range, rate, stuck; ground-effect window near the pad | inflate R near the ground; persistent R0 |
| IBIT | IMU in calibrated range; lidar beams read pad height; GPS origin calibrated | no-arm |

**Guidance and control** (CoDR 27–28, PDR 82):

| Check | Method | Response |
|---|---|---|
| Loss-of-control bound | tilt error beyond the limit for a set time | R1, then R2 |
| Attitude envelope | tilt error or pitch/yaw rate beyond the hard limit (roll rate is the roll-rate monitor, R2) | R4 under the R4 permission above, otherwise R2 |
| Limit cycle | sign-change count on the error rate | R1, then R2 |
| Thrust margin | throttle above 90% of maximum for longer than T | R2 |
| Roll rate | roll-rate bound exceeded | R2 |
| Reference integrity | checksum of the selected trajectory; t=0 reference matches the pad state | no-arm |

**Actuators** (CoDR 28, PDR 83):

| Check | Method | Response |
|---|---|---|
| Gimbal saturation | duty at the limit over a window | high duty R1, sustained R2 |
| Gimbal tracking | command vs. feedback residual | stuck or slow R2, hardover R4 |
| RCS allocation | residual of the allocation | sustained nonzero R1 |
| RCS fan health | delivered vs. commanded thrust (stands in for ESC RPM) | R0: zero that fan |
| Stuck thrust | from the T/W table | below 0.9: hold attitude, dither, land; at 1: dither, slow descent; above 1.1: dither, then R4 |

**Estimator** (CoDR 23–26, PDR 81; `fdir_nav_monitors`), with NAVIGATION on. FDIR gets the filter's 20-state estimate and covariance and checks each new sensor sample against the estimate with the nav team's own measurement models and 99 % gates (NIS = rᵀS⁻¹r, S = HPHᵀ + R):

| Row | Check | Response |
|---|---|---|
| 33 | exactly one sensor group (baro, gyro, GPS, lidar, differential GPS) fails its gate | R0 |
| 34 | two or more fail together | R1, R3 if it persists 2 s |
| 35 | CUSUM of the mean normalized NIS above 2, per second | R1 |
| 36 | main vs shadow estimator | inhibited: there is no second estimator yet |
| 37 | covariance or state not finite, a variance ≤ 0 or beyond its bound, quaternion off unit length, speed or rate beyond physical limits | R3 |

The filter gates on its prior; these checks run on the posterior, so a sample the filter accepted passes and one it rejected, or one that pulled the estimate away, fails. With NAVIGATION `none` the rows never fail. With NAVIGATION `gmekf` the G&C monitors also run on the estimate (what the flight computer knows), not on truth. Still open: sensor masking into the filter and Q/R inflation (the filter has no inputs for them). That limit shows: with NAVIGATION `gmekf`, a stuck gyro latches (row 5) and FDIR flies the abort, but the filter keeps fusing the stuck gyro, the abort flies on a corrupted attitude, and the vehicle crashes. On truth the same fault aborts and lands.

**Limit cycle (row 20) on an estimate:** the check counts sign changes of the tilt-error rate and now also needs the tilt to swing at least 0.5° about its 0.2 s average. Flying on the GMEKF's estimate the vehicle jitters about 0.1°, which gave 22 sign changes in 2 s and a false abort at 0.9 s; real limit cycles swing degrees.

## Responses

| Level | In the sim |
|---|---|
| R0 Mask | Zero a failed RCS fan in the allocation. Sensor masks go out on the mask bus; they act once NAVIGATION consumes them. |
| R1 Degrade | Switch the controller to a low-bandwidth gain set (`IN.gnc.degraded`). Q/R inflation waits for NAVIGATION. |
| R2 Abort to land | Capture the state. Command a vertical velocity ramp from the current v_z to the landing sink rate, zero horizontal velocity, and hold attitude. Control switches to a fixed-gain regulator (hover LQR at the captured mass). Exits only on confirmed touchdown. |
| R3 Nav-degraded abort | The R2 abort-to-land with no horizontal (N, E) position or velocity feedback, since that comes from the distrusted estimate: it descends on attitude, rates and height (`fdir_command_select`). The PDR's shadow vertical estimator would replace the estimate's height. |
| R4 Cutoff & safing | Throttle valves closed (thrust 0), gimbal centered, RCS off, first fault logged, final health frame recorded. |

**Touchdown confirmation** (PDR 85) needs two independent signals with N-of-M persistence:
- lidar median below the contact threshold
- an accelerometer impulse
- v_z near zero
- thrust cut commanded

It is armed only in TERMINAL DESCENT or an abort. One signal alone takes no action. On a missed confirmation after the cut, the actuators are safed and the run is flagged.

## Outputs

FDIR logs the active mode, response level, latched-fault bitfield, first-fault ID and time, per-fault states and counters, and every response command. That covers the PDR's telemetry requirement and the CoDR's downlink health frame. Monte Carlo summaries add the first fault and the final response level.

## Status

| Piece | State |
|---|---|
| `FDIR_SensorFault_Subsystem` (#31) | Its monitor, persistence and arbitration code becomes part of the executive above. Its ladder numbering and its four phases change to the PDR's, and `FDIR_Tests` is updated to match. |
| `hopper_6dof_NED_fault_injection.slx` | Contains no fault injection: it is a stale v2 copy. Retire it. |
| Navigation sensor models (`Src_sim`) | Not used: the shared `Sensors` block feeds the filter, so FDIR and the filter see the same measurements. |
| Navigation filter | The nav team's `Hopper_GMEKF` runs in the harness (NAVIGATION axis, docs/VARIANTS.md). The CTV filter can't be referenced yet (variable-step solver). |
| Open (PDR) | Numeric thresholds and persistence times from Monte Carlo; thrust cut at 0.1 m (PDR) vs. the sim's 0.05 m (`IN.mission.engine_cutoff_alt`); the FMECA catalog for the non-GNC rows. |
