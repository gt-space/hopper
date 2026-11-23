# CTV hop logs

`ctv_1m_hop.csv`, `ctv_2m_hop.csv` and `ctv_4m_hop.csv` (#47) are 100 Hz logs of the CTV flying 1, 2 and 4 m hops in the harness. Each row is the vehicle's truth state and its thrust at time `t`. No hopper code reads them.

| Column | Unit | Meaning |
|---|---|---|
| `t` | s | from the start of the flight |
| `pos_n`, `pos_e`, `pos_d` | m | the sim's position state (the controller's x(1:3)), NED, 0 on the pad. **`pos_d` is positive down**: the 4 m apex is `pos_d` ≈ −4.0 |
| `vel_n`, `vel_e`, `vel_d` | m/s | NED velocity, `vel_d` positive down |
| `p`, `q`, `r` | rad/s | body rates (body x up the thrust axis) |
| `q0`…`q3` | – | attitude quaternion, scalar first, NED to body; vertical is (√2/2, 0, √2/2, 0) |
| `thrust_N` | N | EDF thrust the plant flew, in its 0.05 N throttle steps; 0 before liftoff and after the touchdown cutoff |

The headers were `t,X,Y,Z,Vx,Vy,Vz,P,Q,R,q0,q1,q2,q3,T` until 2026-11-22. They were renamed because a lower-casing reader sees `t` twice, and `Z` and `T` read as a height (positive up) and a thrust they are not.

## What they are not

- **Not the controller's reference.** The LQR's gain schedule, reference state and nominal command come from `LinerizationMaster` / `src/gnc/ctvTuningGains.m` along the profile in `inputs/ctv_inputs.m` (4.5 m apex, 6 s up, 2 s hover, 6 s down, 2 s terminal, plus run-out). These hops have no gain tables, and their row times don't match any `tgrid`. Export the flown schedule with `fsw_bridge/export_plant_gains`.
- **Not a command history.** There are no TVC or roll-torque columns.

Use them as a closed-loop outcome check: fly the flight software against the generated plant (`fsw_bridge/build_plant_c`) on the same hop and compare the trajectory envelope.
