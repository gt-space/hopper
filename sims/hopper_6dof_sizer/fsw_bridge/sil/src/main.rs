//! hopper_sil: luna's Rust LQR controller flying the Simulink-generated
//! hopper environment, in lockstep, one 1 ms step at a time.
//!
//! ```text
//! cargo run --release -- [--gains DIR] [--out FILE] [--luna-k2-sign] [--t-max S]
//! ```
//!
//! Writes a CSV (with a header row) of one row per step: time, true state,
//! command, truth logs, then every sensor measurement (see `CSV_HEADER`).
//!
//! The controller still flies on the true state; the sensor readings are
//! logged for the estimator that will sit between them.

// luna's lqr.rs imports `common::comm::ctv`; point `common` at this crate.
extern crate self as common;

/// luna `common::comm::ctv`, compiled straight from the luna checkout.
pub mod comm {
    pub mod ctv {
        include!(concat!(env!("LUNA_DIR"), "/common/src/comm/ctv.rs"));
    }
}

/// Mirrors luna `flight2/src/control.rs`; `lqr` is luna's file unchanged.
mod control {
    use crate::comm::ctv::{ControlState, ControlVector};

    pub trait Controller {
        /// Perform one step of the control algorithm
        fn step(&mut self, state: ControlState) -> ControlVector;
    }

    #[allow(dead_code, unused_imports)]
    pub mod lqr {
        include!(concat!(env!("LUNA_DIR"), "/flight2/src/control/lqr.rs"));
    }
}

mod env;
mod gains;

use std::{
    env as std_env,
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
    process::ExitCode,
    time::{Duration, Instant},
};

use comm::ctv::{ControlState, Quaternion, Vector3};
use control::Controller;
use env::Environment;

const CSV_HEADER: &str = "t,\
pos_n,pos_e,pos_d,vel_n,vel_e,vel_d,p,q,r,q0,q1,q2,q3,\
u_thrust,u_tvc_pitch,u_tvc_yaw,u_rcs,\
thrust,z,ox_mass,fuel_mass,\
acc_x,acc_y,acc_z,gyro_x_dps,gyro_y_dps,gyro_z_dps,\
mag_x,mag_y,mag_z,\
baro_pa,baro_temp_c,\
gps_lat_deg,gps_lon_deg,gps_alt_m,gps_vn,gps_ve,gps_vd,gps_fix,gps_sats,\
lidar_1,lidar_2,lidar_3,lidar_4";

struct Args {
    gains: PathBuf,
    out: PathBuf,
    negate_k2: bool,
    t_max: f64,
}

fn parse_args() -> Result<Args, String> {
    let build = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../build");
    let mut args = Args {
        gains: build.join("gains"),
        out: build.join("sil_run.csv"),
        // Simulink's -K2, the team's chosen sign (luna's lqr.rs has +K2)
        negate_k2: true,
        t_max: 60.0,
    };
    let mut it = std_env::args().skip(1);
    while let Some(arg) = it.next() {
        match arg.as_str() {
            "--gains" => args.gains = it.next().ok_or("--gains needs a folder")?.into(),
            "--out" => args.out = it.next().ok_or("--out needs a file")?.into(),
            "--luna-k2-sign" => args.negate_k2 = false,
            "--t-max" => {
                args.t_max = it
                    .next()
                    .ok_or("--t-max needs seconds")?
                    .parse()
                    .map_err(|e| format!("--t-max: {e}"))?
            }
            other => return Err(format!("unknown argument {other}")),
        }
    }
    Ok(args)
}

/// Environment truth state → the controller's input.
fn control_state(t: f64, x: &[f64; 13]) -> ControlState {
    ControlState {
        time: Duration::from_secs_f64(t),
        position: Vector3::new(x[0], x[1], x[2]),
        velocity: Vector3::new(x[3], x[4], x[5]),
        body_rate: Vector3::new(x[6], x[7], x[8]),
        attitude: Quaternion::new(x[9], x[10], x[11], x[12]),
    }
}

fn main() -> ExitCode {
    let args = match parse_args() {
        Ok(a) => a,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::from(2);
        }
    };

    let mut controller = match gains::load_lqr(&args.gains, args.negate_k2) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("error loading gains (run export_gains in MATLAB): {e}");
            return ExitCode::FAILURE;
        }
    };
    let file = match File::create(&args.out) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("error: {}: {e}", args.out.display());
            return ExitCode::FAILURE;
        }
    };
    let mut csv = BufWriter::new(file);
    if let Err(e) = writeln!(csv, "{CSV_HEADER}") {
        eprintln!("error writing {}: {e}", args.out.display());
        return ExitCode::FAILURE;
    }

    println!(
        "K2 sign: {}",
        if args.negate_k2 { "Simulink (-K2)" } else { "luna as written (+K2)" }
    );

    let wall = Instant::now();
    let mut world = Environment::new();
    let mut steps = 0u64;
    let mut max_alt = f64::NEG_INFINITY;
    let mut gps_fixes = 0u64;
    let mut last_gps_lat = f64::NAN;

    loop {
        world.output();
        if let Some(err) = world.error() {
            eprintln!("environment error at t = {:.3}: {err}", world.time());
            return ExitCode::FAILURE;
        }
        let t = world.time();
        if world.stop_requested() || t > args.t_max {
            break;
        }

        let x = world.x_true();
        let u = controller.step(control_state(t, &x));
        let u = [u.thrust, u.tvc_pitch, u.tvc_yaw, u.rcs_torque];
        world.set_u_cmd(u);

        let truth = world.truth();
        max_alt = max_alt.max(-truth.z);
        let s = world.sensors();
        if s.gps.has_fix && s.gps.latitude_deg != last_gps_lat {
            gps_fixes += 1;
            last_gps_lat = s.gps.latitude_deg;
        }
        let gps = [
            s.gps.latitude_deg,
            s.gps.longitude_deg,
            s.gps.altitude_m,
            s.gps.north_mps,
            s.gps.east_mps,
            s.gps.down_mps,
            f64::from(u8::from(s.gps.has_fix)),
            f64::from(s.gps.num_satellites),
        ];
        let row = std::iter::once(t)
            .chain(x)
            .chain(u)
            .chain([truth.thrust, truth.z, truth.ox_mass, truth.fuel_mass])
            .chain(s.imu.accelerometer)
            .chain(s.imu.gyroscope)
            .chain(s.magnetometer)
            .chain([s.barometer.pressure, s.barometer.temperature])
            .chain(gps)
            .chain(s.lidar.map(|r| r.unwrap_or(-1.0)))
            .map(|v| format!("{v:.17e}"))
            .collect::<Vec<_>>()
            .join(",");
        if let Err(e) = writeln!(csv, "{row}") {
            eprintln!("error writing {}: {e}", args.out.display());
            return ExitCode::FAILURE;
        }

        world.update();
        steps += 1;
    }

    let elapsed = wall.elapsed().as_secs_f64();
    let t_end = world.time();
    println!(
        "flew {steps} steps to t = {t_end:.3} s{} | max altitude {max_alt:.3} m | {gps_fixes} GPS fixes",
        if world.stop_requested() { " (touchdown stop)" } else { "" }
    );
    println!(
        "wall time {elapsed:.2} s ({:.1}x real time), wrote {}",
        t_end / elapsed,
        args.out.display()
    );
    ExitCode::SUCCESS
}
