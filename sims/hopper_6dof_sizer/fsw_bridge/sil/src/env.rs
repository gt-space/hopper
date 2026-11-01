//! Safe wrapper over the Simulink-generated hopper environment (see
//! `csrc/env_shim.c`).

use std::{
    ffi::CStr,
    os::raw::{c_char, c_int},
    sync::atomic::{AtomicBool, Ordering},
};

extern "C" {
    fn env_initialize();
    fn env_output();
    fn env_update();
    fn env_terminate();
    fn env_time() -> f64;
    fn env_stop_requested() -> c_int;
    fn env_error() -> *const c_char;
    fn env_set_u_cmd(u: *const f64);
    fn env_get_x_true(x: *mut f64);
    fn env_get_truth(out: *mut f64);
    fn env_get_sensors(
        imu: *mut f64,
        mag: *mut f64,
        baro: *mut f64,
        gps: *mut f64,
        lidar: *mut f64,
    );
}

/// IMU sample, same fields and units as luna `fc_sensors::Imu`.
#[derive(Debug, Clone, Copy)]
pub struct Imu {
    /// m/s^2, sensor axes
    pub accelerometer: [f64; 3],
    /// deg/s, sensor axes
    pub gyroscope: [f64; 3],
}

/// Barometer sample, same fields and units as luna `fc_sensors::Barometer`.
#[derive(Debug, Clone, Copy)]
pub struct Barometer {
    /// Pa
    pub pressure: f64,
    /// degC
    pub temperature: f64,
}

/// GPS solution, same fields and units as luna `comm::GpsState`.
#[derive(Debug, Clone, Copy)]
pub struct Gps {
    pub latitude_deg: f64,
    pub longitude_deg: f64,
    pub altitude_m: f64,
    pub north_mps: f64,
    pub east_mps: f64,
    pub down_mps: f64,
    pub has_fix: bool,
    pub num_satellites: u8,
}

/// Latest output of every sensor model. Each sensor updates at its own rate
/// (IMU 1 kHz, magnetometer/barometer/LiDAR 100 Hz, GPS 5 Hz) and holds its
/// value in between.
#[derive(Debug, Clone, Copy)]
pub struct Sensors {
    pub imu: Imu,
    /// Gauss, sensor axes (luna `fc_sensors::Magnetometer`)
    pub magnetometer: [f64; 3],
    pub barometer: Barometer,
    pub gps: Gps,
    /// Range along each beam, m; `None` when there is no valid return.
    pub lidar: [Option<f64>; 4],
}

/// The generated code keeps all model state in globals, so only one
/// environment can exist per process.
static IN_USE: AtomicBool = AtomicBool::new(false);

/// Truth signals logged alongside the state (same as the Monte Carlo logs).
#[derive(Debug, Clone, Copy)]
pub struct Truth {
    /// Delivered thrust after the actuator model, N.
    pub thrust: f64,
    /// Down position logged by the model, m (negative is up).
    pub z: f64,
    pub ox_mass: f64,
    pub fuel_mass: f64,
}

/// The simulated world: actuators, propulsion, slosh, wind, 6DOF dynamics.
///
/// One cycle is: [`output`](Self::output) to compute the state at the
/// current time, read it, [`set_u_cmd`](Self::set_u_cmd), then
/// [`update`](Self::update) to advance one fixed step (1 ms). The command
/// takes effect on the following step, as on the vehicle.
pub struct Environment {
    _not_send: std::marker::PhantomData<*const ()>,
}

impl Environment {
    /// Initializes the model at t = 0.
    ///
    /// # Panics
    /// If another `Environment` is alive in this process.
    pub fn new() -> Self {
        assert!(
            !IN_USE.swap(true, Ordering::SeqCst),
            "only one Environment per process (generated code uses globals)"
        );
        unsafe { env_initialize() };
        Environment {
            _not_send: std::marker::PhantomData,
        }
    }

    /// Computes the outputs for the current time.
    pub fn output(&mut self) {
        unsafe { env_output() }
    }

    /// Applies the last command and advances one fixed step.
    pub fn update(&mut self) {
        unsafe { env_update() }
    }

    /// Simulation time, s.
    pub fn time(&self) -> f64 {
        unsafe { env_time() }
    }

    /// The model's own stop condition (touchdown) has fired.
    pub fn stop_requested(&self) -> bool {
        unsafe { env_stop_requested() != 0 }
    }

    /// Error reported by the generated code, if any.
    pub fn error(&self) -> Option<String> {
        let msg = unsafe { env_error() };
        if msg.is_null() {
            None
        } else {
            Some(unsafe { CStr::from_ptr(msg) }.to_string_lossy().into_owned())
        }
    }

    /// Sets the command `[thrust, tvc_pitch, tvc_yaw, rcs]` for the next update.
    pub fn set_u_cmd(&mut self, u: [f64; 4]) {
        unsafe { env_set_u_cmd(u.as_ptr()) }
    }

    /// True state `[pos NED (3), vel (3), body rates (3), quaternion q0..q3 (4)]`.
    pub fn x_true(&self) -> [f64; 13] {
        let mut x = [0.0; 13];
        unsafe { env_get_x_true(x.as_mut_ptr()) };
        x
    }

    /// Sensor measurements at the current time.
    pub fn sensors(&self) -> Sensors {
        let (mut imu, mut mag, mut baro, mut gps, mut lidar) =
            ([0.0; 6], [0.0; 3], [0.0; 2], [0.0; 8], [0.0; 4]);
        unsafe {
            env_get_sensors(
                imu.as_mut_ptr(),
                mag.as_mut_ptr(),
                baro.as_mut_ptr(),
                gps.as_mut_ptr(),
                lidar.as_mut_ptr(),
            )
        };
        Sensors {
            imu: Imu {
                accelerometer: [imu[0], imu[1], imu[2]],
                gyroscope: [imu[3], imu[4], imu[5]],
            },
            magnetometer: mag,
            barometer: Barometer {
                pressure: baro[0],
                temperature: baro[1],
            },
            gps: Gps {
                latitude_deg: gps[0],
                longitude_deg: gps[1],
                altitude_m: gps[2],
                north_mps: gps[3],
                east_mps: gps[4],
                down_mps: gps[5],
                has_fix: gps[6] != 0.0,
                num_satellites: gps[7] as u8,
            },
            lidar: lidar.map(|r| (r >= 0.0).then_some(r)),
        }
    }

    pub fn truth(&self) -> Truth {
        let mut t = [0.0; 4];
        unsafe { env_get_truth(t.as_mut_ptr()) };
        Truth {
            thrust: t[0],
            z: t[1],
            ox_mass: t[2],
            fuel_mass: t[3],
        }
    }
}

impl Drop for Environment {
    fn drop(&mut self) {
        unsafe { env_terminate() };
        IN_USE.store(false, Ordering::SeqCst);
    }
}
