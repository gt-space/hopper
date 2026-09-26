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
