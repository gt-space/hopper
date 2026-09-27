/*
 * env_shim.c - flat C API over the generated hopper_env code.
 *
 * Rust calls these instead of touching hopper_env_U / hopper_env_Y /
 * hopper_env_M directly, so the Rust side does not depend on the generated
 * struct layouts (which change whenever the model's ports change).
 */
#include "hopper_env.h"

void env_initialize(void) { hopper_env_initialize(); }

/* Compute outputs for the current time; read them with the getters below. */
void env_output(void) { hopper_env_output(); }

/* Apply the command set with env_set_u_cmd and advance one fixed step. */
void env_update(void) { hopper_env_update(); }

void env_terminate(void) { hopper_env_terminate(); }

double env_time(void) { return rtmGetT(hopper_env_M); }

int env_stop_requested(void) { return rtmGetStopRequested(hopper_env_M) ? 1 : 0; }

const char *env_error(void) { return rtmGetErrorStatus(hopper_env_M); }

/* u = [thrust, tvc_pitch, tvc_yaw, rcs] */
void env_set_u_cmd(const double u[4])
{
  for (int i = 0; i < 4; i++) {
    hopper_env_U.u_cmd[i] = u[i];
  }
}

/* x = [pos NED (3), vel (3), body rates (3), quaternion q0..q3 (4)] */
void env_get_x_true(double x[13])
{
  for (int i = 0; i < 13; i++) {
    x[i] = hopper_env_Y.x_true[i];
  }
}

/* Sensor measurements (Environment/Sensors), in luna units:
 * imu   = [accel xyz m/s^2, gyro xyz deg/s]
 * mag   = field xyz, Gauss
 * baro  = [pressure Pa, temperature degC]
 * gps   = [lat deg, lon deg, alt m, vN, vE, vD m/s, has_fix, num_sats]
 * lidar = range per sensor, m (-1 = no valid return) */
static void copy(double *dst, const double *src, int n)
{
  for (int i = 0; i < n; i++) {
    dst[i] = src[i];
  }
}

void env_get_sensors(double imu[6], double mag[3], double baro[2], double gps[8],
                     double lidar[4])
{
  copy(imu, hopper_env_Y.imu, 6);
  copy(mag, hopper_env_Y.mag, 3);
  copy(baro, hopper_env_Y.baro, 2);
  copy(gps, hopper_env_Y.gps, 8);
  copy(lidar, hopper_env_Y.lidar, 4);
}

/* Truth signals the Monte Carlo logs: thrust, z, ox_mass, fuel_mass */
void env_get_truth(double out[4])
{
  out[0] = hopper_env_Y.thrust;
  out[1] = hopper_env_Y.z;
  out[2] = hopper_env_Y.ox_mass;
  out[3] = hopper_env_Y.fuel_mass;
}
