/*
 * replay_main.c - open-loop check of the generated hopper_env C code.
 *
 * Reads a command history (t, thrust, tvc_pitch, tvc_yaw, rcs per line) and
 * runs one 1 ms step per row the way a flight-software loop would:
 * hopper_env_output() gives the outputs at time t, the row's command is
 * applied, and hopper_env_update() advances to t + 1 ms. Writes
 * (time, x_true[13], thrust, z, imu[6], mag[3], baro[2], gps[8], lidar[4])
 * per step for comparison with the Simulink run that produced the commands.
 *
 * usage: replay <commands.csv> <outputs.csv>
 */
#include <stdio.h>
#include <stdlib.h>

#include "hopper_env.h"

int main(int argc, char **argv)
{
  if (argc != 3) {
    fprintf(stderr, "usage: %s <commands.csv> <outputs.csv>\n", argv[0]);
    return 2;
  }

  FILE *in = fopen(argv[1], "r");
  FILE *out = fopen(argv[2], "w");
  if (in == NULL || out == NULL) {
    fprintf(stderr, "cannot open input or output file\n");
    return 2;
  }

  hopper_env_initialize();

  double t, u[4];
  long steps = 0;
  while (fscanf(in, "%lf,%lf,%lf,%lf,%lf", &t, &u[0], &u[1], &u[2], &u[3]) == 5) {
    hopper_env_output();
    if (rtmGetErrorStatus(hopper_env_M) != NULL || rtmGetStopRequested(hopper_env_M)) {
      break;
    }

    fprintf(out, "%.17g", rtmGetT(hopper_env_M));
    for (int i = 0; i < 13; i++) {
      fprintf(out, ",%.17g", hopper_env_Y.x_true[i]);
    }
    fprintf(out, ",%.17g,%.17g", hopper_env_Y.thrust, hopper_env_Y.z);
    const double *sensors[] = {hopper_env_Y.imu, hopper_env_Y.mag, hopper_env_Y.baro,
                               hopper_env_Y.gps, hopper_env_Y.lidar};
    const int widths[] = {6, 3, 2, 8, 4};
    for (int s = 0; s < 5; s++) {
      for (int i = 0; i < widths[s]; i++) {
        fprintf(out, ",%.17g", sensors[s][i]);
      }
    }
    fputc('\n', out);

    for (int i = 0; i < 4; i++) {
      hopper_env_U.u_cmd[i] = u[i];
    }
    hopper_env_update();
    steps++;
  }

  hopper_env_terminate();
  fclose(in);
  fclose(out);

  const char *err = rtmGetErrorStatus(hopper_env_M);
  printf("replayed %ld steps%s%s\n", steps,
         rtmGetStopRequested(hopper_env_M) ? " (model requested stop)" : "",
         err != NULL ? err : "");
  return 0;
}
