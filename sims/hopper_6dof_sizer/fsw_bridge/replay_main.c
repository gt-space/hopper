/*
 * replay_main.c - open-loop check of the generated hopper_env C code.
 *
 * Reads a command history (t, thrust, tvc_pitch, tvc_yaw, rcs per line),
 * feeds one row per 1 ms step into hopper_env_step(), and writes the
 * resulting outputs (time, x_true[13], thrust, z) so they can be compared
 * against the Simulink run that produced the commands.
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
    if (rtmGetErrorStatus(hopper_env_M) != NULL || rtmGetStopRequested(hopper_env_M)) {
      break;
    }

    for (int i = 0; i < 4; i++) {
      hopper_env_U.u_cmd[i] = u[i];
    }

    /* The step applies u_cmd, integrates one fixed step, and leaves
     * hopper_env_Y describing the end of the step, so label the row with the
     * model time after the call (t + 1 ms), not the command time. */
    hopper_env_step();

    fprintf(out, "%.17g", rtmGetT(hopper_env_M));
    for (int i = 0; i < 13; i++) {
      fprintf(out, ",%.17g", hopper_env_Y.x_true[i]);
    }
    fprintf(out, ",%.17g,%.17g\n", hopper_env_Y.thrust, hopper_env_Y.z);
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
