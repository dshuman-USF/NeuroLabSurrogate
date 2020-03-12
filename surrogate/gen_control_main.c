#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gen_control.h"
#include "read_array.h"

int
main (int argc, char **argv)
{
  double g;
  int *spiketime = 0, spikecount;
  int *control, control_spikecount;
  int n;
  bool use_seed = false;
  unsigned long seed = false;

  if (argc < 3) {
    printf ("usage: %s unit_id gamma [seed]\n", argv[0]);
    return 0;
  }
  
  spikecount = read_array ("all", argv[1], &spiketime);
  g = atof (argv[2]);
  if (argc > 3) {
    use_seed = true;
    seed = strtoul (argv[3], NULL, 0);
  }
  control = gen_control (spiketime, spikecount, g, &control_spikecount, seed, use_seed);
  
  for (n = 0; n < control_spikecount; n++)
    printf ("%d\n", control[n]);

  return 0;
}
