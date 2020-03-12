#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "gammafit_search.h"
#include "read_array.h"
#include "unpin.h"
#include "gen_control.h"

int
main (int argc, char **argv)
{
  int *spiketime = 0;
  int *lodif = 0;
  int spikecount, n;
  int *control, control_spikecount;
  double g;
  bool use_seed = false;
  unsigned long seed = 0;
  
  if (argc < 2) {
    printf ("usage: %s unit_id\n", argv[0]);
    return 0;
  }
  if (argc > 2) {
    use_seed = true;
    seed = strtoul (argv[2], NULL, 0);
  }

  spikecount = read_array ("all", argv[1], &spiketime);
  lodif = unpin (spiketime, spikecount);
  g = gammafit_search (spiketime, lodif, spikecount);

  control = gen_control (spiketime, spikecount, g, &control_spikecount, seed, use_seed);

  for (n = 0; n < control_spikecount; n++)
    printf ("%d\n", control[n]);

  return 0;
}
