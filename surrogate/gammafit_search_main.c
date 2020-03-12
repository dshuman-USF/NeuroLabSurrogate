#include <stdio.h>
#include "gammafit_search.h"
#include "read_array.h"

int
main (int argc, char **argv)
{
  int *spiketime = 0;
  int *lodif = 0;
  int spikecount;

  if (argc < 2) {
    printf ("usage: %s unit_id\n", argv[0]);
    return 0;
  }

  spikecount = read_array ("all", argv[1], &spiketime);

#ifndef NOPIN
  if (read_array ("lodif", argv[1], &lodif) != spikecount) {
    fprintf (stderr, "lodif%s must have the same number of entries as all%s\n", argv[1], argv[1]);
    return 1;
  }
#endif
    
  printf ("%g\n", gammafit_search (spiketime, lodif, spikecount));

  return 0;
}
