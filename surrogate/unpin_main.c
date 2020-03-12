#include <stdio.h>
#include <stdlib.h>
#include "unpin.h"

int
main(void)
{
#ifdef NOPIN
  return 0;
#endif
  
  static int *spiketimes, *lodif;
  int alloc = 0, t;
  int n = 0, spikecount;

  while (scanf ("%d", &t) == 1) {
    if (n + 1 >= alloc)
      spiketimes = realloc (spiketimes, (alloc += 16384) * sizeof *spiketimes);
    spiketimes[n] = t;
    n++;
  }
  lodif = unpin (spiketimes, spikecount = n);

  for (n = 0; n < spikecount; n++)
    printf ("%d\n", lodif[n]);

  return 0;
}
