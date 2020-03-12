#define _GNU_SOURCE
#if HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "gammafit_search.h"
#include "unpin.h"
#include "gen_control.h"
#include "edt.h"
#include "tools.h"
#if !HAVE_ASPRINTF && HAVE_CONFIG_H
int asprintf(char **buffer, char *fmt, ...);
#endif

bool use_smoothed_rate;
        
static inline double
get_gamma (int *spiketimes, int spikecount, Pvoid_t *t_est)
{
  int *lodif = 0;
  double g;

  lodif = unpin (spiketimes, spikecount);
  g = gammafit_search (spiketimes, lodif, spikecount);
  if (use_smoothed_rate) *t_est = get_t_est ();
  free (lodif);
  return g;
}

enum {START, INCLUDE, EXCLUDE, REPLACE, KEEP, COUNT, SEED};

int
main (int argc, char **argv)
{
  Edt *edt;
  bool include[1024];
  memset (include, true, 1024);

  edt = read_edt (argv[1], include);
  nodups (edt);

  for (int n = 0; n < edt->digital_count; n++) {
    SpikeTrain st = {.T = edt->digital[n]->spike, .N = edt->digital[n]->spike_count};
    printf ("%d: %g\n", edt->digital[n]->id, miura_gamma (st));
  }

  return 0;
}
