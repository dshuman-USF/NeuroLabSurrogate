#define _GNU_SOURCE
#if HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>
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

#define TICKS_PER_SEC 10000

static double
rate (SpikeTrain st)
{
  return (st.N - 1) / ((double)(st.T[st.N - 1] - st.T[0]) / TICKS_PER_SEC);
}

int
main (int argc, char **argv)
{
  Edt *edt;
  bool include[1024];
  memset (include, true, 1024);

  edt = read_edt (argv[1], include);
  nodups (edt);

  double basegamma, mod1gamma, mod2gamma;
  double baserate , mod1rate , mod2rate ;
  printf ("%s\n", argv[1]);
  for (int dn = 0; dn < edt->digital_count; dn++) {
    int spike_count = edt->digital[dn]->spike_count;
    int *base = malloc (spike_count * sizeof (int)); int basecnt = 0;
    int *mod1 = malloc (spike_count * sizeof (int)); int mod1cnt = 0;
    int *mod2 = malloc (spike_count * sizeof (int)); int mod2cnt = 0;
    for (int n = 0; n < spike_count; n++) {
      int t = edt->digital[dn]->spike[n];
      int rep = t / (2 * TICKS_PER_SEC);
      int off = t % (2 * TICKS_PER_SEC);
      /* 0 350 650 1350 1650 2000 */
      /*     mod      mod         */

      if (off >= 3500 && off <= 6500) {
        assert (off >= 3500 && off <= 6500);
        mod1[mod1cnt++] = rep * 3000 + off - 3500;
      }
      else if (off >= 13500 && off <= 16500)
        mod2[mod2cnt++] = rep * 3000 + off;
      else if (off < 3500) {
        base[basecnt++] = rep * 14000 + off;
      }
      else if (off < 13500)
        base[basecnt++] = rep * 14000 + off - 3000;
      else
        base[basecnt++] = rep * 14000 + off - 6000;
      if (mod1[mod1cnt - 1] == mod1[mod1cnt - 2]) mod1cnt--;
      if (mod2[mod2cnt - 1] == mod2[mod2cnt - 2]) mod2cnt--;
      if (base[basecnt - 1] == base[basecnt - 2]) basecnt--;
    }
    SpikeTrain st;
    st.T = base; st.N = basecnt; basegamma = miura_gamma (st); baserate = rate (st);
    st.T = mod1; st.N = mod1cnt; mod1gamma = miura_gamma (st); mod1rate = rate (st);
    st.T = mod2; st.N = mod2cnt; mod2gamma = miura_gamma (st); mod2rate = rate (st);
    printf ("rate: %2.0f %2.0f %2.0f %2.0f %2.0f\n", baserate , mod1rate , baserate , mod2rate , baserate );
    printf ("gamm: %2.0f %2.0f %2.0f %2.0f %2.0f\n", basegamma, mod1gamma, basegamma, mod2gamma, basegamma);
    free (base);
    free (mod1);
    free (mod2);
  }
  

  return 0;
}
