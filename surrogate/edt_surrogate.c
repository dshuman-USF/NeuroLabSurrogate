#define _GNU_SOURCE
#if HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <error.h>
#include <errno.h>
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

static void
parse_args (int argc, char **argv, bool *include, bool *replace, int *count, unsigned long *seed, bool *use_seed)
{
  int argn;
  unsigned long n;
  int state = START;
  
  for (argn = 0; argn < argc; argn++) {
    if (strcmp (argv[argn], "include") == 0) {
      memset (include, false, 1024);
      state = INCLUDE;
    }
    else if (strcmp (argv[argn], "exclude") == 0)
      state = EXCLUDE;
    else if (strcmp (argv[argn], "replace") == 0) {
      memset (replace, false, 1024);
      state = REPLACE;
    }
    else if (strcmp (argv[argn], "keep") == 0)
      state = KEEP;
    else if (strcmp (argv[argn], "count") == 0)
      state = COUNT;
    else if (strcmp (argv[argn], "seed") == 0)
      state = SEED;
    else {
      n = strtoul (argv[argn], NULL, 0);
      if (state == COUNT)
        *count = n;
      else if (state == INCLUDE)
        include[n] = true;
      else if (state == EXCLUDE)
        include[n] = false;
      else if (state == REPLACE)
        replace[n] = true;
      else if (state == KEEP)
        replace[n] = false;
      else if (state == SEED) {
        *use_seed = true;
        *seed = n;
      }
    }
  }
}

int
main (int argc, char **argv)
{
  printf ("%s\n", PACKAGE_STRING);
  
  Edt *edt, sur;
  bool include[1024], replace[1024];
  int n, surrogate_count = 1, suridx;
  double *g;
  bool use_seed = false;
  unsigned long seed = 0;

  if (argc < 2) {
    printf (
            "\nusage: %s whatever.edt [OPTIONS]\n"
            "\n"
            "Generates N files (in the current directory) named whatever_01.edt,\n"
            "whatever_02.edt, etc., with surrogate spike trains.  Will overwrite\n"
            "any existing files by those names.\n"
            "\n"
            "OPTIONS:\n"
            "include N ... : channels N ... from whatever.edt will be included in\n"
            "                the output.  All others will be excluded.\n"
            "exclude N ... : channels N ... from whatever.edt will be excluded from\n"
            "                the output.  All others will be included.\n"
            "keep N ...    : included digital channels N ... from whatever.edt will\n"
            "                be included as-is in the output file.  All other\n"
            "                included digital channels will be replaced with\n"
            "                surrogate spike trains.\n"
            "replace N ... : included digital channels N ... from whatever.edt will\n"
            "                be replaced with surrogate spike trains.  All other\n"
            "                included channels will be included as-is in the output\n"
            "                file.\n"
            "count N       : N is the number of output files.  Default is 1.\n"
            "seed N        : N is used as the seed of the random number generator.\n"
            "\n"
            "* If neither \"include\" nor \"exclude\" is specified, all channels will\n"
            "  be included.  Don't use both.\n"
            "* If neither \"keep\" nor \"replace\" is specified, all included digital\n"
            "  channels will be replaced.  Don't use both.\n"
            "* To specify an analog channel, add 1000 to the analog ID.\n"
            "* Included analog channels are always included as-is in the output files.\n"
            "* The input file can be .bdt or .edt.  Output files will be the same type.\n"
            "  If the input is .gdt, the outputs files will be shN.rdt\n"
            "\n"
            "EXAMPLE:\n"
            "\n"
            "%s whatever.edt include 101 97 98 1003 replace 101 count 20 seed 37\n"
            , argv[0], argv[0]);
    return 0;
  }

  memset (include, true, 1024);
  memset (replace, true, 1024);
  parse_args (argc - 2, argv + 2, include, replace, &surrogate_count, &seed, &use_seed);

  printf ("reading %s\n", argv[1]);
  edt = read_edt (argv[1], include);
  nodups (edt);

  memcpy (&sur, edt, sizeof sur);
  sur.digital = malloc (edt->digital_count * sizeof *sur.digital);
  for (n = 0; n < edt->digital_count; n++)
    if (replace[edt->digital[n]->id]){
      sur.digital[n] = malloc (sizeof (Digital));
      sur.digital[n]->id = edt->digital[n]->id;
    }
    else
      sur.digital[n] = edt->digital[n];

  char *basename, *extension, *format;

  if (edt->is_gdt) {
    basename = "sh";
    extension = "rdt";
    format = "%s%d.%s";
  }
  else {
    basename = strdup (argv[1]);
    extension = strrchr (basename, '.');
    char *p = strrchr (basename, '/');
    if (p) basename = p + 1;
    if (extension) extension[0] = 0, extension++;
    else           extension = "";
    format = "%s_%02d.%s";
  }

  if (1) { /* Pauluis and Baker surrogate */
    double *g = malloc (edt->digital_count * sizeof *g);
    RateFunc *ri = malloc (edt->digital_count * sizeof *ri);
    
    set_ticks_per_sec (strcmp (edt->format, "%5d%10d\n") == 0 ? 10000 : 2000);
    for (int n = 0; n < edt->digital_count; n++)
      if (replace[edt->digital[n]->id]) {
        SpikeTrain st = {.T = edt->digital[n]->spike, .N = edt->digital[n]->spike_count};
        printf ("getting firing rate for %d", edt->digital[n]->id); fflush (stdout);
        g[n] = miura_gamma (st);
        printf (" (gamma %g) (%d/%d)\n", g[n], n + 1, edt->digital_count);
        ri[n] = recip_isi (st, g[n]);
        if (0)
        {
          FILE *f = fopen ("rate", "w");
          for (int i = 0; i < ri[n].N; i++)
            fprintf (f, "%g\n", ri[n].r[i]);
          fclose (f);
          exit (0);
        }
      }
    printf ("generating surrogates and writing .edt files\n");
    for (int suridx = 0; suridx < surrogate_count; suridx++) {
      printf ("%d of %d\n", suridx + 1, surrogate_count);
      fflush(stdout);
      char *filename;
      if (asprintf (&filename, format, basename, suridx + 1, extension) == -1) exit (1);
      for (int n = 0; n < edt->digital_count; n++)
        if (replace[edt->digital[n]->id]) {
          SpikeTrain sst = spikes_from_rate (ri[n], lrint (g[n] < 1 ? 1.0 : g[n]));
          sur.digital[n]->spike = sst.T;
          sur.digital[n]->spike_count = sst.N;
          if (0) {
            printf ("index %d, surrogate spike count: %d\n", n, sst.N);
            FILE *f = fopen ("times", "w");
            for (int i = 0; i < sst.N; i++)
              fprintf (f, "%d\n", sst.T[i]);
            fclose (f);
          }
        }
      write_edt (filename, &sur);
      free (filename);
      for (n = 0; n < edt->digital_count; n++)
        if (replace[edt->digital[n]->id])
          free (sur.digital[n]->spike);
    }
    for (n = 0; n < edt->digital_count; n++)
      if (replace[edt->digital[n]->id]) {
        free (sur.digital[n]);
        free (ri[n].r);
      }
    free (ri);
    free (sur.digital);
    free (g);
    exit (0);
  }

  use_smoothed_rate = getenv ("SURROGATE_USE_SMOOTHED_RATE");
  Pvoid_t *t_est = 0;
  if (use_smoothed_rate)
    t_est = malloc (edt->digital_count * sizeof *t_est);

  g = malloc (edt->digital_count * sizeof *g);
  for (n = 0; n < edt->digital_count; n++)
    if (replace[edt->digital[n]->id]) {
      printf ("getting gamma for %d\n", edt->digital[n]->id);
      if (use_smoothed_rate) set_preserve_t_est (true);
      g[n] = get_gamma (edt->digital[n]->spike, edt->digital[n]->spike_count, &t_est[n]);
      //      g[n] = 0x1.d33ddfb71f784p+0;
      printf ("  %3d: %5.2f     \n", edt->digital[n]->id, g[n]);
    }

  set_use_smoothed_rate (use_smoothed_rate);
  for (suridx = 0; suridx < surrogate_count; suridx++) {
    printf ("generating surrogate %d\n", suridx + 1);
    for (n = 0; n < edt->digital_count; n++)
      if (replace[edt->digital[n]->id]) {
        if (use_smoothed_rate) put_t_est (t_est[n]);
        sur.digital[n]->spike = gen_control (edt->digital[n]->spike,
                                            edt->digital[n]->spike_count,
                                            g[n],
                                            &sur.digital[n]->spike_count,
                                            seed,
                                            use_seed);
        if (use_smoothed_rate) clear_t_est ();
      }
    char *filename;
    if (asprintf (&filename, format, basename, suridx + 1, extension) == -1) exit (1);
    write_edt (filename, &sur);
    free (filename);
    for (n = 0; n < edt->digital_count; n++)
      if (replace[edt->digital[n]->id])
        free (sur.digital[n]->spike);
  }

  return 0;
}
