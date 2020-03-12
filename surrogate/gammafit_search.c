/*
Copyright 2005-2020 Kendall F. Morris

This file is part of the Xanalysis software suite.

    The Xanalysis software suite is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public
    License as published by the Free Software Foundation, either
    version 3 of the License, or (at your option) any later version.

    The suite is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with the suite.  If not, see <https://www.gnu.org/licenses/>.
*/


#ifndef NOPIN
#define _GNU_SOURCE
#if HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <float.h>
#include <stdbool.h>
#include <time.h>
#include <Judy.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <unistd.h>
#include "gammafit_search.h"
#include "fpu_check.h"

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

typedef union {Word_t w; double d;} Wrd_Dbl;
typedef union {Word_t w; float d;} Wrd_Flt;

typedef struct
{
  double g;
  double d;
} SearchResult;

static int *spiketime;
static int spikecount;
static int *lodif;

static Pvoid_t result = (Pvoid_t) NULL;
static Pvoid_t t_est = (Pvoid_t) NULL;
static Pvoid_t isi = (Pvoid_t) NULL;

double d_min = DBL_MAX;
double g_min;
int s_min;

static void
add_result (int segcnt, double d, double g)
{
  PWord_t  PValue;                          // pointer to return value
  SearchResult *r = malloc (sizeof *r);
  r->g = g;
  r->d = d;
  JLI (PValue, result, segcnt);
  *PValue = (Word_t) r;
  if (d < d_min)
    d_min = d, g_min = g, s_min = segcnt;
}

static inline void
add_t_est (int spiketime, int y)
{
  PWord_t  PValue;                          // pointer to return value
  Word_t   Index;         // array index
  Index = spiketime;
  JLI (PValue, t_est, Index);
  *PValue = y;
}

static double
interp1 (int spiketime)
{
  PWord_t  PValue;
  Word_t   Index;
  int t0, t1;
  double y0, y1;

  Index = spiketime;
  JLL (PValue, t_est, Index);
  t0 = Index;
  y0 = *PValue + .5;


  Index = spiketime;
  JLF (PValue, t_est, Index);
  t1 = Index;
  y1 = *PValue + .5;
  
  if (t1 != t0)
    y0 += (double)(spiketime - t0) / (t1 - t0) * (y1 - y0);
  return y0;
}

static int current_segcnt;

#if 0
static void
write_t_est (void)
{
  Word_t   Index;                     // array index
  Word_t * PValue;                    // pointer to array element value
  FILE *f;

  f = fopen ("t_est", "w");

  Index = 0;
  JLF (PValue, t_est, Index);
  while (PValue != NULL) {
    fprintf (f, "%10ld %10ld\n", Index, *PValue);
    JLN (PValue, t_est, Index);
  }
  fclose (f);
}
#endif

static gsl_error_handler_t *default_gsl_err_handler;
static int cdf_gamma_P_error_count;

static void
gsl_cdf_gamma_P_err (const char * reason, const char * file, int line, int gsl_errno)
{
  if (gsl_errno == GSL_EMAXITER) {
    if (cdf_gamma_P_error_count == 0)
      fprintf (stderr, "gsl_cdf_gamma_P: %s\n", gsl_strerror (gsl_errno));
    cdf_gamma_P_error_count++;
    return;
  }
  default_gsl_err_handler (reason, file, line, gsl_errno);
}

static void
goto_segcnt (int target_segcnt)
{
  Word_t   Index;         // array index
  int      Rc_int;                          // return code - integer

  while (current_segcnt < target_segcnt) {
    int spikeidx = spikecount - 1 - (current_segcnt - 1);
    add_t_est (spiketime[lodif[spikeidx]-2], lodif[spikeidx]-2);
    current_segcnt++;
  }

  while (current_segcnt > target_segcnt) {
    int spikeidx = spikecount - (current_segcnt - 1);
    Index = spiketime[lodif[spikeidx]-2];
    if (Index == 0) {
      fprintf (stderr, "%s line%d: fatal error: removing spiketime 0, spikeidx %d, lodif %d\n",
               __FILE__, __LINE__, spikeidx, lodif[spikeidx]);
      exit (1);
    }
    JLD (Rc_int, t_est, Index);
    current_segcnt--;
  }
}

static double
unpin_gammafit (int target_segcnt, double *gp)
{
  PWord_t  PValue;
  Word_t   Index;         // array index
  int n, N;
  double prev_y, sum = 0, g, scale, max_d;
#if SIZEOF_LONG == 8
  Wrd_Dbl i;
#else
  Wrd_Flt i;
#endif
  Word_t Rc_word;

  //printf ("current: %d, target: %d, spikecount: %d\n", current_segcnt, target_segcnt, spikecount);

  goto_segcnt (target_segcnt);
  
  prev_y = interp1 (spiketime[0]);
  JLFA(Rc_word, isi);
  double sumsq = 0;
  for (n = 1; n < spikecount; n++) {
    double y = interp1 (spiketime[n]);
    double isi_val = y - prev_y;
    //printf ("%g\n", isi_val);
    prev_y = y;
    sum += isi_val;
    sumsq += (isi_val - 1) * (isi_val - 1);
    if (isi_val <= 0) {
      fprintf (stderr, "%s line %d: isi_val: %g, y: %g, prev_y: %g, n: %d, spiketime[n]: %d\n",
               __FILE__, __LINE__, isi_val, y, prev_y, n, spiketime[n]);
      //write_t_est ();
      exit (1);
    }
    i.d = isi_val;
    JLI (PValue, isi, i.w);
    (*PValue)++;
  }
  N = spikecount - 1;
  
  g = sumsq * N == 0 ? 1e6 : N / sumsq;
  scale = 1 / g;

  max_d = 0;

  Index = 0;
  n = 0;
  JLF(PValue, isi, Index);
  cdf_gamma_P_error_count = 0;
  default_gsl_err_handler = gsl_set_error_handler (gsl_cdf_gamma_P_err);
  while (PValue != NULL) {
    double actual_hi = (double)(n + *PValue) / N;
    double actual_lo = (double) n / N;
    double cdf, d;
    i.w = Index;
    //printf ("%g %g %g\n", i.d, g, scale);
    cdf = gsl_cdf_gamma_P (i.d, g, scale);
    d = MAX (actual_hi - cdf, cdf - actual_lo);
    if (d > max_d)
      max_d = d;
    n += *PValue;
    JLN(PValue, isi, Index);
  }
  gsl_set_error_handler (default_gsl_err_handler);

  *gp = g;
  return max_d;
}

static double
get_result (int segcnt)
{
  double d, g;
  d = unpin_gammafit (segcnt, &g);
  add_result (segcnt, d, g);
  //printf ("%6d %9.7f %9.5f\n", segcnt, d, g);
  return d;
}

static time_t last_time;

static bool
get_next_segcnt (int *segcntp)
{
  Word_t   Index, prev_Index;         // array index
  Word_t * PValue;                    // pointer to array element value
  double threshold_delta = .00001;
  double threshold = d_min + threshold_delta;
  SearchResult *r, *prev_r;
  int max_gap = 0;
  int first = 0, gapsum = 0;

  Index = 0;
  JLF (PValue, result, Index);
  r = (SearchResult *) *PValue;
  prev_Index = Index;
  JLN (PValue, result, Index);
  while (PValue != NULL) {
    prev_r = r;
    r = (SearchResult *) *PValue;
    if (prev_r->d <= threshold || r->d <= threshold) {
      int gap = Index - prev_Index;
      gapsum += gap - 1;
      if (first == 0)
        first = prev_Index;
      //printf ("gap %7ld to %7ld (%d), do %d next\n", prev_Index, Index, gap, *segcntp);
      if (gap > max_gap) {
        max_gap = gap;
        *segcntp = (Index + prev_Index) / 2;
        if (*segcntp > spikecount + 1) {
          printf ("%s line %d: BUG: segcnt too large\n", __FILE__, __LINE__);
          printf ("Prev_Index %ld segcnt %d Index %ld spikecount %d\n",
                  prev_Index, *segcntp, Index, spikecount);
          
          exit (1);
        }
      }
    }
    prev_Index = Index;
    JLN(PValue, result, Index);
  }
  time_t now;
  if ((now = time (0)) > last_time) {
    Word_t Rc_word;
    JLC(Rc_word, result, 0, -1);
    last_time = now;
  }

  return max_gap > 1;
}

static void
write_best (void)
{
  goto_segcnt (s_min);

  PWord_t  PValue;
  Word_t   Index = 0;
  FILE *f = fopen ("rate.tmp", "w");

  double tick = atof (getenv ("EDT_SURROGATE_RATE"));

  Index = 0;
  JLF(PValue, t_est, Index);
  while (PValue != NULL) {
    fprintf (f, "%g %g\n", Index * tick, *PValue + .5);
    JLN(PValue, t_est, Index);
  }
  fclose (f);
}

void
plot_best (int *spiketimes, int spikecount)
{
  goto_segcnt (s_min);

  PWord_t  PValue;
  Word_t   Index = 0;


  int n;
  double max_d = 0;
  int max_d_loc = 0;
  for (n = 0; n < spikecount; n++) {
    double y = interp1 (spiketimes[n]);
    double d = fabs (y - (n + .5));
    if (d > max_d) {
      max_d = d;
      max_d_loc = spiketimes[n];
    }
  }
  max_d += .5;
  
#define FNAME  "plot.scr"
  FILE *f = fopen (FNAME, "w");

  fprintf (f, "1;\n");

  fprintf (f, "sx = [\n");
  for (n = 0; n < spikecount; n++)
    fprintf (f, "%d\n%d\n", spiketimes[n], spiketimes[n]);
  fprintf (f, "];\n");

  fprintf (f, "sy = [\n");
  for (n = 0; n < spikecount; n++)
    fprintf (f, "%d\n%d\n", n, n + 1);
  fprintf (f, "];\n");

  fprintf (f, "x = [\n");
  Index = 0;
  JLF(PValue, t_est, Index);
  while (PValue != NULL) {
    fprintf (f, "%lu\n", Index);
    JLN(PValue, t_est, Index);
  }
  fprintf (f, "];\n");

  fprintf (f, "y = [\n");
  Index = 0;
  JLF(PValue, t_est, Index);
  while (PValue != NULL) {
    fprintf (f, "%g\n", *PValue + .5);
    JLN(PValue, t_est, Index);
  }
  fprintf (f, "];\n");
  
  fprintf (f, "d = %g\n",  max_d);

  fprintf (f,
           "hold off\n"
           "plot (sx, sy + d);\n"
           "hold on\n"
           "plot (sx, sy - d);\n"
           "plot (x, y);\n"
           //           "input ('Press Enter when ready');\n"
           );

  fclose (f);
  
  //  system ("octave " FNAME );
  Word_t Rc_word;
  JLC (Rc_word, t_est, 0, -1);
  printf ("%lu points, d_max = %g at %d\n", Rc_word, max_d, max_d_loc);
  fflush (stdout);
}

static bool preserve_t_test = false;

void
set_preserve_t_est (bool val)
{
  preserve_t_test = val;
}
 
Pvoid_t
get_t_est (void)
{
  goto_segcnt (s_min);
  Pvoid_t retval = t_est;
  t_est = NULL;
  return retval;
}

double
gammafit_search (int *st, int *ld, int sc)
{
  int segcnt = 0;

  fpu_check ();

  d_min = DBL_MAX;

  spiketime = st;
  lodif = ld;
  spikecount = sc;

  Word_t Rc_word;
  JLC(Rc_word, t_est, 0, -1);
  if (Rc_word != 0) {
    fprintf (stderr, "%s line %d: t_est is not empty\n", __FILE__, __LINE__);
    exit (1);
  }
  
  add_t_est (0, 0);
  add_t_est (spiketime[spikecount - 1] + 1, spikecount);
  current_segcnt = 1;

  add_result (0, DBL_MAX, 0);
  if (spikecount < 3) {
    get_result (1);
    get_result (2);
  }
  else {
    get_result (spikecount / 3);
    get_result (2 * spikecount / 3);
  }
  add_result (spikecount + 2, DBL_MAX, 0);

  last_time = time(0);
  while (get_next_segcnt (&segcnt))
    get_result (segcnt);

  //  printf ("%s line %d\n", __FILE__, __LINE__);
  if (getenv("EDT_SURROGATE_PLOT")) {
    printf ("%s line %d\n", __FILE__, __LINE__);
    plot_best (st, sc);
  }

  if (getenv("EDT_SURROGATE_RATE")) {
    printf ("%s line %d: saving rate data\n", __FILE__, __LINE__);
    write_best ();
  }

  if (!preserve_t_test)
    JLFA(Rc_word, t_est);
  JLFA(Rc_word, result);
  JLFA(Rc_word, isi);

  return g_min;
}
#endif
