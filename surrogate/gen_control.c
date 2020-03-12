/*
Copyright 2004-2020 Kendall F. Morris

This file is part of the Scope software suite.

    The Scope software suite is free software: you can redistribute
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

#define _GNU_SOURCE
#include <config.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <stdbool.h>
#include <Judy.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <unistd.h>
#include "gen_control.h"
#include "fpu_check.h"

#if SIZEOF_LONG == 8
typedef union {Word_t w; double f;} Wrd_Flt_t;
typedef union {PWord_t w; double *f; Pvoid_t v;} P_Wrd_Flt_t;
#else
typedef union {Word_t w; float f;} Wrd_Flt_t;
typedef union {PWord_t w; float *f; Pvoid_t v;} P_Wrd_Flt_t;
#endif

typedef union {Word_t w; double d;} Wrd_Dbl;
typedef union {Word_t w; float d;} Wrd_Flt;

static Pvoid_t cumspkcnt = (Pvoid_t) NULL;
static double *val;
static int val_idx;

void
insert (double index, double value)
{
  static int val_alloc;
  PWord_t  PValue;
#if SIZEOF_LONG == 8
  Wrd_Dbl i;
#else
  Wrd_Flt i;
#endif

  i.d = index;
  JLI (PValue, cumspkcnt, i.w);

  if (val_idx == val_alloc)
    val = realloc (val, (val_alloc += 1 << 19) * sizeof *val);

  val[val_idx] = value;
  *PValue = val_idx++;
}

double
interp1 (double y)
{
#if SIZEOF_LONG == 8
  Wrd_Dbl i;
#else
  Wrd_Flt i;
#endif
  PWord_t  PValue;
  double t0, t1;
  double y0, y1;

  i.d = y;
  JLL (PValue, cumspkcnt, i.w);
  y0 = i.d;
  t0 = val[*PValue];

  i.d = y;                      /* i.w gets overwritten by JLL */
  JLF (PValue, cumspkcnt, i.w);
  y1 = i.d;
  t1 = val[*PValue];

  if (y1 != y0)
    t0 += (y - y0) / (y1 - y0) * (t1 - t0);
  return t0;
}

int *
gen_control_from_rate (double end, double g, int *control_spikecount, unsigned long seed, bool use_seed)
{
  int alloc;
  int *control = malloc ((alloc = *control_spikecount) * sizeof *control);
  static gsl_rng *rng;
  double scale = 1 / g;

  if (rng == NULL) {
    rng = gsl_rng_alloc (gsl_rng_mt19937); 
    if (!use_seed) seed = time(0);
    gsl_rng_set (rng, seed);
    printf ("seed = %ld\n", seed);
  }
  
  double ysum = 0;
  int n = 0;
  while (true) {
    ysum += gsl_ran_gamma (rng, g, scale);
    if (ysum > end) {
      control = realloc (control, (*control_spikecount = n) * sizeof *control);
      return control;
    }
    if (n >= alloc)
      control = realloc (control, (alloc *= 2) * sizeof *control);
    control[n++] = floor (interp1 (ysum) + .5);
    if (n > 1 && control[n - 1] < control[n - 2]) {
      fprintf (stderr, "%s line %d: BUG\n", __FILE__, __LINE__);
      exit (1);
    }
  }
}

void
init_control (void)
{
  Word_t Rc_word;
  fpu_check ();
  JLFA (Rc_word, cumspkcnt);
  val_idx = 0;
}

static Pvoid_t t_est;
static bool use_smoothed_rate = false;

void
put_t_est (Pvoid_t t_est_val)
{
  t_est = t_est_val;
}

void
clear_t_est (void)
{
  Word_t Rc_word;
  JLFA (Rc_word, t_est);
}

void
zero_t_est (void)
{
  t_est = NULL;
}

void
set_use_smoothed_rate (bool val)
{
  use_smoothed_rate = val;
}

int *
gen_control (int *spiketime, int spikecount, double g, int *control_spikecount, unsigned long seed, bool use_seed)
{
  double cumspikecount;
  double t[4] = {0};
  double y[4] = {0};
  int tidx = 0, n;

  init_control ();

# define T(n) t[(tidx + n) % 4]
# define Y(n) y[(tidx + n) % 4]

  *control_spikecount = spikecount;

  if (use_smoothed_rate) {
    printf ("using smoothed rate\n");
    P_Wrd_Flt_t PValue;
    Wrd_Flt_t Index = {0};
    JLF(PValue.v, t_est, Index.w);
    double last = 0;
    while (PValue.v != NULL) {
#ifdef NOPIN
      insert (*PValue.f, Index.f);
      last = *PValue.f;
#else
      insert (*PValue.w + .5, Index.w);
      last = *PValue.w + .5;
#endif
      JLN(PValue.v, t_est, Index.w);
    }
    return gen_control_from_rate (last, g, control_spikecount, seed, use_seed);
  }

  cumspikecount = .5;
  for (n = 0; n < spikecount; n++) {
    T(0) = spiketime[n]; Y(0) = cumspikecount; tidx++;
    if (!use_smoothed_rate)
      insert (cumspikecount, spiketime[n]);
    if (10 * (T(1) - T(0)) < T(2) - T(1) && 10 * (T(3) - T(2)) < T(2) - T(1) && n >= 3) {
      double dy = .5 - 1e-6;
      insert (Y(1) + dy, T(1) + (T(1) - T(0)) * dy);
      insert (Y(2) - dy, T(2) - (T(3) - T(2)) * dy);
    }
    else if (n == 1)
      insert (0, T(2) - (T(3) - T(2)) / 2);
    cumspikecount++;
  }
  insert (Y(3) + .5, T(3) + (T(3) - T(2)) / 2);
  
  return gen_control_from_rate (Y(3)+.5, g, control_spikecount, seed, use_seed);
}

