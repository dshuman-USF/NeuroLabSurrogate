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

#define _GNU_SOURCE
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include <malloc.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_roots.h>
#include "tools.h"
#define MEM do {struct mallinfo m = mallinfo (); printf ("%s line %d: %d bytes malloc'ed (%d)\n", __FILE__, __LINE__, m.uordblks, m.hblkhd);} while (0)

static gsl_rng *rng;
static int ticks_per_sec;

static SpikeTrain
surrogate (RateFunc rf, int mult)
{
  int alloc = rf.N / 100 * mult;
  SpikeTrain st;
  st.T = malloc (alloc * sizeof *st.T);
  st.N = 0;
  for (int n = 0; n < rf.N; n++) {
    if (rf.r[n] == 0)
      continue;
    for (int i = 0; i < mult; i++) {
      double r = gsl_rng_uniform (rng);
      if (r <= rf.r[n]) {
        if (st.N >= alloc) st.T = realloc (st.T, (alloc *= 2) * sizeof *st.T);
        st.T[st.N++] = n;
      }
    }
  }

  /* realloc doesn't always shrink allocations! */
  /* st.T = realloc (st.T, st.N * sizeof *st.T); */
  int size = st.N * sizeof *st.T;
  int *tmp = malloc (size);
  memcpy (tmp, st.T, size);
  free (st.T);
  st.T = tmp;

  return st;
}

static SpikeTrain
decimate (SpikeTrain st, int mult)
{
  if (mult == 1)
    return st;
  int to = 0;
  for (int from = gsl_rng_uniform (rng) * mult; from < st.N; from += mult)
    st.T[to++] = st.T[from];
  st.T = realloc (st.T, to * sizeof *st.T);
  st.N = to;
  return st;
}

static void
init_rng (void)
{
  if (rng)
    return;
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
}

void
seed_rng (unsigned long seed)
{
  if (!rng)
    init_rng ();
  gsl_rng_set (rng, seed);
}

SpikeTrain
spikes_from_rate (RateFunc rf, int mult)
{
  if (!rng) init_rng ();
  SpikeTrain spikes = surrogate (rf, mult);
  spikes = decimate (spikes, mult);
  return spikes;
}

static RateFunc
get_normal_kernel (double s)
{
  int count = (int) (9 * s);
  double *k = malloc (count * sizeof *k);
  double ss2 = 2 * s * s;
  double ssqrt2pi = s * sqrt (2 * 3.14159265358979323846);
  int n;
  double sum = 0;
  double lastarea = 0;
  for (n = 0; n < count; n++) {
    k[n] = exp (-(n * n) / ss2) / ssqrt2pi;
    sum += k[n];
    double area = 2 * sum - k[0];
    if (area == lastarea)
      break;
    lastarea = area;
  }
  k = realloc (k, n * sizeof *k);
  return (RateFunc){.r = k, .N = n};
}

static void
convolve (RateFunc s, RateFunc k)
{
  //static int callcount;
  int Nfft = s.N + (k.N - 1) * 2;
  
  double *data = (double *) malloc (Nfft * sizeof *data);


  memcpy (data + k.N - 1, s.r, s.N * sizeof *data);
  int lfrom = k.N - 1;
  int lto = lfrom - 1;
  int rto = k.N - 1 + s.N;
  int rfrom = rto - 1;
  for (int n = 0; n < k.N - 1; n++) {
    data[lto--] = data[lfrom++];
    data[rto++] = data[rfrom--];
  }

  double *y = calloc (Nfft, sizeof *y);

  for (int n = 0; n < Nfft; n++) {
    y[n] += data[n] * k.r[0];
    for (int i = 1; i < k.N; i++) {
      double val =  data[n] * k.r[i];
      if (n + i < Nfft)
        y[n + i] += val;
      if (n - i >= 0)
        y[n - i] += val;
    }
  }
  memcpy (s.r, y + k.N - 1, s.N * sizeof *s.r);
  free (data);
  free (y);
}

void
set_ticks_per_sec (double tr)
{
  ticks_per_sec = tr;
}

static void
fill (double *y, int start, int end, double spike)
{
  for (int n = start; n < end; n++)
    y[n] = spike / (end - start);
}

static inline uint64_t 
tick (void)
{
  struct timespec ts;
  clock_gettime (CLOCK_MONOTONIC, &ts);
  return ts.tv_sec * 1000000000 + ts.tv_nsec;
}

RateFunc
recip_isi (SpikeTrain s, double g)
{
  if (0)
  {
    FILE *f = fopen ("t1", "w");
    for (int i = 0; i < s.N; i++)
      fprintf (f, "%d\n", s.T[i]);
    fclose (f);
    printf ("recip_isi: %d spike, g = %g\n", s.N, g);
    exit (0);
  }
  //  time_t start = time (0);
  
  int *T = s.T;
  int N = s.N;

  int tcnt = T[N - 1] + 1;

  double *y = calloc (tcnt, sizeof *y);

  fill (y, T[0], T[1], 1);
  fill (y, T[1], T[2], 1);
  int smooth_start = 0;
  static RateFunc krnl;
  int smooth = ticks_per_sec != 0;
  if (smooth && !krnl.r) krnl = get_normal_kernel (ticks_per_sec / 100.);
  {
    static int been_here;
    FILE *f = fopen ("DONT_LOWPASS_FILTER_SURROGATES", "r");
    if (f) {
      fclose (f);
      smooth = 0;
      if (!been_here) {
        f = fopen ("SURROGATES_GENERATED_WITHOUT_A_LOWPASS_FILTER", "w");
        if (f) fclose (f);
      }
    }
    been_here = 1;
  }
    
  double pval = .95;
  for (int n = 3; n < N - 3; n++) {
    //    int64_t t0 = tick ();
    int mid = T[n] - T[n - 1];
    
    double spike = 1;
    int start = T[n - 1];
    int end = T[n];

    double cdf1 = gsl_cdf_gamma_P (mid, g, (T[n - 2] - T[n - 3]) / g);
    int Ishort = T[n - 1] - T[n - 2];
    double cdf2 = gsl_cdf_gamma_P (mid, g, Ishort / g);
    if (0) {
      printf ("gamma %g, mean %d, long %d, cdf %g\n", g, T[n - 2] - T[n - 3], mid, cdf1);
      printf ("gamma %g, mean %d, long %d, cdf %g\n", g, Ishort, mid, cdf2);
    }
    int first = 0, second = 0;
    if (cdf1 > pval && cdf2 > pval) { /* decrease */
      double Thi = Ishort / 2.0;
      double Fhi = 1.0 / Ishort;
      for (int i = 0; i < Thi; i++)
        y[T[n - 1] + i] = Fhi;
      start += Thi;
      if (Thi - (int)Thi == .5) {
        y[T[n - 1] + (int)(Thi + .5)] = Fhi / 2;
        start++;
      }
      spike -= .5;
      //printf ("decrease at spike %d, %g seconds, smooth from %g\n", n - 1, (double)start/ticks_per_sec, (double)smooth_start/ticks_per_sec);
      first = start;
    }
    
    cdf1 = gsl_cdf_gamma_P (mid, g, (T[n + 2] - T[n + 1]) / g);
    Ishort = T[n + 1] - T[n];
    cdf2 = gsl_cdf_gamma_P (mid, g, Ishort / g);
    if (0) {
      printf ("gamma %g, mean %d, long %d, cdf %g\n", g, T[n + 2] - T[n + 1], mid, cdf1);
      printf ("gamma %g, mean %d, long %d, cdf %g\n", g, Ishort, mid, cdf2);
    }
    if (cdf1 > pval && cdf2 > pval) { /* increase */
      double Thi = Ishort / 2.0;
      double Fhi = 1.0 / Ishort;
      for (int i = 0; i < Thi; i++)
        y[T[n] - 1 - i] = Fhi;
      end -= Thi;
      if (Thi - (int)Thi == .5) {
        y[T[n - 1] - (int)(Thi + .5)] = Fhi / 2;
        end--;
      }
      spike -= .5;
      //printf ("increase at spike %d, %g seconds, smooth from %g\n", n, (double)end/ticks_per_sec, (double)smooth_start/ticks_per_sec);
      second = end;
    }
    fill (y, start, end, spike);
#define CONVOLVE convolve
    if (smooth && (first || second)) {
      CONVOLVE ((RateFunc){.r = y + smooth_start, .N = (first ? first : second) - smooth_start}, krnl);
      smooth_start = second ? second : first;
    }
    //    printf ("%d: %ld ns\n", n, tick () - t0);
  }
  fill (y, T[N - 3], T[N - 2], 1);
  fill (y, T[N - 2], T[N - 1], 1);
  if (smooth) CONVOLVE ((RateFunc){.r = y + smooth_start, .N = tcnt - smooth_start}, krnl);

  //  printf ("recip_isi: %ld seconds\n", time (0) - start);
  return (RateFunc){.r = y, .N = tcnt};
}

static double
get_SI (int *T, int N)
{
  
  double sum = 0;
  for (int i = 1; i < N - 1; i++) {
    int I0 = T[i] - T[i - 1];
    int I1 = T[i + 1] - T[i];
    sum += log (4.0 * I0 * I1 / pow (I0 + I1, 2));
  }
  return -sum / 2 / (N - 2);
}

static double
eqn (double kappa, void *params)
{
  double S_I = *(double *)params;
  double psi2 = gsl_sf_psi (2 * kappa);
  double psi1 = gsl_sf_psi (kappa);
  double retval = -S_I - log (2) + psi2 - psi1;
  //  printf ("kappa %g: %g, %g %g %g\n", kappa, retval, S_I, gsl_sf_psi (2 * kappa),gsl_sf_psi (kappa) );
  return retval;
}

static double
SI_to_gamma (double S_I)
{
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  double x_lo = 0.1, x_hi = 10000;
  gsl_function F;

  //  printf ("S_I: %g, eqn(x_lo): %g, eqn(x_hi): %g\n", S_I, eqn(x_lo,&S_I), eqn(x_hi,&S_I));

  F.function = &eqn;
  F.params = &S_I;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  do {
    iter++;
    gsl_root_fsolver_iterate (s);
    x_lo = gsl_root_fsolver_x_lower (s);
    x_hi = gsl_root_fsolver_x_upper (s);
    status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
  }
  while (status == GSL_CONTINUE && iter < max_iter);
  r = gsl_root_fsolver_root (s);

  gsl_root_fsolver_free (s);

  return r;
}

double
miura_gamma (SpikeTrain st)
{
  if (st.N < 3)
    return 1;
  return SI_to_gamma (get_SI (st.T, st.N));
}
