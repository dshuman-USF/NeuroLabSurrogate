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

#include <stdbool.h>
#include <Judy.h>

#ifdef __cplusplus
#define CDECL extern "C"
#else
#define CDECL
#endif

CDECL double gammafit_search (int *st, int *ld, int sc) ;
CDECL double gammafit (int *spiketime, int spikecount, Pvoid_t esct, double *gp);

CDECL void set_preserve_t_est (bool val);
CDECL Pvoid_t get_t_est (void);
CDECL Pvoid_t cth_t_est (int *st, int sc, double *E_begin, double *E_end, int ecnt);

#ifdef GMP
#include <gmp.h>
#include <gmpxx.h>
#define GETD(x) x.get_d()
#define INIT(x) new(&x) mpq_class
#else
#define mpq_class double
#define GETD(x) x
#define INIT(x) x = 0
#endif

//typedef struct
//{
//  mpq_class x;
//  mpq_class y;
//} Point;

typedef struct
{
  double x;
  double y;
} Point;

CDECL Point * approx (int *spiketimes, int spikecount, mpq_class dq, int *count) ;
CDECL void plot_best (int *spiketimes, int spikecount) ;
extern int count_min;
extern double d_min;
extern double g_min;
extern double d_est_min;
