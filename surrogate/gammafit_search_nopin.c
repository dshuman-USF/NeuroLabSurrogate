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

#ifdef NOPIN
#define _GNU_SOURCE 1
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <fnmatch.h>
#include <time.h>
#include <float.h>
#include <Judy.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <config.h>
#include "util.h"
#include "fpu_check.h"
#include "gammafit_search.h"

#if SIZEOF_LONG == 8
typedef union {Word_t w; double f;} Wrd_Flt_t;
typedef union {PWord_t w; double *f; Pvoid_t v;} P_Wrd_Flt_t;
#else
typedef union {Word_t w; float f;} Wrd_Flt_t;
typedef union {PWord_t w; float *f; Pvoid_t v;} P_Wrd_Flt_t;
#endif


static double
interp (Pvoid_t ary, int spiketime)
{
  P_Wrd_Flt_t PValue;
  Wrd_Flt_t Index;
  double t0, t1;
  double y0, y1;

  Index.f = spiketime;
  JLL (PValue.v, ary, Index.w);
  if (PValue.w == 0) {
    printf ("lookup of spiketime %d failed\n", spiketime);
    exit (1);
  }
  t0 = Index.f;
  y0 = *PValue.f;

  Index.f = spiketime;
  JLF (PValue.v, ary, Index.w);
  t1 = Index.f;
  y1 = *PValue.f;
  
  if (t1 != t0)
    y0 += (double)(spiketime - t0) / (t1 - t0) * (y1 - y0);
  return y0;
}

static gsl_error_handler_t *default_gsl_err_handler;
static int cdf_gamma_P_error_count;

static void
gsl_cdf_gamma_P_err (const char * reason, const char * file, int line, int gsl_errno)
{
  if (gsl_errno == GSL_EMAXITER) {
    if (cdf_gamma_P_error_count == 0 && getenv ("SUPPRESS_GSL_ERRORS") == NULL)
      fprintf (stderr, "gsl_cdf_gamma_P: %s\n", gsl_strerror (gsl_errno));
    cdf_gamma_P_error_count++;
    return;
  }
  default_gsl_err_handler (reason, file, line, gsl_errno);
}

double
gammafit (int *spiketime, int spikecount, Pvoid_t esct, double *gp)
{
  Pvoid_t isi  = (Pvoid_t) NULL;
  int n, N;
  
  double prev_y, sum = 0, sumsq = 0, g, scale, max_d;
  P_Wrd_Flt_t  PValue;
  Word_t Rc_word;
  Wrd_Flt_t i;

  prev_y = interp (esct, spiketime[0]);
  for (n = 1; n < spikecount; n++) {
    double y = interp (esct, spiketime[n]);
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
    i.f = isi_val;
    JLI (PValue.v, isi, i.w);
    (*PValue.w)++;
  }
  N = spikecount - 1;
  
  g = sumsq * N == 0 ? 1e6 : N / sumsq;
  scale = 1 / g;

  max_d = 0;

  Word_t Index = 0;
  n = 0;
  JLF(PValue.v, isi, Index);
  cdf_gamma_P_error_count = 0;
  default_gsl_err_handler = gsl_set_error_handler (gsl_cdf_gamma_P_err);
  while (PValue.v != NULL) {
    double actual_hi = (double)(n + *PValue.w) / N;
    double actual_lo = (double) n / N;
    double cdf, d;
    i.w = Index;
    //printf ("%g %g %g\n", i.f, g, scale);
    cdf = gsl_cdf_gamma_P (i.f, g, scale);
    d = MAX (actual_hi - cdf, cdf - actual_lo);
    if (d > max_d)
      max_d = d;
    n += *PValue.w;
    JLN(PValue.v, isi, Index);
  }
  gsl_set_error_handler (default_gsl_err_handler);
  JLFA(Rc_word, isi);

  *gp = g;
  return max_d;
}

typedef struct
{
  int *top_hull;
  int *bot_hull;
  int top_count;
  int bot_count;
  int low_slope_left;
  int low_slope_right;
  int high_slope_left;
  int high_slope_right;
  int endidx;
} Hulls;

typedef struct 
{
  Hulls hulls;
  mpq_class dqc;
  double sum;
  double sumsq;
  int count;
  int longest;
  int boundary_count;
  double m_max;
  double b_max;
  double d;
  int left;
  int right;
  int pointcount;
  Point *point;
  int spikeidx;
  int spikecount;
  int *spiketimes;
  int intersect_constraint_rows[3];
  bool intersect_constraint;
  bool success;
} Data;

static mpq_class ypch, slope, yo, yl, ya, yb, yc, offset[2];
static int *hull[2];
static int *st;
enum {BOT, TOP};

void
dprinte (char *vn, double v)
{
  printf ("%s = %a\n", vn, GETD(v));
}
#define PA(x) dprinte(#x, GETD(x))

static void
printc (char *vn, mpq_class v)
{
  printf ("%s = %g\n", vn, GETD(v));
}
#define PC(x) printc(#x, x)

static inline mpq_class
ydif (int hc, int c, int ha, int a, int hb, int b)
{
  int xa = st[hull[ha][a]];
  int xb = st[hull[hb][b]];
  int xc = st[hull[hc][c]];

  ya = hull[ha][a] + offset[ha];
  yb = hull[hb][b] + offset[hb];
  yc = hull[hc][c] + offset[hc];

  yl = ya + (yb - ya) / (xb - xa) * (xc - xa);
  return yc - yl;
}

static int chull_count;

static Hulls
chull (int *spiketimes, int startidx, int endidx, mpq_class dq)
{
  static Hulls hulls;
  st = spiketimes;
  int count = endidx - startidx + 1;
  if (count < 2) {
    hulls.endidx = endidx;
    return hulls;
  }
  static int *top_hull;
  static int *bot_hull;
  static int alloc;
  if (alloc < count) {
    TREALLOC (top_hull, count);
    TREALLOC (bot_hull, count);
    alloc = count;
  }
  top_hull[0] = bot_hull[0] = startidx;
  top_hull[1] = bot_hull[1] = startidx + 1;

  static mpq_class mq, dt;
  static bool first_time = true;
  if (first_time) {
    offset[BOT] =  1  - dq;
    offset[TOP] =  dq;
    hull[BOT] = bot_hull;
    hull[TOP] = top_hull;
  }
  
  int spikeidx;
  int top_prev = 1;
  int bot_prev = 1;
#define P(x) printf (#x " = %" "d" "\n", x)
#define PG(x) printf (#x " = %" "g" "\n", GETD(x))
  //#define PA(x) printf (#x " = %" "a" "\n", x)
#define Pf(x,f) printf (#x " = %" #f "\n", x)
  for (spikeidx = startidx + 2; spikeidx <= endidx; spikeidx++) {
    chull_count++;
    int top_prev_save = top_prev;
    int bot_prev_save = bot_prev;
    while (top_prev > 0
           && ((long long)(st[spikeidx] - st[top_hull[top_prev]])
               * (top_hull[top_prev] - top_hull[top_prev - 1])
               >= (st[top_hull[top_prev]] - st[top_hull[top_prev - 1]])
               * (long long)(spikeidx - top_hull[top_prev])
               ))
      top_prev--;
    while (bot_prev > 0
           && ((long long)(st[spikeidx] - st[bot_hull[bot_prev]])
               * (bot_hull[bot_prev] - bot_hull[bot_prev - 1])
               <= (st[bot_hull[bot_prev]] - st[bot_hull[bot_prev - 1]])
               * (long long)(spikeidx - bot_hull[bot_prev])
               ))
      bot_prev--;

    int top_idx, bot_idx;
    bool intersect = false;

    //point is at (st[top_hull[top_idx]], top_hull[top_idx] + d)
    //line is from (st[bot_hull[bot_prev]], bot_hull[bot_prev] + 1 - d)
    //          to (st[spikeidx], spikeidx + 1 - d)
    // intersect if point is below line
    if (st[top_hull[top_prev]] >= st[bot_hull[bot_prev]]) {
      slope = (spikeidx - bot_hull[bot_prev]);
      slope /= (st[spikeidx] - st[bot_hull[bot_prev]]);
      yo = bot_hull[bot_prev] + 1 - dq;
    }
    for (top_idx = top_prev;
         top_idx >= 0 && st[top_hull[top_idx]] >= st[bot_hull[bot_prev]] && !intersect;
         top_idx--) {
      mq = top_hull[top_idx];
      ypch = mq + dq;
      dt = st[top_hull[top_idx]] - st[bot_hull[bot_prev]];
      yl = slope * dt + yo;
      if (ypch - yl < 0)
        intersect = true;
    }

    //point is at (st[bot_hull[bot_idx]], bot_hull[bot_idx] + 1 - d)
    //line is from (st[top_hull[top_prev]], top_hull[top_prev] + d)
    //          to (st[spikeidx], spikeidx + d)
    // intersect fi point is above line
    if (st[bot_hull[bot_prev]] >= st[top_hull[top_prev]] && !intersect) {
      slope = (spikeidx - top_hull[top_prev]);
      slope /= (st[spikeidx] - st[top_hull[top_prev]]);
      yo = top_hull[top_prev] + dq;
    }
    for (bot_idx = bot_prev;
         bot_idx >= 0 && st[bot_hull[bot_idx]] >= st[top_hull[top_prev]] && !intersect;
         bot_idx--) {
      mq = bot_hull[bot_idx] + 1;
      ypch = mq - dq;
      dt = st[bot_hull[bot_idx]] - st[top_hull[top_prev]];
      yl = slope * dt + yo;
      
      if (ypch - yl > 0)
        intersect = true;
    }
    if (intersect) {
      top_prev = top_prev_save;
      bot_prev = bot_prev_save;
      break;
    }
    top_hull[++top_prev] = spikeidx;
    bot_hull[++bot_prev] = spikeidx;
  }
  hulls.top_hull = top_hull;
  hulls.bot_hull = bot_hull;
  hulls.top_count = top_prev + 1;
  hulls.bot_count = bot_prev + 1;
  hulls.endidx = top_hull[top_prev];

  int b = 0;
  int t = top_prev;
  int b0 = b + 1;
  int t0 = 0;
  
  while (b != b0 || t != t0) {
    b0 = b; t0 = t;
    for ( ; spiketimes[bot_hull[b + 1]] < spiketimes[top_hull[t]]; b++)
      if (ydif (BOT, b + 1, BOT, b, TOP, t) <= 0)
        break;
    for ( ; spiketimes[bot_hull[b]] < spiketimes[top_hull[t - 1]]; t--)
      if (ydif (TOP, t - 1, BOT, b, TOP, t) >= 0)
        break;
  }
  hulls.high_slope_left  = bot_hull[b];
  hulls.high_slope_right = top_hull[t];
  
  b = bot_prev;
  t = 0;
  t0 = t + 1;
  while (b != b0 || t != t0) {
    b0 = b; t0 = t;
    for ( ; spiketimes[top_hull[t + 1]] < spiketimes[bot_hull[b]]; t++)
      if (ydif (TOP, t + 1, TOP, t, BOT, b) >= 0)
        break;
    for ( ;  spiketimes[top_hull[t]] < spiketimes[bot_hull[b - 1]]; b--)
      if (ydif (BOT, b - 1, TOP, t, BOT, b) <= 0)
        break;
  }
  hulls.low_slope_left  = top_hull[t];
  hulls.low_slope_right = bot_hull[b];

  return hulls;
}

/*
static void
plot (Data *p, double m_max, double b_max, double m_min, double b_min)
{
  #define FNAME  "plot.scr"
  ofstream f(FNAME, ios::out);
  int n;

  int stir = p->hulls.top_hull[p->hulls.top_count - 1] + 1;
  int stil = stir;
  Point *p0 = &p->point[p->pointcount - 2];
  Point *p1 = &p->point[p->pointcount - 1];
  while (p->spiketimes[stil] >= p0->x)
    stil--;

  f << "1;\n";

  f << "sx = [\n";
  for (n = stil; n <= stir; n++)
    f << p->spiketimes[n] << ",\n" <<p->spiketimes[n] << "\n";
  f << "];\n";

  f << "sy = [\n";
  for (n = stil; n <= stir; n++)
    f << n << ",\n" << n + 1 << "\n";
  f <<"];\n";

  f << "px = [" << p0->x << ", " << p1->x << "];\n";
  f << "py = [" << p0->y << ", " << p1->y << "];\n";
  
  f << "nx = [" << p0->x << ", " <<  p->spiketimes[stir] << "];\n";

  f << "nmaxy = [" <<  m_max*p0->x+b_max << ", " <<  m_max*p->spiketimes[stir]+b_max << "];\n";
  f << "nminy = [" <<  m_min*p0->x+b_min << ", " <<  m_min*p->spiketimes[stir]+b_min << "];\n";

  f <<
  "d = " << p->d << ";\n"
  "hold off\n"
  "plot (sx, sy + d);\n"
  "hold on\n"
  "plot (sx, sy - d);\n"
  "plot (px, py);\n"
  "plot (nx, nmaxy);\n"
  "plot (nx, nminy);\n"
  "input (\"Press Enter when ready\");\n";

  f.close ();

  //  system ("octave " FNAME "&");
  exit (0);
}
*/

static mpq_class yhl, yhr, b_max, m_max, yll, ylr, m_min, b_min, xi, yi;
static mpq_class m0, b0, m, b;
static Point ptmp;

static int back_count;

int
getseg_n (Data *p)
{
  static int passcnt;
  passcnt++;
  int endidx = p->hulls.endidx;
  if (p->spikeidx == -1) printf ("\ngetseg_n: %d to %d\n", p->spikeidx, endidx);
  if (p->spikeidx == p->spikecount - 1 || p->spikeidx == endidx) {
    xi = p->spiketimes[p->spikeidx];
    yi = p->spikeidx + .5;
    m_min = 0;
    m_max = p->spikecount;
  }
  else {
    int xhl = p->spiketimes[p->hulls.high_slope_left];
    int xhr = p->spiketimes[p->hulls.high_slope_right];
    yhl = p->hulls.high_slope_left + 1 - p->dqc;
    yhr = p->hulls.high_slope_right + p->dqc;
    m_max = (yhr - yhl) / (xhr - xhl);
    b_max = yhl - m_max * xhl;

    int xll = p->spiketimes[p->hulls.low_slope_left];
    int xlr = p->spiketimes[p->hulls.low_slope_right];
    yll = p->hulls.low_slope_left + p->dqc;
    ylr = p->hulls.low_slope_right + 1 - p->dqc;
    m_min = (ylr - yll) / (xlr - xll);
    b_min = yll - m_min * xll;
    if (m_min == m_max) {
      xi = xhl;
      yi = yhl;
    }
    else {
      xi = -(b_min-b_max)/(m_min-m_max);
      yi = (b_max*m_min-b_min*m_max)/(m_min-m_max);
    }
    if (m_min < 0) { m_min = 0; b_min = yi; }
    if (m_max <= 0) {
      P(p->hulls.high_slope_left);
      P(p->hulls.high_slope_right);
      PC(p->dqc);
      PC(yhl);
      PC(yhr);
    }
    m_max > 0 || DIE;
  }
  static mpq_class x0, yt, yb, dx, xt, xb, slope, xl, xr, x, y;
  x0 = p->pointcount >= 2 ? p->point[p->pointcount - 2].x : 0;
  if (p->spikeidx == -1) {
    P(p->pointcount);
    PC(p->point[p->pointcount - 2].x);
    PC(x0);
  }
  ((p->spikeidx == 0 && p->pointcount == 0)
   || (p->spikeidx >= 2 && p->pointcount >= 2)) || DIE;
  int n;
  for (n = p->spikeidx - 1; n > 0 && p->spiketimes[n - 1] > x0; n--) {
    back_count++;
    dx = xi - p->spiketimes[n];
    yt = n + p->dqc;
    yb = n + 1 - p->dqc;
    slope = (yi - yt) / dx;
    if (slope > m_min) {
      if (slope > m_max)
        break;
      m_min = slope;
    }
    slope = (yi - yb) / dx;
    if (slope < m_max) {
      if (slope < m_min)
        break;
      m_max = slope;
    }
  }
  static Point *p0, *p1, pm;
  if (p->pointcount < 2) {
    p->pointcount == 0 || DIE;
    TREALLOC (p->point, ++p->pointcount);
    INIT (p->point[p->pointcount - 1].x);
    INIT (p->point[p->pointcount - 1].y);

    p0 = &ptmp;
    p1 = &p->point[0];
    p0->x = -p->spiketimes[p->spikecount - 1];
    p0->y = 0;
    p1->x = p->spiketimes[0];
    p1->y = 0;
    pm.x = -p->spiketimes[p->spikecount - 1];
  }
  else {
    p0 = &p->point[p->pointcount - 2];
    p1 = &p->point[p->pointcount - 1];
    pm.x = p->spiketimes[n];
  }
  m0 = (p1->y - p0->y) / (p1->x - p0->x);
  b0 = p0->y - m0 * p0->x;
  if (p->spikeidx == -1) {PC(m0); PC(b0);}
 

  pm.y = m0 * pm.x + b0;
  static mpq_class ma, mb, tmp;
  ma = (yi - pm.y ) / (xi - pm.x );
  if (xi == p1->x)
    mb = p->spikecount;
  else
    mb = (yi - p1->y) / (xi - p1->x);
  if (ma > mb) {tmp = ma; ma = mb; mb = tmp;}
  if (mb < m_min || ma > m_max) {
    if (p->pointcount == 1)
      p->pointcount = 0;
    endidx > p->spikeidx || DIE;
    return -1;
  }
  if (ma < m_min) ma = m_min;
  if (mb > m_max) mb = m_max;
  if (p->spikeidx == -1) {
    PG(m_min);
    PG(m_max);
  }
  m_min = ma;
  m_max = mb;
  if (m0 == m_min || m0 == m_max) {
    PG(m0);
    PG(b0);
    PA(m0);
    PA(b0);
    PG(m_min);
    PA(m_min);
    PG(m_max);
    PA(m_max);
    PG(xi);
    PG(yi);
    PG(pm.x);
    PG(pm.y);
    PG(p1->x);
    PG(p1->y);
    PG(p0->x);
    PG(p0->y);
    PA(xi);
    PA(yi);
    PA(pm.x);
    PA(pm.y);
    PA(p1->x);
    PA(p1->y);
    PA(p0->x);
    PA(p0->y);
    PG(ma);
    PG(mb);
    PA(ma);
    PA(mb);
    P(p->spikecount);
    P(p->spikeidx);
    P(endidx);
    P(p->spiketimes[p->spikeidx]);
  }
  (m0 != m_min && m0 != m_max) || DIE;
  if (p->spikeidx == -1) {PC(m_min); PC(m_max);}

  //solve([yt=m0*xt+b0,m_min=(yi-yt)/(xi-xt)],[xt,yt]);
  xt = (m_min*xi-yi+b0)/(m_min-m0);
  //  printf ("%d %g %g %g %g %g %a %a \n", passcnt, xt, m_min, xi, yi, b0, m_min, m0);

  yt = (m0*(yi-m_min*xi)-b0*m_min)/(m0-m_min);
  //solve([yb=m0*xb+b0,m_min=(yi-yb)/(xi-xb)],[xb,yb]);
  xb = (m_max*xi-yi+b0)/(m_max-m0);
  yb = (m0*(yi-m_max*xi)-b0*m_max)/(m0-m_max);
  if (xt < xb) {xl = xt; xr = xb; }
  else         {xl = xb; xr = xt; }
  if (xr < pm.x || xl > p1->x) {
    if (p->pointcount == 1)
      p->pointcount = 0;
    if (p->spikeidx == -1) {
      PA(m_min);
      PA(m_max);
      P(p->spiketimes[p->spikeidx - 1]);
      PG(p->spikeidx - 1 + 1 - p->dqc);
      PG(p->spikeidx - 1 + p->dqc);
      P(p->spiketimes[p->spikeidx]);
      PG(p->spikeidx + 1 - p->dqc);
      PG(p->spikeidx + p->dqc);
      PG(p0->x);
      PG(p0->y);
      PG(p1->x);
      PG(p1->y);
      PG(xl);
      PG(xr);
      PG(pm.x);
      PG(p1->x);
      PA(xl);
      PA(xr);
      PA(pm.x);
      PA(p1->x);
    }
    endidx > p->spikeidx || DIE;
    return -1;
  }
  if (xl < pm.x) xl = pm.x;
  if (xr > p1->x) xr = p1->x;

  if (p->spikeidx == -1) {PC(xl); PC(xr);}

  x = (xl + xr) / 2;
  y = m0 * x + b0;
  m = (yi - y) / (xi - x);
  if (m > 1) {
    static mpq_class yl, yr, ml, mr;
    yl = m0 * xl + b0;
    yr = m0 * xr + b0;
    ml = (yi - yl) / (xi - xl);
    mr = (yi - yr) / (xi - xr);
    m = (ml <= 1 || mr <= 1) ? 1 : (ml < mr ? ml : mr);
  }
  b = yi - m * xi;
  if (p->spikeidx == -1) {PC(m0); PC(b);}

  p1->x = (b0-b)/(m-m0);
  p1->y = (b0*m-b*m0)/(m-m0);

  TREALLOC (p->point, ++p->pointcount);
  INIT (p->point[p->pointcount - 1].x);
  INIT (p->point[p->pointcount - 1].y);
  x = (endidx + p->d - b + 1) / m;// intersect the top?
  y = endidx + 1 + p->d;

  if (x < p->spiketimes[endidx]) {
    P(p->spikeidx);
//    P(endidx);
//    PC(m_min);
//    PC(m_max);
//    P(p->spiketimes[p->spikeidx - 1]);
//    P(p->spiketimes[p->spikeidx]);
//    PC(xl);
//    PC(xr);
//    PC(xi);
//    PC(yi);
//    PC(x);
//    PC(m);
//    PC(m0);
//    PC(b0);
//    PC(xt);
//    PC(yt);
//    PC(xb);
//    PC(yb);
  }

  xl = p->spiketimes[endidx];
  if (x < xl) {
    printf ("x: %g, xl: %g\n", x, xl);
    exit (DIE);
  }
  //  x >= xl || DIE;               // die unless right of the left end
  if (endidx + 1 < p->spikecount && x > p->spiketimes[endidx + 1]) {
    xr = p->spiketimes[endidx + 1];
    y = m * xr + b;             // then must intersect the right side
//    if (y > endidx + 2 - p->d) P(p->spikeidx);
//    y <= endidx + 2 - p->d || DIE;  // die unless it goes through the right
//    y >= endidx + 1 - p->d || DIE;  // hand end segment
    x = xr;
  }
  p1 = &p->point[p->pointcount - 1];
  p1->x = x;
  p1->y = y;

  return endidx + 1;
}

static int
findseg0 (Data *p, int spikeidx)
{
  if (spikeidx == -1) printf ("findseg0\n");
  if (spikeidx == -1) P(p->pointcount);
  if (spikeidx == -1) PC(p->point[p->pointcount - 2].x);
  if (spikeidx == -1) PC(p->point[3].x);

  int retval = 0;
  p->spikeidx = spikeidx;
  int endidx = p->spikecount - 1;
  do {
    if (spikeidx == -1 && retval == -1) printf ("getseg failed\n");
    if (spikeidx == -1) printf ("%d %d\n", spikeidx, endidx);
    p->hulls = chull (p->spiketimes, spikeidx, endidx, p->dqc);
    endidx = p->hulls.endidx - 1;
  } while ((retval = getseg_n (p)) == -1);
  if (spikeidx == -1) printf ("getseg succeeded, next spikeidx %d, endidx %d\n", retval, p->hulls.endidx);
  return retval;
}

static int
findseg (Data *p, int spikeidx)
{
  if (spikeidx == -1) printf ("findseg ");
  if (spikeidx == -1) P(p->pointcount);
  if (spikeidx == -1) PC(p->point[p->pointcount - 2].x);
  int retval;
  p->spikeidx = spikeidx;
  int hi = p->spikecount;
  int lo = spikeidx - 1;
  int next_spikeidx = 0;
  int orig_pointcount = p->pointcount;
  Point orig_lastpoint = {0}, new_lastpoint = {0};
  Point new_point[2];
  int added_pointcount = 0;
  
  //printf ("findseg spikeidx %d, spikecount %d\n", spikeidx, p->spikecount);
  if (orig_pointcount > 0){
    orig_lastpoint = p->point[p->pointcount - 1];
  }
  while (hi - lo > 1) {
    int mid = (hi + lo) / 2;
    //printf ("spikeidx: %d, lo: %d, mid: %d, hi: %d\n", spikeidx, lo, mid, hi);
    
    p->hulls = chull (p->spiketimes, spikeidx, mid, p->dqc);

    p->pointcount = orig_pointcount;
    if (orig_pointcount > 0) {
      p->point[p->pointcount - 1] = orig_lastpoint;
    }
    //printf ("mid %d, endidx %d\n", mid, p->hulls.endidx);

    (lo < mid && mid < hi) || DIE;
    if (p->hulls.endidx != mid) {
      p->hulls.endidx < mid || DIE;
      if (hi > p->hulls.endidx + 1)
        hi = p->hulls.endidx + 1;
      mid = p->hulls.endidx;
    }

    if (spikeidx == -1) printf ("%d %d\n", spikeidx, mid);
    //printf ("after: lo %d, mid %d, hi %d\n", lo, mid, hi);

    (lo <= mid && mid < hi) || DIE;
    if ((retval = getseg_n (p)) == -1) {
      //printf ("getseg failed\n");
      hi = mid;
    }
    else {
      if (orig_pointcount > 0)
        new_lastpoint = p->point[orig_pointcount - 1];
      added_pointcount = p->pointcount - orig_pointcount;
      added_pointcount <= 2 || DIE;
      for (int n = 0; n < added_pointcount; n++)
        new_point[n] = p->point[orig_pointcount + n];
      //printf ("getseg succeeded, next spikeidx %d\n", retval);
      lo = mid, next_spikeidx = retval;
    }
  }
  if (orig_pointcount > 0)
    p->point[orig_pointcount - 1] = new_lastpoint;
  p->pointcount = orig_pointcount + added_pointcount;
  TREALLOC (p->point, p->pointcount);
  for (int n = 0; n < added_pointcount; n++)
    p->point[orig_pointcount + n] = new_point[n];
  return next_spikeidx;
}

static bool oldway = false;

Point *
approx (int *spiketimes, int spikecount, mpq_class dq, int *count)
{
  int spikeidx = 0;
  static Data p;
  
  memset (&p, 0, sizeof p);
  p.d = GETD(dq);
  p.dqc = dq;
  p.spikecount = spikecount;
  p.spiketimes = spiketimes;
  
  if (0) findseg0 (&p, spikeidx);
  if (0) findseg (&p, spikeidx);

  time_t lasttime = time (0);
  //  printf ("%s line %d\n", __FILE__, __LINE__);
  while (spikeidx < spikecount) {
    time_t now = time (0);
    if (now - lasttime > 1) {
      //      printf ("spikeidx %d %d %d\n", spikeidx, back_count, chull_count);
      lasttime = now;
    }
    if (spikeidx == -1) {
      P(p.pointcount);
      PC(p.point[p.pointcount - 2].x);
      PC(p.point[p.pointcount - 2].y);
    }
    if (oldway)
      spikeidx = findseg0 (&p, spikeidx);
    else
      spikeidx = findseg (&p, spikeidx);
  }

  *count = p.pointcount;
  return p.point;
}

static double
one_segment_maxd (int *st, int sc)
{
  double m =  (double)(sc - 1) / (st[sc - 1] - st[0]);
  double b = .5 - m * st[0];
  int n;
  double maxd = 0;
  for (n = 0; n < sc; n++) {
    double yl = m * st[n] + b;
    double d = MAX (n + 1 - yl, yl - n);
    if (d > maxd)
      maxd = d;
  }
  return maxd;
}

typedef struct
{
  double g;
  double d;
  int count;
} SearchResult;

static int *spiketime;
static int spikecount;

static Pvoid_t result = (Pvoid_t) NULL;
static Pvoid_t t_est = (Pvoid_t) NULL;
static Pvoid_t isi = (Pvoid_t) NULL;

int count_min;
double d_min;
double g_min;
double d_est_min;

static double threshold_delta;

static bool
get_next_d_est (double *d_estp)
{
  Wrd_Flt_t Index, prev_Index;         // array index
  Word_t * PValue;                    // pointer to array element value
  double threshold = d_min + threshold_delta;
  SearchResult *r, *prev_r;
  double max_gap = 0;
  //  double gap_threshold = .1;
  double gap_threshold = 3;
  
  Index.w = 0;
  JLF (PValue, result, Index.w);
  r = (SearchResult *) *PValue;
  prev_Index = Index;
  JLN (PValue, result, Index.w);
  while (PValue != NULL) {
    prev_r = r;
    r = (SearchResult *) *PValue;
    if ((prev_r->d <= threshold || r->d <= threshold) && Index.f - prev_Index.f > .000001) {
      //      double gap = fabs (r->g - prev_r->g) + (Index.f - prev_Index.f) / 4;
      double gap = (abs (r->count - prev_r->count)
                    + Index.f - prev_Index.f
                    + fabs (r->g - prev_r->g));
      //      printf ("gap: %g\n", gap);
      if (gap > max_gap) {
        max_gap = gap;
        *d_estp = (Index.f + prev_Index.f) / 2;
      }
    }
    prev_Index = Index;
    JLN(PValue, result, Index.w);
  }

  return max_gap > gap_threshold;
}

static bool
add_result (double d_est, double d, double g, int count)
{
  PWord_t  PValue;                          // pointer to return value
  SearchResult *r = malloc (sizeof *r);
  r->g = g;
  r->d = d;
  r->count = count;
  if (d < d_min)
    d_min = d, g_min = g, d_est_min = d_est, count_min = count;
  Wrd_Flt_t Index;
  Index.f = d_est;
  JLI (PValue, result, Index.w);
  bool retval = *PValue == 0;
  *PValue = (Word_t) r;
  return retval;
}

static inline void
insert (Pvoid_t *ary, double x, double y)
{
  Wrd_Flt_t Index;
  P_Wrd_Flt_t PValue;
  
  if (x < 0)
    x = 0;
  Index.f = x;
  JLI (PValue.v, *ary, Index.w);
  *PValue.f = y;
}

static int *
dbl_to_int_times (double *td, int tn)
{
  int *ti = malloc (tn * sizeof *ti);
  for (int n = 0; n < tn; n++)
    ti[n] = floor ((td[n] * 10) + .5);
  return ti;
}

static int *
get_cth (int *st, int sc, int *e_begin, int *e_end, int ecnt)
{
  static int cth[100];

  int max = 0;
  for (int n = 0; n < ecnt; n++) {
    int len = e_end[n] - e_begin[n];
    if (len > max)
      max = len;
  }
  printf ("the longest cycle is %d .1ms ticks\n", max);
  memset (cth, 0, sizeof cth);
  int clen = ecnt ? (e_end[0] - e_begin[0]) : 0;
  int cycle = 0;
  for (int si = 0; si < sc; si++) {
    while (cycle < ecnt && st[si] >= e_end[cycle]) {
      cycle++;
      clen = (e_end[cycle] - e_begin[cycle]);
    }
    if (cycle >= ecnt)
      break;
    if (st[si] < e_begin[cycle])
      continue;
    cth[(st[si] - e_begin[cycle]) * 100LL / clen]++;
  }
  return cth;
}

#define CTH_BINS 100

Pvoid_t
cth_t_est (int *st, int sc, double *E_begin, double *E_end, int ecnt)
{
  int *e_begin = dbl_to_int_times (E_begin, ecnt);
  int *e_end   = dbl_to_int_times (E_end  , ecnt);
  int *cth = get_cth (st, sc, e_begin, e_end, ecnt);
  double y = 0;
  Pvoid_t esct;
  esct = NULL;
  insert (&esct, st[0] - (st[1] - st[0]) / 2, y);
  ecnt >= 0 || DIE;
  sc >= 2 || DIE;
  int stop = ecnt > 0 ?  e_begin[0] : st[sc - 1] + (st[sc - 1] - st[sc - 2]) / 2;
  int sn = 0;
  bool in_cycle = false;
  int cycle = 0;
  double adj = 0;
  double dy = 1e-6;
  int allbins = 0;
  for (int n = 0; n < CTH_BINS; n++)
    allbins += cth[n];
  int count;
  while (sn < sc) {
    for (count = 0; sn < sc && st[sn] < stop; sn++)
      count++;
    if (in_cycle) {
      double cycle_time  = e_end[cycle] - e_begin[cycle];
      int binsum = 0;
      int n;
      double y0 = y;
      for (n = 0; n < CTH_BINS; n++) {
        binsum += cth[n];
        double delta = (double) cth[n] / allbins * count;
        double tmp = 0;
        if (delta == 0) adj += dy, tmp = adj;
        if (delta - adj > dy) tmp = -adj, adj = 0;
        y =  y0 + (double)binsum / allbins * count;
        insert (&esct, e_begin[cycle] + (n + 1.) / CTH_BINS * cycle_time, y + tmp);
      }
      cycle++;
    }
    else {
      double tmp = 0;
      if (count == 0) adj += dy, tmp = adj;
      if (count - adj > dy) tmp = -adj, adj = 0;
      y += count;
      insert (&esct, stop, y + tmp);
    }
    if (sn >= sc)
      break;
    if (cycle >= ecnt) {
      in_cycle = false;
      stop = st[sc - 1] + (st[sc - 1] - st[sc - 2]) / 2;
    }
    else if (stop == e_begin[cycle]) {
      in_cycle = true;
      stop = e_end[cycle];
    }
    else if (stop < e_begin[cycle]) {
      in_cycle = false;
      stop = e_begin[cycle];
    }
    else /* if (stop > e_begin[cycle]) */ exit (DIE);
  }
  free (e_begin);
  free (e_end);
  return esct;
}

static bool
get_result (double d_est)
{
  double d, g;
  int count;
  Point *point;

  point = approx (spiketime, spikecount, d_est, &count);
  Pvoid_t esct;
  esct = NULL;
  int n;
  for (n = 0; n < count; n++)
    insert (&esct, GETD(point[n].x), GETD(point[n].y));
  free (point);

  d = gammafit (spiketime, spikecount, esct, &g);
  Word_t Rc_word;
  JLFA(Rc_word, esct);

  bool retval = add_result (d_est, d, g, count);
  //  printf ("get_result %g: %g %g %d\n", d_est, g, d, count);
  return retval;
}

void
plot_best (int *spiketimes, int spikecount)
{
  int count;
  Point *point;

  if (!getenv("EDT_SURROGATE_PLOT"))
    return;

  point = approx (spiketimes, spikecount, d_est_min, &count);
#define FNAME  "plot"
  FILE *f = fopen (FNAME, "w");
  int n;

  //  <sed 1,/^exit$/d  deleted
  fprintf (f, "set format x '%%9.0f'\n"
           "plot \\\n"
           "'"FNAME"' index 0 using 1:($2-%g) with lines notitle,\\\n"
           "'"FNAME"' index 0 using 1:($2+%g) with lines notitle,\\\n"
           "'"FNAME"' index 1 using 1:2        with lines notitle\n"
           "exit\n",
           GETD (d_est_min), GETD (d_est_min));
  
  for (n = 0; n < spikecount; n++)
    fprintf (f, "%d %d\n%d %d\n", spiketimes[n], n, spiketimes[n], n + 1);

  fprintf (f, "\n\n");

  for (n = 0; n < count; n++)
    fprintf (f, "%.17g %g\n", GETD(point[n].x),  GETD(point[n].y));
  

  fclose (f);
  free (point);
  if (system ("gnuplot -persist < " FNAME ));
  printf ("%d points\n", count);
  fflush (stdout);
}

void
set_preserve_t_est (bool val)
{
}
 
Pvoid_t
get_t_est (void)
{
  int count;
  Point *point = approx (spiketime, spikecount, d_est_min, &count);
  Pvoid_t esct;
  esct = NULL;
  int n;
  for (n = 0; n < count; n++)
    insert (&esct, GETD(point[n].x), GETD(point[n].y));
  free (point);

  return esct;
}

double
gammafit_search (int *st, int *lodif, int sc)
{
  double d_est = 0;

  fpu_check ();

  d_min = DBL_MAX;

  spiketime = st;
  spikecount = sc;

  Word_t Rc_word;
  JLC(Rc_word, t_est, 0, -1);
  if (Rc_word != 0) {
    fprintf (stderr, "%s line %d: t_est is not empty\n", __FILE__, __LINE__);
    exit (1);
  }

  if (0) {
    double osmd = one_segment_maxd (st, sc);
    int split = 100;
    int n;
    for (n = 0; n <= split; n++)
      get_result (.501 + n * (osmd - .501) / split);

    //  get_result (.501);
    //  get_result (.501 + (osmd - .501) / 3);
    //  get_result (.501 + 2 * (osmd - .501) / 3);
    //  get_result (osmd);

    threshold_delta = .03;
    int since_new_min = 0;
    double last_min = d_min;
    int wait_time = 5;
    double reduction_factor = .6;
    printf ("threshold_delta = %8.6f, reduction_factor = %8.6f, wait_time = %3d, split = %d\n", threshold_delta, reduction_factor, wait_time, split);
    while (get_next_d_est (&d_est)) {
      get_result (d_est);
      if (d_min < last_min) {
        since_new_min = 0;
        last_min = d_min;
      }
      else since_new_min++;
      if (since_new_min > wait_time) {
        threshold_delta *= reduction_factor;
        since_new_min = 0;
      }
    }
  }
  else get_result (1);

  JLFA(Rc_word, t_est);
  JLFA(Rc_word, result);
  JLFA(Rc_word, isi);

  plot_best (st, sc);
  
  return g_min;
}
#endif
