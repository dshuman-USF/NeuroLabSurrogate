#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include <time.h>
#include "edt.h"

#define SURCNT 20
#define BINCNT 250
#define TICKS_PER_SEC 10000

/* for repeat pattern spectrum (repat) */
#define WLEN 150
#define MAXPLEX 16
#define MAXREPPAT 20

//#define TERMINAL "set terminal postscript eps color\n"
//#define FILETYPE "eps"

#define TERMINAL "set terminal png\n"
#define FILETYPE "png"
//#define REPLOT "replot\n"
#define REPLOT ""

typedef int HRow[BINCNT];
typedef double HRowD[BINCNT];
typedef int RepRow[MAXREPPAT + 1];
typedef double RepRowD[MAXREPPAT + 1];

typedef struct
{
  double *hist;
  int bincnt;
  double *xval;
  char *title;
  char *color;
  int linetype;
} Plot;

typedef struct
{
  HRowD *hist;
  int bincnt;
  double *xval;
  double min;
  double max;
} Plot2d;

typedef struct
{
  int isi;
  int id_idx;
  int end_spike;
  int start_time;
  int end_time;
  int index;
} Point;

static void
plot_repat (RepRowD *data, char *title, char *fname)
{
  FILE *f = fopen (fname, "w");
  char *psfname;
  asprintf (&psfname, "%s.%s", fname, FILETYPE);
  fprintf (f, 
           "set title \"%s\"\n"
           TERMINAL
           "set output \"%s\"\n"
           "set pm3d map\n"
           "set palette functions 2*gray,1-abs(2*gray-1),2-2*gray\n"
           "set cbrange [-3:3]\n"
           "set pm3d corners2color c1\n"
           "set xrange [2:%d]\n"
           "set yrange [3:%d]\n"
           "set xlabel \"number of repetitions\"\n"
           "set ylabel \"complexity\"\n"
           "set cblabel \"SD\"\n"
           "splot \"< sed -n '1,/^exit/!p' %s\" using ($1+2):($2+3):3 matrix title \"\"\n"
           "set terminal wxt\n"
           REPLOT
           "exit\n"
           ,
           title, psfname, MAXREPPAT + 1, MAXPLEX + 1, fname);
  free (psfname);

  for (int row = 3; row <= MAXPLEX; row++) {
    for (int col = 2; col <= MAXREPPAT; col++)
      fprintf (f, " %g", data[row][col]);
    fprintf (f, "\n");
  }

  fclose (f);
  char *cmd;
  asprintf (&cmd, "gnuplot -persist < %s", fname);
  system (cmd);
  free (cmd);
}

static int qsort_size;

int
compare_mem (const void *ap, const void *bp)
{
  return memcmp (ap, bp, qsort_size);
}

int
compare_points (const void *ap, const void *bp)
{
  const Point *a = (const Point *) ap;
  const Point *b = (const Point *) bp;
  return (a->start_time > b->start_time) - (a->start_time < b->start_time);
}

static inline void
print_this_t (Point *this_t, int count, Edt *edt)
{
  for (int n = 0; n < count; n++) {
    Point *p = &this_t[n];
    Digital *d = edt->digital[p->id_idx];
    printf ("%d %d, ", d->id, p->end_time);
  }
  printf ("\n");
}

static inline void
find_zero (Point *point, int point_count, Point *this_t, int *count_p)
{
  int count = 0;
  for (int n = 0; n < point_count; n++)
    if (point[n].isi == 0) {
      memcpy (&this_t[count++], &point[n], sizeof (Point));
    }
  *count_p = count;
}

static inline int
get_inwin (Point *point, int start, int pcount, int wlen)
{
  int t0 = point[start].start_time;
  int inwin = 1;
  for (int n = 1; n < pcount; n++)
    if (point[start + n].start_time - t0 <= wlen)
      inwin++;
    else
      break;
  return inwin;
}

static inline void
append (int *psv[MAXPLEX + 1], int *psvcnt, int plex, Point *this_t, Edt *edt)
{
  static int callcnt;
  callcnt++;
  if (plex > MAXPLEX)
    return;
  if (psvcnt[plex] % 1000000 == 0)
    psv[plex] = realloc (psv[plex], (psvcnt[plex] + 1000000) * plex * 2 * sizeof (int));

  int *psvp1a = &psv[plex][psvcnt[plex] * 2 * plex];
  int *psvp1b = psvp1a + plex;
  int *psvp2a = psvp1b + plex;
  int *psvp2b = psvp2a + plex;

  for (int n = 0; n < plex; n++) {
    Point *p = &this_t[n];
    Digital *d = edt->digital[p->id_idx];
    *psvp1a++ = d->id;
    *psvp1b++ = p->start_time;
    *psvp2a++ = d->id;
    *psvp2b++ = p->end_time;
  }

  if (0) {
    int *psvp1 = &psv[plex][psvcnt[plex] * 2 * plex];
    int *psvp2 = &psv[plex][(psvcnt[plex] + 1) * 2 * plex];
    if (0) printf ("    %d: ", callcnt);
    for (int n = 0; n < 2 * plex; n++)
      printf (" %d", psvp2[n]);
    printf (" %d\n", psvp2[plex] - psvp1[plex]);
  }

  psvcnt[plex] += 2;
}

typedef struct
{
  int id;
  int time;
} PatEnd;

static PatEnd
get_patend (Point *this_t, int count, Edt *edt)
{
  PatEnd pe;
  Point *p = &this_t[count - 1];
  Digital *d = edt->digital[p->id_idx];
  pe.id = d->id;
  pe.time = p->start_time;
  return pe;
}

static void
diff (int *psv, int psvcnt)
{
  int int_cnt = qsort_size / sizeof (int);
  int plex = int_cnt / 2;
  for (int n = 0; n < psvcnt; n++) {
    int *p = psv + n * int_cnt + plex;
    for (int i = 0; i < plex - 1; i++)
      p[i] = p[i + 1] - p[i];
    p[plex - 1] = 0;
  }
}

static void
countrep (int *psv, int psvcnt, int *repforn)
{
  int int_cnt = qsort_size / sizeof (int);
  int a = 0, b = 0;

  while (a < psvcnt) {
    int *ap = psv + a * int_cnt;
    while (b + 1 < psvcnt && memcmp (ap, psv + (b + 1) * int_cnt, qsort_size) == 0)
      b++;
    int rep = b - a + 1;
    if (rep <= MAXREPPAT)
      repforn[rep]++;
    a = ++b;
  }
}

static int
unique (int *psv, int psvcnt)
{
  int int_cnt = qsort_size / sizeof (int);
  int *to_p = psv;
  for (int from = 1; from < psvcnt; from++) {
    int *from_p = psv + from * int_cnt;
    if (memcmp (from_p, to_p, qsort_size) != 0) {
      to_p += int_cnt;
      if (from_p !=  to_p)
        memcpy (to_p, from_p, qsort_size);
    }
  }
  return (to_p - psv) / int_cnt + 1;
}

static RepRow *
repat (Edt *edt)
{
  int point_count = 0;
  for (int n = 0; n < edt->digital_count; n++)
    point_count += edt->digital[n]->spike_count - 1;

  Point *point = malloc (point_count * sizeof *point);

  int maxt = 0;
  Point *p = point;
  for (int dn = 0; dn < edt->digital_count; dn++) {
    Digital *d = edt->digital[dn];
    if (d->spike[d->spike_count - 1] > maxt)
      maxt = d->spike[d->spike_count - 1];
    for (int sn = 0; sn < d->spike_count - 1; sn++) {
      p->isi = d->spike[sn + 1] - d->spike[sn];
      p->id_idx = dn;
      p->end_spike = sn + 1;
      p->start_time = d->spike[sn];
      p->end_time = d->spike[sn + 1];
      p->index = p - point;
      p++;
    }
  }
  assert (p - point == point_count);
  int *psv[MAXPLEX + 1];
  int psvcnt[MAXPLEX + 1];
  memset (psv, 0, sizeof psv);
  memset (psvcnt, 0, sizeof psvcnt);
  
  Point *this_t = malloc (point_count * sizeof *this_t);

  time_t now, last_time = 0;
  for (int shift = 0; shift <= maxt; shift++) {
    if ((now = time (0)) != last_time) {
      printf ("repat: shift %d of %d\n", shift, maxt);
      last_time = now;
    }
    PatEnd prevpatend = {0};
    int count;
    find_zero (point, point_count, this_t, &count);

    if (0) {
      qsort (this_t, count, sizeof *this_t, compare_points);
      print_this_t (this_t, count, edt);
    }

    if (count >= 3) {
      qsort (this_t, count, sizeof *this_t, compare_points);
      for (int m = 0; m < count - 2; m++) {
        int inwin = get_inwin (this_t, m, count - m, WLEN);
        if (0) printf ("%d %d: inwin is %d\n", m, count, inwin);
        if (inwin >= 3) {
          PatEnd nowpatend = get_patend (this_t + m, inwin, edt);
          if (memcmp (&nowpatend, &prevpatend, sizeof (PatEnd)) != 0) {
            append (psv, psvcnt, inwin, this_t + m, edt);
            prevpatend = nowpatend;
          }
        }
      }
    }
    for (int n = 0; n < count; n++) {
      Point *p = &point[this_t[n].index];
      Digital *d = edt->digital[p->id_idx];
      p->end_spike++;
      if (p->end_spike >= d->spike_count)
        p->isi = INT_MAX;
      else {
        p->end_time = d->spike[p->end_spike];
        p->isi = p->end_time - p->start_time - shift;
      }
    }
    for (int n = 0; n < point_count; n++)
      point[n].isi--;
  }
  free (this_t);
  for (int n = 0; n <= MAXPLEX; n++)
    psv[n] = realloc (psv[n], psvcnt[n] * (2 * n) * sizeof (int));

  if (0) {
    for (int n = 0; n < psvcnt[3]; n += 2) {
      for (int i = 0; i < 3; i++)
        printf (" %d", psv[3][n * 6 + i]);
      for (int i = 0; i < 3; i++)
        printf (" %d", psv[3][(n + 1) * 6 + 3 + i]);
      printf (" %d\n", psv[3][(n + 1) * 6 + 3] - psv[3][n * 6 + 3]);
    }
  }

  RepRow *repforn = calloc (MAXPLEX + 1, sizeof *repforn);

  for (int n = 3; n <= MAXPLEX; n++) {
    if (psvcnt[n] == 0)
      continue;
    qsort_size = n * 2 * sizeof *psv[0];
    qsort (psv[n], psvcnt[n], qsort_size, compare_mem);
    psvcnt[n] = unique (psv[n], psvcnt[n]);
    diff (psv[n], psvcnt[n]);
    qsort (psv[n], psvcnt[n], qsort_size, compare_mem);
    countrep (psv[n], psvcnt[n], &repforn[n][0]);
  }
  for (int n = 0; n <= MAXPLEX; n++)
    free (psv[n]);
  for (int n = 3; n <= MAXPLEX; n++) {
    for (int i = 0; i <= MAXREPPAT; i++)
      printf (" %d", repforn[n][i]);
    printf ("\n");
  }  
  return repforn;
}

static char *yrange = "";

static void
plot (Plot *p, int count, char *title, char *fname)
{
  FILE *f = fopen (fname, "w");
  char *psfname;
  if (asprintf (&psfname, "%s.%s", fname, FILETYPE) == -1) exit (1);
  fprintf (f, 
           "set title \"%s\"\n"
           TERMINAL
           "set output \"%s\"\n"
           "set xrange [%g:%g]\n"
           "%s"
           "set xlabel \"seconds\"\n",
           title, psfname, p[0].xval[0], p[0].xval[p[0].bincnt - 1], yrange);
  free (psfname);
  fprintf (f, "plot \\\n");
  for (int n = 0; n < count; n++)
    fprintf (f, "'%s' index %d using 1:2 with histeps lc rgb \"%s\" lt %d title \"%s\"%s\n",
             fname, n, p[n].color, p[n].linetype, p[n].title, n + 1 < count ? ", \\" : "");
  fprintf (f, 
           "set terminal wxt\n"
           REPLOT
           "exit\n"
           );
  
  for (int histidx = 0; histidx < count; histidx++) {
    for (int n = 0; n < p[histidx].bincnt; n++)
      fprintf (f, "%g %g\n", p[histidx].xval[n], p[histidx].hist[n]);
    fprintf (f, "\n\n");
  }

  fclose (f);
  char *cmd;
  if (asprintf (&cmd, "gnuplot -persist < %s", fname) == -1) exit (1);
  if (system (cmd));
  free (cmd);
}

static void
plot2d (Plot2d *p, int count, char *title, char *fname)
{
  FILE *f = fopen (fname, "w");
  char *psfname;
  if (asprintf (&psfname, "%s.%s", fname, FILETYPE) == -1) exit (1);
  double x0 = p[0].xval[0];
  double xe = p[0].xval[p[0].bincnt - 1];
  double dx = p[0].xval[1] -  p[0].xval[0];
  fprintf (f, 
           "set title \"%s\"\n"
           TERMINAL
           "set output \"%s\"\n"
           "set pm3d map\n"
           "set palette functions 2*gray,1-abs(2*gray-1),2-2*gray\n"
           "set cbrange [%g:%g]\n"
           "set pm3d corners2color c1\n"
           "set xrange [%g:%g]\n"
           "set yrange [%g:%g]\n"
           "set xlabel \"seconds\"\n"
           "set ylabel \"seconds\"\n"
           "set cblabel \"SD\"\n"
           "splot \"< sed -n '1,/^exit/!p' %s\" using ($1*%g):($2*%g):3 matrix title \"\"\n"
           "set terminal wxt\n"
           REPLOT
           "exit\n"
           ,
           title, psfname, p[0].min, p[0].max, x0, xe, x0, xe, fname, dx, dx);
  free (psfname);

  for (int row = 0; row < BINCNT; row++) {
    for (int col = 0; col < BINCNT; col++)
      fprintf (f, " %g", p[0].hist[row][col]);
    fprintf (f, "\n");
  }

  fclose (f);
  char *cmd;
  if (asprintf (&cmd, "gnuplot -persist < %s", fname) == -1) exit (1);
  if (system (cmd));
  free (cmd);
}

static int *
isi (Edt *edt, int *binw_p, int *origin)
{
  int *spike = edt->digital[0]->spike;
  int spcnt = edt->digital[0]->spike_count;
  int *h = calloc (BINCNT, sizeof *h);
  /* 10/sec is .1 mean ISI.  Let's do .3 width */
  int psth_width = .3 * TICKS_PER_SEC;
  int binw = psth_width / BINCNT;
  binw = 20;
  psth_width = BINCNT * binw;
  if (binw_p) *binw_p = binw;
  if (origin) *origin = 0;
  assert (BINCNT * binw == psth_width);
  for (int n = 1; n < spcnt; n++) {
    int i = spike[n] - spike[n - 1];
    int bin = i / binw;
    if (bin < BINCNT)
      h[bin]++;
  }
  return h;
}

static HRow *
jint (Edt *edt, int *binw_p, int *origin)
{
  int *spike = edt->digital[0]->spike;
  int spcnt = edt->digital[0]->spike_count;
  int (*h)[BINCNT] = calloc (BINCNT * BINCNT, sizeof *h);
  /* 10/sec is .1 mean ISI.  Let's do .3 width */
  int psth_width = .3 * TICKS_PER_SEC;
  int binw = psth_width / BINCNT;
  binw = 20;
  psth_width = BINCNT * binw;
  if (binw_p) *binw_p = binw;
  if (origin) *origin = 0;
  assert (BINCNT * binw == psth_width);
  int lastbin = -1;
  for (int n = 1; n < spcnt; n++) {
    int i = spike[n] - spike[n - 1];
    int bin = i / binw;
    if (bin < BINCNT && lastbin >= 0)
      h[lastbin][bin]++;
    lastbin = bin;
  }
  return h;
}

static int *
psth (Edt *edt, int *binw_p, int *origin)
{
  int *spike = edt->digital[0]->spike;
  int spcnt = edt->digital[0]->spike_count;
  int *h = calloc (BINCNT, sizeof *h);
  int psth_width = 2 * TICKS_PER_SEC;
  int binw = psth_width / BINCNT;
  if (binw_p) *binw_p = binw;
  assert (BINCNT * binw == psth_width);
  assert (binw % 10 == 0);
  if (origin) *origin = 0;
  for (int n = 0; n < spcnt; n++) {
    int bin = (spike[n] % psth_width) / binw;
    h[bin]++;
  }
  return h;
}

static HRow *
jpsth (Digital *rdig, Digital *tdig, int *binw_p, int *origin)
{
  int *rspike = rdig->spike;
  int rspcnt = rdig->spike_count;
  int *tspike = tdig->spike;
  int tspcnt = tdig->spike_count;
  HRow *h = calloc (BINCNT, sizeof *h);
  int cyclen = 2 * TICKS_PER_SEC;
  int psth_width = 2 * TICKS_PER_SEC / 8;
  int binw = psth_width / BINCNT;
  if (binw_p) *binw_p = binw;
  assert (BINCNT * binw == psth_width);
  assert (binw % 10 == 0);
  if (origin) *origin = 0;
  int t0 = 0;
  int cyc = 0;
  for (int r = 0; r < rspcnt; r++) {
    if (rspike[r] / cyclen != cyc) {
      cyc = rspike[r] / cyclen;
      while (t0 < tspcnt && tspike[t0] / cyclen < cyc)
        t0++;
    }
    int rbin = (rspike[r] % cyclen) / binw;
    if (rbin < BINCNT)
      for (int t = t0; t < tspcnt && tspike[t] / cyclen == cyc; t++) {
        int tbin = (tspike[t] % cyclen) / binw;
        if (tbin < BINCNT)
          h[rbin][tbin]++;
      }
  }
  return h;
}

static int *
autocorr (Edt *edt, int *binw_p, int *origin)
{
  int *spike = edt->digital[0]->spike;
  int spcnt = edt->digital[0]->spike_count;
  int *h = calloc (BINCNT, sizeof *h);
  int psth_width = 2 * TICKS_PER_SEC / 8;
  int binw = psth_width / BINCNT;
  if (binw % 20 == 0) binw += 10;
  psth_width = BINCNT * binw;
  //  printf ("autocorr binw: %d\n", binw);
  if (binw_p) *binw_p = binw;
  if (origin) *origin = 0;
  int tar0 = 0;
  int tar = 0;
  for (int ref = 0; ref < spcnt; ref++) {
    tar = tar0;
    int bin = 0;
    int rsp = spike[ref];
    while (tar < spcnt && (bin = (spike[tar]-rsp) / binw)) {
      //printf ("loop1: spike-tar: %d bin %d\n", spike[tar]-rsp, bin);
      tar++;
    }
    tar0 = tar;
    while (tar < spcnt && bin < BINCNT) {
      if (spike[tar] != rsp)
        h[bin]++;
      if (++tar == spcnt)
        break;
      bin = (spike[tar]-rsp) / binw;
      //printf ("loop1: spike-tar: %6d  bin %6d\n", spike[tar]-rsp, bin);
    }
  }
  return h;
}

static int *
xcorr (Digital *rdig, Digital *tdig, int *binw_p, int *origin)
{
  int *rspike = rdig->spike;
  int rspcnt = rdig->spike_count;
  int *tspike = tdig->spike;
  int tspcnt = tdig->spike_count;
  int *h = calloc (BINCNT, sizeof *h);
  int psth_width = 2 * TICKS_PER_SEC / 8;
  int binw = psth_width / BINCNT;
  if (binw % 20 == 0) binw += 10;
  psth_width = BINCNT * binw;
  if (binw_p) *binw_p = binw;
  if (origin) *origin = -(binw*BINCNT-binw-1)/2;
  int tar0 = 0;
  int tar = 0;
  for (int ref = 0; ref < rspcnt; ref++) {
    tar = tar0;
    int bin = 0;
    int rsp = rspike[ref];
    while (tar < tspcnt && (bin = floor ((double)((tspike[tar]-rsp)+binw*BINCNT/2-(binw-1)/2-1)/binw)) < 0) {
      //printf ("loop1: spike-tar: %d bin %d\n", spike[tar]-rsp, bin);
      tar++;
    }
    tar0 = tar;
    while (tar < tspcnt && bin < BINCNT) {
      h[bin]++;
      if (++tar == tspcnt)
        break;
      bin = floor ((double)((tspike[tar]-rsp)+binw*BINCNT/2-(binw-1)/2-1)/binw);
      //printf ("loop1: spike-tar: %6d  bin %6d\n", spike[tar]-rsp, bin);
    }
  }
  return h;
}

static void
do_test (Edt *edt[SURCNT + 1], char *datatype, char *surtype, int *(*test)(Edt*,int*,int*), char *tstnm)
{
  int binw, origin;
  int *original = test (edt[0], &binw, &origin);
  double dblorig[BINCNT];
  for (int bin = 0; bin < BINCNT; bin++)
    dblorig[bin] = original[bin];
  free (original);
  int *surr[SURCNT];
  for (int n = 0; n < SURCNT; n++)
    surr[n] = test (edt[n + 1], 0, 0);
  double mean[BINCNT];
  memset (mean, 0, sizeof mean);
  for (int surn = 0; surn < SURCNT; surn++)
    for (int bin = 0; bin < BINCNT; bin++)
      mean[bin] += surr[surn][bin];
  for (int bin = 0; bin < BINCNT; bin++)
    mean[bin] /= SURCNT;
  double sd[BINCNT];
  memset (sd, 0, sizeof sd);
  for (int surn = 0; surn < SURCNT; surn++)
    for (int bin = 0; bin < BINCNT; bin++)
      sd[bin] += pow (surr[surn][bin] - mean[bin], 2);
  for (int surn = 0; surn < SURCNT; surn++)
    free (surr[surn]);
  double upper[BINCNT];
  double lower[BINCNT];
  double xval[BINCNT];
  for (int bin = 0; bin < BINCNT; bin++) {
    sd[bin] = sqrt (sd[bin] / (SURCNT - 1));
    upper[bin] = mean[bin] + 3 * sd[bin];
    lower[bin] = mean[bin] - 3 * sd[bin];
    xval[bin] = (double)(origin + bin * binw) / TICKS_PER_SEC;
  }
  
  Plot plotdata[] = {
    {dblorig, BINCNT, xval, "original" , "black", 1},
    {mean   , BINCNT, xval, "surr mean", "green", 1},
    {upper  , BINCNT, xval, "mean +3SD", "red"  , 2},
    {lower  , BINCNT, xval, "mean -3SD", "red"  , 2},
  };

  char *title, *fname;
  if (asprintf (&title, "%s, %s data, %s surrogate", tstnm, datatype, surtype) == -1) exit (1);
  char *s = strdup (surtype);
  char *space = strchr (s, ' ');
  if (space) *space = 0;
  for (int n = 0; n < strlen (s); n++)
    s[n] = tolower (s[n]);
  if (asprintf (&fname, "%s_%s_%s", tstnm, datatype, s) == -1) exit (1);
  plot (plotdata, sizeof plotdata / sizeof *plotdata, title, fname);
  free (title);
  free (fname);
}

static void
do_x_test_orig (Edt *edt[SURCNT + 1], char *datatype, char *surtype, int *(*xtest)(Digital*,Digital*,int*,int*), char *tstnm)
{
  int binw, origin;
  int *surr[SURCNT];
  for (int n = 0; n < SURCNT; n++)
    surr[n] = xtest (edt[0]->digital[0], edt[n + 1]->digital[0], &binw, &origin);
  double mean[BINCNT];
  memset (mean, 0, sizeof mean);
  for (int surn = 0; surn < SURCNT; surn++)
    for (int bin = 0; bin < BINCNT; bin++)
      mean[bin] += surr[surn][bin];
  for (int bin = 0; bin < BINCNT; bin++)
    mean[bin] /= SURCNT;
  double sd[BINCNT];
  memset (sd, 0, sizeof sd);
  for (int surn = 0; surn < SURCNT; surn++)
    for (int bin = 0; bin < BINCNT; bin++)
      sd[bin] += pow (surr[surn][bin] - mean[bin], 2);
  for (int surn = 0; surn < SURCNT; surn++)
    free (surr[surn]);
  double upper[BINCNT];
  double lower[BINCNT];
  double xval[BINCNT];
  for (int bin = 0; bin < BINCNT; bin++) {
    sd[bin] = sqrt (sd[bin] / (SURCNT - 1));
    upper[bin] = mean[bin] + 3 * sd[bin];
    lower[bin] = mean[bin] - 3 * sd[bin];
    xval[bin] = (double)(origin + bin * binw) / TICKS_PER_SEC;
  }
  
  Plot plotdata[] = {
    {mean   , BINCNT, xval, "surr mean", "green", 1},
    {upper  , BINCNT, xval, "mean +3SD", "red"  , 2},
    {lower  , BINCNT, xval, "mean -3SD", "red"  , 2},
  };

  char *title, *fname;
  if (asprintf (&title, "%s, %s data, %s surrogate", tstnm, datatype, surtype) == -1) exit (1);
  char *s = strdup (surtype);
  char *space = strchr (s, ' ');
  if (space) *space = 0;
  for (int n = 0; n < strlen (s); n++)
    s[n] = tolower (s[n]);
  if (asprintf (&fname, "%s_%s_%s", tstnm, datatype, s) == -1) exit (1);
  plot (plotdata, sizeof plotdata / sizeof *plotdata, title, fname);
  free (title);
  free (fname);
}

static void
do_x_test (Edt *edt[SURCNT + 1], char *datatype, char *surtype, int *(*xtest)(Digital*,Digital*,int*,int*), char *tstnm)
{
  int binw, origin;
  int *original = xtest (edt[0]->digital[0], edt[0]->digital[1], &binw, &origin);
  double dblorig[BINCNT];
  for (int bin = 0; bin < BINCNT; bin++)
    dblorig[bin] = original[bin];
  free (original);

  int *surr[SURCNT];
  for (int n = 0; n < SURCNT; n++)
    surr[n] = xtest (edt[n + 1]->digital[0], edt[n + 1]->digital[1], 0, 0);
  double mean[BINCNT];
  memset (mean, 0, sizeof mean);
  for (int surn = 0; surn < SURCNT; surn++)
    for (int bin = 0; bin < BINCNT; bin++)
      mean[bin] += surr[surn][bin];
  for (int bin = 0; bin < BINCNT; bin++)
    mean[bin] /= SURCNT;
  double sd[BINCNT];
  memset (sd, 0, sizeof sd);
  for (int surn = 0; surn < SURCNT; surn++)
    for (int bin = 0; bin < BINCNT; bin++)
      sd[bin] += pow (surr[surn][bin] - mean[bin], 2);
  for (int surn = 0; surn < SURCNT; surn++)
    free (surr[surn]);
  double upper[BINCNT];
  double lower[BINCNT];
  double xval[BINCNT];
  for (int bin = 0; bin < BINCNT; bin++) {
    sd[bin] = sqrt (sd[bin] / (SURCNT - 1));
    upper[bin] = mean[bin] + 3 * sd[bin];
    lower[bin] = mean[bin] - 3 * sd[bin];
    xval[bin] = (double)(origin + bin * binw) / TICKS_PER_SEC;
  }
  
  Plot plotdata[] = {
    {dblorig, BINCNT, xval, "original" , "black", 1},
    {mean   , BINCNT, xval, "surr mean", "green", 1},
    {upper  , BINCNT, xval, "mean +3SD", "red"  , 2},
    {lower  , BINCNT, xval, "mean -3SD", "red"  , 2},
  };

  char *title, *fname;
  if (asprintf (&title, "%s, %s data, %s surrogate", tstnm, datatype, surtype) == -1) exit (1);
  char *s = strdup (surtype);
  char *space = strchr (s, ' ');
  if (space) *space = 0;
  for (int n = 0; n < strlen (s); n++)
    s[n] = tolower (s[n]);
  if (asprintf (&fname, "%s_%s_%s", tstnm, datatype, s) == -1) exit (1);
  plot (plotdata, sizeof plotdata / sizeof *plotdata, title, fname);
  free (title);
  free (fname);
}

static void
do_test_2d (Edt *edt[SURCNT + 1], char *datatype, char *surtype, HRow *(*test)(Edt*,int*,int*), char *tstnm)
{
  int binw, origin;
  HRow *original = test (edt[0], &binw, &origin);
  double dblorig[BINCNT][BINCNT];
  for (int row = 0; row < BINCNT; row++)
    for (int col = 0; col < BINCNT; col++)
      dblorig[row][col] = original[row][col];
  free (original);
  HRow *surr[SURCNT];
  for (int n = 0; n < SURCNT; n++)
    surr[n] = test (edt[n + 1], 0, 0);
  double mean[BINCNT][BINCNT];
  memset (mean, 0, sizeof mean);
  for (int surn = 0; surn < SURCNT; surn++)
    for (int row = 0; row < BINCNT; row++)
      for (int col = 0; col < BINCNT; col++)
        mean[row][col] += surr[surn][row][col];
  for (int row = 0; row < BINCNT; row++)
    for (int col = 0; col < BINCNT; col++)
      mean[row][col] /= SURCNT;
  double sd[BINCNT][BINCNT];
  memset (sd, 0, sizeof sd);
  for (int surn = 0; surn < SURCNT; surn++)
    for (int row = 0; row < BINCNT; row++)
      for (int col = 0; col < BINCNT; col++)
        sd[row][col] += pow (surr[surn][row][col] - mean[row][col], 2);
  for (int surn = 0; surn < SURCNT; surn++)
    free (surr[surn]);
  double xval[BINCNT];
  for (int row = 0; row < BINCNT; row++) {
    xval[row] = (double)(origin + row * binw) / TICKS_PER_SEC;
    for (int col = 0; col < BINCNT; col++) {
      double s = sqrt (sd[row][col] / (SURCNT - 1));
      if (s == 0) s = .5;
      sd[row][col] = (dblorig[row][col] - mean[row][col]) / s;
    }
  }
  
  Plot2d plotdata[] = {
    {sd     , BINCNT, xval, -3, 3},
  };

  char *title, *fname;
  if (asprintf (&title, "%s, %s data, %s surrogate", tstnm, datatype, surtype) == -1) exit (1);
  char *s = strdup (surtype);
  char *space = strchr (s, ' ');
  if (space) *space = 0;
  for (int n = 0; n < strlen (s); n++)
    s[n] = tolower (s[n]);
  if (asprintf (&fname, "%s_%s_%s", tstnm, datatype, s) == -1) exit (1);
  plot2d (plotdata, sizeof plotdata / sizeof *plotdata, title, fname);
  free (title);
  free (fname);
}

static void
do_x_test_2d (Edt *edt[SURCNT + 1], char *datatype, char *surtype, HRow *(*xtest)(Digital*,Digital*,int*,int*), char *tstnm)
{
  int binw, origin;
  HRow *original = xtest (edt[0]->digital[0], edt[0]->digital[1], &binw, &origin);
  double dblorig[BINCNT][BINCNT];
  double max = 0;
  for (int row = 0; row < BINCNT; row++)
    for (int col = 0; col < BINCNT; col++) {
      dblorig[row][col] = original[row][col];
      if (original[row][col] > max)
        max = original[row][col];
    }
  free (original);
  HRow *surr[SURCNT];
  for (int n = 0; n < SURCNT; n++)
    surr[n] = xtest (edt[n + 1]->digital[0], edt[n + 1]->digital[1], 0, 0);
  double mean[BINCNT][BINCNT];
  memset (mean, 0, sizeof mean);
  for (int surn = 0; surn < SURCNT; surn++)
    for (int row = 0; row < BINCNT; row++)
      for (int col = 0; col < BINCNT; col++)
        mean[row][col] += surr[surn][row][col];
  for (int row = 0; row < BINCNT; row++)
    for (int col = 0; col < BINCNT; col++)
      mean[row][col] /= SURCNT;
  double sd[BINCNT][BINCNT];
  memset (sd, 0, sizeof sd);
  for (int surn = 0; surn < SURCNT; surn++)
    for (int row = 0; row < BINCNT; row++)
      for (int col = 0; col < BINCNT; col++)
        sd[row][col] += pow (surr[surn][row][col] - mean[row][col], 2);
  for (int surn = 0; surn < SURCNT; surn++)
    free (surr[surn]);
  double xval[BINCNT];
  for (int row = 0; row < BINCNT; row++) {
    xval[row] = (double)(origin + row * binw) / TICKS_PER_SEC;
    for (int col = 0; col < BINCNT; col++) {
      double s = sqrt (sd[row][col] / (SURCNT - 1));
      if (s == 0) s = .5;
      sd[row][col] = (dblorig[row][col] - mean[row][col]) / s;
    }
  }
  
  Plot2d plotdata[] = {
    {sd     , BINCNT, xval, -3, 3},
    {dblorig, BINCNT, xval, 0, max},
  };

  char *title, *fname;
  if (asprintf (&title, "%s, %s data, %s surrogate", tstnm, datatype, surtype) == -1) exit (1);
  char *s = strdup (surtype);
  char *space = strchr (s, ' ');
  if (space) *space = 0;
  for (int n = 0; n < strlen (s); n++)
    s[n] = tolower (s[n]);
  if (asprintf (&fname, "%s_%s_%s", tstnm, datatype, s) == -1) exit (1);
  plot2d (plotdata, sizeof plotdata / sizeof *plotdata, title, fname);
  free (title);
  free (fname);
}

static void
do_repat (Edt *edt[SURCNT + 1], char *datatype, char *surtype, char *tstnm)
{
  RepRow *original = repat (edt[0]);
  RepRow *surr[SURCNT];
  for (int n = 0; n < SURCNT; n++)
    surr[n] = repat (edt[n + 1]);
  double mean[MAXPLEX + 1][MAXREPPAT + 1];
  memset (mean, 0, sizeof mean);
  for (int surn = 0; surn < SURCNT; surn++)
    for (int row = 0; row <= MAXPLEX; row++)
      for (int col = 0; col <= MAXREPPAT; col++)
        mean[row][col] += surr[surn][row][col];
  for (int row = 0; row <= MAXPLEX; row++)
    for (int col = 0; col <= MAXREPPAT; col++)
      mean[row][col] /= SURCNT;
  double sd[MAXPLEX + 1][MAXREPPAT + 1];
  memset (sd, 0, sizeof sd);
  for (int surn = 0; surn < SURCNT; surn++)
    for (int row = 0; row <= MAXPLEX; row++)
      for (int col = 0; col <= MAXREPPAT; col++)
        sd[row][col] += pow (surr[surn][row][col] - mean[row][col], 2);
  for (int surn = 0; surn < SURCNT; surn++)
    free (surr[surn]);
  for (int row = 0; row <= MAXPLEX; row++) {
    for (int col = 0; col <= MAXREPPAT; col++) {
      double s = sqrt (sd[row][col] / (SURCNT - 1));
      if (s == 0) s = .5;
      sd[row][col] = (original[row][col] - mean[row][col]) / s;
    }
  }

  char *title, *fname;
  asprintf (&title, "%s, %s data, %s surrogate", tstnm, datatype, surtype);
  char *s = strdup (surtype);
  char *space = strchr (s, ' ');
  if (space) *space = 0;
  for (int n = 0; n < strlen (s); n++)
    s[n] = tolower (s[n]);
  asprintf (&fname, "%s_%s_%s", tstnm, datatype, s);
  plot_repat (sd, title, fname);
  free (title);
  free (fname);
}

static void
do_tests (Edt *edt[SURCNT + 1], char *datatype, char *surtype)
{
  do_test (edt, datatype, surtype, psth, "PSTH");
  do_test (edt, datatype, surtype, isi, "ISI");
  yrange = "set yrange [0:3000]\n";
  do_test (edt, datatype, surtype, autocorr, "AUTOCORR");
  yrange = "set yrange [0:25000]\n";
  do_x_test_orig (edt, datatype, surtype, xcorr, "CCH_v_orig");
  yrange = "";
  do_x_test (edt, datatype, surtype, xcorr, "CCH");
  do_test_2d (edt, datatype, surtype, jint, "JINT");
  do_x_test_2d (edt, datatype, surtype, jpsth, "JPSTH");
  if (0)
    do_repat (edt, datatype, surtype, "REPPAT");
}

int
main (int argc, char **argv)
{
  bool include[1000];
  for (int n = 0; n < 1000; n++) include[n] = true;
  for (int argn = 1; argn < argc; argn++) {
    printf ("%s\n", argv[argn]);
    Edt *edt[SURCNT + 1];
    char *filename;
    if (asprintf (&filename, "%s.edt", argv[argn]) == -1) exit (1);
    edt[0] = read_edt (filename, include);
    free (filename); 
    for (int n = 0; n < SURCNT; n++) {
      if (asprintf (&filename, "%s_%02d.edt", argv[argn], n + 1) == -1) exit (1);
      edt[n + 1] = read_edt (filename, include);
      free (filename); 
    }
    do_tests (edt, argv[argn], "Pauluis and Baker / Miura, Okada, and Amari");
    for (int n = 0; n < SURCNT; n++) {
      free_edt (edt[n + 1]);
      if (asprintf (&filename, "%s_surr_%02d.gdf", argv[argn], n + 1) == -1) exit (1);
      edt[n + 1] = read_edt (filename, include);
      free (filename); 
    }
    do_tests (edt, argv[argn], "Gerstein and Baker");
    for (int n = 0; n < SURCNT + 1; n++)
      free_edt (edt[n]);
  }
  return 0;
}
