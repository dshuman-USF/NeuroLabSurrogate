#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <ctype.h>
#include <limits.h>
#include "edt.h"

#define WLEN 150
#define MAXPLEX 16
#define MAXREPPAT 20

typedef struct
{
  int isi;
  int id_idx;
  int end_spike;
  int start_time;
  int end_time;
  int index;
} Point;
typedef int RepRow[MAXREPPAT + 1];

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

static void
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

  for (int shift = 0; shift <= maxt; shift++) {
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
}

int
main (int argc, char **argv)
{
  bool include[1000];
  for (int n = 0; n < 1000; n++) include[n] = true;

  Edt *edt = read_edt ("rpindat.edt", include);
  repat (edt);

  return 0;
}
