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

#ifndef NOPIN
#if HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <Judy.h>
#include "unpin.h"
#include "fpu_check.h"

struct Point
{
  double x;
  double y;
  double d;
  int idx;
};
typedef struct Point Point;

typedef union {Word_t w; double d;} Wrd_Dbl;
typedef union {Word_t word; Pvoid_t pvoid;} Word_Pvoid;
typedef union {Word_t w; float d;} Wrd_Flt;

static Point *point;
static Pvoid_t xorder = (Pvoid_t) NULL;
static Pvoid_t dorder = (Pvoid_t) NULL;

static inline double
diff (Point *p)
{
  int Rc_int;
  Point *prev = p, *next = p;
  Word_t index;

  index = (Word_t)p; J1P (Rc_int, xorder, index); prev = (Point *)index;
  index = (Word_t)p; J1N (Rc_int, xorder, index); next = (Point *)index;
  return fabs (p->y - (prev->y + (p->x - prev->x) / (next->x - prev->x) * (next->y - prev->y)));
}

static inline void
pdrop (PWord_t  PValue, Word_Pvoid parray, Point *p)
{
  int Rc_int;
  int pcnt;
  Word_t index;
#if SIZEOF_LONG == 8
  Wrd_Dbl i;
#else
  Wrd_Flt i;
#endif

  index = (Word_t)p; J1U (Rc_int, parray.pvoid, index);
  J1C (pcnt, parray.pvoid, 0, -1);
  if (pcnt == 0) {
    i.d = p->d;
    JLD (Rc_int, dorder, i.w);
  }
  else *PValue = parray.word;
  
}

static inline void
ddrop (Point *p)
{
  PWord_t  PValue;
#if SIZEOF_LONG == 8
  Wrd_Dbl i;
#else
  Wrd_Flt i;
#endif
  Word_Pvoid value_parray;
  
  i.d = p->d; JLG (PValue, dorder, i.w);
  value_parray.word = *PValue;
  pdrop (PValue, value_parray, p);
}

static inline void
dadd (Point *p)
{
  int Rc_int;
  PWord_t  PValue;
  Word_Pvoid value_parray;
  Word_t index;
#if SIZEOF_LONG == 8
  Wrd_Dbl i;
#else
  Wrd_Flt i;
#endif

  i.d = diff (p);
  p->d = i.d;
  JLI (PValue, dorder, i.w);
  value_parray.word = *PValue;

  index = (Word_t)p; J1S (Rc_int, value_parray.pvoid, index);
  *PValue = value_parray.word;
}

int *
unpin (int *spiketimes, int spikecount)
{
  int n, count, count0;
  int Rc_int;
  Word_t index;
  static int *lodif;

  fpu_check ();

  point = malloc ((spikecount + 2) * sizeof *point);
  
  n = 0;
  point[n].x = 0;
  point[n].y = 0;
  point[n].idx = n + 1;
  n++;

  for (n = 1; n <= spikecount; n++) {
    point[n].x = spiketimes[n - 1];
    point[n].y = n - .5;
    point[n].idx = n + 1;
  }
  point[n].x = point[n - 1].x + 1;
  point[n].y = n - 1;
  point[n].idx = n + 1;
  count0 = count = n + 1;

  for (n = 0; n < count; n++) {
    index = (Word_t)&point[n];
    J1S (Rc_int, xorder, index);
  }
  for (n = 1; n < count - 1; n++)
    dadd (point + n);
  
  lodif = malloc (spikecount * sizeof *lodif);
  n = 0;
  while (count > 2) {
    Point *p, *prev, *next;
    PWord_t  PValue;
    Word_Pvoid value_parray;

    index = 0; JLF (PValue, dorder, index);
    value_parray.word = *PValue;

    index = 0; J1F (Rc_int, value_parray.pvoid, index);
    
    p = (Point *)index;
    lodif[n++] = p->idx;

    pdrop (PValue, value_parray, p);

    index = (Word_t)p; J1P (Rc_int, xorder, index);
    
    prev = (Point *)index; index = (Word_t)p; J1N (Rc_int, xorder, index); next = (Point *)index;
    if (prev->idx != 1) ddrop (prev);
    if (next->idx != count0) ddrop (next);
    index = (Word_t)p; J1U (Rc_int, xorder, index);
    if (prev->idx != 1) dadd (prev);
    if (next->idx != count0) dadd (next);

    count--;
  }

  Word_t Rc_word;
  J1FA(Rc_word, xorder);

  PWord_t  PValue;
  Word_t Index = 0;
  JLF (PValue, dorder, Index);
  while (PValue != NULL) {
    Pvoid_t p = (Pvoid_t)*PValue;
    J1FA (Rc_word, p);
    JLN (PValue, dorder, Index);
  }
  JLFA(Rc_word, dorder);

  free (point); point = 0;
  return lodif;
}
#endif
