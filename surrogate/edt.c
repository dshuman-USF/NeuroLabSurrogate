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
#if HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdbool.h>
#include <error.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include "edt.h"
#if !HAVE_GETLINE && HAVE_CONFIG_H
ssize_t getline(char **linebuf, size_t *linebufsz, FILE *file) ;
#endif

void
free_edt (Edt *edt)
{
  int n;
  for (n = 0; n < edt->digital_count; n++)
    free (edt->digital[n]);
  free (edt->analog);
  free (edt);
}

Edt *
read_edt (char *filename, bool *include)
{
  size_t read;
  char * line = NULL;
  size_t len = 0;
  FILE *f;
  Edt *edt;
  Digital **digital;
  int n, i, c, timelen, codelen, code;

  // if a large .edt file (long experiment) and using NFS mounts, this
  // runs before NFS notices the file exists.  Try a while before giving up
  int giveup=0;
  while(giveup<45) {
    if ((f = fopen (filename, "r")) == NULL) {
      sleep(5);   // 5 seconds
      giveup += 5;
    }
    else
      break;
  }
  if (giveup >= 45)
    error_at_line (1, errno, __FILE__, __LINE__, "can't open %s for read in %s", filename, getcwd (NULL, 0));

  digital = calloc (1000, sizeof *digital);
  edt = calloc (1, sizeof (Edt));
  char *extension = strrchr (filename, '.');
  edt->is_gdt = strcmp (extension, ".gdt") == 0;
  edt->is_gdf = strcmp (extension, ".gdf") == 0;

  if (getline (&line, &len, f));
  n = 0;
  while (isblank (line[n]))
    n++;
  while (!isblank (line[n]))
    n++;
  codelen = n;
  int time = atoi (line + codelen);
  n = strlen (line);
  while (n > 0 && ((c = line[n - 1]) == '\r' || c == '\n'))
    n--;
  line[n] = 0;
  timelen = n - codelen;
  
  edt->header = strdup (line);
  if (asprintf (&edt->format, "%%%dd%%%dd\n", codelen, timelen) == -1) exit (1);
  line[codelen] = 0;
  code = atoi (line);
  if ((edt->is_gdt && code == 21) || edt->is_gdf) {
    edt->has_header = false;
    edt->gdt21time = time;
  }
  else {
    if ((read = getline (&line, &len, f)) == -1) {
      fprintf (stderr, "Empty input file %s, aborting\n", filename);
      exit (1);
    }
    edt->has_header = true;
  }

  time = 0;
  edt->digital_count = 0;
  while ((read = getline (&line, &len, f)) != -1) {
    int code, id;
    time = atoi (line + codelen);
    line[codelen] = 0;
    code = atoi (line);
    if (edt->is_gdt) {
      if (code == 21) {
        edt->gdt21time = time;
        continue;
      }
      else if (code == 22) {
        edt->gdt22time = time;
        break;
      }
    }
    id = code > 4095 ? 1000 + code / 4096 : code;
    if (!include[id])
      continue;
    if (id < 1000) {
      if (digital[id] == NULL) {
        digital[id] = calloc (1, sizeof *digital[id]);
        digital[id]->id = id;
        edt->digital_count++;
      }
      if (digital[id]->spike_count == digital[id]->spike_alloc)
        digital[id]->spike = realloc (digital[id]->spike, (digital[id]->spike_alloc += 100000) * sizeof (int));
      digital[id]->spike[digital[id]->spike_count++] = time;
      //printf ("id %d, spike %d, time %d\n", id, digital[id]->spike_count - 1, time);
    }
    else {
      int analog_idx;
      if (edt->analog_count == edt->analog_alloc)
        edt->analog = realloc (edt->analog, (edt->analog_alloc += 100000) * sizeof (Analog));
      analog_idx = edt->analog_count++;
      edt->analog[analog_idx].code = code;
      edt->analog[analog_idx].time = time;
    }
  }
  edt->analog = realloc (edt->analog, (edt->analog_alloc = edt->analog_count) * sizeof (Analog));
  edt->digital = malloc (edt->digital_count * sizeof (Digital));
  for (i = n = 0; n < 1000; n++)
    if (digital[n] != NULL)
      edt->digital[i++] = digital[n];
  free (digital);
  fclose (f);
  return edt;
}

static int
compare_time (const void *a, const void *b)
{
  const CodeTime *da = (const CodeTime *) a;
  const CodeTime *db = (const CodeTime *) b;

  return (da->time > db->time) - (da->time < db->time);
}


void
write_edt (char *filename, Edt *edt)
{
  FILE *f;
  int count = edt->analog_count, n, spike_idx, id_idx, sample;
  CodeTime *s;

  if ((f = fopen (filename, "w")) == NULL) {
    fprintf (stderr, "Can't open %s for write, aborting: %s\n", filename, strerror (errno));
    exit (1);
  }
  for (n = 0; n < edt->digital_count; n++)
    count += edt->digital[n]->spike_count;
  s = malloc (count * sizeof *s);
  for (n = id_idx = 0; id_idx < edt->digital_count; id_idx++)
    for (spike_idx = 0; spike_idx < edt->digital[id_idx]->spike_count; spike_idx++) {
      s[n].code = edt->digital[id_idx]->id; 
      s[n].time = edt->digital[id_idx]->spike[spike_idx];
      //printf ("id_idx %d, spike %d, id: %d, time %d\n", id_idx, spike_idx, s[n].code, s[n].time);
      n++;
    }
  for (sample = 0; sample < edt->analog_count; sample++) {
    s[n].code = edt->analog[sample].code; 
    s[n].time = edt->analog[sample].time; 
    n++;
  }
  qsort (s, count, sizeof (CodeTime), compare_time);

  if (edt->has_header) {
    fprintf (f, "%s\n", edt->header);
    fprintf (f, "%s\n", edt->header);
  }
  n = 0;

  if (edt->is_gdt) {
    while (n < count && s[n].time <= edt->gdt21time)
      n++;
    if (n == count)
      error (0, 0, "WARNING: defective gdt file. There is a code 21 is at the end.");
    fprintf (f, edt->format, 21, edt->gdt21time);
  }
  else edt->gdt22time = INT_MAX;

  for (; n < count && s[n].time < edt->gdt22time; n++)
    fprintf (f, edt->format, s[n].code, s[n].time);
  if (edt->is_gdt)
    fprintf (f, edt->format, 22, edt->gdt22time);

  free (s);
  fclose (f);
}


void
nodups (Edt *edt)
{
  int n, from, to;

  for (n = 0; n < edt->digital_count; n++) {
    for (to = from = 1; from < edt->digital[n]->spike_count; from++)
      if (edt->digital[n]->spike[from] != edt->digital[n]->spike[from - 1])
        edt->digital[n]->spike[to++] = edt->digital[n]->spike[from];
    edt->digital[n]->spike_count = to;
  }
}
