#define _GNU_SOURCE
#if HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include "read_array.h"
#if !HAVE_ASPRINTF && HAVE_CONFIG_H
int asprintf(char **buffer, char *fmt, ...);
#endif

int
read_array (char *name, char *unit, int **arrayp)
{
  FILE *f;
  int n, alloc = 0;
  char *filename;
  int t;

  if (asprintf (&filename, "%s%s", name, unit) == -1) exit (1);
  if ((f = fopen (filename, "r")) == NULL) {
    fprintf (stderr, "Can't open %s: %s\n", filename, strerror (errno));
    return 1;
  }
  free (filename);

  n = 0;
  while (fscanf (f, "%d", &t) == 1) {
    if (n + 1 >= alloc)
      (*arrayp) = realloc ((*arrayp), (alloc += 16384) * sizeof *(*arrayp));
    (*arrayp)[n++] = t;
  }
  (*arrayp) = realloc ((*arrayp), n * sizeof *(*arrayp));
  return n;
}

