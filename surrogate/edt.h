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

#ifndef EDT_H
#define EDT_H

typedef struct
{
  int id;
  int spike_count;
  int spike_alloc;
  int *spike;
} Digital;

typedef struct
{
  int code;
  int time;
} Analog;

typedef struct
{
  int digital_count;
  Digital **digital;
  int analog_count;
  int analog_alloc;
  Analog *analog;
  char *header;
  char *format;
  int gdt21time;
  int gdt22time;
  bool is_gdt;
  bool is_gdf;
  bool has_header;
} Edt;

typedef struct
{
  int code;
  int time;
} CodeTime;

void free_edt (Edt *edt);

void write_edt (char *filename, Edt *edt);
Edt * read_edt (char *filename, bool *include);
void nodups (Edt *edt);

#endif
