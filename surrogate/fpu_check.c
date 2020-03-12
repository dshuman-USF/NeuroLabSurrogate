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

#include <stdio.h>
#include <stdlib.h>
#include "fpu_check.h"

void
fpu_check (void)
{
  union {double d; long long i;} d;
  union {float f; int i;} f;

  d.i = 4607802914748677342LL;
  f.i = 1066508979;

  if (d.d == (double)0x123456789abcdeLL/0x10000000000000LL
      && f.f == (double)0x1234566/0x1000000)
    return;
  
  fprintf (stderr, "This program will not work on this machine.\n"
           "This program requires IEEE floating point with the same\n"
           "byte order as integers, and this machine does not have that.\n"
           "Sorry.\n");
}


     
