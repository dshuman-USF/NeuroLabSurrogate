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


typedef struct
{
  int N;
  int *T;
} SpikeTrain;

typedef struct
{
  int N;
  double *r;
} RateFunc;

SpikeTrain spikes_from_rate (RateFunc r, int mult);
RateFunc recip_isi (SpikeTrain s, double g);
double miura_gamma (SpikeTrain st);
void set_ticks_per_sec (double tr);
void seed_rng (unsigned long seed);
