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
int * gen_control (int *spiketime, int spikecount, double g, int *control_spikecount, unsigned long seed, bool use_seed);
void insert (double index, double value);
double interp1 (double y);
int *gen_control_from_rate (double end, double g, int *control_spikecount, unsigned long seed, bool use_seed);

void init_control (void) ;

void put_t_est (Pvoid_t t_est_val);
void clear_t_est (void);
void zero_t_est (void);
void set_use_smoothed_rate (bool val);
