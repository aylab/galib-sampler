/*
Genetic algorithm sampler for zebrafish segmentation
Copyright (C) 2013 Ahmet Ay, Jack Holland, Adriana Sperlea, Sebastian Sangervasi

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
init.hpp contains function declarations for init.cpp.
*/

#include "structs.hpp"

#ifndef INIT_HPP
#define INIT_HPP

char* copy_str(const char*);
void init_terminal();
void free_terminal();
void accept_input_params(int, char**, input_params&);
bool option_set(const char*, const char*, const char*);
void ensure_nonempty(const char*, const char*);
void check_input_params(input_params&);
void add_gradient_index(gradient_index**, int);
void init_verbosity(input_params&);
void create_good_sets_file(input_params&);
void init_sim_args(input_params&);
char** copy_args(char**, int);
void read_ranges(input_params&, input_data&);
void store_pipe(char**, int, int);
void delete_files(input_params&);
void reset_cout(input_params&);

#endif

