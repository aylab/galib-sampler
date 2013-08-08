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
ga.cpp contains function to interact with the genetic algorithm library.
Avoid placing I/O functions here and add them to io.cpp instead.
*/

#include "ga.hpp" // Function declarations

#include "galib.cpp"
#include "io.hpp"

void run_ga (input_params& ip) {
	genotype population[ip.population + 1];
	genotype newpopulation[ip.population + 1];
	for (int i = 0; i < ip.population + 1; i++) {
		population[i].initialize(ip.num_dims);
		newpopulation[i].initialize(ip.num_dims);
	}
	initialize("simple_ga_input.txt", ip, population);
	evaluate(ip, population);
	keep_the_best(ip, population);
	for (int generation = 0; generation < ip.generations; generation++) {
		selector(ip, population, newpopulation);
		crossover(ip, population);
		mutate(ip, population);
		report(generation, ip, population);
		evaluate(ip, population);
		elitist(ip, population);
	}
}

