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
galib.cpp contains the genetic algorithm library taken from http://people.sc.fsu.edu/~jburkardt/cpp_src/simple_ga/simple_ga.html and modified as needed.
*/

#include <cmath> // Needed for sqrt
#include <iomanip> // Needed for setw

#include "io.hpp"

using namespace std;

//
//  Each GENOTYPE is a member of the population, with
//  gene: a string of variables,
//  fitness: the fitness
//  upper: the variable upper bounds,
//  lower: the variable lower bounds,
//  rfitness: the relative fitness,
//  cfitness: the cumulative fitness.
//
struct genotype {
	double* gene;
	double fitness;
	double* upper;
	double* lower;
	double rfitness;
	double cfitness;
	
	genotype () {
		this->gene = NULL;
		this->upper = NULL;
		this->lower = NULL;
	}
	
	void initialize (int num_dims) {
		this->gene = new double[num_dims];
		this->upper = new double[num_dims];
		this->lower = new double[num_dims];
	}
	
	~genotype () {
		delete[] this->gene;
		delete[] this->upper;
		delete[] this->lower;
	}
};

void crossover(input_params&, genotype*);
void elitist(input_params&, genotype*);
void evaluate(input_params&, genotype*);
void initialize(input_params&, genotype*);
void keep_the_best(input_params&, genotype*);
void mutate(input_params&, genotype*);
void r8_swap(double*, double*);
double randval(double, double);
void report(int);
void selector(input_params&, genotype*, genotype*);
void Xover(int, int, input_params&, genotype*);

void crossover (input_params& ip, genotype* population) {
  int mem;
  int one = 0;
  int first = 0;
  double x;

  for ( mem = 0; mem < ip.population; ++mem )
  {
    x = ( rand ( ) % 1000 ) / 1000.0;

    if ( x < ip.prob_crossover )
    {
      ++first;

      if ( first % 2 == 0 )
      {
        Xover (one, mem, ip, population);
      }
      else
      {
        one = mem;
      }

    }
  }
  return;
}

void elitist (input_params& ip, genotype* population) {
  int i;
  double best;
  int best_mem = 0;
  double worst;
  int worst_mem = 0;

  best = population[0].fitness;
  worst = population[0].fitness;

  for ( i = 0; i < ip.population - 1; ++i )
  {
    if ( population[i].fitness > population[i+1].fitness )
    {

      if ( best <= population[i].fitness )
      {
        best = population[i].fitness;
        best_mem = i;
      }

      if ( population[i+1].fitness <= worst )
      {
        worst = population[i+1].fitness;
        worst_mem = i + 1;
      }

    }
    else
    {

      if ( population[i].fitness <= worst )
      {
        worst = population[i].fitness;
        worst_mem = i;
      }

      if ( best <= population[i+1].fitness )
      {
        best = population[i+1].fitness;
        best_mem = i + 1;
      }

    }

  }
// 
//  If the best individual from the new population is better than 
//  the best individual from the previous population, then 
//  copy the best from the new population; else replace the 
//  worst individual from the current population with the 
//  best one from the previous generation                     
//
  if ( best >= population[ip.population].fitness )
  {
    for ( i = 0; i < ip.num_dims; i++ )
    {
      population[ip.population].gene[i] = population[best_mem].gene[i];
    }
    population[ip.population].fitness = population[best_mem].fitness;
  }
  else
  {
    for ( i = 0; i < ip.num_dims; i++ )
    {
      population[worst_mem].gene[i] = population[ip.population].gene[i];
    }
    population[worst_mem].fitness = population[ip.population].fitness;
  } 

  return;

}

void evaluate (input_params& ip, genotype* population) {
  int member;
  int i;
  double x[ip.num_dims];

  for ( member = 0; member < ip.population; member++ )
  {
    for ( i = 0; i < ip.num_dims; i++ )
    {
      x[i] = population[member].gene[i];
    } 
    population[member].fitness = simulate_set(ip, x);
  }
  return;
}

void initialize (input_params& ip, genotype* population) {
  int i;
  int j;
  double lbound;
  double ubound;

  for ( i = 0; i < ip.num_dims; i++ )
  {
    lbound = ip.ranges[i].first;
    ubound = ip.ranges[i].second;

    for ( j = 0; j < ip.population; j++ )
    {
      population[j].fitness = 0;
      population[j].rfitness = 0;
      population[j].cfitness = 0;
      population[j].lower[i] = lbound;
      population[j].upper[i]= ubound;
      population[j].gene[i] = randval ( population[j].lower[i], population[j].upper[i] );
    }
  }

  return;
}

void keep_the_best (input_params& ip, genotype* population) {
  int cur_best;
  int mem;
  int i;

  cur_best = 0;

  for ( mem = 0; mem < ip.population; mem++ )
  {
    if ( population[mem].fitness > population[ip.population].fitness )
    {
      cur_best = mem;
      population[ip.population].fitness = population[mem].fitness;
    }
  }
// 
//  Once the best member in the population is found, copy the genes.
//
  for ( i = 0; i < ip.num_dims; i++ )
  {
    population[ip.population].gene[i] = population[cur_best].gene[i];
  }

  return;
}

void mutate (input_params& ip, genotype* population) {
  double hbound;
  int i;
  int j;
  double lbound;
  double x;

  for ( i = 0; i < ip.population; i++ )
  {
    for ( j = 0; j < ip.num_dims; j++ )
    {
      x = rand ( ) % 1000 / 1000.0;
// 
//  Find the bounds on the variable to be mutated 
//
      if ( x < ip.prob_mutation )
      {
        lbound = population[i].lower[j];
        hbound = population[i].upper[j];  
        population[i].gene[j] = randval ( lbound, hbound );
      }
    }
  }

  return;
}

void r8_swap ( double *x, double *y ) {
  double temp;

  temp = *x;
  *x = *y;
  *y = temp;

  return;
}

double randval ( double low, double high ) {
  double val;

  val = ( ( double ) ( rand() % 1000 ) / 1000.0 ) 
    * ( high - low ) + low;

  return ( val );
}

void report (int generation, input_params& ip, genotype* population) {
  double avg;
  double best_val;
  int i;
  double square_sum;
  double stddev;
  double sum;
  double sum_square;

  sum = 0.0;
  sum_square = 0.0;

  for ( i = 0; i < ip.population; i++ )
  {
    sum = sum + population[i].fitness;
    sum_square = sum_square + population[i].fitness * population[i].fitness;
  }

  avg = sum / ( double ) ip.population;
  square_sum = avg * avg * ip.population;
  stddev = sqrt ( ( sum_square - square_sum ) / ( ip.population - 1 ) );
  best_val = population[ip.population].fitness;

  cout << "  " << setw(8) << generation 
       << " " << best_val 
       << " " << avg 
       << " " << stddev << "\n";

  return;
}

void selector (input_params& ip, genotype* population, genotype* newpopulation) {
  int i;
  int j;
  int mem;
  double p;
  double sum = 0;
//
//  Find total fitness of the population 
//
  for ( mem = 0; mem < ip.population; mem++ )
  {
    sum = sum + population[mem].fitness;
  }
//
//  Calculate the relative fitness.
//
  for ( mem = 0; mem < ip.population; mem++ )
  {
    population[mem].rfitness = population[mem].fitness / sum;
  }
  population[0].cfitness = population[0].rfitness;
// 
//  Calculate the cumulative fitness.
//
  for ( mem = 1; mem < ip.population; mem++ )
  {
    population[mem].cfitness = population[mem-1].cfitness +       
      population[mem].rfitness;
  }
// 
//  Select survivors using cumulative fitness. 
//
  for ( i = 0; i < ip.population; i++ )
  { 
    p = rand() % 1000 / 1000.0;
    if (p < population[0].cfitness)
    {
      newpopulation[i] = population[0];      
    }
    else
    {
      for ( j = 0; j < ip.population; j++ )
      { 
        if ( p >= population[j].cfitness && p < population[j+1].cfitness )
        {
          newpopulation[i] = population[j+1];
        }
      }
    }
  }
// 
//  Once a new population is created, copy it back 
//
  for ( i = 0; i < ip.population; i++ )
  {
    population[i] = newpopulation[i]; 
  }

  return;     
}

void Xover (int one, int two, input_params& ip, genotype* population) {
  int i;
  int point;
// 
//  Select the crossover point.
//
  if ( 1 < ip.num_dims )
  {

    if ( ip.num_dims == 2 )
    {
      point = 1;
    }
    else
    {
      point = ( rand ( ) % ( ip.num_dims - 1 ) ) + 1;
    }

    for ( i = 0; i < point; i++ )
    {
      r8_swap ( &population[one].gene[i], &population[two].gene[i] );
    }

  }
  return;
}

