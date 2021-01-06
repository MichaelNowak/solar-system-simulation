#ifndef ADD_H
#define ADD_H

#include <string>
#include "planet.h"		// adds Planet class (needed for argument of functions)

void countDimensions(std::string filename, int* prows, int* pcolumns, char delimiter);		// count dimensions of .csv data file
void countDimensions(std::string filename, int* prows, int* pcolumns);
double evaluate(int i, double t, double x[4][8], Planet** solar_system, int j);	// evaluate differential equation
void solve(double t, double h, Planet** solar_system, double vars[4][8], const int nVar, int j);		// solve via Runge Kutta; t = last time value; h = time increment
bool is_integer(float k);		// for counting full earth years in simulation
void solveSun(double s[], double x[4][8], Planet** solar_system);		// calculates movement of barycenter /(Sun)
void energy_conservation (Planet** solar_system, double vars[4][8], double stepsize,double dx, double energy[5], double sun_array[2]);  //calculates the energy of the system
#endif
