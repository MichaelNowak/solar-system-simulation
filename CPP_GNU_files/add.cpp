#include "add.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <stdio.h>     // for printf
#include <cmath>      // for pow

const double pi = 4 * atan(1.0);
const double G = 4 * pow(pi,2) * 2.99983 * pow(10,-6);		// au^3/(Earth masses * years^2)
const double sunMass = 332948.6; 		// in Earth masses

void countDimensions(std::string filename, int* prows, int* pcolumns) {
	char delimiter='\t';
	if ( filename.substr(filename.length()-3)=="csv" ) {
		delimiter = ',';
	}
	countDimensions(filename, prows, pcolumns, delimiter);
}

void countDimensions(std::string filename, int* prows, int* pcolumns, char delimiter) {
	int columns = 0;
	int rows = 0;
	std::fstream fs (filename.c_str(), std::fstream::in);
	if ( !fs.is_open()) {	// file not opened
		return;
	}

	// count the columns in the first line
	std::string line;
	getline (fs,line);
	std::stringstream strs (line);
	while (strs.peek() == '#') strs.ignore();
	std::string word;

	while ( getline(strs, word, delimiter) ) columns++;

	// count the rows (excluding first line)
	while ( getline(fs, line) ) rows++;
	*prows = rows;
	*pcolumns = columns;
	return;
}

bool is_integer(float k)
{
	return std::floor(k) == k;
}

/* -------------------------------------------------------------------------------------------- */
void solveSun(double s[], double x[4][8], Planet** solar_system) {
	double totalMass = sunMass;
	s[0] = 0.0;	// xBarycenter
	s[1] = 0.0;	// yBarycenter

	for (int j=0; j<8; j++) {
		totalMass += solar_system[j]->get_mass();
	}

	// calculates coordinates of barycenter
	for (int j=0; j<8; j++) {
		s[0] += (solar_system[j]->get_mass()) * x[0][j]/ totalMass;
		s[1] += (solar_system[j]->get_mass()) * x[1][j]/ totalMass;
	}

}

double evaluate(int i, double t, double x[4][8], Planet** solar_system, int j) {

	double ax = 0.0;
	double ay = 0.0;

	// loop through every planet in x[][i]
	for (int i=0; i<8; i++) {
		if (i == j)		// skip if it is the planet that is currently evaluated
			continue;

		// vectorial summation of the acceleration vectors due to gravitational force
		// x[0][i] is the x-value of the ith planet
		// x[1][i] is the y-value of the ith planet
		ax += -G * (solar_system[i]->get_mass()) * (x[0][j] - x[0][i])/ pow(sqrt(pow(x[0][j]-x[0][i],2)+pow(x[1][j]-x[1][i],2)),3);
		ay += -G * (solar_system[i]->get_mass()) * (x[1][j] - x[1][i])/ pow(sqrt(pow(x[0][j]-x[0][i],2)+pow(x[1][j]-x[1][i],2)),3);
	}

	// add sun mass located at origin
	ax += -G * sunMass * (x[0][j] - 0)/ pow(sqrt(pow(x[0][j] - 0,2)+pow(x[1][j] - 0,2)),3);
	ay += -G * sunMass * (x[1][j] - 0)/ pow(sqrt(pow(x[0][j] - 0,2)+pow(x[1][j] - 0,2)),3);

	switch (i) {
		case 0:  return x[2][j];                          // x' = vx
		case 1:  return x[3][j];							// y' = vy
		case 2:  return ax;
		case 3:  return ay;
		default:
		std::cout << "throw?  problem in evaluate" << std::endl;
		return 0.0;
	}
} // end of evaluate()

// ---
// --- solve via Runge Kutta; t = last time value; h = time increment
// ---

void solve(double t, double h, Planet** solar_system, double vars[4][8], const int nVar, int j) {
	int i;
	double inp[4][8];
	double  k1[nVar];
	double  k2[nVar];
	double  k3[nVar];
	double  k4[nVar];
	for (i=0; i<nVar; i++)
		k1[i] = evaluate(i,t,vars, solar_system, j);       // evaluate at time t
//
	for (i=0; i<nVar; i++)
		inp[i][j] = vars[i][j]+k1[i]*h/2;       // set up input to diffeqs
	for (i=0; i<nVar; i++)
		k2[i] = evaluate(i,t+h/2,inp, solar_system, j);    // evaluate at time t+h/2
//
	for (i=0; i<nVar; i++)
		inp[i][j] = vars[i][j]+k2[i]*h/2;       // set up input to diffeqs
	for (i=0; i<nVar; i++)
		k3[i] = evaluate(i,t+h/2,inp, solar_system, j);    // evaluate at time t+h/2
//
	for (i=0; i<nVar; i++)
		inp[i][j] = vars[i][j]+k3[i]*h;         // set up input to diffeqs
	for (i=0; i<nVar; i++)
		k4[i] = evaluate(i,t+h,inp, solar_system, j);      // evaluate at time t+h
	for (i=0; i<nVar; i++) {
		vars[i][j] = vars[i][j]+(k1[i]+2*k2[i]+2*k3[i]+k4[i])*h/6;
	}
} //end of solve

//function to calculate the energy of the system
void energy_conservation (Planet** solar_system, double vars[4][8], double stepsize, double dx, double energy[5], double sun_array[2]) {
	
	for (int i=0; i<5; i++)
		energy[i] = 0.0;
	
	/*
	energy[0]=E_kin;		//total energy of the whole system
	energy[1]=E_pot;		//kinetic energy
	energy[2]=E_total;		// potential energy
	energy[3]=p_tot;		//momentum of the whole system
	energy[4]=L_tot;		//angular momentum of the whole system
	*/


    for (int i=0; i<8; i++) {
            energy[0]+=0.5*solar_system[i]->get_mass()*(pow(vars[2][i],2)+pow(vars[3][i],2));  //energy[0]=0.5*m*v^2
            energy[3]+=solar_system[i]->get_mass()*sqrt(pow(vars[2][i],2)+pow(vars[3][i],2));             //p=m*v
            energy[4]+=sqrt(pow(vars[0][i],2)+pow(vars[1][i],2))*solar_system[i]->get_mass()*sqrt(pow(vars[2][i],2)+pow(vars[3][i],2));  //L=r*m*v
            for (int j=0; j<8; j++) {
                if (i == j)        // skip if it is the planet that is currently evaluated
                	continue;
                energy[1]+=G*solar_system[i]->get_mass()*solar_system[j]->get_mass()/sqrt(pow(vars[0][j]-vars[0][i],2)+pow(vars[1][j]-vars[1][i],2));
                }
            energy[1]+=G*solar_system[i]->get_mass()*sunMass/sqrt(pow(vars[0][i],2)+pow(vars[1][i],2));   //potential energy of the sun
            }

    energy[0]+=0.5*sunMass*pow((dx/stepsize),2); //E_kin of the sun
    energy[3]+=sunMass*(dx/stepsize);             //momentum of the sun
    energy[4]+=sqrt(pow(sun_array[0], 2)+pow(sun_array[1], 2))*sunMass*(dx/stepsize);          //angular momentum of the sun
    energy[2]=energy[0]-energy[1];

}
