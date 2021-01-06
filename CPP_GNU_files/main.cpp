#include <iostream>      // standard IO
#include <cmath>	       // math
#include <fstream>		// file streams
#include <stdio.h>		// printf
#include "add.h"		// add functions
#include "planet.h"		// add Planet class

using namespace std;

const int       nVar  = 4;                // number of dynamical variables
const int nTimeSteps  = 165000;             // number of time steps
const double       h  = 0.001;             // Delta t
const double pi = 4 * atan(1.0);		// definition of pi
const double G = 4 * pow(pi,2) * 2.99983 * pow(10,-6);		// Newton's gravitational constant au^3/(Earth masses * years^2)
// total mass multiplied with G can be approximated with G_m1_plus_m2 = 4 * pi * pi
//const double G_m1_plus_m2 = 4 * pi * pi;  	// deviation from total mass multiplied with G: 4*pow(pi,2) - 39.4836 = -0.00513704 (0.0130106 %)
const double sunMass = 332948.6; 		// in Earth masses
const double startTime = 623467.001; 		// where barycenter is closest to (0,0); due to distinct angular velocities on a fixed orbit a startTime is sufficient for a certain planetary start configuration (no need to give random position for each planet)

int main()
{
	string fileName = "planet_properties_au.csv";

	ifstream data(fileName);	// filestream for input data
	string head, line, word;
	int columnNr, rowNr, i;

	countDimensions(fileName, &rowNr, &columnNr);
	cout << fileName << " => rows: "<< rowNr << " columns: " << columnNr << endl;

	// static variable initialization
	Planet::dataNr = columnNr;
	Planet::dataName = new string[columnNr];
	int planetNr = rowNr;

	Planet** solar_system = new Planet*[planetNr];
	//Planet* solar_system[rowNr];	// equivalent to the previous one

   // read the table header and store to dataName
	getline(data,head);
	stringstream strs(head);
	if(strs.peek()=='#') { strs.ignore(2); }    // ignore 3 to ignore the spaces after '#', too

	i=0;
	while( getline(strs, word, ',') )
	{
		Planet::dataName[i]=word;
		cout << fileName << " => header: " << Planet::dataName[i] << endl;
		i++;
	}

    // store planets in solar_system array
	i=0;
	while( getline(data,line) )	// store all planets and values
	{
		solar_system[i] = new Planet( line );
		i++;
	}

	// print planets' properties
	cout << endl;
	for (int i=0; i<8; i++)
		solar_system[i]->printPlanet();
	
	
	// not necessary but good for testing

	/*
	for (int i=0; i<8; i++) {
		solar_system[i]->initialize();
		for (int i=0; i<62; i++) { cout << "-"; } cout << endl;
		cout << solar_system[i]->get_name() << "'s " << Planet::dataName[1] << " is: " << solar_system[i]->get_distance() << endl;
		cout << solar_system[i]->get_name() << "'s " << Planet::dataName[2] << " is: " << solar_system[i]->get_mass() << endl;
		cout << solar_system[i]->get_name() << "'s " << Planet::dataName[3] << " is: " << solar_system[i]->get_orbit_period() << endl;
		cout << solar_system[i]->get_name() << "'s angular velocity (2Pi rad/yr) is: " << solar_system[i]->angular_velocity() << endl;
		cout << solar_system[i]->get_name() << "'s orbital velocity is: " << solar_system[i]->get_orbit_velocity() << endl;
	}
	*/

	/* ------------------------------------------------------------------------------- */

	// initialization to access class Planet member functions
	for (int i=0; i<8; i++) {
		solar_system[i]->initialize();
	}

	// Numerical integration via Runge-Kutta

	double vars[nVar][8];                        // the dynamical variables
	double sunCoord[2];						// array for (x,y) values of the sun
	double sunCoord_old[2] = {0.0};

	double totalMass = sunMass;				// starts with Sun mass
	double energy[5];
	double dx;

	for (int j=0; j<8; j++) {
		totalMass += solar_system[j]->get_mass(); 		// adds mass of every planet
	}

	// initial conditions

	for (int j=0; j<8; j++) {
		sunCoord_old[0] += (solar_system[j]->get_mass()) * (solar_system[j]->calculated_distance()) * cos(2*pi* (solar_system[j]->angular_velocity()) * startTime)/ totalMass;
		sunCoord_old[1] += (solar_system[j]->get_mass()) * (solar_system[j]->calculated_distance()) * sin(2*pi* (solar_system[j]->angular_velocity()) * startTime)/ totalMass;
	}

	//sunCoord_old[0] += sunMass * 0 / totalMass;
	//sunCoord_old[1] += sunMass * 0 / totalMass;
	
			//sunCoord_old[0] = xBarycenter;
			//sunCoord_old[1] = yBarycenter;

	double rBarycenter = sqrt(pow(sunCoord_old[0], 2) + pow(sunCoord_old[1], 2));

	cout << endl;
	cout << "Barycenter at (" << sunCoord_old[0] << " A.U."<< "," << sunCoord_old[1] << " A.U.) -> distance from origin: " << rBarycenter << " A.U" << endl;
	cout << endl;

	double varsCopy[4][8];
	double radius[8];
	double orbitVelocity[8];
	double orbitVelocityCopy[8];
	
	ofstream initial;
	initial.open("initial.txt",ios::out);		// writes statistical data in output
	initial.setf(ios::fixed,ios::floatfield);  // for formated output
	initial.precision(10);                      // floating point precision
	
	initial << endl;
	initial << "Barycenter at (" << sunCoord_old[0] << " A.U."<< "," << sunCoord_old[1] << " A.U.) -> distance from origin: " << rBarycenter << " A.U" << endl;
	initial << endl;
	initial << "nTimeSteps: " << nTimeSteps << endl;
	initial << "h: " << h << endl;
	initial << endl;
	
	for (int j=0; j<8; j++) {

		vars[0][j] = (solar_system[j]->calculated_distance()) * cos(2*pi*(solar_system[j]->angular_velocity())* startTime );
		vars[1][j] = (solar_system[j]->calculated_distance()) * sin(2*pi*(solar_system[j]->angular_velocity())* startTime );

		double xDistBarycenter = vars[0][j] - sunCoord_old[0];
		double yDistBarycenter = vars[1][j] - sunCoord_old[1];
		double distBarycenter = sqrt(pow(xDistBarycenter, 2) + pow(yDistBarycenter, 2));
		double phiDistBarycenter = atan2(yDistBarycenter, xDistBarycenter);

		vars[2][j] = sqrt(G * totalMass / distBarycenter) * cos(phiDistBarycenter + pi/2);
		vars[3][j] = sqrt(G * totalMass / distBarycenter) * sin(phiDistBarycenter + pi/2);
		
		varsCopy[0][j] = vars[0][j];
		varsCopy[1][j] = vars[1][j];
		varsCopy[2][j] = vars[2][j];
		varsCopy[3][j] = vars[3][j];
				
		radius[j] = solar_system[j]->calculated_distance();
		orbitVelocity[j] = sqrt(pow(vars[2][j],2) + pow(vars[3][j], 2));
		orbitVelocityCopy[j] = orbitVelocity[j];

		// formatted Output of initial conditions compared to official data
		
		initial << endl;
		for (int i=0; i<50; i++)
			initial << "-";
		initial << endl;
		initial << "Initial conditions: " << solar_system[j]->get_name() << endl;
		initial << endl;
		initial << "Mean orbit velocity from NASA: " << solar_system[j]->get_orbit_velocity() << " A.U/yr" << endl;
		initial << "Calculated orbit velocity: " << orbitVelocity[j] << " A.U/yr" << endl;
		initial << "Absolute deviation of orbit velocity: " << orbitVelocity[j] - solar_system[j]->get_orbit_velocity() << " A.U/yr" << endl;
		initial << "Relative deviation of orbit velocity: " << 1 - orbitVelocity[j]/ solar_system[j]->get_orbit_velocity() << endl;
		initial << endl;
		initial << "Mean orbit period from NASA: " << solar_system[j]->get_orbit_period() << " yr" << endl;
		initial << "Calculated orbit period: " << 2*pi*solar_system[j]->calculated_distance()/orbitVelocity[j] << " yr" << endl;
		initial << "Absolute deviation of orbit period: " << 2*pi*solar_system[j]->calculated_distance()/orbitVelocity[j] - solar_system[j]->get_orbit_period() << " yr" << endl;
		initial << "Relative deviation of orbit period: " << 1 - (2*pi*solar_system[j]->calculated_distance()/orbitVelocity[j])/ solar_system[j]->get_orbit_period() << endl;
		initial << endl;
		initial << "Mean distance to Sun from NASA: " << solar_system[j]->get_distance() << " A.U" << endl;
		initial << "Calculated distance to Sun: " << solar_system[j]->calculated_distance() << " A.U" << endl;
		initial << "Absolute deviation of calculated distance: " << solar_system[j]->calculated_distance() - solar_system[j]->get_distance() << " A.U" << endl;
		initial << "Relative deviation of calculated distance: " << 1 - solar_system[j]->calculated_distance()/ solar_system[j]->get_distance() << endl;
		initial << endl;

	}
	
	initial.close();

	// -----------------------------------------------------------------------------
	/* ############################################################ */
	// For accurate mean values to calculate fluctuations around mean
	/* ############################################################ */
	
	// Simulate once to get mean values of radius and orbital velocity
	
	/*
	for (int iTimes=0;iTimes<nTimeSteps;iTimes++) {
		double t = h*iTimes;

		for (int j=0; j<8; j++) {
			solar_system[j]->initialize();
			solve(t, h, solar_system, varsCopy, nVar, j);		// solves planetary movement via Runge Kutta
			radius[j] += sqrt(pow(varsCopy[0][j],2)+pow(varsCopy[1][j],2));
			orbitVelocity[j] += sqrt(pow(varsCopy[2][j],2)+pow(varsCopy[3][j],2));
		}
		
		cout << "t: " << t << endl;
	}
	*/
	
	// -----------------------------------------------------------------------------
	/* ############################################################ */
	// For faster simulation with calculated initial values to calculate fluctuations around mean (sufficient accurate)
	/* ############################################################ */
	
	for (int i=0; i<8; i++) {
		radius[i] *= nTimeSteps;
		orbitVelocity[i] *= nTimeSteps;
	}
	// -----------------------------------------------------------------------------
	
	
	// -----------------------------------------------------------------------------
	
	ofstream statistics;
	statistics.open("statistics.txt",ios::out);		// writes statistical data in output
	statistics.setf(ios::fixed,ios::floatfield);  // for formated output
	statistics.precision(10);                      // floating point precision
	
	// Print statistical comparisons
	
	statistics << endl;
	statistics << "nTimeSteps: " << nTimeSteps << endl;
	statistics << "h: " << h << endl;
	statistics << endl;
	
	for (int i=0; i<8; i++) {
		statistics << endl;
		for (int i=0; i<50; i++)
			statistics << "-";
		statistics << endl;
		
		statistics << "After simulation via Runge Kutta: " << solar_system[i]->get_name() << endl;
		statistics << endl;
		statistics << "Mean orbit velocity from NASA: " << solar_system[i]->get_orbit_velocity() << " A.U/yr" << endl;
		statistics << "Calculated orbit velocity initial condition: " << orbitVelocityCopy[i] << " A.U/yr" << endl;
		statistics << "Simulated mean orbit velocity: " << orbitVelocity[i]/nTimeSteps << " A.U/yr" << endl;
		statistics << "Absolute deviation of orbit velocity with NASA: " << orbitVelocity[i]/nTimeSteps - solar_system[i]->get_orbit_velocity() << " A.U/yr" << endl;
		statistics << "Relative deviation of orbit velocity with NASA: " << 1 - (orbitVelocity[i]/nTimeSteps)/ solar_system[i]->get_orbit_velocity() << endl;
		statistics << "Absolute deviation of orbit velocity with initial condition: " << orbitVelocity[i]/nTimeSteps - orbitVelocityCopy[i] << " A.U/yr" << endl;
		statistics << "Relative deviation of orbit velocity with initial condition: " << 1 - (orbitVelocity[i]/nTimeSteps)/ orbitVelocityCopy[i] << endl;
		
		statistics << endl;
		
		statistics << "Mean orbit period from NASA: " << solar_system[i]->get_orbit_period() << " yr" << endl;
		statistics << "Calculated orbit period initial condition: " << 2*pi*solar_system[i]->calculated_distance()/orbitVelocityCopy[i] << " yr" << endl;
		statistics << "Simulated mean orbit period: " << (2*pi*radius[i]/nTimeSteps)/(orbitVelocity[i]/nTimeSteps) << " yr" << endl;
		statistics << "Absolute deviation of orbit period with NASA: " << (2*pi*radius[i]/nTimeSteps)/(orbitVelocity[i]/nTimeSteps) - solar_system[i]->get_orbit_period() << " yr" << endl;
		statistics << "Relative deviation of orbit period with NASA: " << 1 - ((2*pi*radius[i]/nTimeSteps)/(orbitVelocity[i]/nTimeSteps))/ solar_system[i]->get_orbit_period() << endl;
		statistics << "Absolute deviation of orbit period with initial condition: " << ((2*pi*radius[i]/nTimeSteps)/(orbitVelocity[i]/nTimeSteps)) - 2*pi*solar_system[i]->calculated_distance()/orbitVelocityCopy[i] << " yr" << endl;
		statistics << "Relative deviation of orbit period with initial condition: " << 1 - ((2*pi*radius[i]/nTimeSteps)/(orbitVelocity[i]/nTimeSteps))/ (2*pi*solar_system[i]->calculated_distance()/orbitVelocityCopy[i]) << endl;
		
		statistics << endl;
		
		statistics << "Mean distance to Sun from NASA: " << solar_system[i]->get_distance() << " A.U" << endl;
		statistics << "Calculated distance initial condition: " << solar_system[i]->calculated_distance() << " A.U" << endl;
		statistics << "Simulated mean distance to Sun: " << radius[i]/nTimeSteps << " A.U" << endl;
		statistics << "Absolute deviation of simulated mean distance with NASA: " << radius[i]/nTimeSteps - solar_system[i]->get_distance() << " A.U" << endl;
		statistics << "Relative deviation of distance with NASA: " << 1 - (radius[i]/nTimeSteps)/ solar_system[i]->get_distance() << endl;
		statistics << "Absolute deviation of distance with initial condition: " << radius[i]/nTimeSteps - solar_system[i]->calculated_distance() << " A.U" << endl;
		statistics << "Relative deviation of distance with initial condition: " << 1 - (radius[i]/nTimeSteps)/ solar_system[i]->calculated_distance() << endl;
		
		statistics << endl;
	}
	
	statistics.close();
	
	// Write data for planets in .txt file

	ofstream Sun, Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptun, Energy_cons;

	Sun.open("SunOutput.txt",ios::out);		// writes data for Sun sepatately
	Sun.setf(ios::fixed,ios::floatfield);  // for formated output
	Sun.precision(10);                      // floating point precision

	Mercury.open("MercuryOutput.txt",ios::out);		// writes data for Mercury sepatately
	Mercury.setf(ios::fixed,ios::floatfield);  // for formated output
	Mercury.precision(10);                      // floating point precision

	Venus.open("VenusOutput.txt",ios::out);		// writes data for Venus sepatately
	Venus.setf(ios::fixed,ios::floatfield);  // for formated output
	Venus.precision(10);                      // floating point precision

	Earth.open("EarthOutput.txt",ios::out);		// writes data for Earth sepatately
	Earth.setf(ios::fixed,ios::floatfield);  // for formated output
	Earth.precision(10);                      // floating point precision

	Mars.open("MarsOutput.txt",ios::out);		// writes data for Mars sepatately
	Mars.setf(ios::fixed,ios::floatfield);  // for formated output
	Mars.precision(10);                      // floating point precision

	Jupiter.open("JupiterOutput.txt",ios::out);		// writes data for Jupiter sepatately
	Jupiter.setf(ios::fixed,ios::floatfield);  // for formated output
	Jupiter.precision(10);                      // floating point precision

	Saturn.open("SaturnOutput.txt",ios::out);		// writes data for Saturn sepatately
	Saturn.setf(ios::fixed,ios::floatfield);  // for formated output
	Saturn.precision(10);                      // floating point precision

	Uranus.open("UranusOutput.txt",ios::out);		// writes data for Uranus sepatately
	Uranus.setf(ios::fixed,ios::floatfield);  // for formated output
	Uranus.precision(10);                      // floating point precision

	Neptun.open("NeptunOutput.txt",ios::out);		// writes data for Neptun sepatately
	Neptun.setf(ios::fixed,ios::floatfield);  // for formated output
	Neptun.precision(10);                      // floating point precision

	Energy_cons.open("Energy_cons.txt",ios::out);		// writes data for Neptun sepatately
	Energy_cons.setf(ios::fixed,ios::floatfield);  // for formated output
	Energy_cons.precision(10);                     // floating point precision
	
	// simualte second time to get all data
	
	int years = 0;
	for (int iTimes=0;iTimes<nTimeSteps;iTimes++) {
		double t = h*iTimes;
		
		cout << "t: " << t << endl;
		
		if (is_integer(t) == true)	// counts earth years in simulation
			++years;

		for (int j=0; j<8; j++) {
			solar_system[j]->initialize();
			solve(t, h, solar_system, vars, nVar, j);		// solves planetary movement via Runge Kutta
		}
		
		/* calculate movement of barycenter with Sun located at origin of the cartesian coordinate system; changing the origin to the barycenter leads to description of the Sun's movement (relative to barycenter);
		planets are actually orbiting the barycenter but since it is very close to the Sun's core it is sufficient to calculate with an orbit around the sun (located at origin) */
		
		solveSun(sunCoord, vars, solar_system);		// calculates position of barycenter (Sun)
		
		dx = sqrt(pow(sunCoord_old[0]-sunCoord[0],2)+pow(sunCoord_old[1]-sunCoord[1],2)); //distance, sun moved in time h
		
		energy_conservation (solar_system, vars, h, dx, energy, sunCoord_old);
		
		for (int i=0; i<2; i++) {
			sunCoord_old[i]=sunCoord[i];
		}

		//if (is_integer(5*t) == false)		// skip writing steps --> now write every 0.02
			//continue;
		
		// write in Ouput file: years, angle, radius, x-coordinate, y-coordinate, deviation from radius mean, deviation from velocity mean
		int j = 0;
		// write in Ouput file // atan2(double,double) considers the angle with respect to the quadrant
		Mercury << t << " " << atan2(vars[1][j], vars[0][j]) << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) << " " << vars[0][j] << " " << vars[1][j] << " " << vars[2][j] << " " << vars[3][j] << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) - radius[j]/nTimeSteps << " " << sqrt(pow(vars[2][j],2)+pow(vars[3][j],2)) - orbitVelocity[j]/nTimeSteps << endl;

		j = 1;
		// write in Ouput file
		Venus << t << " " << atan2(vars[1][j], vars[0][j]) << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) << " " << vars[0][j] << " " << vars[1][j] << " " << vars[2][j] << " " << vars[3][j] << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) - radius[j]/nTimeSteps << " " << sqrt(pow(vars[2][j],2)+pow(vars[3][j],2)) - orbitVelocity[j]/nTimeSteps << endl;

		j = 2;
		// write in Ouput file
		Earth << t << " " << atan2(vars[1][j], vars[0][j]) << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) << " " << vars[0][j] << " " << vars[1][j] << " " << vars[2][j] << " " << vars[3][j] << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) - radius[j]/nTimeSteps << " " << sqrt(pow(vars[2][j],2)+pow(vars[3][j],2)) - orbitVelocity[j]/nTimeSteps << endl;

		j = 3;
		// write in Ouput file
		Mars << t << " " << atan2(vars[1][j], vars[0][j]) << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) << " " << vars[0][j] << " " << vars[1][j] << " " << vars[2][j] << " " << vars[3][j] << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) - radius[j]/nTimeSteps << " " << sqrt(pow(vars[2][j],2)+pow(vars[3][j],2)) - orbitVelocity[j]/nTimeSteps << endl;

		j = 4;
		// write in Ouput file
		Jupiter << t << " " << atan2(vars[1][j], vars[0][j]) << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) << " " << vars[0][j] << " " << vars[1][j] << " " << vars[2][j] << " " << vars[3][j] << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) - radius[j]/nTimeSteps << " " << sqrt(pow(vars[2][j],2)+pow(vars[3][j],2)) - orbitVelocity[j]/nTimeSteps << endl;

		j = 5;
		// write in Ouput file
		Saturn << t << " " << atan2(vars[1][j], vars[0][j]) << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) << " " << vars[0][j] << " " << vars[1][j] << " " << vars[2][j] << " " << vars[3][j] << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) - radius[j]/nTimeSteps << " " << sqrt(pow(vars[2][j],2)+pow(vars[3][j],2)) - orbitVelocity[j]/nTimeSteps << endl;

		j = 6;
		// write in Ouput file
		Uranus << t << " " << atan2(vars[1][j], vars[0][j]) << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) << " " << vars[0][j] << " " << vars[1][j] << " " << vars[2][j] << " " << vars[3][j] << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) - radius[j]/nTimeSteps << " " << sqrt(pow(vars[2][j],2)+pow(vars[3][j],2)) - orbitVelocity[j]/nTimeSteps << endl;

		j = 7;
		// write in Ouput file
		Neptun << t << " " << atan2(vars[1][j], vars[0][j]) << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) << " " << vars[0][j] << " " << vars[1][j] << " " << vars[2][j] << " " << vars[3][j] << " " << sqrt(pow(vars[0][j],2)+pow(vars[1][j],2)) - radius[j]/nTimeSteps << " " << sqrt(pow(vars[2][j],2)+pow(vars[3][j],2)) - orbitVelocity[j]/nTimeSteps << endl;
			
					
		//write in Output file for Sun: years, angle, radius, x-coordinate, y-coordinate
		Sun << t << " " << atan2(sunCoord[1], sunCoord[0]) << " " << sqrt(pow(sunCoord[0],2)+pow(sunCoord[1],2)) << " " << sunCoord[0] << " " << sunCoord[1] << endl;
	                
	    //write in Output file for energy conservation: E_kin, E_pot, E_total, p_tot, L_tot
		Energy_cons << t <<"  "<<energy[0]<<"  "<<energy[1]<<"  "<<energy[2]<<"  "<<energy[3]<<"  "<<energy[4] << endl;
			
	}
	
	// Close Output writing
	Sun.close();
	Mercury.close();
	Venus.close();
	Earth.close();
	Mars.close();
	Jupiter.close();
	Saturn.close();
	Uranus.close();
	Neptun.close();
	Energy_cons.close();

	/* ------------------------------------------------------------------------------- */

	// Call destructor
	// deletes every array that was dynamically created via 'new' ( allocation of new memory )

	for (int i=0; i<8; i++) {
		delete solar_system[i];
	}

	delete[] solar_system;
	delete[] Planet::dataName;

	return 0;
}  // end of main()
