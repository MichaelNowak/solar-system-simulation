#ifndef PLANET_H
#define PLANET_H

#include <string>
#include <sstream>
#include <stdio.h>

static int planet_num = 0;    // initialization of static variable

class Planet {
	public:
		const int planetID = planet_num++;	   // planet ID number
		std::string planetName;                         // name of the planet
		
		// Constructor and Destructor
		Planet() {}		// default constructor
		Planet( std::string val ) {
			this->myData = new std::string[dataNr];
			std::string value;
		    std::stringstream data_stream(val);
		    int i=0;
			while ( getline( data_stream, value, ',' ) and i<dataNr) {
				this->myData[i] = value;
			   	i++;
			}

			this->planetName = findParam("Name");
			printf("Constructor nr. %d called, name is %s\n", planetID, (this->planetName).c_str());
		}
		~Planet() { printf("\nPlanet '%d' is destroyed\n", planetID);}
		
		void printPlanet( void );              // prints selected data of a planet
		
		// setter and getter of any parameter
		std::string findParam( std::string param, int* index=NULL );
		void setParam( std::string key, std::string value );
		
		static int dataNr;                     // number of different data types
		static std::string* dataName;               // static pointer to the an of strings for the data labels
		
		//void initialize(double in_distance, double in_mass, double in_orbit_velocity);
		void initialize();
		std::string get_name(void) {return name;}
		std::string get_symbol(void) {return symbol;}
		double get_distance(void) {return distance;}
		double get_mass(void) {return mass;}
		double get_orbit_period(void) {return orbit_period;}
		double get_orbit_velocity(void) {return orbit_velocity;}
		double angular_velocity(void) {return 1./ orbit_period;}
		double calculated_distance(void) {return GtotalMass / (orbit_velocity * orbit_velocity);}
		
	protected:
		std::string name;
		std::string symbol;
		double distance;		// distance from the sun
		double mass;
		double orbit_period;
		double orbit_velocity;
		double GtotalMass;
		
	private:
		void setProperties( std::string data );     // to initialize all propteries of the planet
		std::string* myData;                        // pointer to an array of strings for storing the data
};

#endif