#include "planet.h"

#include <iostream>
#include <cstdlib>

void Planet::printPlanet() 
{
	int fieldlength = 70;
	std::cout << "planet:  " << planetName << " (ID=" << planetID << ") " << std::endl;
	int linecount = 0;
	//std::cout << dataNr << std::endl;
	for (int i=0; i<dataNr-1; i++) {
		std::string fieldValue = findParam(dataName[i]);
		if (fieldValue!="") {	// only print when not empty
			if (linecount%8==0) { 
				for (int i=0; i<fieldlength; i++) { std::cout << "-"; } std::cout << std::endl;
			}
			printf("  %*s : %*s\n", fieldlength/2-3, dataName[i].c_str(), fieldlength/2-3 , fieldValue.c_str() );
			linecount++;
		}
	}
	for (int i=0; i<fieldlength; i++) { std::cout << "="; }
	std::cout << std::endl << std::endl;
}

// static variables have to be redefined in C++
int Planet::dataNr;
std::string* Planet::dataName;

void Planet::setParam( std::string param, std::string value )
{
	int idx;
	findParam( param, &idx);
	this->myData[idx] = value;
}

std::string Planet::findParam( std::string param, int* index)
{
	for (int i=0; i<dataNr; i++)
	{
		// matching parameter name
		if (dataName[i]==param)
		{
			if ( index!=NULL ) { *index = i; }
			return this->myData[i];
		}
	}
	//oterwise return empty string
	if ( index!=NULL ) { *index = -1; }
	return "";
}

void Planet::initialize()
{
	std::string::size_type sz;     // alias of size_t
	name = planetName;
	symbol = findParam(dataName[5]).c_str();
	orbit_velocity = std::stod (findParam(dataName[4]).c_str(), &sz);
	orbit_period = std::stod (findParam(dataName[3]).c_str(), &sz);
	distance = std::stod (findParam(dataName[1]).c_str(), &sz);
	mass = std::stod (findParam(dataName[2]).c_str(), &sz);
	GtotalMass = 39.4836;
}