#ifndef _SUN_H_
#define _SUN_H_

#include "Clock.hpp"
#include "Sunray.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;

class Sun
{
private:
	//Using the horizontal coordinates system
	//Coordinates of the location for which the sun is computed
	double _latitude;

	//Position of the sun through the day - in degrees
	double _elevation;
	double _azimuth;
	double _sunDeclination;

	double _raysPerSquareMeters;

	Clock _time;
	int _day;
	vector<Sunray*> _listOfSunrays;


public:
	//Construct / Destruct
	Sun(double, Clock, int, double);
	Sun(const Sun&);
	~Sun();

	//-----------------------------------------------------------------------
	void fComputeSunPosition();
	void fAddNewSunray(double, double, double);
	void fUpdateSunrays();

	void set_time(int); //minutes

	//get methods
	double& get_latitude()  { return _latitude; }
	double& get_elevation() { return _elevation; }
	double& get_azimuth()   { return _azimuth; }
	double& get_sunDeclination() { return _sunDeclination; }
	double get_raysPerSquareMeters() const { return _raysPerSquareMeters; }
	Clock& get_time() { return _time; }
	Sunray* get_sunray(unsigned int i) { return _listOfSunrays[i]; }
	vector<Sunray*>& get_listOfSunrays() { return _listOfSunrays; }
	void removeSunray(unsigned int);
	void clearInterceptedSunrays();
	void clearListOfSunrays();

	void fOutputListOfSunrays(string ) const;

};

#endif
