#ifndef _TURBINE_H_
#define _TURBINE_H_

#include "Constants.hpp"
#include "Global.hpp"

struct Turbine
{
	int _ID; //1 to 8
	int _nbOfStages;
	double _minPower; //W
	double _maxPower; //W
	double _inletPressure; //Pa
	double _inletTemperature; //K
	double _inletEnthalpy; //j/kg
	double _basicEfficiency;
	double _outletEnthalpy;
	bool _condensing;
	bool _reheat;

	Turbine(int);
};

#endif
