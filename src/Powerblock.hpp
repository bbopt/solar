#ifndef _POWERBLOCK_H
#define _POWERBLOCK_H

#include "Global.hpp"
#include "Turbine.hpp"
#include <vector>

class Powerblock
{
private:
	Turbine _turbine;

	std::vector<double> _requiredThermalPower;
	std::vector<double> _powerOutput;
	double _Pout;
	double _steamRate;

public:
	Powerblock(int);

	int get_typeOfTurbine(){ return _turbine._ID; }
	double get_hotEnthalpy(){ return _turbine._inletEnthalpy; }
	double& get_coldEnthalpy(){ return _turbine._outletEnthalpy; }
	double get_temperature(){ return _turbine._inletTemperature; }
	double get_pressure(){ return _turbine._inletPressure; }
	double get_powerOfTurbine(){ return _turbine._maxPower; }
	std::vector<double>& get_powerOutput(){ return _powerOutput; }
	double& get_powerOutput(int i){ return _powerOutput[i]; }
	double& get_Pout(){ return _Pout; }
	double& get_turbineOutletEnthalpy(){ return _turbine._outletEnthalpy; }
	double& get_steamRate(){ return _steamRate; }

	std::vector<double>& get_requiredThermalPower(){ return _requiredThermalPower; }
	void adjustPowerData(double&, double&);
	void set_steamRate(double& x){ _steamRate = x; }

	double fComputeRequiredThermalEnergy(double& );
	double fComputeTurbineEnergy(double&);
	
};

#endif
