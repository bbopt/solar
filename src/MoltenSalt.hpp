#ifndef _MOLTEN_SALT_H_
#define _MOLTEN_SALT_H_
#include "Constants.hpp"
#include <cmath>

//Simple class object containing basic thermodynamical properties of molten salt

class MoltenSalt
{
private:
	double _temperature; //in K
	double _enthalpy;    //in J/kg*K ---- Computed merely as internal energy because we assume liquid condition
	double _pressure;    //in kPa
	double _massFlow;    // in kg/s
	double _viscosity;	 //in kg/ms

public:
	MoltenSalt(double, double, double);
	MoltenSalt(double, double);
	MoltenSalt(MoltenSalt&);

	//standard
	double& get_temperature() { return _temperature; }
	double get_enthalpy() { return _enthalpy; }
	double get_pressure() { return _pressure; }
	double get_massFlow() { return _massFlow; }
	double get_viscosity() { return _viscosity; }

	void set_temperature(double);
	void set_enthalpy(double);
	void set_pressure(double pres) { _pressure = pres; }
	void set_massFlow(double masf) { _massFlow = masf; }
	void fComputeViscosity();

	static double fComputeViscosity(double);

	//Methods
	void fModifyEnergy(double);
};

#endif
