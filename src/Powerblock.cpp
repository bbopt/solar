#include "Powerblock.hpp"

Powerblock::Powerblock(int typeOfTurbine)
:
_turbine(typeOfTurbine)
{
	int n = _turbine._nbOfStages;

	//Tuning partial load correction coefficients
	a = A1 + B1*n + C1*pow(n, 2) + D1*pow(n, 3);
	b = A2 + B2*n + C2*pow(n, 2) + D2*pow(n, 3);
	c = A3 + B3*n + C3*pow(n, 2) + D3*pow(n, 3);
	d = A4 + B4*n + C4*pow(n, 2) + D4*pow(n, 3);

	_powerOutput.reserve(24 * 60 * 60);
	_requiredThermalPower.reserve(24 * 60 * 60);
}

double Powerblock::fComputeRequiredThermalEnergy(double& Pout)
{
	double f;
	_Pout = 0.;
	double correctionFactor, efficiency, requiredThermalEnergy;

	if (Pout > 0)
	{

		if (Pout >= _turbine._maxPower)
		{
			_Pout = _turbine._maxPower;
			f = 100;
		}
		else if (Pout < _turbine._minPower)
		{
			_Pout = _turbine._minPower;
			f = _Pout / _turbine._maxPower;
		}
		else {
			_Pout = Pout;
			f = 100 * (_Pout / _turbine._maxPower);
		}

		correctionFactor = exp(a + f*b + f*f*c + f*f*f*d);

		efficiency = _turbine._basicEfficiency * correctionFactor;

		requiredThermalEnergy = _Pout / efficiency;
	}
	else{
		_Pout = 0.;
		requiredThermalEnergy = 0.;
	}

	return requiredThermalEnergy;
}

void Powerblock::adjustPowerData(double& thermalTransfered, double& thermalNeeded)
{
	if (fabs(thermalTransfered - thermalNeeded) > 1000. )
	{
		_powerOutput.pop_back();
		_powerOutput.push_back(0.);
	}
	else
	{
		_powerOutput.push_back(thermalTransfered);
	}
}
