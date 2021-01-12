#include "Turbine.hpp"

Turbine::Turbine(int ID)
:
_ID(ID)
{
	if (_ID == 1) //SST-110
	{
		_nbOfStages = SST110_STAGES;
		_minPower = 1.*SST110_MIN;
		_maxPower = 1.*SST110_MAX;
		_inletPressure = 1.*SST110_PRESSURE;
		_inletTemperature = 1.*SST110_TEMPERATURE;
		_inletEnthalpy = 1.*SST110_ENTHALPY;
		_condensing = SST110_CONDENSE;
		_reheat = SST110_REHEAT;
		_outletEnthalpy = 1.*SST110_OUTLET_ENTHALPY;
	}
	else if (_ID == 2) //SST120
	{
		_nbOfStages = SST120_STAGES;
		_minPower = 1.*SST120_MIN;
		_maxPower = 1.*SST120_MAX;
		_inletPressure = 1.*SST120_PRESSURE;
		_inletTemperature = 1.*SST120_TEMPERATURE;
		_inletEnthalpy = 1.*SST120_ENTHALPY;
		_condensing = SST120_CONDENSE;
		_reheat = SST120_REHEAT;
		_outletEnthalpy = 1.*SST120_OUTLET_ENTHALPY;
	}
	else if (_ID == 3) //SST300
	{
		_nbOfStages = SST300_STAGES;
		_minPower = 1.*SST300_MIN;
		_maxPower = 1.*SST300_MAX;
		_inletPressure = 1.*SST300_PRESSURE;
		_inletTemperature = 1.*SST300_TEMPERATURE;
		_inletEnthalpy = 1.*SST300_ENTHALPY;
		_condensing = SST300_CONDENSE;
		_reheat = SST300_REHEAT;
		_outletEnthalpy = 1.*SST300_OUTLET_ENTHALPY;
	}
	else if (_ID == 4) //SST400
	{
		_nbOfStages = SST400_STAGES;
		_minPower = 1.*SST400_MIN;
		_maxPower = 1.*SST400_MAX;
		_inletPressure = 1.*SST400_PRESSURE;
		_inletTemperature = 1.*SST400_TEMPERATURE;
		_inletEnthalpy = 1.*SST400_ENTHALPY;
		_condensing = SST400_CONDENSE;
		_reheat = SST400_REHEAT;
		_outletEnthalpy = 1.*SST400_OUTLET_ENTHALPY;
	}
	else if (_ID == 5) //SST600
	{
		_nbOfStages = SST600_STAGES;
		_minPower = 1.*SST600_MIN;
		_maxPower = 1.*SST600_MAX;
		_inletPressure = 1.*SST600_PRESSURE;
		_inletTemperature = 1.*SST600_TEMPERATURE;
		_inletEnthalpy = 1.*SST600_ENTHALPY;
		_condensing = SST600_CONDENSE;
		_reheat = SST600_REHEAT;
		_outletEnthalpy = 1.*SST600_OUTLET_ENTHALPY;
	}
	else if (_ID == 6) //SST700
	{
		_nbOfStages = SST700_STAGES;
		_minPower = 1.*SST700_MIN;
		_maxPower = 1.*SST700_MAX;
		_inletPressure = 1.*SST700_PRESSURE;
		_inletTemperature = 1.*SST700_TEMPERATURE;
		_inletEnthalpy = 1.*SST700_ENTHALPY;
		_condensing = SST700_CONDENSE;
		_reheat = SST700_REHEAT;
		_outletEnthalpy = 1.*SST700_OUTLET_ENTHALPY;
	}
	else if (_ID == 7) //SST800
	{
		_nbOfStages = SST800_STAGES;
		_minPower = 1.*SST800_MIN;
		_maxPower = 1.*SST800_MAX;
		_inletPressure = 1.*SST800_PRESSURE;
		_inletTemperature = 1.*SST800_TEMPERATURE;
		_inletEnthalpy = 1.*SST800_ENTHALPY;
		_condensing = SST800_CONDENSE;
		_reheat = SST800_REHEAT;
		_outletEnthalpy = 1.*SST800_OUTLET_ENTHALPY;
	}
	else if (_ID == 8) //SST900
	{
		_nbOfStages = SST900_STAGES;
		_minPower = 1.*SST900_MIN;
		_maxPower = 1.*SST900_MAX;
		_inletPressure = 1.*SST900_PRESSURE;
		_inletTemperature = 1.*SST900_TEMPERATURE;
		_inletEnthalpy = 1.*SST900_ENTHALPY;
		_condensing = SST900_CONDENSE;
		_reheat = SST900_REHEAT;
		_outletEnthalpy = 1.*SST900_OUTLET_ENTHALPY;
	}

	//Tuning basic efficiency coefficients
	double P = _inletPressure/1000.;
	double R_T = _maxPower/1000.;

	if (_condensing)
	{
		a_e = A1_e + B1_e*P + C1_e*pow(P, 2) + D1_e*pow(P, 3);
		b_e = A2_e + B2_e*P + C2_e*pow(P, 2) + D2_e*pow(P, 3);
		c_e = A3_e + B3_e*P + C3_e*pow(P, 2) + D3_e*pow(P, 3);
		d_e = A4_e + B4_e*P + C4_e*pow(P, 2) + D4_e*pow(P, 3);
	}
	else
	{
		a_e = A1_en + B1_en*P + C1_en*pow(P, 2) + D1_en*pow(P, 3);
		b_e = A2_en + B2_en*P + C2_en*pow(P, 2) + D2_en*pow(P, 3);
		c_e = A3_en + B3_en*P + C3_en*pow(P, 2) + D3_en*pow(P, 3);
		d_e = A4_en + B4_en*P + C4_en*pow(P, 2) + D4_en*pow(P, 3);
	}
	_basicEfficiency = exp(a_e + b_e / R_T + c_e / (R_T*R_T) + d_e / (R_T*R_T*R_T))/100.;

}
