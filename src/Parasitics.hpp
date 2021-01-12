#ifndef _PARASITICS_H_
#define _PARASITICS_H_

#include <vector>

struct Parasitics
{
	std::vector<double> _recPressure; //Pa
	std::vector<double> _recFlow; //kg/s
	std::vector<double> _sGenMsPressure;
	std::vector<double> _sGenMsFlow; //kg/s
	std::vector<double> _turbineStRate; //kg/s
	std::vector<double> _sGenStPressure;

	std::vector<double> _vRecPumpPower;
	std::vector<double> _vSGenMsPumpPower;
	std::vector<double> _vSGenStPumpPower;
	std::vector<double> _vCoolerPumpPower;
	std::vector<double> _vCompressorPower;
	std::vector<double> _vHelConsPower;

	//kWh
	double _recPumpNrg;
	double _sGenMsNrg;
	double _sGenStNrg;
	double _coolerNrg;
	double _compressNrg;
	double _helConsNrg;

	Parasitics();
};

#endif