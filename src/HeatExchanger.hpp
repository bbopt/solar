#ifndef _HEAT_EXCHANGER_H_
#define _HEAT_EXCHANGER_H_

#include "MoltenSalt.hpp"
#include "Global.hpp"
#include "Constants.hpp"
#include "Powerblock.hpp"
#include "Turbine.hpp"
#include <stdexcept>
#include <vector>
#include <algorithm>


class HeatExchanger
{
private:
	MoltenSalt* _input;
	MoltenSalt* _output;
	Powerblock* _powerblock;
	double _inletWaterTemperature;
	double _inletWaterPressure;
	double _outletSteamTemperature;
	double _outletSteamPressure;
	int _exchangerModel;
	
	
	double _tubesLength;
	double _tubesDin;
	double _tubesDout;
	double _tubesSpacing;
	double _baffleCut;
	int _nbOfBaffles;
	int _nbOfTubes;
	int _nbOfPassesPerShell;
	int _nbOfShells;

	//corollary attributes
	double _baffleSpacing;
	double _shellWidth;
	double _shellCrossSection;
	double _crossFlowRows;
	double _L_E;
	double _windowAngle;
	double _windowArea_FG;
	double _windowArea_FR;
	double _windowArea_F;
	double _window_Nrows;
	double _totalRows;
	double _window_Ntubes;
	double _window_EqDiameter;
	double _window_EqPerimeter;
	double _longitudinalPitch;
	double _bypassArea;
	double _bundleArea;
	double _bundle_EqDiameter;
	double _nozzlesDiameter;
	double _nozzlesArea;
	double _a, _b, _c, _e;

	//geometric parameters to add
	double fComputeEpsilon();
	double fComputeC1(double, double);
	double fComputeC2();
	double fComputeM(double, double);
	double fComputeWaterEnthalpy(double);

	//Data gathering
	std::vector<double> _heatTransfered;

public:
	HeatExchanger(MoltenSalt*, MoltenSalt*, Powerblock*);
	HeatExchanger(MoltenSalt*, MoltenSalt*, Powerblock*,
		double, double, double, double, double, int, int, int, int);
	~HeatExchanger();

	void fCalculateEnergyTransfered();
	double fComputeRequiredMoltenSaltMassFlow(double, double);
	double fComputeRequiredMoltenSaltMassFlow(double, double, double);
	double fEnergyToPowerBlock(int );
	double computeYieldPressure();
	double computePressureInTubes(double);
	double computePressureInShells();

	std::vector<double>& get_heatTransfered(){ return _heatTransfered; }

	
	int get_exchangerModel() const { return _exchangerModel; };
	double get_tubesSpacing(){ return _tubesSpacing; }
	double get_tubesDin(){ return _tubesDin; }
	double get_tubesDout(){ return _tubesDout; }
	double get_tubesLength() { return _tubesLength; }
	int get_nbOfTubes(){ return _nbOfTubes; }
	int get_nbOfPassesPerShell(){ return _nbOfPassesPerShell; }
	int get_nbOfShells(){ return _nbOfShells; }

};

#endif
