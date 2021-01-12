#ifndef _HTF_CYCLE_H_
#define _HTF_CYCLE_H_

#include "CentralReceiver.hpp"
#include "ThermalStorage.hpp"
#include "HeatExchanger.hpp"
#include "MoltenSalt.hpp"
#include "Constants.hpp"
#include "Powerblock.hpp"
#include "Global.hpp"

#include <cmath>

class HtfCycle
{
private:
	CentralReceiver _centralReceiver;
	ThermalStorage  _hotStorage;
	ThermalStorage  _coldStorage;
	HeatExchanger   _steamGenerator;

	MoltenSalt	_centralReceiverInlet;
	MoltenSalt      _centralReceiverOutlet;
	MoltenSalt      _steamGeneratorInlet;
	MoltenSalt      _steamGeneratorOutlet;

	//Design parameters
	//Hot Storage
	double _centralReceiverOutletDesignTemperature;
	double _storageTankHeight;
	double _storageTankDiameter;
	double _storageInsulationThickness;

	//cold storage is assumed to have the same dimensions as the hot storage

	//Heat Exchanger
	double _steamGeneratorOutletTemperature;
	//geometricParameters

	//Central Receiver
	double _receiverApertureHeight;
	double _receiverApertureWidth;
	double _receiverTubesInsideDiameter;
	double _receiverTubesThickness;
	int _receiverNumberOfTubes;
	int _receiverNumberOfPasses;
	//exchangerModelparameters

	//Simulation parameters
	double _iterationPrecision;
	int    _timeInterval; //minutes

	//Data gathering
	std::vector<double> _steamGenOutletMsRate;
	std::vector<double> _steamGenOutletTemp;
	double _minColdStorageTemp;
	double _minHotStorageTemp;
	double _minSteamGenOutletTemp;
	std::vector<double> _storageHeat;


public:

	HtfCycle(double, double, double, double, double, 
		Powerblock*, double, double, double, double, int, double, int);

	HtfCycle(double, double, double, double, double,
		Powerblock*, double, double, double, double, int, double, int,
		double, double, double, double, double, int, int, int, int);
	~HtfCycle();

	//Standard
	CentralReceiver& get_centralReceiver() { return _centralReceiver; }
	ThermalStorage&  get_hotStorage(){ return _hotStorage; }
	ThermalStorage&  get_coldStorage() { return _coldStorage; }
	HeatExchanger&   get_steamGenerator() { return _steamGenerator; }
	MoltenSalt& get_centralReceiverInlet(){ return _centralReceiverInlet; }
	MoltenSalt& get_centralReceiverOutlet() { return _centralReceiverOutlet; }
	MoltenSalt& get_steamGeneratorInlet() { return _steamGeneratorInlet; }
	MoltenSalt& get_steamGeneratorOutlet() { return _steamGeneratorOutlet; }

	//Data gathering
	std::vector<double>& get_steamGenOutletMsRate(){ return _steamGenOutletMsRate; }
	std::vector<double>& get_steamGenOutletTemp(){ return _steamGenOutletTemp; }
	std::vector<double>& get_steamGenHeatTransfered(){ return _steamGenerator.get_heatTransfered(); }
	double& get_storageHeat(int i){ return _storageHeat[i]; }
	std::vector<double>& get_storageHeatV(){ return _storageHeat; }
	double& get_minColdStorageTemp(){ return _minColdStorageTemp; }
	double& get_minHotStorageTemp(){ return _minHotStorageTemp; }
	double& get_minSteamGenTemp(){ return _minSteamGenOutletTemp; }
	void initiateColdStorage();
	void setStorage(double, double, double);

	void fOperateCycle(int, double, double);

	void fDetermineStatusOfHotStorage();
	void fDetermineStatusOfColdStorage();
};

#endif
