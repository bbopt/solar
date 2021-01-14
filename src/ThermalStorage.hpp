#ifndef _THERMAL_STORAGE_H_
#define _THERMAL_STORAGE_H_

#include "MoltenSalt.hpp"
#include "Constants.hpp"
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <sstream>
using namespace std;

class ThermalStorage
{
	
private:
	double _storedMass;
	double _storedTemperature;
	MoltenSalt* _inputHTF;
	MoltenSalt* _outputHTF;

	//Dimensions parameters in meters
	double _heightOfStorage;
	double _diameterOfStorage;
	double _thicknessOfInsulation;
	
	//height of volume stored
	double _heightOfVolumeStored;

public:
	ThermalStorage(MoltenSalt*, MoltenSalt*, double, double, double);

	double fInitialStorageTemperature(int);
	double fComputeStorageTemperature(int );
	double fComputeStorageLevel();
	double fComputeStorageLevel(double);
	double fComputeStorageMass(int);
	double fInitialStorageMass(int);
	double fComputeEnergyLosses(double, double);
	double fComputeInsideRadiationLosses();
	double fSolveForT(double, double, double, double, double, double);
	double fSolveForT_i(double, double, double, double, double, double);

	void set_storage(double, double);
	void set_storage2(double, double);
	double get_heightOfVolumeStored(){ return _heightOfVolumeStored; }

	double get_storedMass() { return _storedMass; }
	double& get_storedTemperature() { return _storedTemperature; }
	double get_heightOfStorage(){ return _heightOfStorage; }
	double get_diameterOfStorage(){ return _diameterOfStorage; }
	double get_thicknessOfInsulation(){ return _thicknessOfInsulation; }
	MoltenSalt* get_inputHTF() { return _inputHTF; }
	MoltenSalt* get_outputHTF() { return _outputHTF; }
};

#endif
