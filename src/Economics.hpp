#ifndef _ECONOMICS_H
#define _ECONOMICS_H

#include "Global.hpp"

class Economics
{
private:

	//Attributes
	double _hotStorageInsulationThickness;
	double _coldStorageInsulationThickness;
	double _hotStorageHeight;
	double _storageDiameter;
	double _receiverInsulationThickness;
	double _heightOfTower;
	double _heightOfReceiverAperture;
	double _widthOfReceiverAperture;
	int    _receiverNumberOfTubes;
	double _receiverTubesDout;
	double _lengthOfHeliostats;
	double _widthOfHeliostats;
	int    _nbOfHeliostats;
	double _reflectiveArea;
	double _totalMoltenSaltMass;
	double _turbineNominalPowerOutput;
	int    _exchangerModel;
	double _exchangerTubesOutterDiameter;
	double _exchangerTubesLength;
	int    _exchangerNumberOfTubes;
	int    _exchangerTubePassesPerShell;
	int    _exchangerNumberOfShell;

	//Costs
	double _costOfField;
	double _costPerHeliostat;
	double _costOfTower;
	double _costOfStorage;
	double _costOfPowerblock;
	double _costOfSteamGenerator;
	double _costOfReceiver;
	double _totalCost;

public:
	//Heliostats Field only
	Economics();
	//Full plant
	Economics(int, double, double, double, double, double, double, double, double, double, double, double,
		double, double, int, int, int);


	double evaluateCostOfField(); /**/
	double evaluateCostOfHeliostat(); /**/
	double evaluateCostOfTower(); /**/
	double evaluateCostOfStorage(); /**/
	double evaluateCostOfReceiver(); /**/
	double evaluateCostOfPowerblock(); /**/
	double evaluateCostOfSteamGenerator(); /**/
	double evaluateTotalInvestmentCost();

	//Set
	void set_hotStorageInsulationThickness(double x){ _hotStorageInsulationThickness = x; }
	void set_coldStorageInsulationThickness(double x){ _coldStorageInsulationThickness = x; }
	void set_hotStorageHeight(double x){ _hotStorageHeight = x; }
	void set_storageDiameter(double x){ _storageDiameter = x; }
	void set_receiverInsulationThickness(double x){ _receiverInsulationThickness = x; }
	void set_heightOfTower(double x){ _heightOfTower = x; }
	void set_heightOfReceiverAperture(double x){ _heightOfReceiverAperture = x; }
	void set_widthOfReceiverAperture(double x){ _widthOfReceiverAperture = x; }
	void set_receiverNumberOfTubes(int x){ _receiverNumberOfTubes = x; }
	void set_receiverTubesDout(double x) { _receiverTubesDout = x; }
	void set_lengthOfHeliostats(double x){ _lengthOfHeliostats = x; }
	void set_widthOfHeliostats(double x){ _widthOfHeliostats = x; }
	void set_nbOfHeliostats(int x){ _nbOfHeliostats = x; }
	void set_reflectiveArea(double x){ _reflectiveArea = x; }
	void set_totalMoltenSaltMass(double x){ _totalMoltenSaltMass = x; }
	void set_turbineNominalPowerOutput(double x){ _turbineNominalPowerOutput = x; }
	void set_exchangerModel(int x) { _exchangerModel = x; }
	void set_exchangerTubesOutterDiameter(double x){ _exchangerTubesOutterDiameter = x; }
	void set_exchangerTubesLength(double x){ _exchangerTubesLength = x; }
	void set_exchangerNumberOfTubes(int x){ _exchangerNumberOfTubes = x; }
	void set_exchangerTubePassesPerShell(int x){ _exchangerTubePassesPerShell = x; }
	void set_exchangerNumberOfShell(int x){ _exchangerNumberOfShell = x; }

	//get
	double get_costOfField(){ return _costOfField; }
	double get_costOfTower(){ return _costOfTower; }
	double get_costOfStorage(){ return _costOfStorage; }
	double get_costOfReceiver(){ return _costOfReceiver; }
	double get_costOfPowerblock(){ return _costOfPowerblock; }
	double get_costOfSteamGenerator(){ return _costOfSteamGenerator; }

};

#endif
