#ifndef _POWERPLANT_H_
#define _POWERPLANT_H_

#include "HtfCycle.hpp"
#include "MoltenSalt.hpp"
#include "HeliostatField.hpp"
#include "Economics.hpp"
#include "Powerblock.hpp"
#include "Constants.hpp"
#include "Global.hpp"
#include "Clock.hpp"

#include <cmath>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <fstream>

class Powerplant
{
private:

  //Simulation components
  Clock       _time;
  Sun         _sun;
  int         _model_type; // 1:heliostats field; 2: whole plant
  int         _heliostatsFieldModel;
  // std::string _pathToFieldData;  TOTO VIRER

  //Powerplant components
  HtfCycle      * _moltenSaltLoop;
  HeliostatField* _heliostatsField;
  Powerblock    * _powerblock;
  Economics     * _investmentCost;

  //Include day, starting hour, demand info, storage starting conditions, etc.
	
	//Data gathering
	std::vector<double> _sunElevation;
	std::vector<double> _sunAzimuth;
	std::vector<double> _sunEnergyGathered;
	std::vector<double> _heliostatFieldEfficiency;
	std::vector<double> _heliostatFieldPowerOutput;
	std::vector<double> _hotStorageLevel;
	std::vector<double> _hotStorageTemp;
	std::vector<double> _coldStorageLevel;
	std::vector<double> _coldStorageTemp;
	//
	std::vector<double> _receiverOutletFlow;
	std::vector<double> _receiverPumpHead;
	std::vector<double> _pressureShellSide;
	std::vector<double> _pressureTubesSide;
	std::vector<double> _steamRate;
	std::vector<double> _msRateSteamGen;
	std::vector<double> _heliostatsFieldPar;
	//
	std::vector<double> _energyToPowerBlockWatts;
	std::vector<double> _steamGeneratorInletTemperature;
	std::vector<double> _hotStoragePowerLosses;
	std::vector<double> _coldStoragePowerLosses;
	std::vector<double> _powerplantPowerOutput;
	std::vector<double> _demand;
	std::vector<double> _demandCompliance;
	double _reflectiveSurface;
	double _fieldSurface;
	double _costOfHeliostatsField;
	double _totalEnergyConcentrated;
	double _maximumPressureInReceiver;
	double _maximumPressureInExchanger;
	double _yieldPressureReceiver;
	double _yieldPressureExchanger;
	double _overallComplianceToDemand;

public:

	Powerplant(Clock&, Sun&,int,  HeliostatField*, HtfCycle*,Powerblock*, Economics* );

  // TOTO TUTU
  bool set_heliostatFieldPowerOutput_MINCOST_TS   ( void );
  bool set_heliostatFieldPowerOutput_MAXCOMP_HTF1 ( void );

	HtfCycle* get_moltenSaltLoop(){ return _moltenSaltLoop; }
	HeliostatField* get_heliostatField() { return _heliostatsField; }
	Economics* get_investmentCost(){ return _investmentCost; }

	void fSimulatePowerplant();
	void fSimulateHeliostatField();
	double fComputeSteamRate(double& );
	double fComputeThermalEnergy(double&);
	double fComputePressureInExchanger();
	double fComputeParasiticLosses();
	double fComputeParasiticsForPb3();
	double fComputeParasiticsForPb7();
	double fComputeParasiticsForPb9();

	//Get functions
	std::vector<double>& get_sunElevation(){ return _sunElevation; }
	std::vector<double>& get_sunAzimuth(){ return _sunAzimuth; }
	std::vector<double>& get_sunEnergyGathered(){ return _sunEnergyGathered; }
	std::vector<double>& get_heliostatFieldEfficiency(){ return _heliostatFieldEfficiency; }
	std::vector<double>& get_heliostatFieldPowerOutput(){ return _heliostatFieldPowerOutput; }
	std::vector<double>& get_hotStorageLevel(){ return _hotStorageLevel; }
	std::vector<double>& get_hotStorageTemp(){ return _hotStorageTemp; }
	std::vector<double>& get_coldStorageLevel(){ return _coldStorageLevel; }
	std::vector<double>& get_coldStorageTemp(){ return _coldStorageTemp; }
	std::vector<double>& get_receiverSurfaceTemperature(){ return _moltenSaltLoop->get_centralReceiver().get_surfaceTemperature(); }
	std::vector<double>& get_receiverLosses(){ return _moltenSaltLoop->get_centralReceiver().get_losses(); }
	std::vector<double>& get_receiverEfficiency(){ return _moltenSaltLoop->get_centralReceiver().get_efficiency(); }
	std::vector<double>& get_receiverMsRate(){ return _moltenSaltLoop->get_centralReceiver().get_msRate(); }
	std::vector<double>& get_steamGeneratorInletFlow(){ return _moltenSaltLoop->get_steamGenOutletMsRate(); }
	std::vector<double>& get_energyToPowerBlockWatts(){/* return _moltenSaltLoop->get_steamGenHeatTransfered();*/ 
		return _energyToPowerBlockWatts;
	}
	std::vector<double>& get_steamGeneratorInletTemperature(){ return _steamGeneratorInletTemperature; }
	std::vector<double>& get_steamGeneratorOutletTemperature(){ return _moltenSaltLoop->get_steamGenOutletTemp(); }
	std::vector<double>& get_hotStoragePowerLosses(){ return _hotStoragePowerLosses; }
	std::vector<double>& get_coldStoragePowerLosses(){ return _coldStoragePowerLosses; }
	std::vector<double>& get_powerplantPowerOutput(){ return _powerplantPowerOutput; }
	std::vector<double>& get_thermalPowerNeededFromBlock(){ return _powerblock->get_requiredThermalPower(); }

	double get_reflectiveSurface(){ return _reflectiveSurface; }
	double get_fieldSurface(){ return _fieldSurface; }
	double get_costOfHeliostatsField(){ return _costOfHeliostatsField; }
	double get_totalEnergyConcentrated(){ return _totalEnergyConcentrated; }

	double get_costOfHeliostatField(){ return _investmentCost->evaluateCostOfField(); }
	double get_costOfReceiver(){ return _investmentCost->evaluateCostOfReceiver(); }
	double get_costOfTower(){ return _investmentCost->evaluateCostOfTower(); }
	double get_costOfSteamGenerator(){ return _investmentCost->evaluateCostOfSteamGenerator(); }
	double get_costOfPowerblock(){ return _investmentCost->evaluateCostOfPowerblock(); }
	double get_costOfStorage(){ return _investmentCost->evaluateCostOfStorage(); }
	double get_steamTurbineInletTemperature(){ return _powerblock->get_temperature(); }

	double& get_overallComplianceToDemand(){ return _overallComplianceToDemand; }
	double& get_maximumPressureInReceiver(){ return _maximumPressureInReceiver; }
	double& get_yieldPressureInReceiver(){ return _yieldPressureReceiver; }
	double& get_yieldPressureInExchanger() { return _yieldPressureExchanger; }
	double& get_maximumPressureInExchanger() { return _maximumPressureInExchanger; }
	double& get_minColdStorageTemp(){ return _moltenSaltLoop->get_minColdStorageTemp(); }
	double& get_minHotStorageTemp(){ return _moltenSaltLoop->get_minHotStorageTemp(); }
	double& get_minSteamGenTemp(){ return _moltenSaltLoop->get_minSteamGenTemp(); }
	MoltenSalt& get_receiverInlet(){ return _moltenSaltLoop->get_centralReceiverInlet(); }
	MoltenSalt& get_receiverOutlet(){ return _moltenSaltLoop->get_centralReceiverOutlet(); }
	MoltenSalt& get_steamGenInlet(){ return _moltenSaltLoop->get_steamGeneratorInlet(); }
	MoltenSalt& get_steamGenOutlet(){ return _moltenSaltLoop->get_steamGeneratorOutlet(); }

  // Set methods:
  void set_demand(std::vector<double>& demandVector){ _demand = demandVector; }
  void set_heliostatModel(int& x){ _heliostatsFieldModel = x; }
  // void set_pathToFieldData(const std::string x){ _pathToFieldData = x; }  TOTO VIRER
};

#endif
