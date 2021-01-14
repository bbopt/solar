#ifndef _SCENARIO_H_
#define _SCENARIO_H_

#include "Problem.hpp"
#include "Global.hpp"
#include "Powerplant.hpp"

#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <numeric>

/*----------------*/
class Scenario {
/*----------------*/
  
private:

  int      _n; // number of variables
  double * _x; // vector of variables
  
  std::string _problem; // problem id

  // Scenario parameters:
  int _model_type; // 1:heliostats field; 2: whole plant
  int _heliostatsFieldModel;
  int _exchangerModel;
  
  // std::string _pathToFieldData;  // TOTO VIRER

  
	int _storageStartupCondition; //% of hot storage level //[6]

	//1- constant from t_0 to t_f with value X
	//2- Winter say profile
	//3- Summer day profile
	//4- customized profile using external file (not done yet)
	int _demandProfile; //[7]
	double _maximumPowerDemand;//[8]
	int _tStart;//[9]
	int _tEnd;//[10]
	std::vector<double> _demand;

	double _latitude; //[11]
	int _day;//[12]
	//-------------------------------------------------------------------
	
	//Simulation parameters
	//-------------------------------------------------------------------
	int _numberOfTimeIncrements; //[13]
	int _minutesPerTimeIncrement;//[14]
	double _fixedPointsPrecision;//[15] //in percent
	double _raysPerSquareMeters;//[16]

	//Problem-dependant constraints
	//-------------------------------------------------------------------
	double _cBudget; //$ //[17]
	double _cEnergy; //kWh //[18]
	double _cDemandComplianceRatio; //% //[19]
	double _cMaximumDifferenceToDemand; //W //[20]
	double _cFieldEfficiency; // % //[21]
	double _cFieldSurface; //m^2 //[22]
	double _cReflectiveSurface; //m^2 //[23]
	double _cParasitics; // %

	//Variables
	//-------------------------------------------------------------------
	//Heliostats Field
	double _heliostatHeight;//m //[24]
	double _heliostatWidth; //m //[25]
	double _towerHeight; //m //[26]
	double _receiverApertureHeight; //m //[27]
	double _receiverApertureWidth;//m //[28]
	int _maxNumberOfHeliostats; //[29]
	double _fieldMaxAngularSpan; //degrees //[30]
	double _minimumDistanceToTower;//fraction of tower height //[31]
	double _maximumDistanceToTower;//fraction of tower height //[32]

	//Heat Transfer loop
	double _centralReceiverOutletTemperature; //K //[33]
	double _hotStorageHeight; //m //[34]
	double _hotStorageDiameter; //m //[35]
	double _hotStorageInsulThickness; //m //[36]
	double _coldStorageInsulThickness; //m //[37]
	double _coldMoltenSaltMinTemperature; //K //[38]
	int    _receiverNbOfTubes; //[39]
	double _receiverInsulThickness; //m //[40]
	double _receiverTubesInsideDiam;//m //[41]
	double _receiverTubesOutsideDiam;//m //[42]
	
	//steam generator
	double _exchangerTubesSpacing; //m   //[43]
	double _exchangerTubesLength; //m    //[44]
	double _exchangerTubesDin; //m       //[45]
	double _exchangerTubesDout;//m       //[46]
	double _exchangerBaffleCut;    //
	int    _exchangerNbOfBaffles;	   //
	int    _exchangerNbOfTubes;             //[47]
	int    _exchangerNbOfShells;            //[48]
	int    _exchangerNbOfPassesPerShell;    //[49]

	//powerblock
	int _typeOfTurbine;       //[50]

	//-----------------------------------------------------------------------
	double _minReceiverOutletTemp; //for validation. not an input parameter
	void fFillDemandVector();

	Clock      * _time;
	Sun        * _sun;
	Powerplant * _powerplant;

public:
  Scenario  ( const std::string & problem , const std::string & x_file_name );
  ~Scenario ( void ) { delete_x(); }

private:

  void init_maxNrg_H1     ( void ); // #1
  void init_minSurf_H1    ( void ); // #2
  void init_minCost_C1    ( void ); // #3
  void init_minCost_C2    ( void ); // #4
  void init_maxComp_HTF1  ( void ); // #5
  void init_minCost_TS    ( void ); // #6
  void init_maxEff_RE     ( void ); // #7
  void init_maxHF_minCost ( void ); // #8
  void init_maxNrg_minPar ( void ); // #9
  
  bool read               ( const std::string & x_file_name );
  bool read_x             ( const std::string & x_file_name );
  bool read21_ST1_DEM1    ( const std::string & x_file_name );

  bool read_maxNrg_H1     ( const std::string & x_file_name ); // #1
  bool read_minSurf_H1    ( const std::string & x_file_name ); // #2
  bool read_minCost_C1    ( const std::string & x_file_name ); // #3
  bool read_minCost_C2    ( const std::string & x_file_name ); // #4
  bool read_maxComp_HTF1  ( const std::string & x_file_name ); // #5
  bool read_minCost_TS    ( const std::string & x_file_name ); // #6
  bool read_maxEff_RE     ( const std::string & x_file_name ); // #7
  bool read_maxHF_minCost ( const std::string & x_file_name ); // #8
  bool read_maxNrg_minPar ( const std::string & x_file_name ); // #9

	//add function to return instantaneous power requirement depending on scenario.
private:
  bool validateInputValues    ( void ) const;
  bool validate_maxNrg_H1     ( void ) const; // #1
  bool validate_minSurf_H1    ( void ) const; // #2
  bool validate_minCost_C1    ( void ) const; // #3
  bool validate_minCost_C2    ( void ) const; // #4
  bool validate_maxComp_HTF1  ( void ) const; // #5
  bool validate_minCost_TS    ( void ) const; // #6
  bool validate_maxEff_RE     ( void ) const; // #7
  bool validate_maxHF_minCost ( void ) const; // #8
  bool validate_maxNrg_minPar ( void ) const; // #9

	//construct powerplant
private:
  void construct2H1ST1();
  void construct_maxNrg_H1();
  void construct_minSurf_H1();
  void construct_minCost_C1();
  void construct_minCost_C2();
  void construct_maxComp_HTF1();
  void construct_minCost_TS();
  void construct_maxEff_RE();
  void construct_maxHF_minCost();
  void construct_maxNrg_minPar();

public:
  bool simulate ( std::ostream & out , bool & cnt_eval , bool verbose );
private:
  bool simulate_maxNrg_H1     ( std::ostream & out , bool & cnt_eval ); // #1
  bool simulate_minSurf_H1    ( std::ostream & out , bool & cnt_eval ); // #2
  bool simulate_minCost_C1    ( std::ostream & out , bool & cnt_eval ); // #3
  bool simulate_minCost_C2    ( std::ostream & out , bool & cnt_eval ); // #4
  bool simulate_maxComp_HTF1  ( std::ostream & out , bool & cnt_eval ); // #5
  bool simulate_minCost_TS    ( std::ostream & out , bool & cnt_eval ); // #6
  bool simulate_maxEff_RE     ( std::ostream & out , bool & cnt_eval ); // #7
  bool simulate_maxHF_minCost ( std::ostream & out , bool & cnt_eval ); // #8
  bool simulate_maxNrg_minPar ( std::ostream & out , bool & cnt_eval ); // #9

	//Print
public:
  void print() const;   // TOTO: en verbose mode plutot
  void display_x ( std::ostream & out ) const;

private:
        // void printColdStorage();
	// void printHotStorage();
        // void printSteamGeneratorOutlet();
	// void printRequiredThermalPower();
	// void printReceiver();
	// void printHeliostatsFieldLayout();
	// void printHeliostatsFieldPerformance();
	// void printPowerOutput();
	// void printDemand();
	// void printFlows();

private:

  static bool is_int ( const double x ) { return (round(x)==x); }
  static int  round  ( const double x ) { return static_cast<int> ((x < 0.0 ? -std::floor(.5-x) : std::floor(.5+x))); }

  void delete_x ( void );

};

/*------------------*/
class foncteurSum {
/*------------------*/
private:
  double* _sum;
public:
  foncteurSum(double* sum) { _sum = sum; }
  void operator()(double d){ (*_sum) = (*_sum) + d; }
};

/*------------------------*/
struct foncteurPrintPower {
/*------------------------*/
  std::ofstream* _stream;
  foncteurPrintPower(std::ofstream* stream) { _stream = stream; }
  void operator()(double pwr) { *_stream << pwr << std::endl; }
};

#endif
