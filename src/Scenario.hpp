/*-------------------------------------------------------------------------------*/
/*  SOLAR - The solar thermal power plant simulator                              */
/*  https://github.com/bbopt/solar                                               */
/*                                                                               */
/*  Miguel Diago, Sebastien Le Digabel, Mathieu Lemyre-Garneau, Bastien Talgorn  */
/*                                                                               */
/*  Polytechnique Montreal / GERAD                                               */
/*  sebastien.le-digabel@polymtl.ca                                              */
/*                                                                               */
/*  This program is free software: you can redistribute it and/or modify it      */
/*  under the terms of the GNU Lesser General Public License as published by     */
/*  the Free Software Foundation, either version 3 of the License, or (at your   */
/*  option) any later version.                                                   */
/*                                                                               */
/*  This program is distributed in the hope that it will be useful, but WITHOUT  */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  */
/*  for more details.                                                            */
/*                                                                               */
/*  You should have received a copy of the GNU Lesser General Public License     */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.         */
/*                                                                               */
/*-------------------------------------------------------------------------------*/
#ifndef __SCENARIO_H__
#define __SCENARIO_H__

#include "Problem.hpp"
#include "Global.hpp"
#include "Powerplant.hpp"
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <numeric>

/*----------------*/
class Scenario {
/*----------------*/
  
private:
  
  std::string _problem; // problem id
  
  // Scenario parameters:
  int _model_type; // 1:heliostats field; 2: whole plant
  int _heliostatsFieldModel;
  int _exchangerModel;
  
  int _storageStartupCondition; // % of hot storage level

  // 1- constant from t_0 to t_f with value X
  // 2- Winter say profile
  // 3- Summer day profile
  // 4- customized profile using external file (not done yet)
  int    _demandProfile;
  double _maximumPowerDemand;
  int    _tStart;
  int    _tEnd;
  std::vector<double> _demand;

  double _latitude;
  int    _day;

  // Simulation parameters:
  int    _numberOfTimeIncrements;
  int    _minutesPerTimeIncrement;
  double _fixedPointsPrecision; //in percent
  double _raysPerSquareMeters;

  // Problem-dependant constraints:
  double _cBudget;                    // $
  double _cEnergy;                    // kWh
  double _cDemandComplianceRatio;     // %
  double _cMaximumDifferenceToDemand; // W
  double _cFieldEfficiency;           // %
  double _cFieldSurface;              // m^2
  double _cReflectiveSurface;         // m^2
  double _cParasitics;                // %

  // Heliostats Field:
  double _heliostatLength;        // m
  double _heliostatWidth;         // m
  double _towerHeight;            // m
  double _receiverApertureHeight; // m
  double _receiverApertureWidth;  // m
  int    _numberOfHeliostats; 
  double _fieldAngularWidth;      // degrees
  double _minimumDistanceToTower; // fraction of tower height
  double _maximumDistanceToTower; // fraction of tower height
  
  // Heat Transfer loop:
  double _centralReceiverOutletTemperature; // K
  double _hotStorageHeight;                 // m
  double _hotStorageDiameter;               // m
  double _hotStorageInsulThickness;         // m
  double _coldStorageInsulThickness;        // m
  double _coldMoltenSaltMinTemperature;     // K
  int    _receiverNbOfTubes;
  double _receiverInsulThickness;           // m
  double _receiverTubesInsideDiam;          // m
  double _receiverTubesOutsideDiam;         // m
  
  // Steam generator:
  double _exchangerTubesSpacing; // m
  double _exchangerTubesLength;  // m
  double _exchangerTubesDin;     // m
  double _exchangerTubesDout;    // m
  double _exchangerBaffleCut;
  int    _exchangerNbOfBaffles;
  int    _exchangerNbOfTubes;
  int    _exchangerNbOfShells;
  int    _exchangerNbOfPassesPerShell;
  
  // Powerblock:
  int _typeOfTurbine;
  
  // Powerplant:
  Powerplant * _powerplant;

  double _minReceiverOutletTemp; //for validation; not an input parameter
  
private:

  // for the variable fidelity surrogates:
  static int compute_numberOfTimeIncrements ( int min, int max, double fidelity ) {
    return min + static_cast<int>(floor(fidelity*(max-min)));
  }

  static double compute_raysPerSquareMeters ( double min, double max, double fidelity ) {
    return min + fidelity*(max-min);
  }
  
  void init_maxNrg_H1             ( void );            // # 1
  void init_minSurf_H1            ( double fidelity ); // # 2
  void init_minCost_C1            ( double fidelity ); // # 3
  void init_minCost_C2            ( double fidelity ); // # 4
  void init_maxComp_HTF1          ( void );            // # 5
  void init_minCost_TS            ( void );            // # 6
  void init_maxEff_RE             ( double fidelity ); // # 7
  void init_maxHF_minCost         ( double fidelity ); // # 8
  void init_maxNrg_minPar         ( double fidelity ); // # 9
  void init_minCost_unconstrained ( double fidelity ); // #10

  bool set_typeOfTurbine ( int tot );
  
  bool set_x_maxNrg_H1             ( const double * x ); // # 1
  bool set_x_minSurf_H1            ( const double * x ); // # 2
  bool set_x_minCost_C1            ( const double * x ); // # 3
  bool set_x_minCost_C2            ( const double * x ); // # 4
  bool set_x_maxComp_HTF1          ( const double * x ); // # 5
  bool set_x_minCost_TS            ( const double * x ); // # 6
  bool set_x_maxEff_RE             ( const double * x ); // # 7
  bool set_x_maxHF_minCost         ( const double * x ); // # 8
  bool set_x_maxNrg_minPar         ( const double * x ); // # 9
  bool set_x_minCost_unconstrained ( const double * x ); // #10

  bool check_bounds_maxNrg_H1             ( void ) const; // # 1
  bool check_bounds_minSurf_H1            ( void ) const; // # 2
  bool check_bounds_minCost_C1            ( void ) const; // # 3
  bool check_bounds_minCost_C2            ( void ) const; // # 4
  bool check_bounds_maxComp_HTF1          ( void ) const; // # 5
  bool check_bounds_minCost_TS            ( void ) const; // # 6
  bool check_bounds_maxEff_RE             ( void ) const; // # 7
  bool check_bounds_maxHF_minCost         ( void ) const; // # 8
  bool check_bounds_maxNrg_minPar         ( void ) const; // # 9
  bool check_bounds_minCost_unconstrained ( void ) const; // #10
  
  bool check_apriori_constraints_maxNrg_H1             ( void ) const; // # 1
  bool check_apriori_constraints_minSurf_H1            ( void ) const; // # 2
  bool check_apriori_constraints_minCost_C1            ( void ) const; // # 3
  bool check_apriori_constraints_minCost_C2            ( void ) const; // # 4
  bool check_apriori_constraints_maxComp_HTF1          ( void ) const; // # 5
  bool check_apriori_constraints_minCost_TS            ( void ) const; // # 6
  bool check_apriori_constraints_maxEff_RE             ( void ) const; // # 7
  bool check_apriori_constraints_maxHF_minCost         ( void ) const; // # 8
  bool check_apriori_constraints_maxNrg_minPar         ( void ) const; // # 9
  bool check_apriori_constraints_minCost_unconstrained ( void ) const; // #10
  
  void fFillDemandVector ( void );

  void construct_maxNrg_H1             ( bool & cnt_eval ); // # 1
  void construct_minSurf_H1            ( bool & cnt_eval ); // # 2
  void construct_minCost_C1            ( bool & cnt_eval ); // # 3
  void construct_minCost_C2            ( bool & cnt_eval ); // # 4
  void construct_maxComp_HTF1          ( bool & cnt_eval ); // # 5
  void construct_minCost_TS            ( bool & cnt_eval ); // # 6
  void construct_maxEff_RE             ( bool & cnt_eval ); // # 7
  void construct_maxHF_minCost         ( bool & cnt_eval ); // # 8
  void construct_maxNrg_minPar         ( bool & cnt_eval ); // # 9
  void construct_minCost_unconstrained ( bool & cnt_eval ); // #10

  bool simulate_maxNrg_H1             ( double * outputs, bool & cnt_eval ); // # 1
  bool simulate_minSurf_H1            ( double * outputs, bool & cnt_eval ); // # 2
  bool simulate_minCost_C1            ( double * outputs, bool & cnt_eval ); // # 3
  bool simulate_minCost_C2            ( double * outputs, bool & cnt_eval ); // # 4
  bool simulate_maxComp_HTF1          ( double * outputs, bool & cnt_eval ); // # 5
  bool simulate_minCost_TS            ( double * outputs, bool & cnt_eval ); // # 6
  bool simulate_maxEff_RE             ( double * outputs, bool & cnt_eval ); // # 7
  bool simulate_maxHF_minCost         ( double * outputs, bool & cnt_eval ); // # 8
  bool simulate_maxNrg_minPar         ( double * outputs, bool & cnt_eval ); // # 9
  bool simulate_minCost_unconstrained ( double * outputs              ,
					double * intermediate_outputs ,
					bool     low_fid              ,
					bool   & cnt_eval               ); // #10

  void display_x ( const double * , std::ostream & out ) const;
  
public:

  Scenario  ( const std::string & problem, double fidelity );
  ~Scenario ( void );
  
  bool set_x    ( const double * x );
  bool simulate ( double * outputs, double * intermediate_outputs, double fidelity, bool & cnt_eval );
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

#endif
