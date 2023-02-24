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
#include "Scenario.hpp"

/*-------------------------------------*/
/*              constructor            */
/*-------------------------------------*/
Scenario::Scenario ( const std::string & problem , double fidelity ) :
  _problem                          ( problem ) ,
  _model_type                       ( 0       ) ,
  _heliostatsFieldModel             ( 0       ) ,
  _exchangerModel                   ( 0       ) ,
  _storageStartupCondition          ( 0       ) ,
  _demandProfile                    ( 0       ) ,
  _maximumPowerDemand               ( 0.0     ) ,
  _tStart                           ( 0       ) ,
  _tEnd                             ( 0       ) ,
  _latitude                         ( 0.0     ) ,
  _day                              ( 0       ) ,
  _numberOfTimeIncrements           ( 0       ) ,
  _minutesPerTimeIncrement          ( 0       ) ,
  _fixedPointsPrecision             ( 0.0     ) ,
  _raysPerSquareMeters              ( 0.0     ) ,
  _cBudget                          ( 0.0     ) ,
  _cEnergy                          ( 0.0     ) ,
  _cDemandComplianceRatio           ( 0.0     ) ,
  _cMaximumDifferenceToDemand       ( 0.0     ) ,
  _cFieldEfficiency                 ( 0.0     ) ,
  _cFieldSurface                    ( 0.0     ) ,
  _cReflectiveSurface               ( 0.0     ) ,
  _cParasitics                      ( 0.0     ) ,
  _heliostatLength                  ( 0.0     ) ,
  _heliostatWidth                   ( 0.0     ) ,
  _towerHeight                      ( 0.0     ) ,
  _receiverApertureHeight           ( 0.0     ) , 
  _receiverApertureWidth            ( 0.0     ) , 
  _numberOfHeliostats               ( 0       ) , 
  _fieldAngularWidth                ( 0.0     ) ,
  _minimumDistanceToTower           ( 0.0     ) ,
  _maximumDistanceToTower           ( 0.0     ) ,
  _centralReceiverOutletTemperature ( 0.0     ) ,
  _hotStorageHeight                 ( 0.0     ) ,
  _hotStorageDiameter               ( 0.0     ) ,
  _hotStorageInsulThickness         ( 0.0     ) ,
  _coldStorageInsulThickness        ( 0.0     ) ,
  _coldMoltenSaltMinTemperature     ( 0.0     ) ,
  _receiverNbOfTubes                ( 0       ) ,
  _receiverInsulThickness           ( 0.0     ) ,
  _receiverTubesInsideDiam          ( 0.0     ) ,
  _receiverTubesOutsideDiam         ( 0.0     ) ,
  _exchangerTubesSpacing            ( 0.0     ) ,
  _exchangerTubesLength             ( 0.0     ) ,
  _exchangerTubesDin                ( 0.0     ) ,
  _exchangerTubesDout               ( 0.0     ) ,
  _exchangerBaffleCut               ( 0.0     ) ,
  _exchangerNbOfBaffles             ( 0       ) ,
  _exchangerNbOfTubes               ( 0       ) ,
  _exchangerNbOfShells              ( 0       ) ,
  _exchangerNbOfPassesPerShell      ( 0       ) ,
  _typeOfTurbine                    ( 0       ) ,
  _powerplant                       ( NULL    ) ,
  _minReceiverOutletTemp            ( 0.0     )   {

  // Check fidelity:
  // ---------------
  if ( fidelity <= 0.0 || fidelity > 1.0 )
    throw std::invalid_argument ( "fidelity is not valid" );

  if ( fidelity < 1.0 && _problem == "MAXNRG_H1" )
    throw std::invalid_argument ( "fidelity must be 1.0 for Problem MAXNRG_H1 (#1)" );

  if ( fidelity < 1.0 && _problem == "MAXCOMP_HTF1" )
    throw std::invalid_argument ( "fidelity must be 1.0 for Problem MAXCOMP_HTF1 (#5)" );

  if ( fidelity < 1.0 && _problem == "MINCOST_TS" )
    throw std::invalid_argument ( "fidelity must be 1.0 for Problem MINCOST_TS (#6)" );
  
  // Problem #1:
  if ( _problem == "MAXNRG_H1" )
    init_maxNrg_H1();

  // Problem #2:
  if ( _problem == "MINSURF_H1" )
    init_minSurf_H1 ( fidelity );
  
  // Problem #3:
  if ( _problem == "MINCOST_C1" )
    init_minCost_C1 ( fidelity );

  // Problem #4:
  if ( _problem == "MINCOST_C2" )
    init_minCost_C2 ( fidelity );
  
  // Problem #5:
  if ( _problem == "MAXCOMP_HTF1" )
    init_maxComp_HTF1();

  // Problem #6:
  if ( _problem == "MINCOST_TS" )
    init_minCost_TS();

  // Problem #7:
  if ( _problem == "MAXEFF_RE" )
    init_maxEff_RE ( fidelity );
  
  // Problem #8:
  if ( _problem == "MAXHF_MINCOST" )
    init_maxHF_minCost ( fidelity );

  // Problem #9:
  if ( _problem == "MAXNRG_MINPAR" )
    init_maxNrg_minPar ( fidelity );

  // Problem #10:
  if ( _problem == "MINCOST_UNCONSTRAINED" )
    init_minCost_unconstrained ( fidelity );
  
  // Check problem id:
  // -----------------
  if ( problem != "MAXNRG_H1"             &&    // # 1
       problem != "MINSURF_H1"            &&    // # 2
       problem != "MINCOST_C1"            &&    // # 3
       problem != "MINCOST_C2"            &&    // # 4
       problem != "MAXCOMP_HTF1"          &&    // # 5
       problem != "MINCOST_TS"            &&    // # 6
       problem != "MAXEFF_RE"             &&    // # 7
       problem != "MAXHF_MINCOST"         &&    // # 8
       problem != "MAXNRG_MINPAR"         &&    // # 9
       problem != "MINCOST_UNCONSTRAINED"    )  // # 10
    throw std::invalid_argument ( problem + " is not a valid problem ID" );
}

/*-----------------------------------------*/
/*                 destructor              */
/*-----------------------------------------*/
Scenario::~Scenario ( void ) {
  if ( _powerplant )
    delete _powerplant;
}

/*-----------------------------------------*/
/*           general set_x function        */
/*-----------------------------------------*/
bool Scenario::set_x ( const double * x ) {
  if ( _problem == "MAXNRG_H1"             ) { return set_x_maxNrg_H1             ( x ); } // # 1
  if ( _problem == "MINSURF_H1"            ) { return set_x_minSurf_H1            ( x ); } // # 2
  if ( _problem == "MINCOST_C1"            ) { return set_x_minCost_C1            ( x ); } // # 3
  if ( _problem == "MINCOST_C2"            ) { return set_x_minCost_C2            ( x ); } // # 4
  if ( _problem == "MAXCOMP_HTF1"          ) { return set_x_maxComp_HTF1          ( x ); } // # 5
  if ( _problem == "MINCOST_TS"            ) { return set_x_minCost_TS            ( x ); } // # 6
  if ( _problem == "MAXEFF_RE"             ) { return set_x_maxEff_RE             ( x ); } // # 7
  if ( _problem == "MAXHF_MINCOST"         ) { return set_x_maxHF_minCost         ( x ); } // # 8
  if ( _problem == "MAXNRG_MINPAR"         ) { return set_x_maxNrg_minPar         ( x ); } // # 9
  if ( _problem == "MINCOST_UNCONSTRAINED" ) { return set_x_minCost_unconstrained ( x ); } // #10
  return false;
}

/*-----------------------------------------*/
/*    initialize problem maxNrg_H1 (#1)    */
/*-----------------------------------------*/
void Scenario::init_maxNrg_H1 ( void ) {
  _model_type              = 1; // heliostats field
  _heliostatsFieldModel    = 1;
  _exchangerModel          = 1;
  _latitude                = 44.95;
  _day                     = 100; // April 10th
  _numberOfTimeIncrements  = 24;
  _minutesPerTimeIncrement = 60;
  _raysPerSquareMeters     = 0.01;
  _cFieldSurface           = 1.95e6; // 195 hectares
  _cBudget                 = 50e6;
}

/*-----------------------------------------*/
/*  set inputs for problem maxNRG_H1 (#1)  */
/*-----------------------------------------*/
bool Scenario::set_x_maxNrg_H1 ( const double * x ) {

  // check discrete variables:
  // -------------------------
  if ( !is_int(x[5]) )
    throw std::invalid_argument ( "Problem with input: One of the discrete variables has a non-integer value" );

  // assign variables:
  // -----------------
  _heliostatLength        = x[0];
  _heliostatWidth         = x[1];
  _towerHeight            = x[2];
  _receiverApertureHeight = x[3];
  _receiverApertureWidth  = x[4];
  _numberOfHeliostats     = myround(x[5]);
  _fieldAngularWidth      = x[6];
  _minimumDistanceToTower = x[7];
  _maximumDistanceToTower = x[8];

  // check bounds:
  // -------------
  if ( !check_bounds_maxNrg_H1() )
    throw std::invalid_argument ( "Problem with input: One of the inputs is outside its bounds" );
  
  return true;
}

/*-----------------------------------------*/
/*    initialize problem minSurf_H1 (#2)   */
/*-----------------------------------------*/
void Scenario::init_minSurf_H1 ( double fidelity ) {
  
  // scenario parameters:
  _model_type           = 2; // whole plant
  _heliostatsFieldModel = 1;
  _exchangerModel       = 1;
  
  _latitude                = 37.55;
  _day                     = 180; // June 29th
  _storageStartupCondition = 60;
  _demandProfile           = 3;
  _maximumPowerDemand      = 20e6;
  _minutesPerTimeIncrement = 60;

  // variable fidelity surrogate:
  _numberOfTimeIncrements = Scenario::compute_numberOfTimeIncrements (    1,   72, fidelity );
  _raysPerSquareMeters    = Scenario::compute_raysPerSquareMeters    ( 1e-7, 0.01, fidelity );
   
  _cFieldSurface           = 4e6;   // 400 hectares
  _cDemandComplianceRatio  = 100;
  _cBudget                 = 300e6; // 300 millions
  _cParasitics             = 0.18;  // not used

  // design parameters:
  _hotStorageDiameter           = 23.0;
  _hotStorageHeight             = 10.5;
  _hotStorageInsulThickness     = 0.3;
  _coldStorageInsulThickness    = 0.2;
  _coldMoltenSaltMinTemperature = 555.0; // 282 + 273

  // type of turbine:
  set_typeOfTurbine(3);
}

/*------------------------------------------*/
/*  set inputs for problem minSurf_H1 (#2)  */
/*------------------------------------------*/
bool Scenario::set_x_minSurf_H1 ( const double * x ) {

  // check discrete variables:
  // -------------------------
  if ( !is_int(x[5]) || !is_int(x[10]) )
    throw std::invalid_argument ( "Problem with input: One of the discrete variables has a non-integer value" );
 
  // assign variables:
  // -----------------

  // Heliostats Field:
  _heliostatLength        = x[0];
  _heliostatWidth         = x[1];
  _towerHeight            = x[2];
  _receiverApertureHeight = x[3];
  _receiverApertureWidth  = x[4];
  _numberOfHeliostats     = myround(x[5]);
  _fieldAngularWidth      = x[6];
  _minimumDistanceToTower = x[7];
  _maximumDistanceToTower = x[8];

  // Htf cycle:
  _centralReceiverOutletTemperature = x[ 9];
  _receiverNbOfTubes                = myround(x[10]);
  _receiverInsulThickness           = x[11];
  _receiverTubesInsideDiam          = x[12];
  _receiverTubesOutsideDiam         = x[13];

  // check bounds:
  // -------------
  if ( !check_bounds_minSurf_H1() )
    throw std::invalid_argument ( "Problem with input: One of the inputs is outside its bounds" );

  return true;
}
  
/*-----------------------------------------*/
/*    initialize problem minCost_C1 (#3)   */
/*-----------------------------------------*/
void Scenario::init_minCost_C1 ( double fidelity ) {

  // Demand profile: 1, from 3 pm to 9 pm
  // Maximum demand: 10 MW
  // Latitude      : 35 deg
  // Day           : 1
  // Duration      : 48 hours
  // maximum field surface of 800,000 m^2 (80 hectares)

  // Scenario parameters:
  _model_type           = 2; // whole plant
  _heliostatsFieldModel = 1;
  _exchangerModel       = 1;
  
  _latitude                = 35.0;
  _day                     = 270;  // https://www.epochconverter.com/days/2019: Sept. 27th
  _demandProfile           = 1;
  _tStart                  = 900;  // 15*60
  _tEnd                    = 1260; // 21*60;
  _maximumPowerDemand      = 1e7;
  _storageStartupCondition = 0;
  _minutesPerTimeIncrement = 60;
  _fixedPointsPrecision    = 0.001;
  _cFieldSurface           = 800000; // 80 hectares
  _cDemandComplianceRatio  = 100;
  
  // variable fidelity surrogate:
  _numberOfTimeIncrements = Scenario::compute_numberOfTimeIncrements (    1,   24, fidelity );
  _raysPerSquareMeters    = Scenario::compute_raysPerSquareMeters    ( 1e-6, 0.01, fidelity ); 
}
  
/*-----------------------------------------*/
/*  set inputs for problem minCost_C1 (#3) */
/*-----------------------------------------*/
bool Scenario::set_x_minCost_C1 ( const double * x ) {

  // check discrete variables:
  // -------------------------
  if ( !is_int(x[5]) || !is_int(x[15]) || !is_int(x[19]) )
    throw std::invalid_argument ( "Problem with input: One of the discrete variables has a non-integer value" );
 
  // assign variables:
  // -----------------

  // Heliostats Field:
  _heliostatLength        = x[ 0];
  _heliostatWidth         = x[ 1];
  _towerHeight            = x[ 2];
  _receiverApertureHeight = x[ 3];
  _receiverApertureWidth  = x[ 4];
  _numberOfHeliostats     = myround(x[5]);
  _fieldAngularWidth      = x[ 6];
  _minimumDistanceToTower = x[ 7];
  _maximumDistanceToTower = x[ 8];

  // Htf cycle:
  _centralReceiverOutletTemperature = x[ 9];
  _hotStorageHeight                 = x[10];
  _hotStorageDiameter               = x[11];
  _hotStorageInsulThickness         = x[12];
  _coldStorageInsulThickness        = x[13];
  _coldMoltenSaltMinTemperature     = x[14];
  _receiverNbOfTubes                = myround(x[15]);
  _receiverInsulThickness           = x[16];
  _receiverTubesInsideDiam          = x[17];
  _receiverTubesOutsideDiam         = x[18];

  // Powerblock:
  if ( !set_typeOfTurbine ( myround(x[19]) ) )
    throw std::invalid_argument ( "Problem with input: Type of turbine is not in {1, 2, ..., 8}" );

  // check bounds:
  // -------------
  if ( !check_bounds_minCost_C1() )
    throw std::invalid_argument ( "Problem with input: One of the inputs is outside its bounds" );
  
  return true;
}

/*------------------------------------------*/
/*    initialize problem minCost_C2 (#4)    */
/*------------------------------------------*/
void Scenario::init_minCost_C2 ( double fidelity ) {
 
  // Demand profile : 1, from 3 pm to 9 pm
  // Maximum demand : 25 MW
  // Latitude : 35 deg
  // Day : 1
  // Duration : 24 hours
  // maximum field surface of 200 hectares

  // Scenario parameters:
  _model_type           = 2; // whole plant
  _heliostatsFieldModel = 1;
  _exchangerModel       = 2;
  
  _latitude                = 35.0;
  _day                     = 1;
  _demandProfile           = 1;
  _tStart                  = 900;  // 15x60
  _tEnd                    = 1260; // 21x60
  _maximumPowerDemand      = 25e6;
  _storageStartupCondition = 50;
  _minutesPerTimeIncrement = 60;

  _fixedPointsPrecision   = 0.001;
  _cFieldSurface          = 2000000; // 200 hectares
  _cDemandComplianceRatio = 100;
  _cParasitics            = 0.18;
  
  // variable fidelity surrogate:
  _numberOfTimeIncrements = Scenario::compute_numberOfTimeIncrements (    1,   24, fidelity );
  _raysPerSquareMeters    = Scenario::compute_raysPerSquareMeters    ( 1e-8, 0.01, fidelity );
}

/*------------------------------------------*/
/*  set inputs for problem minCost_C2 (#4)  */
/*------------------------------------------*/
bool Scenario::set_x_minCost_C2 ( const double * x ) {

  // check discrete variables:
  // -------------------------
  if ( !is_int(x[ 5]) ||
       !is_int(x[15]) ||
       !is_int(x[24]) ||
       !is_int(x[25]) ||
       !is_int(x[26]) ||
       !is_int(x[27]) ||
       !is_int(x[28])    )
    throw std::invalid_argument ( "Problem with input: One of the discrete variables has a non-integer value" );
 
  // assign variables:
  // -----------------

  // Heliostats Field:
  _heliostatLength        = x[0];
  _heliostatWidth         = x[1];
  _towerHeight            = x[2];
  _receiverApertureHeight = x[3];
  _receiverApertureWidth  = x[4];
  _numberOfHeliostats     = myround(x[5]);
  _fieldAngularWidth      = x[6];
  _minimumDistanceToTower = x[7];
  _maximumDistanceToTower = x[8];

  // Htf cycle:
  _centralReceiverOutletTemperature = x[ 9];
  _hotStorageHeight                 = x[10];
  _hotStorageDiameter               = x[11];
  _hotStorageInsulThickness         = x[12];
  _coldStorageInsulThickness        = x[13];
  _coldMoltenSaltMinTemperature     = x[14];
  _receiverNbOfTubes                = myround(x[15]);
  _receiverInsulThickness           = x[16];
  _receiverTubesInsideDiam          = x[17];
  _receiverTubesOutsideDiam         = x[18];
  _exchangerTubesSpacing            = x[19];
  _exchangerTubesLength             = x[20];
  _exchangerTubesDin                = x[21];
  _exchangerTubesDout               = x[22];
  _exchangerBaffleCut               = x[23];
  _exchangerNbOfBaffles             = myround(x[24]);
  _exchangerNbOfTubes               = myround(x[25]);
  _exchangerNbOfShells              = myround(x[26]);
  _exchangerNbOfPassesPerShell      = myround(x[27]);

  // Powerblock:
  if ( !set_typeOfTurbine ( myround(x[28]) ) )
    throw std::invalid_argument ( "Problem with input: Type of turbine is not in {1, 2, ..., 8}" );

  // check bounds:
  // -------------
  if ( !check_bounds_minCost_C2() )
    throw std::invalid_argument ( "Problem with input: One of the inputs is outside its bounds" );

  return true;
}

/*--------------------------------------------*/
/*    initialize problem maxComp_HTF1X (#5)   */
/*--------------------------------------------*/
void Scenario::init_maxComp_HTF1 ( void ) {
  
  _model_type           = 2; // whole plant
  _heliostatsFieldModel = 2;
  _exchangerModel       = 2;
  
  _latitude           = 37.5581;
  _day                = 30;
  _demandProfile      = 1;
  _tStart             = 0;
  _tEnd               = 1380; // 23x60
  _maximumPowerDemand = 12e6; // 12 MW
 
  _minutesPerTimeIncrement = 60;
  _fixedPointsPrecision    = 0.001;
  _cBudget                 = 100e6;
  _cParasitics             = 0.18;

  _heliostatLength        = 2.1336;
  _heliostatWidth         = 3.048;
  _towerHeight            = 100.0;
  _receiverApertureHeight = 6.0;
  _receiverApertureWidth  = 6.0;
  _numberOfHeliostats     = 3800;
  _fieldAngularWidth      = 89;
  _minimumDistanceToTower = 0.5;
  _maximumDistanceToTower = 10;
 
  // variable fidelity surrogate: disabled

  _numberOfTimeIncrements = 720;
  _raysPerSquareMeters    = 0.01;
}

/*--------------------------------------------*/
/*  set inputs for problem maxComp_HTF1 (#5)  */
/*--------------------------------------------*/
bool Scenario::set_x_maxComp_HTF1 ( const double * x ) {

  // check discrete variables:
  // -------------------------
  if ( !is_int(x[ 6]) ||
       !is_int(x[15]) ||
       !is_int(x[16]) ||
       !is_int(x[17]) ||
       !is_int(x[18]) ||
       !is_int(x[19])    )
    throw std::invalid_argument ( "Problem with input: One of the discrete variables has a non-integer value" );
 
  // assign variables:
  // -----------------

  // htf cycle:
  _centralReceiverOutletTemperature = x[ 0];
  _hotStorageHeight                 = x[ 1];
  _hotStorageDiameter               = x[ 2];
  _hotStorageInsulThickness         = x[ 3];
  _coldStorageInsulThickness        = x[ 4];
  _coldMoltenSaltMinTemperature     = x[ 5];
  _receiverNbOfTubes                = myround(x[6]);
  _receiverInsulThickness           = x[ 7];
  _receiverTubesInsideDiam          = x[ 8];
  _receiverTubesOutsideDiam         = x[ 9];
  _exchangerTubesSpacing            = x[10];
  _exchangerTubesLength             = x[11];
  _exchangerTubesDin                = x[12];
  _exchangerTubesDout               = x[13];
  _exchangerBaffleCut               = x[14];
  _exchangerNbOfBaffles             = myround(x[15]);
  _exchangerNbOfTubes               = myround(x[16]);
  _exchangerNbOfShells              = myround(x[17]);
  _exchangerNbOfPassesPerShell      = myround(x[18]);
	
  // powerblock:
  if ( !set_typeOfTurbine ( myround(x[19]) ) )
    throw std::invalid_argument ( "Problem with input: Type of turbine is not in {1, 2, ..., 8}" );

  // check bounds:
  // -------------
  if ( !check_bounds_maxComp_HTF1() )
    throw std::invalid_argument ( "Problem with input: One of the inputs is outside its bounds" );
  
  return true;
}

/*-----------------------------------------*/
/*    initialize problem minCost_TS (#6)   */
/*-----------------------------------------*/
void Scenario::init_minCost_TS ( void ) {
  
  _model_type           = 2; // whole plant
  _heliostatsFieldModel = 2;
  _exchangerModel       = 1;

  _latitude                     = 30.05;
  _day                          = 1;
  _demandProfile                = 1;
  _tStart                       =    0;  // 0 * 60
  _tEnd                         = 1380; // 23 * 60
  _maximumPowerDemand           = 120e6;
  _storageStartupCondition      = 50;
  _minutesPerTimeIncrement      = 60;
  _fixedPointsPrecision         = 0.001;
  _cDemandComplianceRatio       = 100.0;
  _coldMoltenSaltMinTemperature = 530;

  set_typeOfTurbine(7);
	
  _heliostatLength          = 9;
  _heliostatWidth           = 9;
  _towerHeight              = 250.0;
  _receiverApertureHeight   = 15;
  _receiverApertureWidth    = 20;
  _numberOfHeliostats       = 12232;
  _fieldAngularWidth        = 65;
  _minimumDistanceToTower   = 1;
  _maximumDistanceToTower   = 10.5;
  _receiverInsulThickness   = 0.5;
  _receiverNbOfTubes        = 85;
  _receiverTubesInsideDiam  = 0.033;
  _receiverTubesOutsideDiam = 0.050;

  // variable fidelity surrogate: disabled

  _numberOfTimeIncrements = 24;
  _raysPerSquareMeters    = 0.01;
}

/*------------------------------------------*/
/*  set inputs for problem minCost_TS (#6)  */
/*------------------------------------------*/
bool Scenario::set_x_minCost_TS ( const double * x ) {

  // assign variables:
  // -----------------
  _centralReceiverOutletTemperature = x[0];
  _hotStorageHeight                 = x[1];
  _hotStorageDiameter               = x[2];
  _hotStorageInsulThickness         = x[3];
  _coldStorageInsulThickness        = x[4];

  // check bounds:
  // -------------
  if ( !check_bounds_minCost_TS() )
    throw std::invalid_argument ( "Problem with input: One of the inputs is outside its bounds" );

  return true;
}

/*-----------------------------------------*/
/*    initialize problem maxEff_RE (#7)    */
/*-----------------------------------------*/
void Scenario::init_maxEff_RE ( double fidelity ) {

  _model_type              = 2; // whole plant
  _heliostatsFieldModel    = 1;
  _exchangerModel          = 1;
  _latitude                = 30.05;
  _day                     = 1;
  _demandProfile           = 1;
  _tStart                  =    0; //  0x60
  _tEnd                    = 1380; // 23x60
  _maximumPowerDemand      = 0;
  _minutesPerTimeIncrement = 60;
  _fixedPointsPrecision    = 0.001;
  _cParasitics             = 0.03;
  _cBudget                 = 45e6;

  set_typeOfTurbine(1);
  
  _hotStorageDiameter           = 25;
  _hotStorageHeight             = 30;
  _hotStorageInsulThickness     = 5;
  _coldStorageInsulThickness    = 5;
  _coldMoltenSaltMinTemperature = 550;

  _heliostatLength        = 9;
  _heliostatWidth         = 9;
  _towerHeight            = 175.0;
  _numberOfHeliostats     = 5000;
  _fieldAngularWidth      = 55;
  _minimumDistanceToTower = 1;
  _maximumDistanceToTower = 7;

  // variable fidelity surrogate:
  _numberOfTimeIncrements = Scenario::compute_numberOfTimeIncrements (    8,   24, fidelity );
  _raysPerSquareMeters    = Scenario::compute_raysPerSquareMeters    ( 1e-7, 0.01, fidelity );
}

/*------------------------------------------*/
/*  set inputs for problem maxEff_RE (#7)   */
/*------------------------------------------*/
bool Scenario::set_x_maxEff_RE ( const double * x ) {

  // check discrete variables:
  // -------------------------
  if ( !is_int(x[3]) )
    throw std::invalid_argument ( "Problem with input: One of the discrete variables has a non-integer value" );

  // assign variables:
  // -----------------
  _receiverApertureHeight           = x[0];
  _receiverApertureWidth            = x[1];
  _centralReceiverOutletTemperature = x[2];
  _receiverNbOfTubes                = myround(x[3]);
  _receiverInsulThickness           = x[4];
  _receiverTubesInsideDiam          = x[5];
  _receiverTubesOutsideDiam         = x[6];

  // check bounds:
  // -------------
  if ( !check_bounds_maxEff_RE() )
    throw std::invalid_argument ( "Problem with input: One of the inputs is outside its bounds" );

  return true;
}

/*--------------------------------------------*/
/*    initialize problem maxHF_minCost (#8)   */
/*--------------------------------------------*/
void Scenario::init_maxHF_minCost ( double fidelity ) {
  
  _model_type              = 2; // whole plant
  _heliostatsFieldModel    = 1;
  _exchangerModel          = 1;
  _latitude                = 45.0;
  _day                     = 1;
  _demandProfile           = 1;
  _tStart                  = 0;    //  0x60;
  _tEnd                    = 1380; // 23x60;
  _maximumPowerDemand      = 0;
  _minutesPerTimeIncrement = 60;
  _storageStartupCondition = 0;
  _fixedPointsPrecision    = 0.001;
  _cFieldSurface           = 4e6;
  _cParasitics             = 0.08; 

  set_typeOfTurbine(1);
  
  _centralReceiverOutletTemperature = 950;
  _hotStorageDiameter               = 25;
  _hotStorageHeight                 = 30;
  _hotStorageInsulThickness         = 5;
  _coldStorageInsulThickness        = 5;
  _coldMoltenSaltMinTemperature     = 550;

  // variable fidelity surrogate:
  _numberOfTimeIncrements = Scenario::compute_numberOfTimeIncrements (   12,   24, fidelity );
  _raysPerSquareMeters    = Scenario::compute_raysPerSquareMeters    ( 1e-8, 0.01, fidelity );
}

/*---------------------------------------------*/
/*  set inputs for problem maxHF_minCost (#8)  */
/*---------------------------------------------*/
bool Scenario::set_x_maxHF_minCost ( const double * x ) {

  // check discrete variables:
  // -------------------------
  if ( !is_int(x[5]) || !is_int(x[9]) )
    throw std::invalid_argument ( "Problem with input: One of the discrete variables has a non-integer value" );
 
  // assign variables:
  // -----------------

  // Heliostats Field:
  _heliostatLength        = x[0];
  _heliostatWidth         = x[1];
  _towerHeight            = x[2];
  _receiverApertureHeight = x[3];
  _receiverApertureWidth  = x[4];
  _numberOfHeliostats     = myround(x[5]);
  _fieldAngularWidth      = x[6];
  _minimumDistanceToTower = x[7];
  _maximumDistanceToTower = x[8];
  	
  // Htf cycle:
  _receiverNbOfTubes        = myround(x[9]);
  _receiverInsulThickness   = x[10];
  _receiverTubesInsideDiam  = x[11];
  _receiverTubesOutsideDiam = x[12];

  // check bounds:
  // -------------
  if ( !check_bounds_maxHF_minCost() )
    throw std::invalid_argument ( "Problem with input: One of the inputs is outside its bounds" );

  return true;
}

/*--------------------------------------------*/
/*    initialize problem maxNrg_minPar (#9)   */
/*--------------------------------------------*/
void Scenario::init_maxNrg_minPar ( double fidelity ) {

  // Demand profile: 1, from 3 pm to 9 pm
  // Maximum demand: 250MW
  // Latitude      : 25 deg
  // Day           : 180
  // Duration      : 24 hours
  // maximum field surface of 5M m^2 (500 hectares)
  
  // Scenario parameters:
 _model_type           = 2; // whole plant
 _heliostatsFieldModel = 1;
 _exchangerModel       = 2;

 _latitude      = 25.0;
 _day           = 180;
 _demandProfile = 1;
 _tStart        = 0;    //  0x60;
 _tEnd          = 1380; // 23x60;

 _maximumPowerDemand      = 250e6;
 _minutesPerTimeIncrement = 60;
 _fixedPointsPrecision    = 0.001;
 _cFieldSurface           = 5e6;
 _cParasitics             = 0.2;
 _cBudget                 = 1.2e9;

 // variable fidelity surrogate:
 _numberOfTimeIncrements = Scenario::compute_numberOfTimeIncrements (    6,   24, fidelity );
 _raysPerSquareMeters    = Scenario::compute_raysPerSquareMeters    ( 1e-7, 0.01, fidelity );
}

/*---------------------------------------------*/
/*  set inputs for problem maxNrg_minPar (#9)  */
/*---------------------------------------------*/
bool Scenario::set_x_maxNrg_minPar ( const double * x ) {

  // check discrete variables:
  // -------------------------
  if ( !is_int(x[ 5]) ||
       !is_int(x[15]) ||
       !is_int(x[24]) ||
       !is_int(x[25]) ||
       !is_int(x[26]) ||
       !is_int(x[27]) ||
       !is_int(x[28])    )
    throw std::invalid_argument ( "Problem with input: One of the discrete variables has a non-integer value" );
 
  // assign variables:
  // -----------------

  // Heliostats Field:
  _heliostatLength        = x[0];
  _heliostatWidth         = x[1];
  _towerHeight            = x[2];
  _receiverApertureHeight = x[3];
  _receiverApertureWidth  = x[4];
  _numberOfHeliostats     = myround(x[5]);
  _fieldAngularWidth      = x[6];
  _minimumDistanceToTower = x[7];
  _maximumDistanceToTower = x[8];

  // Htf cycle:
  _centralReceiverOutletTemperature = x[ 9];
  _hotStorageHeight                 = x[10];
  _hotStorageDiameter               = x[11];
  _hotStorageInsulThickness         = x[12];
  _coldStorageInsulThickness        = x[13];
  _coldMoltenSaltMinTemperature     = x[14];
  _receiverNbOfTubes                = myround(x[15]);
  _receiverInsulThickness           = x[16];
  _receiverTubesInsideDiam          = x[17];
  _receiverTubesOutsideDiam         = x[18];
  _exchangerTubesSpacing            = x[19];
  _exchangerTubesLength             = x[20];
  _exchangerTubesDin                = x[21];
  _exchangerTubesDout               = x[22];
  _exchangerBaffleCut               = x[23];
  _exchangerNbOfBaffles             = myround(x[24]);
  _exchangerNbOfTubes               = myround(x[25]);
  _exchangerNbOfShells              = myround(x[26]);
  _exchangerNbOfPassesPerShell      = myround(x[27]);

  // Powerblock:
  if ( !set_typeOfTurbine(myround(x[28])) )
    throw std::invalid_argument ( "Problem with input: Type of turbine is not in {1, 2, ..., 8}" );

  // check bounds:
  // -------------
  if ( !check_bounds_maxNrg_H1() )
    throw std::invalid_argument ( "Problem with input: One of the inputs is outside its bounds" );

  return true;
}

/*-----------------------------------------------------*/
/*    initialize problem mincost_unconstrained (#10)   */
/*-----------------------------------------------------*/
void Scenario::init_minCost_unconstrained ( double fidelity ) {
  
  _model_type           = 2; // whole plant
  _heliostatsFieldModel = 2;
  _exchangerModel       = 1;

  _latitude                     = 30.05;
  _day                          = 1;
  _demandProfile                = 1;
  _tStart                       =    0;  // 0 * 60
  _tEnd                         = 1380; // 23 * 60
  _maximumPowerDemand           = 120e6;
  _storageStartupCondition      = 50;
  _fixedPointsPrecision         = 0.001;
  _cDemandComplianceRatio       = 100.0;
  _coldMoltenSaltMinTemperature = 530;

  set_typeOfTurbine(7);
	
  _heliostatLength          = 9;
  _heliostatWidth           = 9;
  _towerHeight              = 250.0;
  _receiverApertureHeight   = 15;
  _receiverApertureWidth    = 20;
  _numberOfHeliostats       = 12232;
  _fieldAngularWidth        = 65;
  _minimumDistanceToTower   = 1;
  _maximumDistanceToTower   = 10.5;
  _receiverInsulThickness   = 0.5;
  _receiverNbOfTubes        = 85;
  _receiverTubesInsideDiam  = 0.033;
  _receiverTubesOutsideDiam = 0.050;

  // variable fidelity surrogate:
  // real in ]0;1] --> integer in {1,2,...,60}:
  _minutesPerTimeIncrement = myround(fidelity*60);
  if ( _minutesPerTimeIncrement <= 0 )
    _minutesPerTimeIncrement = 1;
  if ( fidelity < 1.0 && _minutesPerTimeIncrement==60 )
    _minutesPerTimeIncrement = 59;

  _numberOfTimeIncrements = 24;
  _raysPerSquareMeters    = 0.01;  
}

/*------------------------------------------------------*/
/*  set inputs for problem mincost_unconstrained (#10)  */
/*------------------------------------------------------*/
bool Scenario::set_x_minCost_unconstrained ( const double * x ) {

  // assign variables:
  // -----------------
  _centralReceiverOutletTemperature = x[0];
  _hotStorageHeight                 = x[1];
  _hotStorageDiameter               = x[2];
  _hotStorageInsulThickness         = x[3];
  _coldStorageInsulThickness        = x[4];

  // check bounds:
  // -------------
  if ( !check_bounds_minCost_TS() )
    throw std::invalid_argument ( "Problem with input: One of the inputs is outside its bounds" );
  
  return true;
}

/*----------------------------------------------*/
/*  functions to launch simulation and output   */
/*  the values according to the problem chosen  */
/*----------------------------------------------*/
bool Scenario::simulate ( double * outputs, double * intermediate_outputs, double fidelity, bool & cnt_eval ) {
 
  // #1:
  if ( _problem == "MAXNRG_H1" )
    return simulate_maxNrg_H1 ( outputs, cnt_eval );
  
  // #2:
  if ( _problem == "MINSURF_H1" )
    return simulate_minSurf_H1 ( outputs, cnt_eval );

  // #3:
  if ( _problem == "MINCOST_C1" )
    return simulate_minCost_C1 ( outputs, cnt_eval );

  // #4:
  if ( _problem == "MINCOST_C2" )
    return simulate_minCost_C2 ( outputs, cnt_eval );

  // #5:
  if ( _problem == "MAXCOMP_HTF1" )
    return simulate_maxComp_HTF1 ( outputs, cnt_eval );

  // #6:
  if ( _problem == "MINCOST_TS" )
    return simulate_minCost_TS ( outputs, cnt_eval );

  // #7:
  if ( _problem == "MAXEFF_RE" )
    return simulate_maxEff_RE ( outputs, cnt_eval );

  // #8:
  if ( _problem == "MAXHF_MINCOST" )
    return simulate_maxHF_minCost ( outputs, cnt_eval );

  // #9:
  if ( _problem == "MAXNRG_MINPAR" )
    return simulate_maxNrg_minPar ( outputs, cnt_eval );

  // #10:
  if ( _problem == "MINCOST_UNCONSTRAINED" )    
    return simulate_minCost_unconstrained ( outputs, intermediate_outputs, fidelity <  1.0, cnt_eval );

  return false;
}

/*---------------------------*/
/*  simulate_maxNRG_H1 (#1)  */
/*---------------------------*/
bool Scenario::simulate_maxNrg_H1 ( double * outputs, bool & cnt_eval ) {

  // x1: _heliostatLength
  // x2: _heliostatWidth
  // x3: _towerHeight
  // x4: _receiverApertureHeight
  // x5: _receiverApertureWidth
  // x6: _numberOfHeliostats (int)
  // x7: _fieldAngularWidth
  // x8: _minimumDistanceToTower
  // x9: _maximumDistanceToTower

  for ( int i = 0 ; i < 6 ; ++i )
    outputs[i] = 1e20;

  cnt_eval = true;
  
  try {

    // check a priori constraints:
    if ( !check_apriori_constraints_maxNrg_H1() ) {
      cnt_eval = false;     
      throw std::invalid_argument ( "one of the a priori constraints is violated" );
    }
   
    // creating required objects:    
    construct_maxNrg_H1 ( cnt_eval );

    // Launching simulation:
    _powerplant->fSimulatePowerplant ( false );
    
    // Objective function: total energy gathered in kWh:
    outputs[0] = -_powerplant->get_totalEnergyConcentrated();
   
    // c1: check if budget is respected:
    outputs[1] =
      _powerplant->get_costOfHeliostatField() +
      _powerplant->get_costOfTower         () +
      _powerplant->get_costOfReceiver      () -
      _cBudget;

    // c2: check total land area:
    // PI*x3*x3 ( x9*x9 - x8*x8 ) * x7/180 <= 1.95e6
    outputs[2] = PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0)) * _fieldAngularWidth / 180.0 - _cFieldSurface;
   
    // Check basic geometric requirements (c3, c4):
    outputs[3] = 2 * _heliostatLength - _towerHeight;               // c3: 2*x1-x3 <= 0
    outputs[4] = _minimumDistanceToTower - _maximumDistanceToTower; // c4: x8 <= x9

    // c5: check that x6 heliostats can fit in the field:
    // c5: x6 <= _powerplant->get_heliostatField()->get_nb_heliostats()
    // _powerplant->get_heliostatField()->get_nb_heliostats() is the number of heliostats generated in the grid
    outputs[5] = _numberOfHeliostats - 1.0*_powerplant->get_heliostatField()->get_nb_heliostats();
    
  }

  catch ( const std::exception & e ) {

    // output a priori constraints:
    
    // c2: check total land area:
    // PI*x3*x3 ( x9*x9 - x8*x8 ) * x7/180 <= 1.95e6
    outputs[2] = PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0)) * _fieldAngularWidth / 180.0 - _cFieldSurface;
    outputs[3] = 2 * _heliostatLength - _towerHeight;               // c3: 2*x1-x3 <= 0
    outputs[4] = _minimumDistanceToTower - _maximumDistanceToTower; // c4: x8 <= x9
    
    throw Simulation_Interruption ( "Simulation could not go through: " + std::string(e.what()) );
  }
  
  return true;
}

/*----------------------------*/
/*  simulate_minSurf_H1 (#2)  */
/*----------------------------*/
bool Scenario::simulate_minSurf_H1 ( double * outputs , bool & cnt_eval ) {

  // x1 : _heliostatLength
  // x2 : _heliostatWidth
  // x3 : _towerHeight
  // x4 : _receiverApertureHeight
  // x5 : _receiverApertureWidth
  // x6 : _numberOfHeliostats (int)
  // x7 : _fieldAngularWidth
  // x8 : _minimumDistanceToTower
  // x9 : _maximumDistanceToTower
  // x10: _centralReceiverOutletTemperature
  // x11: _receiverNbOfTubes (int)
  // x12: _receiverInsulThickness
  // x13: _receiverTubesInsideDiam
  // x14: _receiverTubesOutsideDiam
  
  for ( int i = 0 ; i < 13 ; ++i )
    outputs[i] = 1e20;

  cnt_eval = true;
    
  try {

    // check a priori constraints:
    if ( !check_apriori_constraints_minSurf_H1() ) {
      cnt_eval = false;
      throw std::invalid_argument ( "one of the a priori constraints (possibly hidden) is violated" );
    }
    
    // Creating required objects:
    construct_minSurf_H1 ( cnt_eval );

    // Launching simulation:
    _powerplant->fSimulatePowerplant ( false );

    // Objective function: field surface (m^2):
    outputs[0] = _fieldAngularWidth * (PI / 180.0) *
      (pow(_towerHeight*_maximumDistanceToTower, 2.0) - pow(_towerHeight*_minimumDistanceToTower, 2.0));
    
    // c1: check if surface constraint is respected: f <= max. surface:
    outputs[1] = _fieldAngularWidth * (PI / 180.0) *
      (pow(_towerHeight*_maximumDistanceToTower, 2.0) - pow(_towerHeight*_minimumDistanceToTower, 2.0))
      - _cFieldSurface;

    // c2: check if demand is met:
    outputs[2] = _cDemandComplianceRatio - _powerplant->get_overallComplianceToDemand();
  
    // c3: check if total cost does not exceed limit:
    outputs[3] = _powerplant->get_costOfHeliostatField()
      + _powerplant->get_costOfTower()
      + _powerplant->get_costOfReceiver()
      + _powerplant->get_costOfStorage()
      + _powerplant->get_costOfSteamGenerator()
      + _powerplant->get_costOfPowerblock()
      - _cBudget;

    // Check basic geometric requirements: c4 and c5:
    outputs[4] = 2 * _heliostatLength - _towerHeight;               // c4: 2x1 - x3 <= 0
    outputs[5] = _minimumDistanceToTower - _maximumDistanceToTower; // c5: x8 <= x9

    // c6: check that _numberOfHeliostats (x6) heliostats can fit in the field:
    outputs[6] = _numberOfHeliostats - 1.0*_powerplant->get_heliostatField()->get_nb_heliostats();
   
    // c7: check pressure in tubes:
    outputs[7] = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();

    // Check if molten salt temperature does not drop below the melting point (c8, c9):
    outputs[ 8] = MELTING_POINT - _powerplant->get_minHotStorageTemp ();
    outputs[ 9] = MELTING_POINT - _powerplant->get_minColdStorageTemp();

    // This was the 10th constraint before version 0.5.4. It is now a hidden constraint.
    // (it seems to be always feasible)
    if ( MELTING_POINT - _powerplant->get_minSteamGenTemp() > 0.0 )
      throw Simulation_Interruption ( "Problem with the molten salt temperature" );

    // c10: x13 <= x14:
    outputs[10] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
   
    // c11: check if tubes fit in receiver: x11*x14 - x5 * PI / 2.0 <= 0:
    outputs[11] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;

    // c12:
    outputs[12] = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;        
  }
  catch ( const std::exception & e ) {

    // Old version (bug): Removed in version 0.4.0:
    // outputs[1] = PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0)
    // 		 - pow(_minimumDistanceToTower*_towerHeight, 2.0))
    //              * (2 * _fieldAngularWidth / 360.0) - _cFieldSurface;
    // outputs[ 4] = 2 * _heliostatLength - _towerHeight;
    // outputs[ 5] = _minimumDistanceToTower - _maximumDistanceToTower;
    // outputs[ 6] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    // outputs[11] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;
    // outputs[13] = _minReceiverOutletTemp - _centralReceiverOutletTemperature;
   
    // corrected version:

    // c1: x7 * (PI / 180.0) * (pow(x3*x9, 2.0) - pow(x3*x8, 2.0)) <=  _cFieldSurface:
    outputs[1] = _fieldAngularWidth * (PI / 180.0) *
      (pow(_towerHeight*_maximumDistanceToTower, 2.0) - pow(_towerHeight*_minimumDistanceToTower, 2.0))
      - _cFieldSurface;
    
    outputs[ 4] = 2 * _heliostatLength - _towerHeight;                  //  c4: 2x1 - x3 <= 0
    outputs[ 5] = _minimumDistanceToTower - _maximumDistanceToTower;    //  c5:  x8 <= x9
    outputs[10] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam; // c10: x13 <= x14

    // c11: check if tubes fit in receiver: x11*x14 - x5 * PI / 2.0 <= 0:
    outputs[11] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;
    
    throw Simulation_Interruption ( "Simulation could not go through: " + std::string(e.what()) );
  }
  
  return true;
}

/*----------------------------*/
/*  simulate_minCost_C1 (#3)  */
/*----------------------------*/
bool Scenario::simulate_minCost_C1 ( double * outputs , bool & cnt_eval ) {

  //  x1: _heliostatLength
  //  x2: _heliostatWidth
  //  x3: _towerHeight
  //  x4: _receiverApertureHeight
  //  x5: _receiverApertureWidth
  //  x6: _numberOfHeliostats (int)
  //  x7: _fieldAngularWidth
  //  x8: _minimumDistanceToTower
  //  x9: _maximumDistanceToTower
  // x10: _centralReceiverOutletTemperature
  // x11: _hotStorageHeight
  // x12: _hotStorageDiameter
  // x13: _hotStorageInsulThickness
  // x14: _coldStorageInsulThickness
  // x15: _coldMoltenSaltMinTemperature
  // x16: _receiverNbOfTubes (int)
  // x17: _receiverInsulThickness
  // x18: _receiverTubesInsideDiam
  // x19: _receiverTubesOutsideDiam
  // x20: _typeOfTurbine (int)
  
  for ( int i = 0 ; i < 14 ; ++i )
    outputs[i] = 1e20;

  cnt_eval = true;
  
  try {

    // check a priori constraints:
    if ( !check_apriori_constraints_minCost_C1() ) {
      cnt_eval = false;
      throw std::invalid_argument ( "one of the a priori constraints (possibly hidden) is violated" );
    }
    
    // creating required objects:
    construct_minCost_C1 ( cnt_eval );
  
    // launching simulation:
    _powerplant->fSimulatePowerplant ( false );
    
    // objective function: total investment cost:
    outputs[0] = _powerplant->get_costOfHeliostatField()
      + _powerplant->get_costOfTower()
      + _powerplant->get_costOfReceiver()
      + _powerplant->get_costOfStorage()
      + _powerplant->get_costOfSteamGenerator()
      + _powerplant->get_costOfPowerblock();
    
    // c1: check total land area is below 800000 m^2:
    // PI * x3 * x3 * ( x9*x9 - x8*x8 ) * x7 / 180 <= 800000 :
    outputs[1] = PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0))
      * _fieldAngularWidth / 180.0 - _cFieldSurface;
    
    // c2: check if compliance to demand is 100%:
    outputs[2] = _cDemandComplianceRatio - _powerplant->get_overallComplianceToDemand();

    // c3: check if tower at least twice as high as heliostats: 2*x1 - x3 <= 0:
    outputs[3] = 2 * _heliostatLength - _towerHeight;

    // c4: check Rmin <= Rmax: x8 <= x9:
    outputs[4] = _minimumDistanceToTower - _maximumDistanceToTower;

    // c5: check that _numberOfHeliostats (x6) heliostats can fit in the field:
    outputs[5] = _numberOfHeliostats - 1.0*_powerplant->get_heliostatField()->get_nb_heliostats();
    
    // c6: check maximum pressure in receiver tubes do not exceed yield pressure:
    outputs[6] = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();

    // check molten Salt temperature does not drop below the melting point: c7, c8, and c9:
    outputs[7] = MELTING_POINT - _powerplant->get_minHotStorageTemp();
    outputs[8] = MELTING_POINT - _powerplant->get_minColdStorageTemp();
    outputs[9] = MELTING_POINT - _powerplant->get_minSteamGenTemp();
   
    // c10: check receiver tubes Din < Dout: x18 <= x19:
    outputs[10] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;

    // c11: check if tubes fit in receiver: x16 * x19 - x5 * PI / 2 <= 0:
    outputs[11] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;
		
    // c12: check central receiver outlet is higher than that required by the turbine:
    outputs[12] = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;

    // c13: check if storage is back to initial conditions:
    outputs[13] = 1.0*_storageStartupCondition
      - (_powerplant->get_moltenSaltLoop()->get_hotStorage().get_heightOfVolumeStored() / _hotStorageHeight);
  }
  catch ( const std::exception & e ) {
    
    // output a priori constraints:
    
    outputs[1] = PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0))
      * _fieldAngularWidth / 180.0 - _cFieldSurface;
   
    outputs[3] = 2 * _heliostatLength - _towerHeight;
    outputs[4] = _minimumDistanceToTower - _maximumDistanceToTower;
    
    outputs[10] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    outputs[11] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;
    // outputs[12] = _minReceiverOutletTemp - _centralReceiverOutletTemperature;  removed in version 0.4.0

    throw Simulation_Interruption ( "Simulation could not go through: " + std::string(e.what()) );
  }
  return true;
}

/*----------------------------*/
/*  simulate_minCost_C2 (#4)  */
/*----------------------------*/
bool Scenario::simulate_minCost_C2 ( double * outputs , bool & cnt_eval ) {

  //  x1: _heliostatLength
  //  x2: _heliostatWidth
  //  x3: _towerHeight
  //  x4: _receiverApertureHeight
  //  x5: _receiverApertureWidth
  //  x6: _numberOfHeliostats (int)
  //  x7: _fieldAngularWidth
  //  x8: _minimumDistanceToTower
  //  x9: _maximumDistanceToTower
  // x10: _centralReceiverOutletTemperature
  // x11: _hotStorageHeight
  // x12: _hotStorageDiameter
  // x13: _hotStorageInsulThickness
  // x14: _coldStorageInsulThickness
  // x15: _coldMoltenSaltMinTemperature
  // x16: _receiverNbOfTubes (int)
  // x17: _receiverInsulThickness
  // x18: _receiverTubesInsideDiam
  // x19: _receiverTubesOutsideDiam
  // x20: _exchangerTubesSpacing
  // x21: _exchangerTubesLength
  // x22: _exchangerTubesDin
  // x23: _exchangerTubesDout
  // x24: _exchangerBaffleCut
  // x25: _exchangerNbOfBaffles (int)
  // x26: _exchangerNbOfTubes (int)
  // x27: _exchangerNbOfShells (int)
  // x28: _exchangerNbOfPassesPerShell (int)
  // x29: _typeOfTurbine (int)
  
  for ( int i = 0 ; i < 17 ; ++i )
    outputs[i] = 1e20;

  cnt_eval = true;
  
  try {

    // check a priori constraints:
    if ( !check_apriori_constraints_minCost_C2() ) {
      cnt_eval = false;
      throw std::invalid_argument ( "one of the a priori constraints (possibly hidden) is violated" );
    }
    
    // Creating required objects:
    construct_minCost_C2 ( cnt_eval );
    
    // Launching simulation:
    _powerplant->fSimulatePowerplant ( false );

    // Objective function: total investment cost:
    outputs[0] = _powerplant->get_costOfHeliostatField()
      + _powerplant->get_costOfTower()
      + _powerplant->get_costOfReceiver()
      + _powerplant->get_costOfStorage()
      + _powerplant->get_costOfSteamGenerator()
      + _powerplant->get_costOfPowerblock();

    // c1: total land area is below 200 hectares:
    // PI*x3*x3(x9*x9- x8*x8) * x7 / 180.0 <= 2000000
    outputs[1] = PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0))
      * (_fieldAngularWidth / 180.0) - _cFieldSurface;
    
    // c2: compliance to demand is 100%:
    outputs[2] = _cDemandComplianceRatio - _powerplant->get_overallComplianceToDemand();

    // c3: tower at least twice as high as heliostats: 2x1-x3 <= 0:
    outputs[3] = 2 * _heliostatLength - _towerHeight;

    // c4: Rmin <= Rmax: x8 <= x9:
    outputs[4]= _minimumDistanceToTower - _maximumDistanceToTower;
    
    // c5: Check that _numberOfHeliostats (x6) heliostats can fit in the field:
    outputs[5] = _numberOfHeliostats - 1.0*_powerplant->get_heliostatField()->get_nb_heliostats();
    
    // c6: maximum pressure in receiver tubes do not exceed yield pressure
    outputs[6] = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();

    // Molten Salt temperature does not drop below the melting point: c7, c8, and c9:
    outputs[7] = MELTING_POINT - _powerplant->get_minHotStorageTemp();
    outputs[8] = MELTING_POINT - _powerplant->get_minColdStorageTemp();
    outputs[9] = MELTING_POINT - _powerplant->get_minSteamGenTemp();
    
    // c10: Receiver tubes Din < Dout: x18 <= x19:
    outputs[10] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;

    // c11: Tubes fit in receiver: x16*x19 <= x5*PI/2:
    outputs[11] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;

    // c12: central receiver outlet is higher than that required by the turbine:
    outputs[12] = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;
 
    // c13: Parasitics do not exceed 20% of energy production
    {
      double sum = 1.0;
      for ( unsigned int i = 0; i < _powerplant->get_powerplantPowerOutput().size(); ++i )
	sum += _powerplant->get_powerplantPowerOutput()[i];
      outputs[13] = _powerplant->fComputeParasiticLosses()/sum - _cParasitics;
    }
   
    // c14: x23 <= x20:
    outputs[14] = _exchangerTubesDout - _exchangerTubesSpacing;

    // c15: x22 <= x23:
    outputs[15] = _exchangerTubesDin  - _exchangerTubesDout;

    // c16: Pressure in steam generator tubes does not exceed yield pressure:
    outputs[16] = _powerplant->get_maximumPressureInExchanger() - _powerplant->get_yieldPressureInExchanger();
  }
  catch ( const std::exception & e ) {

    // output a priori constraints:

    // c1: PI*x3*x3(x9*x9- x8*x8) * x7 / 180.0 <= 2000000:
    outputs[1] = PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0))
      * (_fieldAngularWidth / 180.0) - _cFieldSurface;

    // c3: 2x1-x3 <= 0:
    outputs[3] = 2 * _heliostatLength - _towerHeight;

    // c4: x8 <= x9:
    outputs[4] = _minimumDistanceToTower - _maximumDistanceToTower;

    // c10: x18 <= x19:
    outputs[10] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;

    // c11: x16*x19 <= x5*PI/2:
    outputs[11] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;

    // c12: Removed in version 0.4.0:
    // outputs[12] = _minReceiverOutletTemp - _centralReceiverOutletTemperature;

    // c14: x23 <= x20:
    outputs[14] = _exchangerTubesDout - _exchangerTubesSpacing;

    // c15: x22 <= x23:
    outputs[15] = _exchangerTubesDin  - _exchangerTubesDout;
    
    throw Simulation_Interruption ( "Simulation could not go through: " + std::string(e.what()) );
  }
  
  return true;
}

/*------------------------------*/
/*  simulate_maxComp_HTF1 (#5)  */
/*------------------------------*/
bool Scenario::simulate_maxComp_HTF1 ( double * outputs , bool & cnt_eval ) {

  //  x1: _centralReceiverOutletTemperature
  //  x2: _hotStorageHeight
  //  x3: _hotStorageDiameter
  //  x4: _hotStorageInsulThickness
  //  x5: _coldStorageInsulThickness
  //  x6: _coldMoltenSaltMinTemperature
  //  x7: _receiverNbOfTubes (int)
  //  x8: _receiverInsulThickness
  //  x9: _receiverTubesInsideDiam
  // x10: _receiverTubesOutsideDiam
  // x11: _exchangerTubesSpacing
  // x12: _exchangerTubesLength
  // x13: _exchangerTubesDin
  // x14: _exchangerTubesDout
  // x15: _exchangerBaffleCut
  // x16: _exchangerNbOfBaffles (int)
  // x17: _exchangerNbOfTubes (int)
  // x18: _exchangerNbOfShells (int)
  // x19: _exchangerNbOfPassesPerShell (int)
  // x20: _typeOfTurbine (int)
  
  for ( int i = 0 ; i < 13 ; ++i )
    outputs[i] = 1e20;

  cnt_eval = true;
  
  try {

    // check a priori constraints:
    if ( !check_apriori_constraints_maxComp_HTF1() ) {
      cnt_eval = false;
      throw std::invalid_argument ( "one of the a priori constraints (possibly hidden) is violated" );
    }

    // create required objects:
    construct_maxComp_HTF1 ( cnt_eval );

    // launch simulation:
    _powerplant->fSimulatePowerplant ( false );

    // objective function: time for which the demand is met:
    outputs[0] = - _powerplant->get_overallComplianceToDemand();

    // c1: cost constraint:
    outputs[1] =  _powerplant->get_costOfReceiver()
      + _powerplant->get_costOfStorage()
      + _powerplant->get_costOfSteamGenerator()
      + _powerplant->get_costOfPowerblock()
      - _cBudget;
    
    // c2: pressure in receiver tubes does not exceed yield pressure:
    outputs[2] = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();

    // c3, c4, c5: Molten Salt temperature does not drop below the melting point:
    outputs[3] = MELTING_POINT - _powerplant->get_minHotStorageTemp();
    outputs[4] = MELTING_POINT - _powerplant->get_minColdStorageTemp();
    outputs[5] = MELTING_POINT - _powerplant->get_minSteamGenTemp();
    
    // c6: Receiver tubes Din <= Dout: x9 <= x10:
    outputs[6] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    // c7: Tubes fit in receiver: x7*x10 <= 6*PI/2 (_receiverApertureWidth is fixed to 6):
    outputs[7] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;
    
    // c8: central receiver outlet is higher than that required by the turbine:
    outputs[8] = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;
    
    // c9: parasitics do not exceed limit:
    {
      double sum = 1;
      foncteurSum somme(&sum);
      std::for_each(_powerplant->get_powerplantPowerOutput().begin(),
		    _powerplant->get_powerplantPowerOutput().end(),
		    somme);
      outputs[9] = _powerplant->fComputeParasiticLosses() / sum - _cParasitics;
    }
    
    // c10 and c11: spacing and tubes dimensions in steam gen.
    outputs[10] = _exchangerTubesDout - _exchangerTubesSpacing; // x14 <= x11
    outputs[11] = _exchangerTubesDin  - _exchangerTubesDout;    // x13 <= x14
    
    // c12: pressure in steam gen. tubes do not exceed yield
    outputs[12] = _powerplant->get_maximumPressureInExchanger() - _powerplant->get_yieldPressureInExchanger();
  }
  catch ( const std::exception & e ) {
    
    // output a priori constraints:

    outputs[6] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    outputs[7] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;
    
    // c8: central receiver outlet is higher than that required by the turbine: Removed in version 0.4.0
    // outputs[8] = _minReceiverOutletTemp - _centralReceiverOutletTemperature;
    
    outputs[10] = _exchangerTubesDout - _exchangerTubesSpacing;
    outputs[11] = _exchangerTubesDin  - _exchangerTubesDout;
   
    throw Simulation_Interruption ( "Simulation could not go through: " + std::string(e.what()) );
  }
  return true;
}

/*----------------------------*/
/*  simulate_minCost_TS (#6)  */
/*----------------------------*/
bool Scenario::simulate_minCost_TS ( double * outputs , bool & cnt_eval ) {

  // x1: _centralReceiverOutletTemperature
  // x2: _hotStorageHeight
  // x3: _hotStorageDiameter
  // x4: _hotStorageInsulThickness
  // x5: _coldStorageInsulThickness
  
  for ( int i = 0 ; i < 7 ; ++i )
    outputs[i] = 1e20;

  cnt_eval = true;
  
  try {

    // check a priori constraints:
    if ( !check_apriori_constraints_minCost_TS() ) {
      cnt_eval = false;
      throw std::invalid_argument ( "one of the hidden constraints is violated (detected a priori)" );
    }
    
    // create required objects:
    construct_minCost_TS ( cnt_eval );
    
    // launch simulation:
    _powerplant->fSimulatePowerplant ( false );
   
    // objective function: cost of storage:
    outputs[0] = _powerplant->get_costOfStorage();

    // c1: demand met:
    outputs[1] = _cDemandComplianceRatio - _powerplant->get_overallComplianceToDemand();

    // c2: pressure in receiver tubes does not exceed yield pressure:
    outputs[2] = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();
    
    // c3 and c4: molten Salt temperature does not drop below the melting point:
    outputs[3] = MELTING_POINT - _powerplant->get_minHotStorageTemp ();
    outputs[4] = MELTING_POINT - _powerplant->get_minColdStorageTemp();

    // c5: central receiver outlet is higher than that required by the turbine:
    outputs[5] = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;

    // c6: storage is back to initial conditions:
    outputs[6] = 1.0*_storageStartupCondition
      - 100*(_powerplant->get_moltenSaltLoop()->get_hotStorage().get_heightOfVolumeStored() /_hotStorageHeight);

  }
  catch ( const std::exception & e ) {
    
    // c5: central receiver outlet is higher than that required by the turbine: Removed in version 0.4.0
    // outputs[5] = _minReceiverOutletTemp - _centralReceiverOutletTemperature;
    
    throw Simulation_Interruption ( "Simulation could not go through: " + std::string(e.what()) );
  }
  
  return true;
}

/*-----------------------------*/
/*   simulate_maxEff_Re (#7)   */
/*-----------------------------*/
bool Scenario::simulate_maxEff_RE ( double * outputs , bool & cnt_eval ) {

  // x1: _receiverApertureHeight
  // x2: _receiverApertureWidth
  // x3: _centralReceiverOutletTemperature
  // x4: _receiverNbOfTubes (int)
  // x5: _receiverInsulThickness
  // x6: _receiverTubesInsideDiam
  // x7: _receiverTubesOutsideDiam
  
  for ( int i = 0 ; i < 7 ; ++i )
    outputs[i] = 1e20;

  cnt_eval = true;
  
  try {

    // check a priori constraints:
    if ( !check_apriori_constraints_maxEff_RE() ) {
      cnt_eval = false;
      throw std::invalid_argument ( "one of the a priori constraints (possibly hidden) is violated" );
    }
    
    // create required objects:
    construct_maxEff_RE ( cnt_eval );

    // launch simulation:
    _powerplant->fSimulatePowerplant ( false );

    // objective function: absorbed energy:
    double Q_abs = _powerplant->get_moltenSaltLoop()->get_hotStorage().get_storedMass()
      * (_centralReceiverOutletTemperature - _coldMoltenSaltMinTemperature)
      * HEAT_CAPACITY;
    outputs[0] = -Q_abs*1e-9;  

    // c1: budget is respected:
    outputs[1] = _powerplant->get_costOfReceiver() - _cBudget;
    
    // c2: maximum pressure in tubes <= yield:
    outputs[2] = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();
    
    // c3: receiver inner tubes diameter <= outer diameter: x6 <= x7:
    outputs[3] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    // c4: central receiver outlet is higher than that required by the turbine:
    outputs[4] = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;
    
    // c5: tubes fit in receiver: x4*x7 <= x2*PI/2:
    outputs[5] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;
    
    // c6: work to drive receiver pump does not exceed 5% of the absorbed energy:
    outputs[6] = ( Q_abs > 0.01 ) ? _powerplant->fComputeParasiticsForPb7() / Q_abs - _cParasitics : 1.0;
  }
  catch ( const std::exception & e ) {

    // ouput a priori constraints:
   
    outputs[3] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    // c4: central receiver outlet is higher than that required by the turbine: Removed in version 0.4.0
    // outputs[4] = _minReceiverOutletTemp - _centralReceiverOutletTemperature;
    
    outputs[5] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;
    
    throw Simulation_Interruption ( "Simulation could not go through: " + std::string(e.what()) );
  }
  return true;
}

/*----------------------------------*/
/*    simulate_maxHF_minCost (#8)   */
/*----------------------------------*/
bool Scenario::simulate_maxHF_minCost ( double * outputs , bool & cnt_eval ) {

  //  x1: _heliostatLength
  //  x2: _heliostatWidth
  //  x3: _towerHeight
  //  x4: _receiverApertureHeight
  //  x5: _receiverApertureWidth
  //  x6: _numberOfHeliostats (int)
  //  x7: _fieldAngularWidth
  //  x8: _minimumDistanceToTower
  //  x9: _maximumDistanceToTower
  // x10: _receiverNbOfTubes (int)
  // x11: _receiverInsulThickness
  // x12: _receiverTubesInsideDiam
  // x13: _receiverTubesOutsideDiam
  
  for ( int i = 0 ; i < 11 ; ++i )
    outputs[i] = 1e20;

  cnt_eval = true;
  
  try {

    // check a priori constraints:
    if ( !check_apriori_constraints_maxHF_minCost() ) {
      cnt_eval = false;
      throw std::invalid_argument ( "one of the a priori constraints is violated" );
    }
    
    // create required objects:
    construct_maxHF_minCost ( cnt_eval );
    
    // launch simulation:
    _powerplant->fSimulatePowerplant ( false );

    // objective function #1: absorbed energy:
    double Q_abs = _powerplant->get_moltenSaltLoop()->get_hotStorage().get_storedMass()
      * (_centralReceiverOutletTemperature - _coldMoltenSaltMinTemperature)
      * HEAT_CAPACITY;
    outputs[0] = -Q_abs;

    // objective function #2: total investment cost:
    outputs[1] = _powerplant->get_costOfHeliostatField()
      + _powerplant->get_costOfTower()
      + _powerplant->get_costOfReceiver();
    
    // c1: total land area:
    // PI*x3*x3*( x9*x9 - x8*x8 ) * x7 / 180 <= 4e6:
    outputs[2] = PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0))
      * (_fieldAngularWidth / 180.0) - _cFieldSurface;

    // c2: tower at least twice as high as heliostats: 2x1-x3 <= 0:
    outputs[3] = 2 * _heliostatLength - _towerHeight;
    
    // c3: Rmin < Rmax: x8 <= x9:
    outputs[4] = _minimumDistanceToTower - _maximumDistanceToTower;
   
    // c4: Check that _numberOfHeliostats (x6) heliostats can fit in the field:
    outputs[5] = _numberOfHeliostats - 1.0*_powerplant->get_heliostatField()->get_nb_heliostats();
    
    // c5: maximum pressure in receiver tubes do not exceed yield pressure:
    outputs[6] = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();
    
    // c6: Receiver tubes Din < Dout: x12 <= x13:
    outputs[7] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    // c7: Tubes fit in receiver: x10*x13 <= x5*PI/2:
    outputs[8] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;
    
    // c8: minimal energy production satisfied:
    outputs[9] = 1.44e+12 - Q_abs; // 1.44e+12 = 400e6 * 3600
    
    // c9: work to drive receiver pump does not exceed 8% of the absorbed energy:
    if ( Q_abs > 1.0 )
      outputs[10] = (_powerplant->fComputeParasiticsForPb9() / Q_abs) - _cParasitics;
  }
  catch ( const std::exception & e ) {

    // output a priori constraints:
    
    outputs[2] = PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0))
      * (_fieldAngularWidth / 180.0) - _cFieldSurface;

    outputs[3] = 2 * _heliostatLength - _towerHeight;
    
    outputs[4] = _minimumDistanceToTower - _maximumDistanceToTower;
    
    outputs[7] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    outputs[8] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;
    
    throw Simulation_Interruption ( "Simulation could not go through: " + std::string(e.what()) );
  }
  return true;
}

/*----------------------------------*/
/*    simulate_maxNrg_minPar (#9)   */
/*----------------------------------*/
bool Scenario::simulate_maxNrg_minPar ( double * outputs , bool & cnt_eval ) {

  //  x1: _heliostatLength
  //  x2: _heliostatWidth
  //  x3: _towerHeight
  //  x4: _receiverApertureHeight
  //  x5: _receiverApertureWidth
  //  x6: _numberOfHeliostats (int)
  //  x7: _fieldAngularWidth
  //  x8: _minimumDistanceToTower
  //  x9: _maximumDistanceToTower
  // x10: _centralReceiverOutletTemperature
  // x11: _hotStorageHeight
  // x12: _hotStorageDiameter
  // x13: _hotStorageInsulThickness
  // x14: _coldStorageInsulThickness
  // x15: _coldMoltenSaltMinTemperature
  // x16: _receiverNbOfTubes (int)
  // x17: _receiverInsulThickness
  // x18: _receiverTubesInsideDiam
  // x19: _receiverTubesOutsideDiam
  // x20: _exchangerTubesSpacing
  // x21: _exchangerTubesLength
  // x22: _exchangerTubesDin
  // x23: _exchangerTubesDout
  // x24: _exchangerBaffleCut
  // x25: _exchangerNbOfBaffles (int)
  // x26: _exchangerNbOfTubes (int)
  // x27: _exchangerNbOfShells (int)
  // x28: _exchangerNbOfPassesPerShell (int)
  // x29: _typeOfTurbine (int)
  
  for ( int i = 0 ; i < 19 ; ++i )
    outputs[i] = 1e20;

  cnt_eval = true;
  
  try {

    // check a priori constraints:
    if ( !check_apriori_constraints_maxNrg_H1() ) {
      cnt_eval = false;
      throw std::invalid_argument ( "one of the a priori constraints (possibly hidden) is violated" );
    }
    
    // create required objects:
    construct_maxNrg_minPar ( cnt_eval );

    // launch simulation:
    _powerplant->fSimulatePowerplant ( false );
    
    // objective #1: power output (MWe)
    double sum = 1.0;
    foncteurSum somme(&sum);
    std::for_each ( _powerplant->get_powerplantPowerOutput().begin(),
		    _powerplant->get_powerplantPowerOutput().end(),
		    somme );
    outputs[0] = -sum;

    // objective #2: parasitic losses:
    outputs[1] = _powerplant->fComputeParasiticLosses();
    
    // c1: total investment cost:
    outputs[2] = _powerplant->get_costOfHeliostatField()
      + _powerplant->get_costOfTower()
      + _powerplant->get_costOfReceiver()
      + _powerplant->get_costOfStorage()
      + _powerplant->get_costOfSteamGenerator()
      + _powerplant->get_costOfPowerblock()
      - _cBudget;

    // c2:
    outputs[3] = 4.32e+11 - sum; // 4.32e+11 = 3600.0 * 120.0e6
    
    // c3: total land area: PI*x3*x3*( x9*x9 - x8*x8 ) * x7 / 180 <= 5e6:
    outputs[4] = PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0))
      * (_fieldAngularWidth / 180.0) - _cFieldSurface;
    
    // c4: tower at least twice as high as heliostats: 2x1-x3 <= 0:
    outputs[5] = 2 * _heliostatLength - _towerHeight;
    
    // c5: Rmin < Rmax: x8 <= x9:
    outputs[6] = _minimumDistanceToTower - _maximumDistanceToTower;
    
    // c6: Check that _numberOfHeliostats (x6) heliostats can fit in the field:
    outputs[7] = _numberOfHeliostats - 1.0*_powerplant->get_heliostatField()->get_nb_heliostats();
    
    // c7: maximum pressure in receiver tubes do not exceed yield pressure:
    outputs[8] = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();
        
    // Molten Salt temperature does not drop below the melting point: c8, c9, and c10:
    outputs[ 9] = MELTING_POINT - _powerplant->get_minHotStorageTemp ();
    outputs[10] = MELTING_POINT - _powerplant->get_minColdStorageTemp();
    outputs[11] = MELTING_POINT - _powerplant->get_minSteamGenTemp   ();
    
    // c11: Receiver tubes Din < Dout: x18 <= x19:
    outputs[12] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    // c12: Tubes fit in receiver: x16*x19 <= x5*PI/2:
    outputs[13] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;
    
    // c13: central receiver outlet is higher than that required by the turbine:
    outputs[14] = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;

    // c14: Ratio parasitics vs power output <= 20%:
    if ( sum > 1.0 )
      outputs[15] = outputs[1]/sum - _cParasitics;

    // c15: x23 <= x20:
    outputs[16] = _exchangerTubesDout - _exchangerTubesSpacing;

    // c16: x22 <= x23:
    outputs[17] = _exchangerTubesDin  - _exchangerTubesDout;

    // c17: Pressure in steam generator tubes <= yield pressure:
    outputs[18] = _powerplant->get_maximumPressureInExchanger() - _powerplant->get_yieldPressureInExchanger();
  }
  catch ( const std::exception & e ) {

    // output a priori constraints:
    
    // total land area is below 2000k m^2:
    outputs[4] = PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0))
      * (_fieldAngularWidth / 180.0) - _cFieldSurface;
    
    // tower at least twice as high as heliostats:
    outputs[5] = 2 * _heliostatLength - _towerHeight;
    
    // Rmin <= Rmax:
    outputs[6] = _minimumDistanceToTower - _maximumDistanceToTower;
    
    // Receiver tubes Din <= Dout:
    outputs[12] = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    // Tubes fit in receiver:
    outputs[13] = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0;
    
    // c13: central receiver outlet is higher than that required by the turbine: Removed in version 0.4.0
    // outputs[14] = _minReceiverOutletTemp - _centralReceiverOutletTemperature;   
    
    outputs[16] = _exchangerTubesDout - _exchangerTubesSpacing;
    outputs[17] = _exchangerTubesDin  - _exchangerTubesDout;
    
    throw Simulation_Interruption ( "Simulation could not go through: " + std::string(e.what()) );
  }
  return true;
}

/*----------------------------------------*/
/*  simulate_minCost_unconstrained (#10)  */
/*----------------------------------------*/
bool Scenario::simulate_minCost_unconstrained ( double * outputs              ,
						double * intermediate_outputs ,
						bool     low_fid              ,
						bool   & cnt_eval               ) {
 
  // x1: _centralReceiverOutletTemperature
  // x2: _hotStorageHeight
  // x3: _hotStorageDiameter
  // x4: _hotStorageInsulThickness
  // x5: _coldStorageInsulThickness

  if ( !intermediate_outputs ) {
    cnt_eval = false;
    return false;
  }
  
  // reset output:
  outputs[0] = 1e20;
  for ( int i = 0 ; i < 7 ; ++i )
    intermediate_outputs[i] = 1e20;

  cnt_eval = true;
  
  try {

    // check a priori constraints:
    if ( !check_apriori_constraints_minCost_unconstrained() ) {
      cnt_eval = false;
      throw std::invalid_argument ( "one of the hidden constraints is violated (detected a priori)" );
    }

    // create required objects:
    construct_minCost_unconstrained ( cnt_eval );
   
    // launch simulation:
    _powerplant->fSimulatePowerplant ( low_fid );
    
    // objective function: cost of storage:
    intermediate_outputs[0] = _powerplant->get_costOfStorage();
			    
    // c1: demand met:
    double c1 = intermediate_outputs[1] = _cDemandComplianceRatio - _powerplant->get_overallComplianceToDemand();
    if ( c1 < 0.0 )
      c1 = 0.0;

    // c2: pressure in receiver tubes does not exceed yield pressure:
    double c2 = intermediate_outputs[2] = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();
    if ( c2 < 0.0 )
      c2 = 0.0;
    
    // c3 and c4: molten Salt temperature does not drop below the melting point:
    double c3 = intermediate_outputs[3] = MELTING_POINT - _powerplant->get_minHotStorageTemp ();
    if ( c3 < 0.0 )
      c3 = 0.0;
    
    double c4 = intermediate_outputs[4] = MELTING_POINT - _powerplant->get_minColdStorageTemp();
    if ( c4 < 0.0 )
      c4 = 0.0;
    
    // c5: central receiver outlet is higher than that required by the turbine:
    double c5 = intermediate_outputs[5] = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;
    if ( c5 < 0.0 )
      c5 = 0.0;
    
    // c6: storage is back to initial conditions:
    double c6 = intermediate_outputs[6] = 1.0*_storageStartupCondition
      - 100*(_powerplant->get_moltenSaltLoop()->get_hotStorage().get_heightOfVolumeStored() /_hotStorageHeight);
    if ( c6 < 0.0 )
      c6 = 0.0;
      
    // Aggregate the original objective (cost) and penalties on the constraints violations (+scaling):
    outputs[0] =    intermediate_outputs[0]*1e-6 +
      0.5 * ( pow ( c1     , 2.0 ) +
	      pow ( c2*2e-6, 2.0 ) +
	      pow ( c3     , 2.0 ) +
	      pow ( c4     , 2.0 ) +
	      pow ( c5     , 2.0 ) +
	      pow ( c6     , 2.0 )   );
  }
  catch ( const std::exception & e ) {
    throw Simulation_Interruption ( "Simulation could not go through: " + std::string(e.what()) );
  }
  
  return true;
}

/*-------------------------------------------------*/
/*      check bounds for problem maxNrg_H1 (#1)    */
/*-------------------------------------------------*/
bool Scenario::check_bounds_maxNrg_H1 ( void ) const {

  if ( _heliostatLength < 1 || _heliostatLength > 40 )
    return false;

  if ( _heliostatWidth < 1 || _heliostatWidth > 40 )
    return false;

  if ( _towerHeight < 20 || _towerHeight > 250 )
    return false;

  if ( _receiverApertureHeight < 1 || _receiverApertureHeight > 30 )
    return false;

  if ( _receiverApertureWidth < 1 || _receiverApertureWidth > 30 )
    return false;

  if ( _numberOfHeliostats < 1 )
    return false;

  if ( _fieldAngularWidth < 1 || _fieldAngularWidth > 89 )
    return false;
  
  if ( _minimumDistanceToTower < 0 || _minimumDistanceToTower > 20 )
    return false;

  if ( _maximumDistanceToTower < 1 || _maximumDistanceToTower > 20 )
    return false;

  return true;
}

/*-------------------------------------------------*/
/*    check a priori constraints for instance #1   */
/*-------------------------------------------------*/
bool Scenario::check_apriori_constraints_maxNrg_H1 ( void ) const {
  
  if ( _towerHeight < 2 * _heliostatLength )
    return false;
    
  if ( _maximumDistanceToTower < _minimumDistanceToTower )
    return false;
  
  if ( PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0)) *
       (_fieldAngularWidth / 180.0) > _cFieldSurface )
    return false;
  
  return true;
}

/*-------------------------------------------------*/
/*     check bounds for problem minSurf_H1 (#2)    */
/*-------------------------------------------------*/
bool Scenario::check_bounds_minSurf_H1 ( void ) const {

  if ( _heliostatLength < 1.0 || _heliostatLength > 40 )
    return false;

  if ( _heliostatWidth < 1.0 || _heliostatWidth > 40 )
    return false;

  if ( _towerHeight < 20 || _towerHeight > 250 )
    return false;

  if ( _receiverApertureHeight < 1 || _receiverApertureHeight > 30 )
    return false;

  if ( _receiverApertureWidth < 1 || _receiverApertureWidth > 30 )
    return false;

  if ( _numberOfHeliostats < 1)
    return false;

  if ( _fieldAngularWidth < 1 || _fieldAngularWidth > 89 )
    return false;

  if ( _minimumDistanceToTower < 0 || _minimumDistanceToTower > 20 )
    return false;
	
  if ( _maximumDistanceToTower < 1 || _maximumDistanceToTower > 20 )
    return false;
 
  if ( _centralReceiverOutletTemperature > 995 )
    return false;

  if ( _receiverNbOfTubes < 1 || _receiverNbOfTubes > 9424 )
    return false;

  if ( _receiverInsulThickness < 0.01 || _receiverInsulThickness > 5 )
    return false;

  if ( _receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1 )
    return false;

  if ( _receiverTubesOutsideDiam < 0.005 || _receiverTubesOutsideDiam > 0.1 )
    return false;

  return true;
}

/*-------------------------------------------------*/
/*    check a priori constraints for instance #2   */
/*-------------------------------------------------*/
bool Scenario::check_apriori_constraints_minSurf_H1 ( void ) const {
  
  if ( _centralReceiverOutletTemperature < _minReceiverOutletTemp ) // hidden constraint
    return false;
  
  if ( _towerHeight < 2 * _heliostatLength)
    return false;
  
  if ( _maximumDistanceToTower < _minimumDistanceToTower )
    return false;

  if ( PI*(pow(_maximumDistanceToTower * _towerHeight, 2.0) - pow(_minimumDistanceToTower * _towerHeight, 2.0)) *
       (_fieldAngularWidth / 180.0) > _cFieldSurface )
    return false;  // hidden constraint
  
  if ( _receiverTubesOutsideDiam < _receiverTubesInsideDiam )
    return false;

  if ( _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * PI / 2.0 > 0.0 )
    return false;
  
  return true;
}

/*-------------------------------------------------*/
/*      check bounds for problem minCost_C1 (#3)   */
/*-------------------------------------------------*/
bool Scenario::check_bounds_minCost_C1 ( void ) const {
 
  if ( _heliostatLength < 1 || _heliostatLength > 40 )
    return false;
    
  if ( _heliostatWidth < 1 || _heliostatWidth > 40 )
    return false;
    
  if ( _towerHeight < 20 || _towerHeight > 250 )
    return false;
  
  if ( _receiverApertureHeight < 1 || _receiverApertureHeight > 30 )
    return false;

  if ( _receiverApertureWidth < 1 || _receiverApertureWidth > 30 )
    return false;
  
  if ( _numberOfHeliostats < 1 )
    return false;

  if ( _fieldAngularWidth < 1 || _fieldAngularWidth > 89 )
    return false;

  if ( _minimumDistanceToTower < 0 || _minimumDistanceToTower > 20 )
    return false;

  if ( _maximumDistanceToTower < 1 || _maximumDistanceToTower > 20 )
    return false;

  if ( _centralReceiverOutletTemperature > 995 )
    return false;
  
  if ( _hotStorageHeight < 1 || _hotStorageHeight > 50 )
    return false;

  if ( _hotStorageDiameter < 1 || _hotStorageDiameter > 30 )
    return false;
      
  if ( _hotStorageInsulThickness < 0.01 || _hotStorageInsulThickness > 5 )
    return false;
    
  if ( _coldStorageInsulThickness < 0.01 || _coldStorageInsulThickness > 5 )
    return false;
  
  if ( _coldMoltenSaltMinTemperature < MELTING_POINT || _coldMoltenSaltMinTemperature > 650 )
    return false;
  
  if ( _receiverNbOfTubes < 1 || _receiverNbOfTubes > 9424 )
    return false;

  if ( _receiverInsulThickness < 0.01 || _receiverInsulThickness > 5 )
    return false;
  
  if ( _receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1 )
    return false;
  
  if ( _receiverTubesOutsideDiam < 0.005 || _receiverTubesOutsideDiam > 0.1 )
    return false;

  if ( _typeOfTurbine < 1 || _typeOfTurbine > 8 )
    return false;

  return true;
}

/*-------------------------------------------------*/
/*    check a priori constraints for instance #3   */
/*-------------------------------------------------*/
bool Scenario::check_apriori_constraints_minCost_C1 ( void ) const {

  if ( _centralReceiverOutletTemperature < _minReceiverOutletTemp ) // hidden constraint
    return false;

  if ( _receiverNbOfTubes * _receiverTubesOutsideDiam > PI*_receiverApertureWidth / 2.0 )
    return false;

  if ( _receiverTubesOutsideDiam < _receiverTubesInsideDiam )
    return false;

  if ( _towerHeight < 2 * _heliostatLength )
    return false;
    
  if ( _maximumDistanceToTower < _minimumDistanceToTower )
    return false;
  
  if ( PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0)) *
      (_fieldAngularWidth / 180.0) > _cFieldSurface )
    return false;

  return true;
}

/*-------------------------------------------------*/
/*      check bounds for problem minCost_C2 (#4)   */
/*-------------------------------------------------*/
bool Scenario::check_bounds_minCost_C2 ( void ) const {

  if ( _heliostatLength < 1 || _heliostatLength > 40 )
    return false;

  if ( _heliostatWidth < 1 || _heliostatWidth > 40 )
    return false;

  if ( _towerHeight < 20 || _towerHeight > 250 )
    return false;
   
  if ( _receiverApertureHeight < 1 || _receiverApertureHeight > 30 )
    return false;
  
  if ( _receiverApertureWidth < 1 || _receiverApertureWidth > 30 )
    return false;
  
  if ( _numberOfHeliostats < 1 )
    return false;
  
  if ( _fieldAngularWidth < 1 || _fieldAngularWidth > 89 )
    return false;
  
  if ( _minimumDistanceToTower < 0 || _minimumDistanceToTower > 20 )
    return false;
  
  if ( _maximumDistanceToTower < 1 || _maximumDistanceToTower > 20 )
    return false;
  
  if ( _centralReceiverOutletTemperature > 995 )
    return false;
 
  if ( _hotStorageHeight < 1 || _hotStorageHeight > 50 )
    return false;
 
  if ( _hotStorageDiameter < 1 || _hotStorageDiameter > 30 )
    return false;
 
  if ( _hotStorageInsulThickness < 0.01 || _hotStorageInsulThickness > 5 )
    return false;
  
  if ( _coldStorageInsulThickness < 0.01 || _coldStorageInsulThickness > 5 )
    return false;
  
  if ( _coldMoltenSaltMinTemperature < MELTING_POINT || _coldMoltenSaltMinTemperature > 650 )
    return false;
  
  if ( _receiverNbOfTubes < 1 || _receiverNbOfTubes > 7853 )
    return false;
  
  if ( _receiverInsulThickness < 0.01 || _receiverInsulThickness > 5 )
    return false;

  if ( _receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1 )
    return false;

  if ( _receiverTubesOutsideDiam < 0.006 || _receiverTubesOutsideDiam > 0.1 )
    return false;

  if ( _exchangerTubesSpacing > 0.3 )
    return false;

  if ( _exchangerTubesLength < 0.5 || _receiverTubesOutsideDiam > 10 )
    return false;
  
  if ( _exchangerTubesDin < 0.005 || _exchangerTubesDin > 0.1 )
    return false;

  if ( _exchangerTubesDout < 0.006 || _exchangerTubesDout > 0.1 )
    return false;

  if ( _exchangerNbOfBaffles < 2 )
    return false;

  if ( _exchangerBaffleCut < 0.15 || _exchangerBaffleCut > 0.4 )
    return false;

  if ( _exchangerNbOfTubes < 1 )
    return false;

  if ( _exchangerNbOfShells < 1 || _exchangerNbOfShells > 10 )
    return false;

  if ( _exchangerNbOfPassesPerShell < 1 || _exchangerNbOfPassesPerShell > 9 )
    return false;

  if ( _typeOfTurbine < 1 || _typeOfTurbine >  8 )
    return false;

   return true;
}

/*-------------------------------------------------*/
/*    check a priori constraints for instance #4   */
/*-------------------------------------------------*/
bool Scenario::check_apriori_constraints_minCost_C2 ( void ) const {
  
  if ( _towerHeight < 2 * _heliostatLength )
    return false;
  
  if ( _maximumDistanceToTower < _minimumDistanceToTower )
    return false;
  
  if ( PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0)) *
       (_fieldAngularWidth / 180.0) > _cFieldSurface )
    return false;

  if ( _centralReceiverOutletTemperature < _minReceiverOutletTemp ) // hidden constraint
    return false;
 
  if ( _receiverNbOfTubes * _receiverTubesOutsideDiam > PI*_receiverApertureWidth / 2.0 )
    return false;
  
  if ( _receiverTubesOutsideDiam < _receiverTubesInsideDiam )
    return false;

  if ( _exchangerTubesSpacing < _exchangerTubesDout )
    return false;

  if ( _exchangerTubesDout < _exchangerTubesDin )
    return false;

   return true;
}

/*-------------------------------------------------*/
/*     check bounds for problem maxComp_HTF1 (#5)  */
/*-------------------------------------------------*/
bool Scenario::check_bounds_maxComp_HTF1 ( void ) const {
 
  if ( _centralReceiverOutletTemperature > 995 )
    return false;
  
  if ( _hotStorageHeight < 1 || _hotStorageHeight > 30 )
    return false;
  
  if ( _hotStorageDiameter < 1 || _hotStorageDiameter > 30 )
    return false;
  
  if ( _hotStorageInsulThickness < 0.01 || _hotStorageInsulThickness > 2 )
    return false;
  
  if ( _coldStorageInsulThickness < 0.01 || _coldStorageInsulThickness > 2 )
    return false;
  
  if ( _coldMoltenSaltMinTemperature < MELTING_POINT || _coldMoltenSaltMinTemperature > 650 )
    return false;
  
  if ( _receiverNbOfTubes < 1 || _receiverNbOfTubes > 1884 )
    return false;
  
  if ( _receiverInsulThickness < 0.1 || _receiverInsulThickness > 2 )
    return false;
  
  if ( _receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1 )
    return false;
  
  if ( _receiverTubesOutsideDiam < 0.005 || _receiverTubesOutsideDiam > 0.1 )
    return false;

  if ( _exchangerTubesSpacing > 0.2 )
    return false;
 
  if ( _exchangerTubesLength < 0.5 || _exchangerTubesLength > 10 )
    return false;
  
  if ( _exchangerTubesDin < 0.005 || _exchangerTubesDin > 0.1 )
    return false;
  
  if ( _exchangerTubesDout < 0.006 || _exchangerTubesDout > 0.1 )
    return false;
 
  if ( _exchangerNbOfBaffles < 2 )
    return false;
  
  if ( _exchangerBaffleCut < 0.15 || _exchangerBaffleCut > 0.4 )
    return false;
  
  if ( _exchangerNbOfTubes < 1 )
    return false;
  
  if ( _exchangerNbOfShells < 1 || _exchangerNbOfShells > 10 )
    return false;
  
  if ( _exchangerNbOfPassesPerShell < 1 || _exchangerNbOfPassesPerShell > 9 )
    return false;
  
  if ( _typeOfTurbine < 1 || _typeOfTurbine >  8 )
    return false;

   return true;
}

/*-------------------------------------------------*/
/*    check a priori constraints for instance #5   */
/*-------------------------------------------------*/
bool Scenario::check_apriori_constraints_maxComp_HTF1 ( void ) const {

  if ( _centralReceiverOutletTemperature < _minReceiverOutletTemp ) // hidden constraint
    return false;
  
  if ( _receiverTubesOutsideDiam < _receiverTubesInsideDiam )
    return false;
  
  if ( _exchangerTubesSpacing < _exchangerTubesDout )
    return false;
  
  if ( _exchangerTubesDout < _exchangerTubesDin )
    return false;

  if ( _receiverNbOfTubes * _receiverTubesOutsideDiam > _receiverApertureWidth * PI / 2.0 )
    return false;
  
  return true;
}

/*-------------------------------------------------*/
/*      check bounds for problem minCost_TS (#6)   */
/*-------------------------------------------------*/
bool Scenario::check_bounds_minCost_TS ( void ) const {
 
  if ( _centralReceiverOutletTemperature > 995 )
    return false;
  
  if ( _hotStorageHeight < 2 || _hotStorageHeight > 50 )
    return false;
  
  if ( _hotStorageDiameter < 2 || _hotStorageDiameter > 30 )
    return false;
  
  if ( _hotStorageInsulThickness < 0.01 || _hotStorageInsulThickness > 5 )
    return false;
  
  if ( _coldStorageInsulThickness < 0.01 || _coldStorageInsulThickness > 5 )
    return false;

  return true;
}

/*-------------------------------------------------*/
/*    check a priori constraints for instance #6   */
/*-------------------------------------------------*/
bool Scenario::check_apriori_constraints_minCost_TS ( void ) const {
  if ( _centralReceiverOutletTemperature < _minReceiverOutletTemp ) // hidden constraint
    return false;
  return true;
}

/*-------------------------------------------------*/
/*       check bounds for problem maxEff_RE (#7)   */
/*-------------------------------------------------*/
bool Scenario::check_bounds_maxEff_RE ( void ) const {

  if ( _receiverApertureHeight < 1 || _receiverApertureHeight > 30 )
    return false;

  if ( _receiverApertureWidth < 1 || _receiverApertureWidth > 30 )
    return false;

  if ( _centralReceiverOutletTemperature > 995 )
    return false;
 
  if ( _receiverNbOfTubes < 1 || _receiverNbOfTubes > 8567 )
    return false;
  
  if ( _receiverInsulThickness < 0.01 || _receiverInsulThickness > 5.0 )
    return false;
  
  if ( _receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1 )
    return false;
  
  if ( _receiverTubesOutsideDiam < 0.0055 || _receiverTubesOutsideDiam > 0.1 )
    return false;

  return true;
}

/*-------------------------------------------------*/
/*    check a priori constraints for instance #7   */
/*-------------------------------------------------*/
bool Scenario::check_apriori_constraints_maxEff_RE ( void ) const {

  if ( _centralReceiverOutletTemperature < _minReceiverOutletTemp ) // hidden constraint
    return false;
    
  if ( _receiverNbOfTubes * _receiverTubesOutsideDiam > PI*_receiverApertureWidth / 2.0 )
    return false;
  
  if ( _receiverTubesOutsideDiam < _receiverTubesInsideDiam )
    return false;

   return true;
}

/*-------------------------------------------------*/
/*    check bounds for problem maxHF_minCost (#8)  */
/*-------------------------------------------------*/
bool Scenario::check_bounds_maxHF_minCost ( void ) const {

  if ( _heliostatLength < 1 || _heliostatLength > 40 )
    return false;

  if ( _heliostatWidth < 1 || _heliostatWidth > 40 )
    return false;
  
  if ( _towerHeight < 20 || _towerHeight > 250 )
    return false;
   
  if ( _receiverApertureHeight < 1 || _receiverApertureHeight > 30 )
    return false;
  
  if ( _receiverApertureWidth < 1 || _receiverApertureWidth > 30 )
    return false;
  
  if ( _numberOfHeliostats < 1 )
    return false;
  
  if ( _fieldAngularWidth < 1 || _fieldAngularWidth > 89 )
    return false;
  
  if ( _minimumDistanceToTower < 0 || _minimumDistanceToTower > 20 )
    return false;
  
  if ( _maximumDistanceToTower < 1 || _maximumDistanceToTower > 20 )
    return false;
   
  if ( _receiverNbOfTubes < 1 || _receiverNbOfTubes > 7853 )
    return false;
  
  if ( _receiverInsulThickness < 0.01 || _receiverInsulThickness > 5 )
    return false;
  
  if ( _receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1 )
    return false;
  
  if ( _receiverTubesOutsideDiam < 0.006 || _receiverTubesOutsideDiam > 0.1 )
    return false;

  return true;
}

/*-------------------------------------------------*/
/*    check a priori constraints for instance #8   */
/*-------------------------------------------------*/
bool Scenario::check_apriori_constraints_maxHF_minCost ( void ) const {

  if ( _maximumDistanceToTower < _minimumDistanceToTower )
    return false;
  
  if ( PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0)) *
       (_fieldAngularWidth / 180.0) > _cFieldSurface )
    return false;

  if ( _towerHeight < 2 * _heliostatLength )
    return false;

  if ( _receiverTubesOutsideDiam < _receiverTubesInsideDiam )
    return false;

  if ( _receiverNbOfTubes*_receiverTubesOutsideDiam > _receiverApertureWidth * PI / 2.0 )
    return false;
  
  return true;
}

/*-------------------------------------------------*/
/*    check bounds for problem maxNrg_minPar (#9)  */
/*-------------------------------------------------*/
bool Scenario::check_bounds_maxNrg_minPar ( void ) const {

  if ( _heliostatLength < 1 || _heliostatLength > 40 )
    return false;

  if ( _heliostatWidth  < 1 || _heliostatWidth > 40 )
    return false;

  if ( _towerHeight < 20 || _towerHeight > 250 )
    return false;
  
  if ( _receiverApertureHeight < 1 || _receiverApertureHeight > 30 )
    return false;
  
  if ( _receiverApertureWidth  < 1 || _receiverApertureWidth  > 30 )
    return false;
  
  if ( _fieldAngularWidth < 1 || _fieldAngularWidth > 89 )
    return false;

  if ( _minimumDistanceToTower < 0 || _minimumDistanceToTower > 20 )
    return false;

  if ( _maximumDistanceToTower < 1 || _maximumDistanceToTower > 20 )
    return false;

  if ( _numberOfHeliostats < 1 )
    return false;

  if ( _centralReceiverOutletTemperature > 995 )
    return false;
 
  if ( _hotStorageHeight < 1 || _hotStorageHeight > 50 )
    return false;

  if ( _hotStorageDiameter < 1 || _hotStorageDiameter > 30 )
    return false;

  if ( _hotStorageInsulThickness < 0.01 || _hotStorageInsulThickness > 5 )
    return false;

  if ( _coldStorageInsulThickness < 0.01 || _coldStorageInsulThickness > 5 )
    return false;

  if ( _coldMoltenSaltMinTemperature < MELTING_POINT || _coldMoltenSaltMinTemperature > 650 )
    return false;

  if ( _receiverNbOfTubes < 1 || _receiverNbOfTubes > 7853 )
    return false;

  if ( _receiverInsulThickness < 0.01 || _receiverInsulThickness > 5 )
    return false;

  if ( _receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1 )
    return false;

  if ( _receiverTubesOutsideDiam < 0.006 || _receiverTubesOutsideDiam > 0.1 )
    return false;
  
  if ( _receiverTubesOutsideDiam < _receiverTubesInsideDiam + 0.0001 )
    return false;

  if ( _exchangerTubesSpacing < 0.007 || _exchangerTubesSpacing > 0.2 )
    return false;

  if ( _exchangerTubesLength < 0.5 || _receiverTubesOutsideDiam > 10 )
    return false;

  if ( _exchangerTubesDin < 0.005 || _exchangerTubesDin > 0.1 )
    return false;
  
  if ( _exchangerTubesDout < 0.006 || _exchangerTubesDout > 0.1 )
    return false;
 
  if ( _exchangerNbOfBaffles < 2 )
    return false;
  
  if ( _exchangerBaffleCut < 0.15 || _exchangerBaffleCut > 0.4 )
    return false;

  if ( _exchangerNbOfTubes < 1 )
    return false;

  if ( _exchangerNbOfShells < 1 || _exchangerNbOfShells > 10 )
    return false;
  
  if ( _exchangerNbOfPassesPerShell < 1 || _exchangerNbOfPassesPerShell > 9 )
    return false;

  if ( _typeOfTurbine < 1 || _typeOfTurbine >  8 )
    return false;

   return true;
}

/*-------------------------------------------------*/
/*    check a priori constraints for instance #9   */
/*-------------------------------------------------*/
bool Scenario::check_apriori_constraints_maxNrg_minPar ( void ) const {

  if ( _centralReceiverOutletTemperature < _minReceiverOutletTemp ) // hidden constraint
    return false;

  if ( _receiverNbOfTubes * _receiverTubesOutsideDiam > PI*_receiverApertureWidth / 2.0 )
    return false;
  
  if ( _towerHeight < 2 * _heliostatLength )
    return false;
    
  if ( _maximumDistanceToTower < _minimumDistanceToTower )
    return false;
  
  if ( PI*(pow(_maximumDistanceToTower*_towerHeight, 2.0) - pow(_minimumDistanceToTower*_towerHeight, 2.0))
       *(_fieldAngularWidth / 180.0) > _cFieldSurface )
    return false;
 
  if ( _exchangerTubesDout < _exchangerTubesDin )
    return false;
  
  if ( _exchangerTubesSpacing < _exchangerTubesDout )
    return false;

  if ( _receiverTubesInsideDiam > _receiverTubesOutsideDiam )
    return false;
  
  return true;
}

/*-------------------------------------------------------------*/
/*         validate problem mincost_unconstrained (#10)        */
/*-------------------------------------------------------------*/
bool Scenario::check_bounds_minCost_unconstrained ( void ) const {
 
  if ( _centralReceiverOutletTemperature > 995 )
    return false;
  
  if ( _hotStorageHeight < 2 || _hotStorageHeight > 50 )
    return false;
  
  if ( _hotStorageDiameter < 2 || _hotStorageDiameter > 30 )
    return false;
  
  if ( _hotStorageInsulThickness < 0.01 || _hotStorageInsulThickness > 5 )
    return false;
  
  if ( _coldStorageInsulThickness < 0.01 || _coldStorageInsulThickness > 5 )
    return false;

  return true;
}

/*--------------------------------------------------*/
/*    check a priori constraints for instance #10   */
/*--------------------------------------------------*/
bool Scenario::check_apriori_constraints_minCost_unconstrained ( void ) const {

  if ( _centralReceiverOutletTemperature < _minReceiverOutletTemp ) // hidden constraint
    return false;

  return true;
}

/*-----------------------------------------------*/
/*      constructing model components (1/10)     */
/*-----------------------------------------------*/
void Scenario::construct_maxNrg_H1 ( bool & cnt_eval ) {

  if ( _powerplant ) {
    delete _powerplant;
    _powerplant = NULL;
  }
  
  Time_Manager time ( _numberOfTimeIncrements , 0 , _minutesPerTimeIncrement );
  Sun          sun  ( _latitude, time, _day, _raysPerSquareMeters );

  HeliostatField * field     = NULL;
  Economics      * economics = NULL;
  
  try {
    
    field = new HeliostatField ( _numberOfHeliostats     ,
				 _heliostatLength        ,
				 _heliostatWidth         ,
				 _towerHeight            ,
				 _receiverApertureHeight ,
				 _receiverApertureWidth  ,
				 _minimumDistanceToTower ,
				 _maximumDistanceToTower ,
				 _fieldAngularWidth      ,
				 sun                       );

    economics = new Economics;
    economics->set_heightOfTower            ( _towerHeight            );
    economics->set_heightOfReceiverAperture ( _receiverApertureHeight );
    economics->set_widthOfReceiverAperture  ( _receiverApertureWidth  );
    economics->set_lengthOfHeliostats       ( _heliostatLength        );
    economics->set_widthOfHeliostats        ( _heliostatWidth         );
    economics->set_exchangerModel           ( _exchangerModel         );
  }
  catch ( const std::exception & e ) {   
    if ( field     ) delete field;
    if ( economics ) delete economics;
    cnt_eval = false;
    throw Simulation_Interruption ( e.what() );
  }
    
  _powerplant = new Powerplant ( time        ,
				 _model_type ,
				 field       ,
				 NULL        ,
				 NULL        ,
				 economics     );
  
  _powerplant->set_heliostatModel ( _heliostatsFieldModel );
}

/*-----------------------------------------------*/
/*      constructing model components (2/10)     */
/*-----------------------------------------------*/
void Scenario::construct_minSurf_H1 ( bool & cnt_eval ) {

  if ( _powerplant ) {
    delete _powerplant;
    _powerplant = NULL;
  }
  
  _model_type = 2; // whole plant

  Time_Manager time ( _numberOfTimeIncrements, 0, _minutesPerTimeIncrement );
  Sun          sun  ( _latitude, time, _day, _raysPerSquareMeters );
  fFillDemandVector();

  HeliostatField * field      = NULL;
  HtfCycle       * htfCycle   = NULL;
  Powerblock     * powerblock = NULL;
  Economics      * economics  = NULL;

  try {
    
    // powerblock model:
    powerblock = new Powerblock ( _typeOfTurbine );

    // constructing heliostats field:
    field = new HeliostatField ( _numberOfHeliostats     ,
				 _heliostatLength        ,
				 _heliostatWidth         ,
				 _towerHeight            ,
				 _receiverApertureHeight ,
				 _receiverApertureWidth  ,
				 _minimumDistanceToTower ,
				 _maximumDistanceToTower ,
				 _fieldAngularWidth      ,
				 sun                       );

    // constructing Htf Cycle with desired steam generator model:
    htfCycle = new HtfCycle ( _centralReceiverOutletTemperature ,
			      _hotStorageHeight                 ,
			      _hotStorageDiameter               ,
			      _hotStorageInsulThickness         ,
			      _coldMoltenSaltMinTemperature     ,
			      powerblock                        ,
			      _receiverApertureHeight           ,
			      _receiverApertureWidth            ,
			      _receiverTubesInsideDiam          ,
			      _receiverTubesOutsideDiam         ,
			      _receiverNbOfTubes                ,
			      _minutesPerTimeIncrement            );

    htfCycle->setStorage ( _storageStartupCondition               ,
			   0.98*_centralReceiverOutletTemperature ,
			   _coldMoltenSaltMinTemperature            );

    // investment cost model:
    economics = new Economics;
    economics->set_heightOfTower                  ( _towerHeight                     );
    economics->set_heightOfReceiverAperture       ( _receiverApertureHeight          );
    economics->set_widthOfReceiverAperture        ( _receiverApertureWidth           );
    economics->set_receiverNumberOfTubes          ( _receiverNbOfTubes               );
    economics->set_receiverTubesDout              ( _receiverTubesOutsideDiam        );
    economics->set_lengthOfHeliostats             ( _heliostatLength                 );
    economics->set_widthOfHeliostats              ( _heliostatWidth                  );
    economics->set_hotStorageHeight               ( _hotStorageHeight                );
    economics->set_storageDiameter                ( _hotStorageDiameter              );
    economics->set_hotStorageInsulationThickness  ( _hotStorageInsulThickness        );
    economics->set_coldStorageInsulationThickness ( _coldStorageInsulThickness       );
    economics->set_receiverInsulationThickness    ( _receiverInsulThickness          );
    economics->set_turbineNominalPowerOutput      ( powerblock->get_powerOfTurbine() );
    economics->set_exchangerModel                 ( _exchangerModel                  );
  }
  catch ( const std::exception & e ) {
    if ( field      ) delete field;
    if ( htfCycle   ) delete htfCycle;
    if ( powerblock ) delete powerblock;
    if ( economics  ) delete economics;
    cnt_eval = false;
    throw Simulation_Interruption ( e.what() );
  }
    
  _powerplant = new Powerplant ( time        ,
				 _model_type ,
				 field       ,
				 htfCycle    ,
				 powerblock  ,
				 economics     );
  _powerplant->set_demand         ( _demand               );
  _powerplant->set_heliostatModel ( _heliostatsFieldModel );
}

/*-----------------------------------------------*/
/*      constructing model components (3/10)     */
/*-----------------------------------------------*/
void Scenario::construct_minCost_C1 ( bool & cnt_eval ) {

  if ( _powerplant ) {
    delete _powerplant;
    _powerplant = NULL;
  }
  
  _model_type = 2; // whole plant

  Time_Manager time ( _numberOfTimeIncrements, 0, _minutesPerTimeIncrement) ;
  Sun          sun  ( _latitude, time, _day, _raysPerSquareMeters );
  fFillDemandVector();

  HeliostatField * field      = NULL;
  HtfCycle       * htfCycle   = NULL;
  Powerblock     * powerblock = NULL;
  Economics      * economics  = NULL;

  try {
    
    // powerblock model:
    powerblock = new Powerblock ( _typeOfTurbine );

    // constructing heliostats field:
    field = new HeliostatField ( _numberOfHeliostats     ,
				 _heliostatLength        ,
				 _heliostatWidth         ,
				 _towerHeight            ,
				 _receiverApertureHeight ,
				 _receiverApertureWidth  ,
				 _minimumDistanceToTower ,
				 _maximumDistanceToTower ,
				 _fieldAngularWidth      ,
				 sun                       );

    // constructing Htf Cycle with desired steam generator model:
    htfCycle = new HtfCycle ( _centralReceiverOutletTemperature ,
			      _hotStorageHeight                 ,
			      _hotStorageDiameter               ,
			      _hotStorageInsulThickness         ,
			      _coldMoltenSaltMinTemperature     ,
			      powerblock                        ,
			      _receiverApertureHeight           ,
			      _receiverApertureWidth            ,
			      _receiverTubesInsideDiam          ,
			      _receiverTubesOutsideDiam         ,
			      _receiverNbOfTubes                ,
			      _minutesPerTimeIncrement 	          );
    htfCycle->initiateColdStorage();
  
    // investment cost model:
    economics = new Economics;
    economics->set_heightOfTower                  ( _towerHeight                     );
    economics->set_heightOfReceiverAperture       ( _receiverApertureHeight          );
    economics->set_widthOfReceiverAperture        ( _receiverApertureWidth           );
    economics->set_receiverNumberOfTubes          ( _receiverNbOfTubes               );
    economics->set_receiverTubesDout              ( _receiverTubesOutsideDiam        );
    economics->set_lengthOfHeliostats             ( _heliostatLength                 );
    economics->set_widthOfHeliostats              ( _heliostatWidth                  );
    economics->set_hotStorageHeight               ( _hotStorageHeight                );
    economics->set_storageDiameter                ( _hotStorageDiameter              );
    economics->set_hotStorageInsulationThickness  ( _hotStorageInsulThickness        );
    economics->set_coldStorageInsulationThickness ( _coldStorageInsulThickness       );
    economics->set_receiverInsulationThickness    ( _receiverInsulThickness          );
    economics->set_turbineNominalPowerOutput      ( powerblock->get_powerOfTurbine() );
    economics->set_exchangerModel                 ( _exchangerModel                  );
  }
  catch ( const std::exception & e ) {
    if ( field      ) delete field;
    if ( htfCycle   ) delete htfCycle;
    if ( powerblock ) delete powerblock;
    if ( economics  ) delete economics;
    cnt_eval = false;
    throw Simulation_Interruption ( e.what() );
  }

  _powerplant = new Powerplant ( time        ,
				 _model_type ,
				 field       ,
				 htfCycle    ,
				 powerblock  ,
				 economics     );
  
  _powerplant->set_demand              ( _demand               );
  _powerplant->set_heliostatModel      ( _heliostatsFieldModel );
}

/*-----------------------------------------------*/
/*      constructing model components (4/10)     */
/*-----------------------------------------------*/
void Scenario::construct_minCost_C2 ( bool & cnt_eval ) {

  if ( _powerplant ) {
    delete _powerplant;
    _powerplant = NULL;
  }
  
  _model_type = 2; // whole plant
	
  Time_Manager time ( _numberOfTimeIncrements, 0, _minutesPerTimeIncrement );
  Sun          sun  ( _latitude, time, _day, _raysPerSquareMeters );
  fFillDemandVector();
	
  HeliostatField * field      = NULL;
  HtfCycle       * htfCycle   = NULL;
  Powerblock     * powerblock = NULL;
  Economics      * economics  = NULL;

  try {
    
    // powerblock model:
    powerblock = new Powerblock ( _typeOfTurbine );

    // constructing heliostats field:
    field = new HeliostatField ( _numberOfHeliostats     ,
				 _heliostatLength        ,
				 _heliostatWidth         ,
				 _towerHeight            ,
				 _receiverApertureHeight ,
				 _receiverApertureWidth  ,
				 _minimumDistanceToTower ,
				 _maximumDistanceToTower ,
				 _fieldAngularWidth      ,
				 sun                       );
  
    // constructing Htf Cycle with desired steam generator model:
    htfCycle = new HtfCycle ( _centralReceiverOutletTemperature ,
			      _hotStorageHeight                 ,
			      _hotStorageDiameter               ,
			      _hotStorageInsulThickness         ,
			      _coldMoltenSaltMinTemperature     ,
			      powerblock                        ,
			      _receiverApertureHeight           ,
			      _receiverApertureWidth            ,
			      _receiverTubesInsideDiam          ,
			      _receiverTubesOutsideDiam         ,
			      _receiverNbOfTubes                ,
			      _minutesPerTimeIncrement          ,
			      _exchangerTubesDin                , 
			      _exchangerTubesDout               , 
			      _exchangerTubesLength             ,
			      _exchangerTubesSpacing            ,
			      _exchangerBaffleCut               ,
			      _exchangerNbOfBaffles             ,
			      _exchangerNbOfTubes               ,
			      _exchangerNbOfPassesPerShell      ,
			      _exchangerNbOfShells                );
    
    htfCycle->initiateColdStorage();
    
    // investment cost model:
    economics = new Economics;
    economics->set_heightOfTower                  ( _towerHeight                     );
    economics->set_heightOfReceiverAperture       ( _receiverApertureHeight          );
    economics->set_widthOfReceiverAperture        ( _receiverApertureWidth           );
    economics->set_receiverNumberOfTubes          ( _receiverNbOfTubes               );
    economics->set_receiverTubesDout              ( _receiverTubesOutsideDiam        );
    economics->set_lengthOfHeliostats             ( _heliostatLength                 );
    economics->set_widthOfHeliostats              ( _heliostatWidth                  );
    economics->set_hotStorageHeight               ( _hotStorageHeight                );
    economics->set_storageDiameter                ( _hotStorageDiameter              );
    economics->set_hotStorageInsulationThickness  ( _hotStorageInsulThickness        );
    economics->set_coldStorageInsulationThickness ( _coldStorageInsulThickness       );
    economics->set_receiverInsulationThickness    ( _receiverInsulThickness          );
    economics->set_turbineNominalPowerOutput      ( powerblock->get_powerOfTurbine() );
    economics->set_exchangerModel                 ( _exchangerModel                  );
  }
  catch ( const std::exception & e ) {
    if ( field      ) delete field;
    if ( htfCycle   ) delete htfCycle;
    if ( powerblock ) delete powerblock;
    if ( economics  ) delete economics;
    cnt_eval = false;
    throw Simulation_Interruption ( e.what() );
  }
    
  _powerplant = new Powerplant ( time        ,
				 _model_type ,
				 field       ,
				 htfCycle    ,
				 powerblock  ,
				 economics     );

  _powerplant->set_demand         ( _demand               );
  _powerplant->set_heliostatModel ( _heliostatsFieldModel );
}

/*-----------------------------------------------*/
/*      constructing model components (5/10)     */
/*-----------------------------------------------*/
void Scenario::construct_maxComp_HTF1 ( bool & cnt_eval ) {

  if ( _powerplant ) {
    delete _powerplant;
    _powerplant = NULL;
  }
  
  _model_type = 2; // whole plant

  Time_Manager time ( _numberOfTimeIncrements, 0, _minutesPerTimeIncrement );
  Sun          sun  ( _latitude, time, _day, _raysPerSquareMeters        );
  fFillDemandVector();

  HeliostatField * field      = NULL;
  HtfCycle       * htfCycle   = NULL;
  Powerblock     * powerblock = NULL;
  Economics      * economics  = NULL;

  try {
    
    // powerblock model:
    powerblock = new Powerblock ( _typeOfTurbine );

    // constructing heliostats field:
    field = new HeliostatField ( _numberOfHeliostats     ,
				 _heliostatLength        ,
				 _heliostatWidth         ,
				 _towerHeight            ,
				 _receiverApertureHeight ,
				 _receiverApertureWidth  ,
				 _minimumDistanceToTower ,
				 _maximumDistanceToTower ,
				 _fieldAngularWidth      ,
				 sun                       );

    // constructing Htf Cycle with desired steam generator model:
    htfCycle = new HtfCycle ( _centralReceiverOutletTemperature ,
			      _hotStorageHeight                 ,
			      _hotStorageDiameter               ,
			      _hotStorageInsulThickness         ,
			      _coldMoltenSaltMinTemperature     ,
			      powerblock                        ,
			      _receiverApertureHeight           ,
			      _receiverApertureWidth            ,
			      _receiverTubesInsideDiam          ,
			      _receiverTubesOutsideDiam         ,
			      _receiverNbOfTubes                ,
			      _minutesPerTimeIncrement          ,
			      _exchangerTubesDin                ,
			      _exchangerTubesDout               ,
			      _exchangerTubesLength             ,
			      _exchangerTubesSpacing            ,
			      _exchangerBaffleCut               ,
			      _exchangerNbOfBaffles             ,
			      _exchangerNbOfTubes               ,
			      _exchangerNbOfPassesPerShell      ,
			      _exchangerNbOfShells 	        );
    
    htfCycle->initiateColdStorage();

    // investment cost model:
    economics = new Economics;
    economics->set_heightOfTower                  ( _towerHeight                     );
    economics->set_heightOfReceiverAperture       ( _receiverApertureHeight          );
    economics->set_widthOfReceiverAperture        ( _receiverApertureWidth           );
    economics->set_receiverNumberOfTubes          ( _receiverNbOfTubes               );
    economics->set_receiverTubesDout              ( _receiverTubesOutsideDiam        );
    economics->set_lengthOfHeliostats             ( _heliostatLength                 );
    economics->set_widthOfHeliostats              ( _heliostatWidth                  );
    economics->set_hotStorageHeight               ( _hotStorageHeight                );
    economics->set_storageDiameter                ( _hotStorageDiameter              );
    economics->set_hotStorageInsulationThickness  ( _hotStorageInsulThickness        );
    economics->set_coldStorageInsulationThickness ( _coldStorageInsulThickness       );
    economics->set_receiverInsulationThickness    ( _receiverInsulThickness          );
    economics->set_turbineNominalPowerOutput      ( powerblock->get_powerOfTurbine() );
    economics->set_exchangerModel                 ( _exchangerModel                  );
  }
  catch ( const std::exception & e ) {
    if ( field      ) delete field;
    if ( htfCycle   ) delete htfCycle;
    if ( powerblock ) delete powerblock;
    if ( economics  ) delete economics;
    cnt_eval = false;
    throw Simulation_Interruption ( e.what() );
  }
       
  _powerplant = new Powerplant ( time        ,
				 _model_type ,
				 field       ,
				 htfCycle    ,
				 powerblock  ,
				 economics     );

  _powerplant->set_demand         ( _demand               );
  _powerplant->set_heliostatModel ( _heliostatsFieldModel );
		
  if  ( !_powerplant->set_heliostatFieldPowerOutput_MAXCOMP_HTF1() ) {
    delete _powerplant;
    _powerplant = NULL;
    cnt_eval = false;
    throw Simulation_Interruption ( "error in the construction of the problem: Powerplant initialization" );
  }
}

/*-----------------------------------------------*/
/*      constructing model components (6/10)     */
/*-----------------------------------------------*/
void Scenario::construct_minCost_TS ( bool & cnt_eval ) {

  if ( _powerplant ) {
    delete _powerplant;
    _powerplant = NULL;
  }
  
  _model_type = 2; // whole plant

  Time_Manager time ( _numberOfTimeIncrements, 0, _minutesPerTimeIncrement);
  Sun          sun  ( _latitude, time, _day, _raysPerSquareMeters);
  fFillDemandVector();

  HtfCycle       * htfCycle   = NULL;
  Powerblock     * powerblock = NULL;
  Economics      * economics  = NULL;

  try {
    
    // powerblock model:
    powerblock = new Powerblock ( _typeOfTurbine );

    // constructing Htf Cycle with desired steam generator model:
    htfCycle = new HtfCycle ( _centralReceiverOutletTemperature ,
			      _hotStorageHeight                 ,
			      _hotStorageDiameter               ,
			      _hotStorageInsulThickness         ,
			      _coldMoltenSaltMinTemperature     ,
			      powerblock                        ,
			      _receiverApertureHeight           ,
			      _receiverApertureWidth            ,
			      _receiverTubesInsideDiam          ,
			      _receiverTubesOutsideDiam         ,
			      _receiverNbOfTubes                ,
			      _minutesPerTimeIncrement            );
  
    htfCycle->setStorage ( _storageStartupCondition               ,
			   0.98*_centralReceiverOutletTemperature ,
			   _coldMoltenSaltMinTemperature            );

    // investment cost model:
    economics = new Economics;
    economics->set_heightOfTower                  ( _towerHeight                     );
    economics->set_heightOfReceiverAperture       ( _receiverApertureHeight          );
    economics->set_widthOfReceiverAperture        ( _receiverApertureWidth           );
    economics->set_receiverNumberOfTubes          ( _receiverNbOfTubes               );
    economics->set_receiverTubesDout              ( _receiverTubesOutsideDiam        );
    economics->set_lengthOfHeliostats             ( _heliostatLength                 );
    economics->set_widthOfHeliostats              ( _heliostatWidth                  );
    economics->set_hotStorageHeight               ( _hotStorageHeight                );
    economics->set_storageDiameter                ( _hotStorageDiameter              );
    economics->set_hotStorageInsulationThickness  ( _hotStorageInsulThickness        );
    economics->set_coldStorageInsulationThickness ( _coldStorageInsulThickness       );
    economics->set_receiverInsulationThickness    ( _receiverInsulThickness          );
    economics->set_turbineNominalPowerOutput      ( powerblock->get_powerOfTurbine() );
    economics->set_exchangerModel                 ( _exchangerModel                  );
  }
  catch ( const std::exception & e ) {
    if ( htfCycle   ) delete htfCycle;
    if ( powerblock ) delete powerblock;
    if ( economics  ) delete economics;
    cnt_eval = false;
    throw Simulation_Interruption ( e.what() );
  }
    
  _powerplant = new Powerplant ( time        ,
				 _model_type ,
				 NULL        ,
				 htfCycle    ,
				 powerblock  ,
				 economics     );

  // setting demand vector:
  _powerplant->set_demand         ( _demand               );
  _powerplant->set_heliostatModel ( _heliostatsFieldModel );

  if  ( !_powerplant->set_heliostatFieldPowerOutput_MINCOST_TS() ) {
    delete _powerplant;
    _powerplant = NULL;
    cnt_eval = false;
    throw Simulation_Interruption ( "error in the construction of the problem" );
  }
}

/*-----------------------------------------------*/
/*      constructing model components (7/10)     */
/*-----------------------------------------------*/
void Scenario::construct_maxEff_RE ( bool & cnt_eval ) {

  if ( _powerplant ) {
    delete _powerplant;
    _powerplant = NULL;
  }
  
  _model_type = 2; // whole plant

  Time_Manager time ( _numberOfTimeIncrements, 0, _minutesPerTimeIncrement);
  Sun          sun  ( _latitude, time, _day, _raysPerSquareMeters);
  fFillDemandVector();

  HeliostatField * field      = NULL;
  HtfCycle       * htfCycle   = NULL;
  Powerblock     * powerblock = NULL;
  Economics      * economics  = NULL;

  try {
    
    // powerblock model:
    powerblock = new Powerblock ( _typeOfTurbine );

    // constructing heliostats field:
    field = new HeliostatField ( _numberOfHeliostats     ,
				 _heliostatLength        ,
				 _heliostatWidth         ,
				 _towerHeight            ,
				 _receiverApertureHeight ,
				 _receiverApertureWidth  ,
				 _minimumDistanceToTower ,
				 _maximumDistanceToTower ,
				 _fieldAngularWidth      ,
				 sun                       );

    // constructing Htf Cycle with desired steam generator model:
    htfCycle = new HtfCycle ( _centralReceiverOutletTemperature ,
			      _hotStorageHeight                 ,
			      _hotStorageDiameter               ,
			      _hotStorageInsulThickness         ,
			      _coldMoltenSaltMinTemperature     ,
			      powerblock                        ,
			      _receiverApertureHeight           ,
			      _receiverApertureWidth            ,
			      _receiverTubesInsideDiam          ,
			      _receiverTubesOutsideDiam         ,
			      _receiverNbOfTubes                ,
			      _minutesPerTimeIncrement            );
    htfCycle->initiateColdStorage();

    // investment cost model:
    economics = new Economics;
    economics->set_heightOfTower                  ( _towerHeight                     );
    economics->set_heightOfReceiverAperture       ( _receiverApertureHeight          );
    economics->set_widthOfReceiverAperture        ( _receiverApertureWidth           );
    economics->set_receiverNumberOfTubes          ( _receiverNbOfTubes               );
    economics->set_receiverTubesDout              ( _receiverTubesOutsideDiam        );
    economics->set_lengthOfHeliostats             ( _heliostatLength                 );
    economics->set_widthOfHeliostats              ( _heliostatWidth                  );
    economics->set_hotStorageHeight               ( _hotStorageHeight                );
    economics->set_storageDiameter                ( _hotStorageDiameter              );
    economics->set_hotStorageInsulationThickness  ( _hotStorageInsulThickness        );
    economics->set_coldStorageInsulationThickness ( _coldStorageInsulThickness       );
    economics->set_receiverInsulationThickness    ( _receiverInsulThickness          );
    economics->set_turbineNominalPowerOutput      ( powerblock->get_powerOfTurbine() );
    economics->set_exchangerModel                 ( _exchangerModel                  );
  }
  catch ( const std::exception & e ) {
    if ( field      ) delete field;
    if ( htfCycle   ) delete htfCycle;
    if ( powerblock ) delete powerblock;
    if ( economics  ) delete economics;
    cnt_eval = false;
    throw Simulation_Interruption ( e.what() );
  }
    
  _powerplant = new Powerplant ( time      ,
				 _model_type ,
				 field       ,
				 htfCycle    ,
				 powerblock  ,
				 economics     );
  
  // setting demand vector:
  _powerplant->set_demand         ( _demand               );
  _powerplant->set_heliostatModel ( _heliostatsFieldModel );
}

/*-----------------------------------------------*/
/*      constructing model components (8/10)     */
/*-----------------------------------------------*/
void Scenario::construct_maxHF_minCost ( bool & cnt_eval ) {

  if ( _powerplant ) {
    delete _powerplant;
    _powerplant = NULL;
  }
  
  _model_type = 2; // whole plant

  Time_Manager time ( _numberOfTimeIncrements, 0, _minutesPerTimeIncrement);
  Sun          sun  ( _latitude, time, _day, _raysPerSquareMeters);
  fFillDemandVector();

  HeliostatField * field      = NULL;
  HtfCycle       * htfCycle   = NULL;
  Powerblock     * powerblock = NULL;
  Economics      * economics  = NULL;

  try {
    
    // powerblock model:
    powerblock = new Powerblock ( _typeOfTurbine );
    
    // constructing heliostats field:
    field = new HeliostatField ( _numberOfHeliostats     ,
				 _heliostatLength        ,
				 _heliostatWidth         ,
				 _towerHeight            ,
				 _receiverApertureHeight ,
				 _receiverApertureWidth  ,
				 _minimumDistanceToTower ,
				 _maximumDistanceToTower ,
				 _fieldAngularWidth      ,
				 sun                       );
    
    // constructing Htf Cycle with desired steam generator model:
    htfCycle = new HtfCycle ( _centralReceiverOutletTemperature ,
			      _hotStorageHeight                 ,
			      _hotStorageDiameter               ,
			      _hotStorageInsulThickness         ,
			      _coldMoltenSaltMinTemperature     ,
			      powerblock                        ,
			      _receiverApertureHeight           ,
			      _receiverApertureWidth            ,
			      _receiverTubesInsideDiam          ,
			      _receiverTubesOutsideDiam         ,
			      _receiverNbOfTubes                ,
			      _minutesPerTimeIncrement 	        );
    htfCycle->initiateColdStorage();
    
    // investment cost model:
    economics = new Economics;
    economics->set_heightOfTower                  ( _towerHeight                     );
    economics->set_heightOfReceiverAperture       ( _receiverApertureHeight          );
    economics->set_widthOfReceiverAperture        ( _receiverApertureWidth           );
    economics->set_receiverNumberOfTubes          ( _receiverNbOfTubes               );
    economics->set_receiverTubesDout              ( _receiverTubesOutsideDiam        );
    economics->set_lengthOfHeliostats             ( _heliostatLength                 );
    economics->set_widthOfHeliostats              ( _heliostatWidth                  );
    economics->set_hotStorageHeight               ( _hotStorageHeight                );
    economics->set_storageDiameter                ( _hotStorageDiameter              );
    economics->set_hotStorageInsulationThickness  ( _hotStorageInsulThickness        );
    economics->set_coldStorageInsulationThickness ( _coldStorageInsulThickness       );
    economics->set_receiverInsulationThickness    ( _receiverInsulThickness          );
    economics->set_turbineNominalPowerOutput      ( powerblock->get_powerOfTurbine() );
    economics->set_exchangerModel                 ( _exchangerModel                  );
  }
  catch ( const std::exception & e ) {
    if ( field      ) delete field;
    if ( htfCycle   ) delete htfCycle;
    if ( powerblock ) delete powerblock;
    if ( economics  ) delete economics;
    cnt_eval = false;
    throw Simulation_Interruption ( e.what() );
  }

  _powerplant = new Powerplant ( time        ,
				 _model_type ,
				 field       ,
				 htfCycle    ,
				 powerblock  ,
				 economics     );
  
  // setting demand vector:
  _powerplant->set_demand         ( _demand               );
  _powerplant->set_heliostatModel ( _heliostatsFieldModel );
}

/*-----------------------------------------------*/
/*      constructing model components (9/10)     */
/*-----------------------------------------------*/
void Scenario::construct_maxNrg_minPar ( bool & cnt_eval ) {

  if ( _powerplant ) {
    delete _powerplant;
    _powerplant = NULL;
  }
  
  _model_type = 2; // whole plant

  Time_Manager time ( _numberOfTimeIncrements, 0, _minutesPerTimeIncrement);
  Sun          sun  ( _latitude, time, _day, _raysPerSquareMeters);
  fFillDemandVector();

  HeliostatField * field      = NULL;
  HtfCycle       * htfCycle   = NULL;
  Powerblock     * powerblock = NULL;
  Economics      * economics  = NULL;

  try {

    // powerblock model:
    powerblock = new Powerblock ( _typeOfTurbine );

    // constructing heliostats field:
    field = new HeliostatField ( _numberOfHeliostats     ,
				 _heliostatLength        ,
				 _heliostatWidth         ,
				 _towerHeight            ,
				 _receiverApertureHeight ,
				 _receiverApertureWidth  ,
				 _minimumDistanceToTower ,
				 _maximumDistanceToTower ,
				 _fieldAngularWidth      ,
				 sun                       );
    
    // constructing Htf Cycle with desired steam generator model:
    htfCycle = new HtfCycle ( _centralReceiverOutletTemperature ,
			      _hotStorageHeight                 ,
			      _hotStorageDiameter               ,
			      _hotStorageInsulThickness         ,
			      _coldMoltenSaltMinTemperature     ,
			      powerblock                        ,
			      _receiverApertureHeight           ,
			      _receiverApertureWidth            ,
			      _receiverTubesInsideDiam          ,
			      _receiverTubesOutsideDiam         ,
			      _receiverNbOfTubes                ,
			      _minutesPerTimeIncrement          ,
			      _exchangerTubesDin                ,
			      _exchangerTubesDout               ,
			      _exchangerTubesLength             ,
			      _exchangerTubesSpacing            ,
			      _exchangerBaffleCut               ,
			      _exchangerNbOfBaffles             ,
			      _exchangerNbOfTubes               ,
			      _exchangerNbOfPassesPerShell      ,
			      _exchangerNbOfShells                );

    htfCycle->initiateColdStorage();
    
    // investment cost model:
    economics = new Economics;
    economics->set_heightOfTower                  ( _towerHeight                     );
    economics->set_heightOfReceiverAperture       ( _receiverApertureHeight          );
    economics->set_widthOfReceiverAperture        ( _receiverApertureWidth           );
    economics->set_receiverNumberOfTubes          ( _receiverNbOfTubes               );
    economics->set_receiverTubesDout              ( _receiverTubesOutsideDiam        );
    economics->set_lengthOfHeliostats             ( _heliostatLength                 );
    economics->set_widthOfHeliostats              ( _heliostatWidth                  );
    economics->set_hotStorageHeight               ( _hotStorageHeight                );
    economics->set_storageDiameter                ( _hotStorageDiameter              );
    economics->set_hotStorageInsulationThickness  ( _hotStorageInsulThickness        );
    economics->set_coldStorageInsulationThickness ( _coldStorageInsulThickness       );
    economics->set_receiverInsulationThickness    ( _receiverInsulThickness          );
    economics->set_turbineNominalPowerOutput      ( powerblock->get_powerOfTurbine() );
    economics->set_exchangerModel                 ( _exchangerModel                  );
  }
  catch ( const std::exception & e ) {
    if ( field      ) delete field;
    if ( htfCycle   ) delete htfCycle;
    if ( powerblock ) delete powerblock;
    if ( economics  ) delete economics;
    cnt_eval = false;
    throw Simulation_Interruption ( e.what() );
  }
   
  _powerplant = new Powerplant ( time        ,
				 _model_type ,
				 field       ,
				 htfCycle    ,
				 powerblock  ,
				 economics     );
    
  // setting demand vector:
  _powerplant->set_demand         ( _demand               );
  _powerplant->set_heliostatModel ( _heliostatsFieldModel );
}

/*-----------------------------------------------*/
/*      constructing model components (10/10)    */
/*-----------------------------------------------*/
void Scenario::construct_minCost_unconstrained ( bool & cnt_eval ) {

  if ( _powerplant ) {
    delete _powerplant;
    _powerplant = NULL;
  }
  
  _model_type = 2; // whole plant

  Time_Manager time ( _numberOfTimeIncrements, 0, _minutesPerTimeIncrement);
  Sun          sun  ( _latitude, time, _day, _raysPerSquareMeters);
  fFillDemandVector();

  HtfCycle       * htfCycle   = NULL;
  Powerblock     * powerblock = NULL;
  Economics      * economics  = NULL;

  try {
    
    // powerblock model:
    powerblock = new Powerblock ( _typeOfTurbine );

    // constructing Htf Cycle with desired steam generator model:
    htfCycle = new HtfCycle ( _centralReceiverOutletTemperature ,
			      _hotStorageHeight                 ,
			      _hotStorageDiameter               ,
			      _hotStorageInsulThickness         ,
			      _coldMoltenSaltMinTemperature     ,
			      powerblock                        ,
			      _receiverApertureHeight           ,
			      _receiverApertureWidth            ,
			      _receiverTubesInsideDiam          ,
			      _receiverTubesOutsideDiam         ,
			      _receiverNbOfTubes                ,
			      _minutesPerTimeIncrement            );
  
    htfCycle->setStorage ( _storageStartupCondition               ,
			   0.98*_centralReceiverOutletTemperature ,
			   _coldMoltenSaltMinTemperature            );

    // investment cost model:
    economics = new Economics;
    economics->set_heightOfTower                  ( _towerHeight                     );
    economics->set_heightOfReceiverAperture       ( _receiverApertureHeight          );
    economics->set_widthOfReceiverAperture        ( _receiverApertureWidth           );
    economics->set_receiverNumberOfTubes          ( _receiverNbOfTubes               );
    economics->set_receiverTubesDout              ( _receiverTubesOutsideDiam        );
    economics->set_lengthOfHeliostats             ( _heliostatLength                 );
    economics->set_widthOfHeliostats              ( _heliostatWidth                  );
    economics->set_hotStorageHeight               ( _hotStorageHeight                );
    economics->set_storageDiameter                ( _hotStorageDiameter              );
    economics->set_hotStorageInsulationThickness  ( _hotStorageInsulThickness        );
    economics->set_coldStorageInsulationThickness ( _coldStorageInsulThickness       );
    economics->set_receiverInsulationThickness    ( _receiverInsulThickness          );
    economics->set_turbineNominalPowerOutput      ( powerblock->get_powerOfTurbine() );
    economics->set_exchangerModel                 ( _exchangerModel                  );
  }
  catch ( const std::exception & e ) {
    if ( htfCycle   ) delete htfCycle;
    if ( powerblock ) delete powerblock;
    if ( economics  ) delete economics;
    cnt_eval = false;
    throw Simulation_Interruption ( e.what() );
  }
    
  _powerplant = new Powerplant ( time        ,
				 _model_type ,
				 NULL        ,
				 htfCycle    ,
				 powerblock  ,
				 economics     );

  // setting demand vector:
  _powerplant->set_demand         ( _demand               );
  _powerplant->set_heliostatModel ( _heliostatsFieldModel );

  if  ( !_powerplant->set_heliostatFieldPowerOutput_MINCOST_TS() ) {
    delete _powerplant;
    _powerplant = NULL;
    cnt_eval    = false;
    throw Simulation_Interruption ( "error in the construction of the problem" );
  }
}

/*-------------------------------------------*/
void Scenario::fFillDemandVector ( void ) {
/*-------------------------------------------*/
  std::vector<int>    timeVector;
  std::vector<int>    basicTime;
  std::vector<double> demandVector;
  timeVector.reserve   ( _numberOfTimeIncrements );
  basicTime.reserve    ( _numberOfTimeIncrements );
  demandVector.reserve ( _numberOfTimeIncrements );

  int time = 0;

  for ( int i = 0; i < 24; ++i )
    basicTime.push_back ( i * 60 );

  for ( int i = 0; i < _numberOfTimeIncrements; ++i )
    timeVector.push_back ( i*_minutesPerTimeIncrement );

  if ( _demandProfile == 1 ) {
    for ( int i = 0; i < _numberOfTimeIncrements; ++i ) {
      time = timeVector[i] % 1440;
      if ( time < _tStart || time > _tEnd )
	_demand.push_back(0.0);
      else
	_demand.push_back ( _maximumPowerDemand );
    }
  }

  else if ( _demandProfile == 2 ) {
    for ( int i = 0; i < 24; ++i )
      demandVector.push_back(demandProfile_W[i]);
    for ( int i = 0; i < _numberOfTimeIncrements; ++i )
      _demand.push_back(kernelSmoothing(basicTime,demandVector,timeVector[i] % 1440) * _maximumPowerDemand);
  }
  
  if ( _demandProfile == 3 ) {
    for ( int i = 0; i < 24; ++i )
      demandVector.push_back(demandProfile_S[i]);
    for (int i = 0; i < _numberOfTimeIncrements; ++i)
      _demand.push_back(kernelSmoothing(basicTime, demandVector, timeVector[i]%1440) * _maximumPowerDemand);
  }
}

/*-----------------------------------------------------*/
/*  set the type of turbine and associated parameters  */
/*-----------------------------------------------------*/
bool Scenario::set_typeOfTurbine ( int tot ) {

  _minReceiverOutletTemp = 0.0;
  _typeOfTurbine         = 0;
  
  if ( tot < 1 || tot > 8 )
    return false;

  _typeOfTurbine = tot;
  
  switch ( _typeOfTurbine ) {
  case 1:
    _minReceiverOutletTemp = SST110_TEMPERATURE;
    break;
  case 2:
    _minReceiverOutletTemp = SST120_TEMPERATURE;
    break;
  case 3:
    _minReceiverOutletTemp = SST300_TEMPERATURE;
    break;
  case 4:
    _minReceiverOutletTemp = SST400_TEMPERATURE;
    break;
  case 5:
    _minReceiverOutletTemp = SST600_TEMPERATURE;
    break;
  case 6:
    _minReceiverOutletTemp = SST700_TEMPERATURE;
    break;
  case 7:
    _minReceiverOutletTemp = SST800_TEMPERATURE;
    break;
  case 8:
    _minReceiverOutletTemp = SST900_TEMPERATURE;
    break;
  }

  return true;
}
