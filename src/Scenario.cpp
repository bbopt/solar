#include "Scenario.hpp"

/*-------------------------------------*/
/*              constructor            */
/*-------------------------------------*/
Scenario::Scenario ( const std::string & problem , const std::string & x_file_name ) :
  _n                       ( 0       ) ,
  _x                       ( NULL    ) ,
  _problem                 ( problem ) ,
  _model_type              ( 0       ) ,
  _heliostatsFieldModel    ( 0       ) ,
  _exchangerModel          ( 0       ) ,
  _storageStartupCondition ( 0       )   {   // TOTO continuer les init comme cela


  // initialiser toutes les variables ici: TOTO
  // TOTO: Initialiser ces variables inutiles:  Dans le constructeur
  
  // // Htf cycle:
  // _centralReceiverOutletTemperature = _x[ 9];
  // _hotStorageHeight                 = _x[10];
  // _hotStorageDiameter               = _x[11];
  // _hotStorageInsulThickness         = _x[12];
  // _coldStorageInsulThickness        = _x[13];
  // _coldMoltenSaltMinTemperature     = _x[14];
  // _receiverNbOfTubes                = _x[15];
  // _receiverInsulThickness           = _x[16];
  // _receiverTubesInsideDiam          = _x[17];
  // _receiverTubesOutsideDiam         = _x[18];

  // // Powerblock:
  // _typeOfTurbine                    = _x[19];

  //    _cBudget = 0;
  
  // Problem #1:
  if ( _problem == "MAXNRG_H1"  || _problem == "MAXNRG_H1S" )
    init_maxNrg_H1();

  // Problem #2:
  if ( _problem == "MINSURF_H1" || _problem == " MINSURF_H1S" )
    init_minSurf_H1();
  
  // Problem #3:
  if ( _problem == "MINCOST_C1" || _problem == "MINCOST_C1S" )
    init_minCost_C1();

  // Problem #4:
  if ( _problem == "MINCOST_C2" || _problem == "MINCOST_C2S" )
    init_minCost_C2();
  
  // Problem #5:
  if ( _problem == "MAXCOMP_HTF1" || _problem == "MAXCOMP_HTF1S" )
    init_maxComp_HTF1();

  // Problem #6:
  if ( _problem == "MINCOST_TS" || _problem == "MINCOST_TSS" )
    init_minCost_TS();

  // Problem #7:
  if ( _problem == "MAXEFF_RE" || _problem == "MAXEFF_RES" )
    init_maxEff_RE();
  
  // Problem #8:
  if ( _problem == "MAXHF_MINCOST" || _problem == "MAXHF_MINCOSTS" )
    init_maxHF_minCost();

  // Problem #9:
  if ( _problem == "MAXNRG_MINPAR" || _problem == "MAXNRG_MINPARS" )
    init_maxNrg_minPar();
     
  
  // Check problem id:
  // -----------------
  if ( problem != "MAXNRG_H1"     && problem != "MAXNRG_H1S"     &&      // #1
       problem != "MINSURF_H1"    && problem != "MINSURF_H1S"    &&      // #2
       problem != "MINCOST_C1"    && problem != "MINCOST_C1S"    &&      // #3
       problem != "MINCOST_C2"    && problem != "MINCOST_C2S"    &&      // #4
       problem != "MAXCOMP_HTF1"  && problem != "MAXCOMP_HTF1S"  &&      // #5
       problem != "MINCOST_TS"    && problem != "MINCOST_TSS"    &&      // #6
       problem != "MAXEFF_RE"     && problem != "MAXEFF_RES"     &&      // #7
       problem != "MAXHF_MINCOST" && problem != "MAXHF_MINCOSTS" &&      // #8
       problem != "MAXNRG_MINPAR" && problem != "MAXNRG_MINPARS"    ) {  // #9
    delete_x();
    throw invalid_argument ( problem + " is not a valid problem ID" );
  }

  // Read x file:
  // ------------
  if ( !read(x_file_name) ) {
    delete_x();
    throw invalid_argument ( "Problem with input file " + x_file_name );
  }
}

/*-----------------------------------------*/
/*            general read function        */
/*-----------------------------------------*/
bool Scenario::read ( const std::string & x_file_name ) {
  if ( _problem == "MAXNRG_H1"     || _problem == "MAXNRG_H1S"     ) { return read_maxNrg_H1     ( x_file_name ); } // #1
  if ( _problem == "MINSURF_H1"    || _problem == "MINSURF_H1S"    ) { return read_minSurf_H1    ( x_file_name ); } // #2
  if ( _problem == "MINCOST_C1"    || _problem == "MINCOST_C1S"    ) { return read_minCost_C1    ( x_file_name ); } // #3
  if ( _problem == "MINCOST_C2"    || _problem == "MINCOST_C2S"    ) { return read_minCost_C2    ( x_file_name ); } // #4
  if ( _problem == "MAXCOMP_HTF1"  || _problem == "MAXCOMP_HTF1S"  ) { return read_maxComp_HTF1  ( x_file_name ); } // #5
  if ( _problem == "MINCOST_TS"    || _problem == "MINCOST_TSS"    ) { return read_minCost_TS    ( x_file_name ); } // #6
  if ( _problem == "MAXEFF_RE"     || _problem == "MAXEFF_RES"     ) { return read_maxEff_RE     ( x_file_name ); } // #7
  if ( _problem == "MAXHF_MINCOST" || _problem == "MAXHF_MINCOSTS" ) { return read_maxHF_minCost ( x_file_name ); } // #8
  if ( _problem == "MAXNRG_MINPAR" || _problem == "MAXNRG_MINPARS" ) { return read_maxNrg_minPar ( x_file_name ); } // #9
  return false;
}


// TOTO: la boite noire doit checker que les variables sont dans les bornes et que les entiers sont bien entiers.


/*-----------------------------------------*/
/*    initialize problem maxNrg_H1 (#1)    */
/*-----------------------------------------*/
void Scenario::init_maxNrg_H1 ( void ) {

  // Scenario parameters:
  _model_type           = 1; // heliostats field
  _heliostatsFieldModel = 1;
  _exchangerModel       = 1;
  
  _latitude = 44.95;
  _day = 100; // April 10th
  _numberOfTimeIncrements = 24;
  _minutesPerTimeIncrement = 60;
  _raysPerSquareMeters     = ( _problem == "MAXNRG_H1S" ) ? 0.0005 : 0.01;
  _cFieldSurface = 1.95e6; // 195 hectares
  _cBudget = 50e6;
  
  // vector of variables:
  delete_x();
  _n = 9;
  _x = new double[9];
}

/*------------------------------------------*/
/*  read inputs for problem maxNRG_H1 (#1)  */
/*------------------------------------------*/
bool Scenario::read_maxNrg_H1 ( const std::string & x_file_name ) {

  // read file and fill _x:
  if ( !read_x ( x_file_name ) ) {
    delete_x();
    return false;
  }

  // check discrete variables:
  // -------------------------
  if ( !Scenario::is_int(_x[5]) ) {
    throw invalid_argument ( "Problem with input file \"" + x_file_name + "\": One of the discrete variables has a non-integer value" );
    return false;
  }
  
  // assign variables:
  // -----------------

  // Heliostats Field:
  _heliostatHeight        = _x[0];
  _heliostatWidth         = _x[1];
  _towerHeight            = _x[2];
  _receiverApertureHeight = _x[3];
  _receiverApertureWidth  = _x[4];
  _maxNumberOfHeliostats  = Scenario::round(_x[5]);
  _fieldMaxAngularSpan    = _x[6];
  _minimumDistanceToTower = _x[7];
  _maximumDistanceToTower = _x[8];
 
  return true;
}

/*-----------------------------------------*/
/*    initialize problem minSurf_H1 (#2)   */
/*-----------------------------------------*/
void Scenario::init_minSurf_H1 ( void ) {

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
  _numberOfTimeIncrements  = ( _problem == "MINSURF_H1S" ) ? 24 : 72; 
  _raysPerSquareMeters     = ( _problem == "MINSURF_H1S" ) ? 0.01 : 0.01; // TOTO: C'etaient les memes valeurs... changer pour le surrogate ?
  _cFieldSurface           = 4e6; // 400 hectares
  _cDemandComplianceRatio  = 100;
  _cBudget                 = 300e6; // 300 millions
  _cParasitics             = 0.18;

  // design parameters:
  _hotStorageDiameter = 23;
  _hotStorageHeight = 10.5;
  _hotStorageInsulThickness = 0.3;
  _coldStorageInsulThickness = 0.2;
  _coldMoltenSaltMinTemperature = 282 + 273;
  _typeOfTurbine = 3;

  // vector of variables:
  delete_x();
  _n = 14;
  _x = new double[14];
}

/*-------------------------------------------*/
/*  read inputs for problem minSurf_H1 (#2)  */
/*-------------------------------------------*/
bool Scenario::read_minSurf_H1 ( const std::string & x_file_name ) {

  // read file and fill _x:
  if ( !read_x ( x_file_name ) ) {
    delete_x();
    return false;
  }

  // check discrete variables:
  // -------------------------
  if ( !Scenario::is_int(_x[5]) || !Scenario::is_int(_x[10]) ) {
    throw invalid_argument ( "Problem with input file \"" + x_file_name + "\": One of the discrete variables has a non-integer value" );
    return false;
  }
  
  // assign variables:
  // -----------------

  // Heliostats Field:
  _heliostatHeight        = _x[0];
  _heliostatWidth         = _x[1];
  _towerHeight            = _x[2];
  _receiverApertureHeight = _x[3];
  _receiverApertureWidth  = _x[4];
  _maxNumberOfHeliostats  = Scenario::round(_x[5]);
  _fieldMaxAngularSpan    = _x[6];
  _minimumDistanceToTower = _x[7];
  _maximumDistanceToTower = _x[8];

  // Htf cycle:
  _centralReceiverOutletTemperature = _x[ 9];
  _receiverNbOfTubes                = Scenario::round(_x[10]);
  _receiverInsulThickness           = _x[11];
  _receiverTubesInsideDiam          = _x[12];
  _receiverTubesOutsideDiam         = _x[13];

  return true;
}
  
/*-----------------------------------------*/
/*    initialize problem minCost_C1 (#3)   */
/*-----------------------------------------*/
void Scenario::init_minCost_C1 ( void ) {

  // Demand profile: 1, from 3 pm to 9 pm
  // Maximum demand: 10 MW
  // Latitude      : 35 deg
  // Day           : 1
  // Duration      : 48 hours
  // maximum field surface of 800 000 m^2

  // parameters:
  // -----------

  // Scenario parameters:
  _model_type           = 2; // whole plant
  _heliostatsFieldModel = 1;
  _exchangerModel       = 1;
  
  _latitude                = 35.;
  _day                     = 270;   // https://www.epochconverter.com/days/2019: Sept. 27th
  _demandProfile           = 1;
  _tStart                  =  900; // 15*60
  _tEnd                    = 1260; // 21*60;
  _maximumPowerDemand      = 1e7;
  _numberOfTimeIncrements  = 24;
  _storageStartupCondition = 0;
  _minutesPerTimeIncrement = 60;
  _fixedPointsPrecision    = 0.001;
  _cFieldSurface           = 800000;
  _cDemandComplianceRatio  = 100;
  _raysPerSquareMeters     = ( _problem == "MINCOST_C1S" ) ? 0.0005 : 0.01;
  
  // vector of variables:
  delete_x();
  _n = 20;
  _x = new double[20];
}
  
/*------------------------------------------*/
/*  read inputs for problem minCost_C1 (#3) */
/*------------------------------------------*/
bool Scenario::read_minCost_C1 ( const std::string & x_file_name ) {

  // read file and fill _x:
  if ( !read_x ( x_file_name ) ) {
    delete_x();
    return false;
  }

  // check discrete variables:
  // -------------------------
  if ( !Scenario::is_int(_x[5]) || !Scenario::is_int(_x[15]) || !Scenario::is_int(_x[19]) ) {
    throw invalid_argument ( "Problem with input file \"" + x_file_name + "\": One of the discrete variables has a non-integer value" );
    return false;
  }
  
  // assign variables:
  // -----------------

  // Heliostats Field:
  _heliostatHeight                  = _x[ 0];
  _heliostatWidth                   = _x[ 1];
  _towerHeight                      = _x[ 2];
  _receiverApertureHeight           = _x[ 3];
  _receiverApertureWidth            = _x[ 4];
  _maxNumberOfHeliostats            = Scenario::round(_x[5]);
  _fieldMaxAngularSpan              = _x[ 6];
  _minimumDistanceToTower           = _x[ 7];
  _maximumDistanceToTower           = _x[ 8];

  // Htf cycle:
  _centralReceiverOutletTemperature = _x[ 9];
  _hotStorageHeight                 = _x[10];
  _hotStorageDiameter               = _x[11];
  _hotStorageInsulThickness         = _x[12];
  _coldStorageInsulThickness        = _x[13];
  _coldMoltenSaltMinTemperature     = _x[14];
  _receiverNbOfTubes                = Scenario::round(_x[15]);
  _receiverInsulThickness           = _x[16];
  _receiverTubesInsideDiam          = _x[17];
  _receiverTubesOutsideDiam         = _x[18];

  // Powerblock:
  _typeOfTurbine                    = Scenario::round(_x[19]);
 
  return true;
}

/*------------------------------------------*/
/*    initialize problem minCost_C2 (#4)    */
/*------------------------------------------*/
void Scenario::init_minCost_C2 ( void ) {

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

  
  _latitude = 35.;
  _day = 1;
  _demandProfile = 1;
  _tStart = 15 * 60;
  _tEnd = 21 * 60;
  _maximumPowerDemand = 25e6;
  _storageStartupCondition = 50;
  _numberOfTimeIncrements = 24;
  _minutesPerTimeIncrement = 60;
  _raysPerSquareMeters = 0.01;
  if (_problem == "MINCOST_C2S"){ _raysPerSquareMeters = 0.0005; }
  _fixedPointsPrecision = 0.001;
  _cFieldSurface = 2000000.; // 200 hectares
  _cDemandComplianceRatio = 100;
  _cParasitics = 0.18;
  
  // vector of variables:
  delete_x();
  _n = 29;
  _x = new double[29];
}

/*-------------------------------------------*/
/*  read inputs for problem minCost_C2 (#4)  */
/*-------------------------------------------*/
bool Scenario::read_minCost_C2 ( const std::string & x_file_name ) {

  // read file and fill _x:
  if ( !read_x ( x_file_name ) ) {
    delete_x();
    return false;
  }

  // check discrete variables:
  // -------------------------
  if ( !Scenario::is_int(_x[ 5]) ||
       !Scenario::is_int(_x[15]) ||
       !Scenario::is_int(_x[24]) ||
       !Scenario::is_int(_x[25]) ||
       !Scenario::is_int(_x[26]) ||
       !Scenario::is_int(_x[27]) ||
       !Scenario::is_int(_x[28])    ) {
    throw invalid_argument ( "Problem with input file \"" + x_file_name + "\": One of the discrete variables has a non-integer value" );
    return false;
  }
  
  // assign variables:
  // -----------------

  // Heliostats Field:
  _heliostatHeight        = _x[0];
  _heliostatWidth         = _x[1];
  _towerHeight            = _x[2];
  _receiverApertureHeight = _x[3];
  _receiverApertureWidth  = _x[4];
  _maxNumberOfHeliostats  = Scenario::round(_x[5]);
  _fieldMaxAngularSpan    = _x[6];
  _minimumDistanceToTower = _x[7];
  _maximumDistanceToTower = _x[8];

  // Htf cycle:
  _centralReceiverOutletTemperature = _x[ 9];
  _hotStorageHeight                 = _x[10];
  _hotStorageDiameter               = _x[11];
  _hotStorageInsulThickness         = _x[12];
  _coldStorageInsulThickness        = _x[13];
  _coldMoltenSaltMinTemperature     = _x[14];
  _receiverNbOfTubes                = Scenario::round(_x[15]);
  _receiverInsulThickness           = _x[16];
  _receiverTubesInsideDiam          = _x[17];
  _receiverTubesOutsideDiam         = _x[18];
  _exchangerTubesSpacing            = _x[19];
  _exchangerTubesLength             = _x[20];
  _exchangerTubesDin                = _x[21];
  _exchangerTubesDout               = _x[22];
  _exchangerBaffleCut               = _x[23];
  _exchangerNbOfBaffles             = Scenario::round(_x[24]);
  _exchangerNbOfTubes               = Scenario::round(_x[25]);
  _exchangerNbOfShells              = Scenario::round(_x[26]);
  _exchangerNbOfPassesPerShell      = Scenario::round(_x[27]);

  // Powerblock:
  _typeOfTurbine = Scenario::round(_x[28]);

  return true;
}

/*--------------------------------------------*/
/*    initialize problem maxComp_HTF1X (#5)   */
/*--------------------------------------------*/
void Scenario::init_maxComp_HTF1 ( void ) {

  // Scenario parameters:
  _model_type = 2; // whole plant
  _heliostatsFieldModel = 2;
  _exchangerModel = 2;
  
  _latitude = 37.5581;
  _day = 30;
  _demandProfile = 1;
  _tStart = 0 * 60;
  _tEnd = 23 * 60;
  _maximumPowerDemand = 12e6;  // 12 MW

  _numberOfTimeIncrements = 720; // 24*30
  if (_problem == "MAXCOMP_HTF1S"){ _numberOfTimeIncrements = 72; }
  
  _minutesPerTimeIncrement = 60;
  _fixedPointsPrecision = 0.001;
  _raysPerSquareMeters = 0.01;
  _cBudget = 100e6;
  _cParasitics = 0.18;

  _heliostatHeight = 2.1336;
  _heliostatWidth = 3.048;
  _towerHeight = 100.;
  _receiverApertureHeight = 6.;
  _receiverApertureWidth = 6.;
  _maxNumberOfHeliostats = 3800;
  _fieldMaxAngularSpan = 89;
  _minimumDistanceToTower = 0.5;
  _maximumDistanceToTower = 10;

  // vector of variables:
  delete_x();
  _n = 20;
  _x = new double[20];
}

/*---------------------------------------------*/
/*  read inputs for problem maxComp_HTF1 (#5)  */
/*---------------------------------------------*/
bool Scenario::read_maxComp_HTF1 ( const std::string & x_file_name ) {

  // read file and fill _x:
  if ( !read_x ( x_file_name ) ) {
    delete_x();
    return false;
  }

  // check discrete variables:
  // -------------------------
  if ( !Scenario::is_int(_x[ 6]) ||
       !Scenario::is_int(_x[15]) ||
       !Scenario::is_int(_x[16]) ||
       !Scenario::is_int(_x[17]) ||
       !Scenario::is_int(_x[18]) ||
       !Scenario::is_int(_x[19])    ) {
    throw invalid_argument ( "Problem with input file \"" + x_file_name + "\": One of the discrete variables has a non-integer value" );
    return false;
  }
  
  // assign variables:
  // -----------------

  // htf cycle:
  _centralReceiverOutletTemperature = _x[ 0];
  _hotStorageHeight                 = _x[ 1];
  _hotStorageDiameter               = _x[ 2];
  _hotStorageInsulThickness         = _x[ 3];
  _coldStorageInsulThickness        = _x[ 4];
  _coldMoltenSaltMinTemperature     = _x[ 5];
  _receiverNbOfTubes                = Scenario::round(_x[ 6]);
  _receiverInsulThickness           = _x[ 7];
  _receiverTubesInsideDiam          = _x[ 8];
  _receiverTubesOutsideDiam         = _x[ 9];
  _exchangerTubesSpacing            = _x[10];
  _exchangerTubesLength             = _x[11];
  _exchangerTubesDin                = _x[12];
  _exchangerTubesDout               = _x[13];
  _exchangerBaffleCut               = _x[14];
  _exchangerNbOfBaffles             = Scenario::round(_x[15]);
  _exchangerNbOfTubes               = Scenario::round(_x[16]);
  _exchangerNbOfShells              = Scenario::round(_x[17]);
  _exchangerNbOfPassesPerShell      = Scenario::round(_x[18]);
	
  // powerblock:
  _typeOfTurbine = Scenario::round(_x[19]);

  return true;
}

/*-----------------------------------------*/
/*    initialize problem minCost_TS (#6)   */
/*-----------------------------------------*/
void Scenario::init_minCost_TS ( void ) {

  // Scenario parameters:
  _model_type = 2; // whole plant
  _heliostatsFieldModel = 2;
  _exchangerModel = 1;

  _latitude = 30.05;
  _day = 1;
  _demandProfile = 1;
  _tStart =    0;  // 0 * 60
  _tEnd   = 1380; // 23 * 60
  _maximumPowerDemand = 120e6;
  _storageStartupCondition = 50;
  _numberOfTimeIncrements = 24;
  _minutesPerTimeIncrement = 60;
  _fixedPointsPrecision = 0.001;
  _raysPerSquareMeters = 0.01;
  if (_problem == "MINCOST_TSS"){ 
    _raysPerSquareMeters = 0.0005;
  }
  _cDemandComplianceRatio = 100.;
  _coldMoltenSaltMinTemperature = 530;
  _typeOfTurbine = 7;
	
  _heliostatHeight = 9;
  _heliostatWidth = 9;
  _towerHeight = 250.;
  _receiverApertureHeight = 15;
  _receiverApertureWidth = 20;
  _maxNumberOfHeliostats = 12232;
  _fieldMaxAngularSpan = 65;
  _minimumDistanceToTower = 1;
  _maximumDistanceToTower = 10.5;
  _receiverInsulThickness = 0.5;
  _receiverNbOfTubes = 85;
  _receiverTubesInsideDiam = 0.033;
  _receiverTubesOutsideDiam = 0.050;

  // vector of variables:
  delete_x();
  _n = 5;
  _x = new double[5];
}

/*-------------------------------------------*/
/*  read inputs for problem minCost_TS (#6)  */
/*-------------------------------------------*/
bool Scenario::read_minCost_TS ( const std::string & x_file_name ) {

  // read file and fill _x:
  if ( !read_x ( x_file_name ) ) {
    delete_x();
    return false;
  }
  
  // assign variables:
  // -----------------

  // htf cycle:
  _centralReceiverOutletTemperature = _x[0];
  _hotStorageHeight                 = _x[1];
  _hotStorageDiameter               = _x[2];
  _hotStorageInsulThickness         = _x[3];
  _coldStorageInsulThickness        = _x[4];

  return true;
}


/*-----------------------------------------*/
/*    initialize problem maxEff_RE (#7)    */
/*-----------------------------------------*/
void Scenario::init_maxEff_RE ( void ) {

  // Scenario parameters:
  _model_type = 2; // whole plant
  _heliostatsFieldModel = 1;
  _exchangerModel = 1;

  _latitude = 30.05;
  _day = 1;
  _demandProfile = 1;
  _tStart =    0; // 0*60
  _tEnd   = 1380; // 23*60
  _maximumPowerDemand = 0;
  _numberOfTimeIncrements = 24;
  _minutesPerTimeIncrement = 60;
  _raysPerSquareMeters = 0.01;
  if (_problem == "MAXEFF_RES"){ _raysPerSquareMeters = 0.0005; }
  _fixedPointsPrecision = 0.001;
  _cParasitics = 0.03;
  _cBudget = 45e6;

  _typeOfTurbine = 1;
  _hotStorageDiameter = 25;
  _hotStorageHeight = 30;
  _hotStorageInsulThickness = 5;
  _coldStorageInsulThickness = 5;
  _coldMoltenSaltMinTemperature = 550;

  _heliostatHeight = 9;
  _heliostatWidth = 9;
  _towerHeight = 175.;
  _maxNumberOfHeliostats = 5000;
  _fieldMaxAngularSpan = 55;
  _minimumDistanceToTower = 1;
  _maximumDistanceToTower = 7;

  
  // vector of variables:
  delete_x();
  _n = 7;
  _x = new double[7];
}

/*-------------------------------------------*/
/*  read inputs for problem maxEff_RE (#7)   */
/*-------------------------------------------*/
bool Scenario::read_maxEff_RE ( const std::string & x_file_name ) {

  // read file and fill _x:
  if ( !read_x ( x_file_name ) ) {
    delete_x();
    return false;
  }

  // check discrete variables:
  // -------------------------
  if ( !Scenario::is_int(_x[3]) ) {
    throw invalid_argument ( "Problem with input file \"" + x_file_name + "\": One of the discrete variables has a non-integer value" );
    return false;
  }
  
  // assign variables:
  // -----------------

  // Receiver design parameters:
  _receiverApertureHeight           = _x[0];
  _receiverApertureWidth            = _x[1];
  _centralReceiverOutletTemperature = _x[2];
  _receiverNbOfTubes                = Scenario::round(_x[3]);
  _receiverInsulThickness           = _x[4];
  _receiverTubesInsideDiam          = _x[5];
  _receiverTubesOutsideDiam         = _x[6];

  return true;
}


/*--------------------------------------------*/
/*    initialize problem maxHF_minCost (#8)   */
/*--------------------------------------------*/
void Scenario::init_maxHF_minCost ( void ) {

  // Scenario parameters:
  _model_type = 2; // whole plant
  _heliostatsFieldModel = 1;
  _exchangerModel = 1;
  
  _latitude = 45.;
  _day = 1;
  _demandProfile = 1;
  _tStart =    0; //  0*60;
  _tEnd   = 1380; // 23*60;
  _maximumPowerDemand = 0;
  _numberOfTimeIncrements = 24;
  _minutesPerTimeIncrement = 60;
  _raysPerSquareMeters = 0.01;
  _storageStartupCondition = 0;
  if (_problem == "MAXHF_MINCOSTS"){ _raysPerSquareMeters = 0.0005; }
  _fixedPointsPrecision = 0.001;
  _cFieldSurface = 4e6;
  _cParasitics = 0.08;
  
  _typeOfTurbine = 1;
  _centralReceiverOutletTemperature = 950;
  _hotStorageDiameter = 25;
  _hotStorageHeight = 30;
  _hotStorageInsulThickness = 5;
  _coldStorageInsulThickness = 5;
  _coldMoltenSaltMinTemperature = 550;
    
  // vector of variables:
  delete_x();
  _n = 13;
  _x = new double[13];
}

/*----------------------------------------------*/
/*  read inputs for problem maxHF_minCost (#8)  */
/*----------------------------------------------*/
bool Scenario::read_maxHF_minCost ( const std::string & x_file_name ) {

  // read file and fill _x:
  if ( !read_x ( x_file_name ) ) {
    delete_x();
    return false;
  }

    // check discrete variables:
  // -------------------------
  if ( !Scenario::is_int(_x[5]) || !Scenario::is_int(_x[9]) ) {
    throw invalid_argument ( "Problem with input file \"" + x_file_name + "\": One of the discrete variables has a non-integer value" );
    return false;
  }
  
  // assign variables:
  // -----------------

  // Heliostats Field:
  _heliostatHeight        = _x[0];
  _heliostatWidth         = _x[1];
  _towerHeight            = _x[2];
  _receiverApertureHeight = _x[3];
  _receiverApertureWidth  = _x[4];
  _maxNumberOfHeliostats  = Scenario::round(_x[5]);
  _fieldMaxAngularSpan    = _x[6];
  _minimumDistanceToTower = _x[7];
  _maximumDistanceToTower = _x[8];
  	
  // Htf cycle:
  _receiverNbOfTubes        = Scenario::round(_x[ 9]);
  _receiverInsulThickness   = _x[10];
  _receiverTubesInsideDiam  = _x[11];
  _receiverTubesOutsideDiam = _x[12];

  return true;
}


/*--------------------------------------------*/
/*    initialize problem maxNrg_minPar (#9)   */
/*--------------------------------------------*/
void Scenario::init_maxNrg_minPar ( void ) {

  // Demand profile : 1, from 3 pm to 9 pm
  // Maximum demand : 250MW
  // Latitude : 25 deg
  // Day : 180
  // Duration : 24 hours
  // maximum field surface of 5M m^2

  // Scenario parameters:
 _model_type = 2; // whole plant
 _heliostatsFieldModel = 1;
 _exchangerModel = 2;

 _latitude = 25.;
 _day = 180;
 _demandProfile = 1;
 _tStart =    0; //  0*60;
 _tEnd   = 1380; // 23*60;
 _maximumPowerDemand = 250e6;
 _numberOfTimeIncrements = 24;
 _minutesPerTimeIncrement = 60;
 _raysPerSquareMeters = 0.01;
 if (_problem == "MAXNRG_MINPARS"){ _raysPerSquareMeters = 0.0005; }
 _fixedPointsPrecision = 0.001;
 _cFieldSurface = 5e6;
 _cParasitics = 0.2;
 _cBudget = 1.2e9;

  // vector of variables:
  delete_x();
  _n = 29;
  _x = new double[29];
}

/*----------------------------------------------*/
/*  read inputs for problem maxNrg_minPar (#9)  */
/*----------------------------------------------*/
bool Scenario::read_maxNrg_minPar ( const std::string & x_file_name ) {

  // read file and fill _x:
  if ( !read_x ( x_file_name ) ) {
    delete_x();
    return false;
  }

  // check discrete variables:
  // -------------------------
  if ( !Scenario::is_int(_x[ 5]) ||
       !Scenario::is_int(_x[15]) ||
       !Scenario::is_int(_x[24]) ||
       !Scenario::is_int(_x[25]) ||
       !Scenario::is_int(_x[26]) ||
       !Scenario::is_int(_x[27]) ||
       !Scenario::is_int(_x[28])    ) {
    throw invalid_argument ( "Problem with input file \"" + x_file_name + "\": One of the discrete variables has a non-integer value" );
    return false;
  }
  
  // assign variables:
  // -----------------

  // Heliostats Field:
  _heliostatHeight        = _x[0];
  _heliostatWidth         = _x[1];
  _towerHeight            = _x[2];
  _receiverApertureHeight = _x[3];
  _receiverApertureWidth  = _x[4];
  _maxNumberOfHeliostats  = Scenario::round(_x[5]);
  _fieldMaxAngularSpan    = _x[6];
  _minimumDistanceToTower = _x[7];
  _maximumDistanceToTower = _x[8];

  // Htf cycle:
  _centralReceiverOutletTemperature = _x[ 9];
  _hotStorageHeight                 = _x[10];
  _hotStorageDiameter               = _x[11];
  _hotStorageInsulThickness         = _x[12];
  _coldStorageInsulThickness        = _x[13];
  _coldMoltenSaltMinTemperature     = _x[14];
  _receiverNbOfTubes                = Scenario::round(_x[15]);
  _receiverInsulThickness           = _x[16];
  _receiverTubesInsideDiam          = _x[17];
  _receiverTubesOutsideDiam         = _x[18];
  _exchangerTubesSpacing            = _x[19];
  _exchangerTubesLength             = _x[20];
  _exchangerTubesDin                = _x[21];
  _exchangerTubesDout               = _x[22];
  _exchangerBaffleCut               = _x[23];
  _exchangerNbOfBaffles             = Scenario::round(_x[24]);
  _exchangerNbOfTubes               = Scenario::round(_x[25]);
  _exchangerNbOfShells              = Scenario::round(_x[26]);
  _exchangerNbOfPassesPerShell      = Scenario::round(_x[27]);

  // Powerblock:
  _typeOfTurbine = Scenario::round(_x[28]);

  return true;
}

/*----------------------------------------------*/
/*  functions to launch simulation and output   */
/*  the values according to the problem chosen  */
/*----------------------------------------------*/
bool Scenario::simulate ( std::ostream & out , bool & cnt_eval , bool verbose ) {

  bool simulation_completed = false;
  cnt_eval                  = false;
  
  // #1:
  if ( _problem == "MAXNRG_H1" || _problem == "MAXNRG_H1S" ) {
    if ( verbose ) {
      out << "BEGIN SIMULATE MAXNRG_H1 for x=( ";
      display_x(out);
      out <<  ") ..." << std::endl;
    }
    simulation_completed = simulate_maxNrg_H1 ( out , cnt_eval );
    if ( verbose )
      out << "... END SIMULATE" << std::endl;
    return simulation_completed;
  }

  // #2:
  if ( _problem == "MINSURF_H1" || _problem == "MINSURF_H1S" ) {
    if ( verbose ) {
      out << "BEGIN SIMULATE MINSURF_H1 for x=( ";
      display_x(out);
      out <<  ") ..." << std::endl;
    }
    simulation_completed = simulate_minSurf_H1 ( out , cnt_eval );
    if ( verbose )
      out << "... END SIMULATE" << std::endl;
    return simulation_completed;
  }

  // #3:
  if ( _problem == "MINCOST_C1" || _problem == "MINCOST_C1S" ) {  // TOTO: VERSION S ??
    if ( verbose ) {
      out << "BEGIN SIMULATE MINCOST_C1 for x=( ";
      display_x(out);
      out <<  ") ..." << std::endl;
    }
    simulation_completed = simulate_minCost_C1 ( out , cnt_eval );
    if ( verbose )
      out << "... END SIMULATE" << std::endl;
    return simulation_completed;
  }

  // #4:
  if ( _problem == "MINCOST_C2" || _problem == "MINCOST_C2S" ) {
    if ( verbose ) {
      out << "BEGIN SIMULATE MINCOST_C2 for x=( ";
      display_x(out);
      out <<  ") ..." << std::endl;
    }
    simulation_completed = simulate_minCost_C2 ( out , cnt_eval );
    if ( verbose )
      out << "... END SIMULATE" << std::endl;
    return simulation_completed;
  }

  // #5:
  if ( _problem == "MAXCOMP_HTF1" || _problem == "MAXCOMP_HTF1S" ) {
    if ( verbose ) {
      out << "BEGIN SIMULATE MAXCOMP_HTF1 for x=( ";
      display_x(out);
      out <<  ") ..." << std::endl;
    }
    simulation_completed = simulate_maxComp_HTF1 ( out , cnt_eval );
    if ( verbose )
      out << "... END SIMULATE" << std::endl;
    return simulation_completed;
  }

  // #6:
  if ( _problem == "MINCOST_TS" || _problem == "MINCOST_TSS" ) {
    if ( verbose ) {
      out << "BEGIN SIMULATE MINCOST_TS for x=( ";
      display_x(out);
      out <<  ") ..." << std::endl;
    }
    simulation_completed = simulate_minCost_TS ( out , cnt_eval );
    if ( verbose )
      out << "... END SIMULATE" << std::endl;
    return simulation_completed;
  }

  // #7:
  if ( _problem == "MAXEFF_RE" || _problem == "MAXEFF_RES" ) {
    if ( verbose ) {
      out << "BEGIN SIMULATE MAXEFF_RE for x=( ";
      display_x(out);
      out <<  ") ..." << std::endl;
    }
    simulation_completed = simulate_maxEff_RE ( out , cnt_eval );
    if ( verbose )
      out << "... END SIMULATE" << std::endl;
    return simulation_completed;
  }

  // #8:
  if ( _problem == "MAXHF_MINCOST" || _problem == "MAXHF_MINCOSTS" ) {
    if ( verbose ) {
      out << "BEGIN SIMULATE MAXHF_MINCOST for x=( ";
      display_x(out);
      out <<  ") ..." << std::endl;
    }
    simulation_completed = simulate_maxHF_minCost ( out , cnt_eval );
    if ( verbose )
      out << "... END SIMULATE" << std::endl;
    return simulation_completed;
  }

  // #9:
  if ( _problem == "MAXNRG_MINPAR" || _problem == "MAXNRG_MINPARS" ) {
    if ( verbose ) {
      out << "BEGIN SIMULATE MAXNRG_MINPAR for x=( ";
      display_x(out);
      out <<  ") ..." << std::endl;
    }
    simulation_completed = simulate_maxNrg_minPar ( out , cnt_eval );
    if ( verbose )
      out << "... END SIMULATE" << std::endl;
    return simulation_completed;
  }

  return false;
}

/*---------------------------*/
/*  simulate_maxNRG_H1 (#1)  */
/*---------------------------*/
bool Scenario::simulate_maxNrg_H1 ( std::ostream & out , bool & cnt_eval ) {

  cnt_eval = false;
  
  double f1, c1, c2, c3, c4, c5;
  f1 = c1 = c2 = c3 = c4 = c5 = 1e20;
  
  try {
		
    // Verifying bounds and a priori constraints:
    cnt_eval = validateInputValues();
    
    // Creating required objects:
    construct_maxNrg_H1();

    // Launching simulation:
    _powerplant->fSimulatePowerplant();
  
    // Objective function: total energy gathered in kWh:
    f1 = -_powerplant->get_totalEnergyConcentrated();  
   
    // Check if budget is respected:
    c1 = _powerplant->get_costOfHeliostatField()
      + _powerplant->get_costOfTower()
      + _powerplant->get_costOfReceiver()
      - _cBudget;

    // Check if total land area is below 1M m^2:
    c2 = M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
      *(2 * _fieldMaxAngularSpan / 360.) - _cFieldSurface;
    
    // Check basic geometric requirements:
    c3 = 2 * _heliostatHeight - _towerHeight;
    c4 = _minimumDistanceToTower - _maximumDistanceToTower;
    c5 = _maxNumberOfHeliostats - _powerplant->get_heliostatField()->get_listOfHeliostats().size();

    out << f1 << " " << c1 << " " << c2 << " " << c3 << " " << c4 << " " << c5 << std::endl;
  }
  catch ( ... ) {
    
    c2 = M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
      *(2 * _fieldMaxAngularSpan / 360.) - _cFieldSurface;
    
    c3 = 2 * _heliostatHeight - _towerHeight;
    c4 = _minimumDistanceToTower - _maximumDistanceToTower;
    
    out << f1 << " " << c1 << " " << c2 << " " << c3 << " " << c4 << " " << c5 << std::endl;
    
    throw logic_error ( "simulation could not go through" );
  }
  
  return true;
}

/*----------------------------*/
/*  simulate_minSurf_H1 (#2)  */
/*----------------------------*/
bool Scenario::simulate_minSurf_H1 ( std::ostream & out , bool & cnt_eval ) {

  cnt_eval = false;

  double f1 = 1e20, c1 = 1e20, c2 = 1e20,
    c3 = 1e20, c4 = 1e20, c5 = 1e20, c6 = 1e20, c7 = 1e20, 
    c8 = 1e20, c9 = 1e20, c10 = 1e20, c11 = 1e20, c12 = 1e20, 
    c13 = 1e20;

  try{

    // Verifying a priori constraints:
    cnt_eval = validateInputValues();

    // Creating required objects:
    construct_minSurf_H1();

    // Launching simulation:
    _powerplant->fSimulatePowerplant();

    // Objective function: field surface (m^2):
    f1 = _powerplant->get_fieldSurface();

    // Check if surface constraint is respected:
    c1 = _powerplant->get_fieldSurface() - _cFieldSurface;

    // Check if demand is met:
    c2 = _cDemandComplianceRatio - _powerplant->get_overallComplianceToDemand();

    // Check if total cost does not exceed limit:
    c3 = _powerplant->get_costOfHeliostatField()
      + _powerplant->get_costOfTower()
      + _powerplant->get_costOfReceiver()
      + _powerplant->get_costOfStorage()
      + _powerplant->get_costOfSteamGenerator()
      + _powerplant->get_costOfPowerblock()
      - _cBudget;

    // Check basic geometric requirements:
    c4 = 2 * _heliostatHeight - _towerHeight;
    c5 = _minimumDistanceToTower - _maximumDistanceToTower;
		
    // Check number of heliostats requested fits in the field:
    c6 = _maxNumberOfHeliostats - _powerplant->get_heliostatField()->get_listOfHeliostats().size();

    // Check pressure in tubes:
    c7 = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();

    // Check if molten salt temperature does not drop below the melting point:
    c8  = MELTING_POINT - _powerplant->get_minHotStorageTemp();
    c9  = MELTING_POINT - _powerplant->get_minColdStorageTemp();
    c10 = MELTING_POINT - _powerplant->get_minSteamGenTemp();

    c11 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;

    // Check if tubes fit in receiver:
    c12 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;
    
    c13 = (_powerplant->get_steamTurbineInletTemperature()) - _centralReceiverOutletTemperature;
    
    out << f1 << " " << c1 << " " << c2 << " "
	<< c3 << " " << c4 << " " << c5 << " " << c6 << " " << c7 << " "
	<< c8 << " " << c9 << " " << c10 << " " << c11 << " "
	<< c12 << " " << c13
	<< std::endl;
  }
  catch (...) {
    c1 = M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
      *(2 * _fieldMaxAngularSpan / 360.) - _cFieldSurface;

    c4 = 2 * _heliostatHeight - _towerHeight;
    c5 = _minimumDistanceToTower - _maximumDistanceToTower;

    c6 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    c11 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;
    c13 = _minReceiverOutletTemp - _centralReceiverOutletTemperature;

    out << f1 << " " << c1 << " " << c2 << " "
	<< c3 << " " << c4 << " " << c5 << " " << c6 << " " << c7 << " "
	<< c8 << " " << c9 << " " << c10 << " " << c11 << " "
	<< c12 << " " << c13
	<< std::endl;
    
    throw logic_error ( "Simulation could not go through" );
  }
  
  return true;
}

/*----------------------------*/
/*  simulate_minCost_C1 (#3)  */
/*----------------------------*/
bool Scenario::simulate_minCost_C1 ( std::ostream & out , bool & cnt_eval ) {

  cnt_eval = false;

  double f1 = 1e20, c1 = 1e20, c2 = 1e20,
    c3 = 1e20, c4 = 1e20, c5 = 1e20, c6 = 1e20,
    c7 = 1e20, c8 = 1e20, c9 = 1e20, c10 = 1e20, c11 = 1e20, 
    c12 = 1e20, c13 = 1e20;

  try{

    // verifying a priori constraints:
    cnt_eval = validateInputValues();

    // creating required objects:
    construct_minCost_C1();

    // launching simulation:
    _powerplant->fSimulatePowerplant();

    // objective function: total investment cost:
    f1 = _powerplant->get_costOfHeliostatField()
      + _powerplant->get_costOfTower()
      + _powerplant->get_costOfReceiver()
      + _powerplant->get_costOfStorage()
      + _powerplant->get_costOfSteamGenerator()
      + _powerplant->get_costOfPowerblock();

    // check total land area is below 700k m^2:
    c1 = M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
      *(2 * _fieldMaxAngularSpan / 360.) - _cFieldSurface;

    // check if compliance to demand is 100%:
    c2 = _cDemandComplianceRatio - _powerplant->get_overallComplianceToDemand();

    // check if tower at least twice as high as heliostats:
    c3 = 2 * _heliostatHeight - _towerHeight;

    // check Rmin < Rmax:
    c4 = _minimumDistanceToTower - _maximumDistanceToTower;

    // check number of heliostats fits in the field:
    c5 = _maxNumberOfHeliostats - _powerplant->get_heliostatField()->get_listOfHeliostats().size();

    // check maximum pressure in receiver tubes do not exceed yield pressure:
    c6 = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();

    // check molten Salt temperature does not drop below the melting point:
    c7 = MELTING_POINT - _powerplant->get_minHotStorageTemp();
    c8 = MELTING_POINT - _powerplant->get_minColdStorageTemp();
    c9 = MELTING_POINT - _powerplant->get_minSteamGenTemp();

    // check receiver tubes Din < Dout:
    c10 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;

    // check if tubes fit in receiver:
    c11 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;
		
    // check central receiver outlet is higher than that required by the turbine:
    c12 = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;

    // check if storage is back to initial conditions:
    double storageFinalConditions = 0.;
    storageFinalConditions = _powerplant->get_moltenSaltLoop()->get_hotStorage().get_heightOfVolumeStored()
      / _hotStorageHeight;
    c13 = 1.*_storageStartupCondition - storageFinalConditions;

    out << f1 << " " << c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6 << " " << c7 << " " << 
      c8 << " " << c9 << " " << c10 << " " << c11 << " " << c12 << " " << c13
	<< std::endl;
  }
  catch (...){
    c1 = M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
      *(2 * _fieldMaxAngularSpan / 360.) - _cFieldSurface;
    c3 = 2 * _heliostatHeight - _towerHeight;
    c4 = _minimumDistanceToTower - _maximumDistanceToTower;
    
    c10 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    c11 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;
    c12 =  _minReceiverOutletTemp - _centralReceiverOutletTemperature;
    
    out << f1 << " " << c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6 << " " << c7 << " " << 
      c8 << " " << c9 << " " << c10 << " " << c11 << " " << c12 << " " << c13
	<< std::endl;
    
    throw logic_error ( "Simulation could not go through" );
  }
  return true;
}

/*----------------------------*/
/*  simulate_minCost_C2 (#4)  */
/*----------------------------*/
bool Scenario::simulate_minCost_C2 ( std::ostream & out , bool & cnt_eval ) {

  cnt_eval = false;

  double f1 = 1e20, c1 = 1e20, c2 = 1e20,
    c3 = 1e20, c4 = 1e20, c5 = 1e20, c6 = 1e20,
    c7 = 1e20, c8 = 1e20, c9 = 1e20, c10 = 1e20, c11 = 1e20, c12 = 1e20,
    c13 = 1e20, c14 = 1e20, c15 = 1e20, c16 = 1e20;

  try {

    // Verifying a priori constraints:
    cnt_eval = validateInputValues();
  
    // Creating required objects:
    construct_minCost_C2();

    
    // Launching simulation
    _powerplant->fSimulatePowerplant();

    // Objective function : total investment cost
    f1 = _powerplant->get_costOfHeliostatField()
      + _powerplant->get_costOfTower()
      + _powerplant->get_costOfReceiver()
      + _powerplant->get_costOfStorage()
      + _powerplant->get_costOfSteamGenerator()
      + _powerplant->get_costOfPowerblock();

    // Verify: total land area is below 700k m^2
    c1 = M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
      *(2 * _fieldMaxAngularSpan / 360.) - _cFieldSurface;

    // Verify: compliance to demand is 100%
    c2 = _cDemandComplianceRatio - _powerplant->get_overallComplianceToDemand();

    // Verify: tower at least twice as high as heliostats
    c3 = 2 * _heliostatHeight - _towerHeight;

    // Verify: Rmin < Rmax
    c4 = _minimumDistanceToTower - _maximumDistanceToTower;
    
    // Verify: number of heliostats fits in the field
    c5 = _maxNumberOfHeliostats - _powerplant->get_heliostatField()->get_listOfHeliostats().size();
    
    // Verify: maximum pressure in receiver tubes do not exceed yield pressure
    c6 = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();

    // Verify: Molten Salt temperature does not drop below the melting point
    c7 = MELTING_POINT - _powerplant->get_minHotStorageTemp();
    c8 = MELTING_POINT - _powerplant->get_minColdStorageTemp();
    c9 = MELTING_POINT - _powerplant->get_minSteamGenTemp();
    
    // Verify: Receiver tubes Din < Dout
    c10 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;

    // Verify: Tubes fit in receiver
    c11 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;

    // Verify: central receiver outlet is higher than that required by the turbine
    c12 = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;

    double sum = 1;
    for (unsigned int i = 0; i < _powerplant->get_powerplantPowerOutput().size(); ++i)
      {
	sum += _powerplant->get_powerplantPowerOutput()[i];
      }
    c13 = _powerplant->fComputeParasiticLosses()/sum - _cParasitics;
    c14 = _exchangerTubesDout - _exchangerTubesSpacing;
    c15 = _exchangerTubesDin - _exchangerTubesDout;
    c16 = _powerplant->get_maximumPressureInExchanger() - _powerplant->get_yieldPressureInExchanger();
    
    out << f1 << " " << c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6 << " " <<
      c7 << " " << c8 << " " << c9 << " " << c10 << " " << c11 << " " << 
      c12 << " " << c13 << " " << c14 << " " << c15 << " " << c16
	<< std::endl;
  }
  catch (...){
    c1 = M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
      *(2 * _fieldMaxAngularSpan / 360.) - _cFieldSurface;
    c3 = 2 * _heliostatHeight - _towerHeight;
    c4 = _minimumDistanceToTower - _maximumDistanceToTower;
    
    c10 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    c11 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;
    c12 = _minReceiverOutletTemp - _centralReceiverOutletTemperature;
    
    c14 = _exchangerTubesDout - _exchangerTubesSpacing;
    c15 = _exchangerTubesDin - _exchangerTubesDout;
    out << f1 << " " << c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6 << " " <<
      c7 << " " << c8 << " " << c9 << " " << c10 << " " << c11 << " " <<
      c12 << " " << c13 << " " << c14 << " " << c15 << " " << c16
	<< std::endl;
    
    throw logic_error ( "Simulation could not go through" );
  }
  
  return true;
}

/*------------------------------*/
/*  simulate_maxComp_HTF1 (#5)  */
/*------------------------------*/
bool Scenario::simulate_maxComp_HTF1 ( std::ostream & out , bool & cnt_eval ) {

  cnt_eval = false;

  double
    f1 = 1e20,  c1 = 1e20,  c2 = 1e20,  c3 = 1e20,
    c4 = 1e20,  c5 = 1e20,  c6 = 1e20,  c7 = 1e20,
    c8 = 1e20,  c9 = 1e20, c10 = 1e20, c11 = 1e20, c12 = 1e20;

  try{

    //Verifying a priori constraints
    cnt_eval = validateInputValues();

    //Creating required objects
    construct_maxComp_HTF1();

    //Launching simulation
    _powerplant->fSimulatePowerplant();

    //Objective function : time for which the demand is met
    f1 = - _powerplant->get_overallComplianceToDemand();

    //Verify : cost constraint
    double totalCost = _powerplant->get_costOfReceiver()
      + _powerplant->get_costOfStorage()
      + _powerplant->get_costOfSteamGenerator()
      + _powerplant->get_costOfPowerblock();
    
    c1 =  totalCost - _cBudget;

    //Verify : pressure in receiver tubes doesn't exceed yield pressure
    c2 = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();

    //Verify : Molten Salt temperature does not drop below the melting point
    c3 = MELTING_POINT - _powerplant->get_minHotStorageTemp();

    c4 = MELTING_POINT - _powerplant->get_minColdStorageTemp();

    c5 = MELTING_POINT - _powerplant->get_minSteamGenTemp();
    
    //Verify : Receiver tubes Din < Dout
    c6 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    //Verify : Tubes fit in receiver
    c7 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;
    
    //Verify : central receiver outlet is higher than that required by the turbine
    c8 = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;
    
    //Verify : parasitics do not exceed limit
    double sum = 1;
    foncteurSum somme(&sum);
    std::for_each(_powerplant->get_powerplantPowerOutput().begin(),
		  _powerplant->get_powerplantPowerOutput().end(),
		  somme);
    c9 = _powerplant->fComputeParasiticLosses() / sum - _cParasitics;
    
    //Verify : spacing and tubes dimensions in steam gen.
    c10 = _exchangerTubesDout - _exchangerTubesSpacing;
    c11 = _exchangerTubesDin - _exchangerTubesDout;
    
    //Verify : pressure in steam gen. tubes do not exceed yield
    c12 = _powerplant->get_maximumPressureInExchanger() - _powerplant->get_yieldPressureInExchanger();
    
    out << f1 << " " << c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6 << " " << c7 << " " << 
      c8 << " " << c9 << " " << c10 << " " << c11 << " " << c12
	<< std::endl;
  }
  catch (...){
    
    //Verify : Receiver tubes Din < Dout
    c6 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    //Verify : Tubes fit in receiver
    c7 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;
    
    //Verify : central receiver outlet is higher than that required by the turbine
    c8 = _minReceiverOutletTemp - _centralReceiverOutletTemperature;
    
    //Verify : spacing and tubes dimensions in steam gen.
    c10 = _exchangerTubesDout - _exchangerTubesSpacing;
    c11 = _exchangerTubesDin - _exchangerTubesDout;
    
    out << f1 << " " << c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6 << " " << c7 << " " <<
      c8 << " " << c9 << " " << c10 << " " << c11 << " " << c12
	<< std::endl;
    
    throw logic_error ( "Simulation could not go through" );
  }
  return true;
}


/*----------------------------*/
/*  simulate_minCost_TS (#6)  */
/*----------------------------*/
bool Scenario::simulate_minCost_TS ( std::ostream & out , bool & cnt_eval ) {

  cnt_eval = false;
  
  double
    f1 = 1e20, c1 = 1e20, c2 = 1e20,
    c3 = 1e20, c4 = 1e20, c5 = 1e20, c6 = 1e20;
	
  try{

    //Verifying a priori constraints
    cnt_eval = validateInputValues();

    //Creating required objects
    construct_minCost_TS();

    //Launching simulation
    _powerplant->fSimulatePowerplant();

    //Objective function : cost of storage
    f1 = _powerplant->get_costOfStorage();

    //Verify : demand met 100%
    c1 = _cDemandComplianceRatio - _powerplant->get_overallComplianceToDemand();

    //Verify : pressure in receiver tubes doesn't exceed yield pressure
    c2 = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();
    
    //Verify : Molten Salt temperature does not drop below the melting point
    c3 = MELTING_POINT - _powerplant->get_minHotStorageTemp();

    c4 = MELTING_POINT - _powerplant->get_minColdStorageTemp();

    //Verify : central receiver outlet is higher than that required by the turbine
    c5 = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;

    //Verify: storage is back to initial conditions
    double storageFinalConditions = 0.;
    storageFinalConditions = _powerplant->get_moltenSaltLoop()->get_hotStorage().get_heightOfVolumeStored()
      /_hotStorageHeight;
    c6 = 1.*_storageStartupCondition - storageFinalConditions*100;
        
    out << f1 << " " << c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6
	<< std::endl;
  }
  catch (...) {
    
    //Verify : central receiver outlet is higher than that required by the turbine
    c5 = _minReceiverOutletTemp - _centralReceiverOutletTemperature;

    out << f1 << " " << c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6
	<< std::endl;
    
    throw logic_error ( "Simulation could not go through" );
  }
  
  return true;
}

/*-----------------------------*/
/*   simulate_maxEff_Re (#7)   */
/*-----------------------------*/
bool Scenario::simulate_maxEff_RE ( std::ostream & out , bool & cnt_eval ) {

  cnt_eval = false;

  double f1 = 1e20, c1 = 1e20, c2 = 1e20, c3 = 1e20, c4 = 1e20, c5 = 1e20, c6 = 1e20;

  try{

    //Verifying a priori constraints
    cnt_eval = validateInputValues();

    //Creating required objects
    construct_maxEff_RE();

    //Launching simulation
    _powerplant->fSimulatePowerplant();

    //Objective function : absorbed energy
    double Q_abs = _powerplant->get_moltenSaltLoop()->get_hotStorage().get_storedMass()
      *(_centralReceiverOutletTemperature - _coldMoltenSaltMinTemperature)
      *HEAT_CAPACITY;
    f1 = -Q_abs*1e-9;

    //Verify : budget is respected
    c1 = _powerplant->get_costOfReceiver() - _cBudget;
    
    //Verify : maximum pressure in tubes < yield
    c2 = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();
    
    //Verify : receiver inner tubes diameter < outer diameter
    c3 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    //Verify : central receiver outlet is higher than that required by the turbine
    c4 = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;
    
    //Verify : tubes fit in receiver
    c5 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;
    
    //Verify : work to drive receiver pump does not exceed 5% of the absorbed energy
    if (Q_abs > 0.01){
      c6 = _powerplant->fComputeParasiticsForPb7() / Q_abs - _cParasitics;
    }
    else
      c6 = 1.0;

    out << f1 << " " << c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6
	<< std::endl;
  }
  catch (...){
    //Verify : tower at least twice as high as heliostats
    c3 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    //Verify : central receiver outlet is higher than that required by the turbine
    c4 = _minReceiverOutletTemp - _centralReceiverOutletTemperature;
    
    //Verify : tubes fit in receiver
    c5 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;
    
    out << f1 << " " << c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6
	<< std::endl;
    
    throw logic_error ( "Simulation could not go through" );
  }
  return true;
}

/*----------------------------------*/
/*    simulate_maxHF_minCost (#8)   */
/*----------------------------------*/
bool Scenario::simulate_maxHF_minCost ( std::ostream & out , bool & cnt_eval ) {

  cnt_eval = false;

  double
    f1 = 1e20, f2 = 1e20, c1 = 1e20, c2 = 1e20,
    c3 = 1e20, c4 = 1e20, c5 = 1e20, c6 = 1e20,
    c7 = 1e20, c8 = 1e20, c9 = 1e20;
  
  try{
	  
    //Verifying a priori constraints
    cnt_eval = validateInputValues();
    
    //Creating required objects
    construct_maxHF_minCost();
    
    //Launching simulation
    _powerplant->fSimulatePowerplant();

    //Objective function : absorbed energy
    double Q_abs = _powerplant->get_moltenSaltLoop()->get_hotStorage().get_storedMass()
      *(_centralReceiverOutletTemperature - _coldMoltenSaltMinTemperature)
      *HEAT_CAPACITY;
    f1 = -Q_abs;

    //Objective function : total investment cost
    f2 = _powerplant->get_costOfHeliostatField()
      + _powerplant->get_costOfTower()
      + _powerplant->get_costOfReceiver();
    
    //Verify : total land area
    c1 = M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
      *(2 * _fieldMaxAngularSpan / 360.) - _cFieldSurface;

    //Verify : tower at least twice as high as heliostats
    c2 = 2 * _heliostatHeight - _towerHeight;
    
    //Verify : Rmin < Rmax
    c3 = _minimumDistanceToTower - _maximumDistanceToTower;
    
    //Verify : number of heliostats fits in the field
    c4 = _maxNumberOfHeliostats - _powerplant->get_heliostatField()->get_listOfHeliostats().size();
    
    //Verify : maximum pressure in receiver tubes do not exceed yield pressure
    c5 = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();
    
    //Verify : Receiver tubes Din < Dout
    c6 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    //Verify : Tubes fit in receiver
    c7 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;
    
    //Verify : minimal energy production satisfied
    c8 = ((400e6)*3600.) - Q_abs;
    
    //Verify : work to drive receiver pump does not exceed 8% of the absorbed energy
    double parasitics = 0;
    parasitics = _powerplant->fComputeParasiticsForPb9();
    if (Q_abs > 1.){
      c9 = (parasitics / Q_abs) - _cParasitics;
    }
    
    out << f1 << " " << f2 << " " <<  c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6 << " " <<
      c7 << " " << c8 << " " << c9
	<< std::endl;
  }
  catch (...){
    
    //Verify : total land area
    c1 = M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
      *(2 * _fieldMaxAngularSpan / 360.) - _cFieldSurface;
    
    //Verify : tower at least twice as high as heliostats
    c2 = 2 * _heliostatHeight - _towerHeight;
    
    //Verify : Rmin < Rmax
    c3 = _minimumDistanceToTower - _maximumDistanceToTower;
    
    //Verify : Receiver tubes Din < Dout
    c6 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    //Verify : Tubes fit in receiver
    c7 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;
    
    out << f1 << " " << f2 << " " << c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6 << " " <<
      c7 << " " << c8 << " " << c9
	<< std::endl;
    
    throw logic_error ( "Simulation could not go through" );
  }
  return true;
}

/*----------------------------------*/
/*    simulate_maxNrg_minPar (#9)   */
/*----------------------------------*/
bool Scenario::simulate_maxNrg_minPar ( std::ostream & out , bool & cnt_eval ) {

  cnt_eval = false;

  double
    f1  = 1e20,  f2 = 1e20,  c1 = 1e20,  c2 = 1e20,
    c3  = 1e20,  c4 = 1e20,  c5 = 1e20,  c6 = 1e20,
    c7  = 1e20,  c8 = 1e20,  c9 = 1e20, c10 = 1e20,
    c11 = 1e20, c12 = 1e20, c13 = 1e20, c14 = 1e20,
    c15 = 1e20, c16 = 1e20, c17 = 1e20;

  try{

    //Verifying a priori constraints
    cnt_eval = validateInputValues();

    //Creating required objects
    construct_maxNrg_minPar();

    //Launching simulation
    _powerplant->fSimulatePowerplant();
    
    //Objective : power output (MWe)
    double sum = 1.;
    foncteurSum somme(&sum);
    std::for_each(_powerplant->get_powerplantPowerOutput().begin(),
		  _powerplant->get_powerplantPowerOutput().end(),
		  somme);

    f1 = -sum;

    //Objective : parasitic losses
    f2 = _powerplant->fComputeParasiticLosses();
    
    //Verify : total investment cost
    c1 = _powerplant->get_costOfHeliostatField()
      + _powerplant->get_costOfTower()
      + _powerplant->get_costOfReceiver()
      + _powerplant->get_costOfStorage()
      + _powerplant->get_costOfSteamGenerator()
      + _powerplant->get_costOfPowerblock()
      - _cBudget;
    
    c2 = 3600. * 120.e6 - sum;
    
    //Verify : total land area is below 2000k m^2
    c3 = M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
      *(2 * _fieldMaxAngularSpan / 360.) - _cFieldSurface;
    
    //Verify : tower at least twice as high as heliostats
    c4 = 2 * _heliostatHeight - _towerHeight;
    
    //Verify : Rmin < Rmax
    c5 = _minimumDistanceToTower - _maximumDistanceToTower;
    
    //Verify : number of heliostats fits in the field
    c6 = _maxNumberOfHeliostats - _powerplant->get_heliostatField()->get_listOfHeliostats().size();
    
    //Verify : maximum pressure in receiver tubes do not exceed yield pressure
    c7 = _powerplant->get_maximumPressureInReceiver() - _powerplant->get_yieldPressureInReceiver();
        
    //Verify : Molten Salt temperature does not drop below the melting point
    c8 = MELTING_POINT - _powerplant->get_minHotStorageTemp();
    
    c9 = MELTING_POINT - _powerplant->get_minColdStorageTemp();
    
    c10 = MELTING_POINT - _powerplant->get_minSteamGenTemp();
    
    //Verify : Receiver tubes Din < Dout
    c11 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    //Verify : Tubes fit in receiver
    c12 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;
    
    //Verify : central receiver outlet is higher than that required by the turbine
    c13 = _powerplant->get_steamTurbineInletTemperature() - _centralReceiverOutletTemperature;
    
    if (sum > 1.){
      c14 = f2/sum - 0.20;
    }
    c15 = _exchangerTubesDout - _exchangerTubesSpacing;
    c16 = _exchangerTubesDin - _exchangerTubesDout;
    c17 = _powerplant->get_maximumPressureInExchanger() - _powerplant->get_yieldPressureInExchanger();
        
    out << f1 << " " << f2 << " " << c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6 << " " <<
      c7 << " " << c8 << " " << c9 << " " << c10 << " " << c11 << " " <<
      c12 << " " << c13 << " " << c14 << " " << c15 << " " << c16 << " " << c17
	<< std::endl;
  }
  catch (...){
    //Verify : total land area is below 2000k m^2
    c3 = M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
      *(2 * _fieldMaxAngularSpan / 360.) - _cFieldSurface;
    
    //Verify : tower at least twice as high as heliostats
    c4 = 2 * _heliostatHeight - _towerHeight;
    
    //Verify : Rmin < Rmax
    c5 = _minimumDistanceToTower - _maximumDistanceToTower;
    
    //Verify : Receiver tubes Din < Dout
    c11 = _receiverTubesInsideDiam - _receiverTubesOutsideDiam;
    
    //Verify : Tubes fit in receiver
    c12 = _receiverNbOfTubes*_receiverTubesOutsideDiam - _receiverApertureWidth * M_PI / 2.;
    
    //Verify : central receiver outlet is higher than that required by the turbine
    c13 = _minReceiverOutletTemp - _centralReceiverOutletTemperature;
    
    c15 = _exchangerTubesDout - _exchangerTubesSpacing;
    c16 = _exchangerTubesDin - _exchangerTubesDout;
    
    out << f1 << " " << f2 << " " << c1 << " " << c2 << " " <<
      c3 << " " << c4 << " " << c5 << " " << c6 << " " <<
      c7 << " " << c8 << " " << c9 << " " << c10 << " " << c11 << " " <<
      c12 << " " << c13 << " " << c14 << " " << c15 << " " << c16 << " " << c17
	<< std::endl;
    
    throw logic_error ( "Simulation could not go through" );
  }
  return true;
}

/*-----------------------------------------------------------------*/
/*  validation functions (check bounds and a priori constraints)   */
/*  return true if all is ok and the evaluation should be counted  */
/*-----------------------------------------------------------------*/
bool Scenario::validateInputValues ( void ) const {
  if ( _problem == "MAXNRG_H1"     || _problem == "MAXNRG_H1S"     ) { return validate_maxNrg_H1();     } // #1
  if ( _problem == "MINSURF_H1"    || _problem == "MINSURF_H1S"    ) { return validate_minSurf_H1();    } // #2
  if ( _problem == "MINCOST_C1"    || _problem == "MINCOST_C1S"    ) { return validate_minCost_C1();    } // #3
  if ( _problem == "MINCOST_C2"    || _problem == "MINCOST_C2S"    ) { return validate_minCost_C2();    } // #4
  if ( _problem == "MAXCOMP_HTF1"  || _problem == "MAXCOMP_HTF1S"  ) { return validate_maxComp_HTF1();  } // #5
  if ( _problem == "MINCOST_TS"    || _problem == "MINCOST_TSS"    ) { return validate_minCost_TS();    } // #6
  if ( _problem == "MAXEFF_RE"     || _problem == "MAXEFF_RES"     ) { return validate_maxEff_RE();     } // #7
  if ( _problem == "MAXHF_MINCOST" || _problem == "MAXHF_MINCOSTS" ) { return validate_maxHF_minCost(); } // #8
  if ( _problem == "MAXNRG_MINPAR" || _problem == "MAXNRG_MINPARS" ) { return validate_maxNrg_minPar(); } // #9

  return false;
}

/*-------------------------------------------------*/
/*          validate problem maxNrg_H1 (#1)        */
/*-------------------------------------------------*/
bool Scenario::validate_maxNrg_H1 ( void ) const {

  // Bounds:
  // -------
  if ( _heliostatHeight < 1 || _heliostatHeight > 40 )
    throw invalid_argument ( "detected unacceptable value for heliostat height" );

  if ( _heliostatWidth < 1 || _heliostatWidth > 40 )
    throw invalid_argument ( "detected unacceptable value for heliostat width" );

  if ( _towerHeight < 20 || _towerHeight > 250. )
    throw invalid_argument ( "detected unacceptable value for tower height" );

  if ( _receiverApertureHeight < 1 || _receiverApertureHeight > 30 )
    throw invalid_argument ( "detected unacceptable value for receiver aperture height" );

  if ( _receiverApertureWidth < 1 || _receiverApertureWidth > 30 )
    throw invalid_argument ( "detected unacceptable value for receiver aperture width" );

  if ( _maxNumberOfHeliostats < 1 )
    throw invalid_argument ( "detected unacceptable value for number of heliostats" );

  if ( _fieldMaxAngularSpan < 1 || _fieldMaxAngularSpan > 89 )
    throw invalid_argument ( "detected unacceptable value for field max angular span" );
  
  if ( _minimumDistanceToTower < 0 || _minimumDistanceToTower > 20 )
    throw invalid_argument ( "detected unacceptable value for minimum distance to tower" );

  if ( _maximumDistanceToTower < 0 || _maximumDistanceToTower > 20 )
    throw invalid_argument ( "detected unacceptable value for maximum distance to tower" );

  // A priori constraints:
  // ---------------------
  if ( _towerHeight < 2 * _heliostatHeight )
    throw invalid_argument ( "detected unacceptable value for tower height" );
    
  if ( _maximumDistanceToTower <= _minimumDistanceToTower )
    throw invalid_argument ( "detected unacceptable value (1)" );
  
  if ( M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
       *(2 * _fieldMaxAngularSpan / 360.) > _cFieldSurface )
    throw invalid_argument ( "detected unacceptable value (2)" );

  return true;
}

/*-------------------------------------------------*/
/*         validate problem minSurf_H1 (#2)        */
/*-------------------------------------------------*/
bool Scenario::validate_minSurf_H1 ( void ) const {

  invalid_argument inputOutOfRange ( "detected unacceptable value for one or more input parameters" );

  if ( _heliostatHeight < 1. || _heliostatHeight > 40 )
    throw inputOutOfRange;

  if ( _heliostatWidth < 1. || _heliostatWidth > 40 )
    throw inputOutOfRange;

  if ( _towerHeight < 20 || _towerHeight > 250 )
    throw inputOutOfRange;

  if ( _receiverApertureHeight < 1 || _receiverApertureHeight > 30 )
    throw inputOutOfRange;

  if ( _receiverApertureWidth < 1 || _receiverApertureWidth > 30 )
    throw inputOutOfRange;

  if ( _maxNumberOfHeliostats < 1)
    throw inputOutOfRange;

  if ( _fieldMaxAngularSpan < 1 || _fieldMaxAngularSpan > 89 )
    throw inputOutOfRange;

  if ( _minimumDistanceToTower < 0 || _minimumDistanceToTower > 20 )
    throw inputOutOfRange;
	
  if ( _maximumDistanceToTower < 1 || _maximumDistanceToTower > 20 )
    throw inputOutOfRange;

  if ( _towerHeight < 2 * _heliostatHeight)
    throw inputOutOfRange;
  
  if ( _maximumDistanceToTower <= _minimumDistanceToTower )
    throw inputOutOfRange;
  
  if ( M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
       *(2 * _fieldMaxAngularSpan / 360.) > _cFieldSurface )
    throw inputOutOfRange;

  double minReceiverOutletTemp = 0.;
  switch ( _typeOfTurbine ) {
  case 1:
    minReceiverOutletTemp = SST110_TEMPERATURE;
    break;
  case 2:
    minReceiverOutletTemp = SST120_TEMPERATURE;
    break;
  case 3:
    minReceiverOutletTemp = SST300_TEMPERATURE;
    break;
  case 4:
    minReceiverOutletTemp = SST400_TEMPERATURE;
    break;
  case 5:
    minReceiverOutletTemp = SST600_TEMPERATURE;
    break;
  case 6:
    minReceiverOutletTemp = SST700_TEMPERATURE;
    break;
  case 7:
    minReceiverOutletTemp = SST800_TEMPERATURE;
    break;
  case 8:
    minReceiverOutletTemp = SST900_TEMPERATURE;
    break;
  }
  
  if ( _centralReceiverOutletTemperature < minReceiverOutletTemp || _centralReceiverOutletTemperature > 995 )
    throw inputOutOfRange;

  if ( _receiverNbOfTubes < 1 )
    throw inputOutOfRange;

  if ( _receiverInsulThickness < 0.01 || _receiverInsulThickness > 5 )
    throw inputOutOfRange;

  if ( _receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1 )
    throw inputOutOfRange;

  if ( _receiverTubesOutsideDiam < 0.005 || _receiverTubesOutsideDiam > 0.1 )
    throw inputOutOfRange;
	
  if ( _receiverTubesOutsideDiam < _receiverTubesInsideDiam + 0.0005 )
    throw inputOutOfRange;

   return true;
}

/*-------------------------------------------------*/
/*          validate problem minCost_C1 (#3)       */
/*-------------------------------------------------*/
bool Scenario::validate_minCost_C1 ( void ) const {
 
  invalid_argument inputOutOfRange ( "detected unacceptable value for one or more input parameters" );

  // check bounds and other a priori constraints:
  if (_heliostatHeight < 1 || _heliostatHeight > 40)
    throw inputOutOfRange;
    
  if (_heliostatWidth < 1 || _heliostatWidth > 40)
    throw inputOutOfRange;
    
  if (_towerHeight < 20 || _towerHeight > 250)
    throw inputOutOfRange;
  
  if (_receiverApertureHeight < 1 || _receiverApertureHeight > 30)
    throw inputOutOfRange;

  if (_receiverApertureWidth < 1 || _receiverApertureWidth > 30)
    throw inputOutOfRange;
  
  if (_maxNumberOfHeliostats < 1)
    throw inputOutOfRange;

  if (_fieldMaxAngularSpan < 1 || _fieldMaxAngularSpan > 89)
    throw inputOutOfRange;

  if (_minimumDistanceToTower < 0 || _minimumDistanceToTower > 20)
    throw inputOutOfRange;

  if (_maximumDistanceToTower < 1 || _maximumDistanceToTower > 20)
    throw inputOutOfRange;

  if (_towerHeight < 2 * _heliostatHeight)
    throw inputOutOfRange;
    
  if (_maximumDistanceToTower <= _minimumDistanceToTower)
    throw inputOutOfRange;
  
  if (M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
      *(2 * _fieldMaxAngularSpan / 360.) > _cFieldSurface)
    throw inputOutOfRange;
      
  double minReceiverOutletTemp = 0.;
  switch (_typeOfTurbine) {
  case 1:
    minReceiverOutletTemp = SST110_TEMPERATURE;
    break;
  case 2:
    minReceiverOutletTemp = SST120_TEMPERATURE;
    break;
  case 3:
    minReceiverOutletTemp = SST300_TEMPERATURE;
    break;
  case 4:
    minReceiverOutletTemp = SST400_TEMPERATURE;
    break;
  case 5:
    minReceiverOutletTemp = SST600_TEMPERATURE;
    break;
  case 6:
    minReceiverOutletTemp = SST700_TEMPERATURE;
    break;
  case 7:
    minReceiverOutletTemp = SST800_TEMPERATURE;
    break;
  case 8:
    minReceiverOutletTemp = SST900_TEMPERATURE;
    break;
  }
  
  if (_centralReceiverOutletTemperature < minReceiverOutletTemp || _centralReceiverOutletTemperature > 995)
    throw inputOutOfRange;
  
  if (_receiverNbOfTubes * _receiverTubesOutsideDiam > M_PI*_receiverApertureWidth / 2.)
    throw inputOutOfRange;
  
  if (_hotStorageHeight < 1 || _hotStorageHeight > 50)
    throw inputOutOfRange;

  if (_hotStorageDiameter < 1 || _hotStorageDiameter > 30)
    throw inputOutOfRange;
      
  if (_hotStorageInsulThickness < 0.01 || _hotStorageInsulThickness > 5)
    throw inputOutOfRange;
    
  if (_coldStorageInsulThickness < 0.01 || _coldStorageInsulThickness > 5)
    throw inputOutOfRange;
  
  if (_coldMoltenSaltMinTemperature < MELTING_POINT || _coldMoltenSaltMinTemperature > 650)
    throw inputOutOfRange;
  
  if (_receiverNbOfTubes < 1)
    throw inputOutOfRange;

  if (_receiverInsulThickness < 0.01 || _receiverInsulThickness > 5)
    throw inputOutOfRange;
  
  if (_receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1)
    throw inputOutOfRange;
  
  if (_receiverTubesOutsideDiam < 0.005 || _receiverTubesOutsideDiam > 0.1)
    throw inputOutOfRange;

  if (_receiverTubesOutsideDiam < _receiverTubesInsideDiam + 0.0005)
    throw inputOutOfRange;

  if ( _typeOfTurbine < 1 || _typeOfTurbine > 8 )
    throw inputOutOfRange;

  return true;
}

/*-------------------------------------------------*/
/*          validate problem minCost_C2 (#4)       */
/*-------------------------------------------------*/
bool Scenario::validate_minCost_C2 ( void ) const {

  invalid_argument inputOutOfRange ( "detected unacceptable value for one or more input parameters" );

  if (_heliostatHeight < 1 || _heliostatHeight > 40)
    throw inputOutOfRange;

  if (_heliostatWidth < 1 || _heliostatWidth > 40)
    throw inputOutOfRange;

  if (_towerHeight < 20 || _towerHeight > 250)
    throw inputOutOfRange;
  
  if (_towerHeight < 2 * _heliostatHeight)
    throw inputOutOfRange;
  
  if (_receiverApertureHeight < 1 || _receiverApertureHeight > 30)
    throw inputOutOfRange;
  
  if (_receiverApertureWidth < 1 || _receiverApertureWidth > 30)
    throw inputOutOfRange;
  
  if (_maxNumberOfHeliostats < 1)
    throw inputOutOfRange;
  
  if (_fieldMaxAngularSpan < 1 || _fieldMaxAngularSpan > 89)
    throw inputOutOfRange;
  
  if (_minimumDistanceToTower < 0 || _minimumDistanceToTower > 20)
    throw inputOutOfRange;
  
  if (_maximumDistanceToTower < 1 || _maximumDistanceToTower > 20)
    throw inputOutOfRange;
  
  if (_maximumDistanceToTower <= _minimumDistanceToTower)
    throw inputOutOfRange;
  
  if ( M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
      *(2 * _fieldMaxAngularSpan / 360.) > _cFieldSurface )
    throw inputOutOfRange;
  
  double minReceiverOutletTemp = 0.;
  switch ( _typeOfTurbine ) {
  case 1:
    minReceiverOutletTemp = SST110_TEMPERATURE;
    break;
  case 2:
    minReceiverOutletTemp = SST120_TEMPERATURE;
    break;
  case 3:
    minReceiverOutletTemp = SST300_TEMPERATURE;
    break;
  case 4:
    minReceiverOutletTemp = SST400_TEMPERATURE;
    break;
  case 5:
    minReceiverOutletTemp = SST600_TEMPERATURE;
    break;
  case 6:
    minReceiverOutletTemp = SST700_TEMPERATURE;
    break;
  case 7:
    minReceiverOutletTemp = SST800_TEMPERATURE;
    break;
  case 8:
    minReceiverOutletTemp = SST900_TEMPERATURE;
    break;
  }

  if ( _centralReceiverOutletTemperature < minReceiverOutletTemp || _centralReceiverOutletTemperature > 995 )
    throw inputOutOfRange;

  //         x[15]              x[18]                               x[4]
  if (_receiverNbOfTubes * _receiverTubesOutsideDiam > M_PI*_receiverApertureWidth / 2.) {
    throw inputOutOfRange;
  }
 
  if (_hotStorageHeight < 1 || _hotStorageHeight > 50)
    throw inputOutOfRange;
 
  if (_hotStorageDiameter < 1 || _hotStorageDiameter > 30)
    throw inputOutOfRange;
 
  if (_hotStorageInsulThickness < 0.01 || _hotStorageInsulThickness > 5)
    throw inputOutOfRange;
  
  if (_coldStorageInsulThickness < 0.01 || _coldStorageInsulThickness > 5)
    throw inputOutOfRange;
  
  if (_coldMoltenSaltMinTemperature < MELTING_POINT || _coldMoltenSaltMinTemperature > 650)
    throw inputOutOfRange;
  
  if (_receiverNbOfTubes < 1)
    throw inputOutOfRange;
  
  if (_receiverInsulThickness < 0.01 || _receiverInsulThickness > 5)
    throw inputOutOfRange;

  if (_receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1)
    throw inputOutOfRange;

  if (_receiverTubesOutsideDiam < 0.006 || _receiverTubesOutsideDiam > 0.1)
    throw inputOutOfRange;

  if (_receiverTubesOutsideDiam < _receiverTubesInsideDiam + 0.0001)
    throw inputOutOfRange;

  if (_exchangerTubesSpacing <= _exchangerTubesDout || _exchangerTubesSpacing > 0.3)
    throw inputOutOfRange;

  if (_exchangerTubesLength < 0.5 || _receiverTubesOutsideDiam > 10)
    throw inputOutOfRange;
  
  if (_exchangerTubesDin < 0.005 || _exchangerTubesDin > 0.1)
    throw inputOutOfRange;

  if (_exchangerTubesDout < 0.006 || _exchangerTubesDout > 0.1)
    throw inputOutOfRange;

  if (_exchangerTubesDout < _exchangerTubesDin + 0.0005)
    throw inputOutOfRange;

  if (_exchangerNbOfBaffles < 2)
    throw inputOutOfRange;

  if (_exchangerBaffleCut < 0.15 || _exchangerBaffleCut > 0.4)
    throw inputOutOfRange;

  if (_exchangerNbOfTubes < 1)
    throw inputOutOfRange;

  if (_exchangerNbOfShells < 1 || _exchangerNbOfShells > 10)
    throw inputOutOfRange;

  if (_exchangerNbOfPassesPerShell < 1 || _exchangerNbOfPassesPerShell > 9)
    throw inputOutOfRange;

  if (_typeOfTurbine < 1 || _typeOfTurbine >  8)
    throw inputOutOfRange;

   return true;
}

/*-------------------------------------------------*/
/*         validate problem maxComp_HTF1 (#5)      */
/*-------------------------------------------------*/
bool Scenario::validate_maxComp_HTF1 ( void ) const {

  invalid_argument inputOutOfRange ( "detected unacceptable value for one or more input parameters" );

  double minReceiverOutletTemp = 0.;
  switch ( _typeOfTurbine ) {
  case 1:
    minReceiverOutletTemp = SST110_TEMPERATURE;
    break;
  case 2:
    minReceiverOutletTemp = SST120_TEMPERATURE;
    break;
  case 3:
    minReceiverOutletTemp = SST300_TEMPERATURE;
    break;
  case 4:
    minReceiverOutletTemp = SST400_TEMPERATURE;
    break;
  case 5:
    minReceiverOutletTemp = SST600_TEMPERATURE;
    break;
  case 6:
    minReceiverOutletTemp = SST700_TEMPERATURE;
    break;
  case 7:
    minReceiverOutletTemp = SST800_TEMPERATURE;
    break;
  case 8:
    minReceiverOutletTemp = SST900_TEMPERATURE;
    break;
  }
  
  if ( _centralReceiverOutletTemperature < minReceiverOutletTemp || _centralReceiverOutletTemperature > 995 )
    throw inputOutOfRange;
  
  if (_hotStorageHeight < 1 || _hotStorageHeight > 30)
    throw inputOutOfRange;
  
  if (_hotStorageDiameter < 1 || _hotStorageDiameter > 30)
    throw inputOutOfRange;
  
  if (_hotStorageInsulThickness < 0.01 || _hotStorageInsulThickness > 2)
    throw inputOutOfRange;
  
  if (_coldStorageInsulThickness < 0.01 || _coldStorageInsulThickness > 2)
    throw inputOutOfRange;
  
  if (_coldMoltenSaltMinTemperature < MELTING_POINT || _coldMoltenSaltMinTemperature > 650)
    throw inputOutOfRange;
  
  if (_receiverNbOfTubes < 1)
    throw inputOutOfRange;
  
  if (_receiverInsulThickness < 0.1 || _receiverInsulThickness > 2)
    throw inputOutOfRange;
  
  if (_receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1)
    throw inputOutOfRange;
  
  if (_receiverTubesOutsideDiam < 0.005 || _receiverTubesOutsideDiam > 0.1)
    throw inputOutOfRange;
  
  if (_receiverTubesOutsideDiam < _receiverTubesInsideDiam + 0.0005)
    throw inputOutOfRange;
  
  if (_exchangerTubesSpacing <= _exchangerTubesDout || _exchangerTubesSpacing > 0.2)
    throw inputOutOfRange;
  
  if (_exchangerTubesLength < 0.5 || _exchangerTubesLength > 10)
    throw inputOutOfRange;
  
  if (_exchangerTubesDin < 0.005 || _exchangerTubesDin > 0.1)
    throw inputOutOfRange;
  
  if (_exchangerTubesDout < 0.006 || _exchangerTubesDout > 0.1)
    throw inputOutOfRange;
  
  if (_exchangerTubesDout < _exchangerTubesDin + 0.0005)
    throw inputOutOfRange;
  
  if (_exchangerNbOfBaffles < 2)
    throw inputOutOfRange;
  
  if (_exchangerBaffleCut < 0.15 || _exchangerBaffleCut > 0.4)
    throw inputOutOfRange;
  
  if (_exchangerNbOfTubes < 1)
    throw inputOutOfRange;
  
  if (_exchangerNbOfShells < 1 || _exchangerNbOfShells > 10)
    throw inputOutOfRange;
  
  if (_exchangerNbOfPassesPerShell < 1 || _exchangerNbOfPassesPerShell > 9)
    throw inputOutOfRange;
  
  if (_typeOfTurbine < 1 || _typeOfTurbine >  8)
    throw inputOutOfRange;

   return true;
}

/*-------------------------------------------------*/
/*         validate problem minCost_TS (#6)        */
/*-------------------------------------------------*/
bool Scenario::validate_minCost_TS ( void ) const {

  invalid_argument inputOutOfRange ( "detected unacceptable value for one or more input parameters" );

  double minReceiverOutletTemp = 0.;
  switch ( _typeOfTurbine ) {
  case 1:
    minReceiverOutletTemp = SST110_TEMPERATURE;
    break;
  case 2:
    minReceiverOutletTemp = SST120_TEMPERATURE;
    break;
  case 3:
    minReceiverOutletTemp = SST300_TEMPERATURE;
    break;
  case 4:
    minReceiverOutletTemp = SST400_TEMPERATURE;
    break;
  case 5:
    minReceiverOutletTemp = SST600_TEMPERATURE;
    break;
  case 6:
    minReceiverOutletTemp = SST700_TEMPERATURE;
    break;
  case 7:
    minReceiverOutletTemp = SST800_TEMPERATURE;
    break;
  case 8:
    minReceiverOutletTemp = SST900_TEMPERATURE;
    break;
  }
  
  if ( _centralReceiverOutletTemperature < minReceiverOutletTemp || _centralReceiverOutletTemperature > 995 )
    throw inputOutOfRange;
  
  if (_hotStorageHeight < 2 || _hotStorageHeight > 50)
    throw inputOutOfRange;
  
  if (_hotStorageDiameter < 2 || _hotStorageDiameter > 30)
    throw inputOutOfRange;
  
  if (_hotStorageInsulThickness < 0.01 || _hotStorageInsulThickness > 5)
    throw inputOutOfRange;
  
  if (_coldStorageInsulThickness < 0.01 || _coldStorageInsulThickness > 5)
    throw inputOutOfRange;

  return true;
}

/*-------------------------------------------------*/
/*         validate problem maxEff_RE (#7)         */
/*-------------------------------------------------*/
bool Scenario::validate_maxEff_RE ( void ) const {

  invalid_argument inputOutOfRange("Detected unacceptable value for one or more input parameters");

  if (_receiverApertureHeight < 1 || _receiverApertureHeight > 30)
    throw inputOutOfRange;

  if (_receiverApertureWidth < 1 || _receiverApertureWidth > 30)
    throw inputOutOfRange;

  double minReceiverOutletTemp = 0.0;
  switch ( _typeOfTurbine ) {
  case 1:
    minReceiverOutletTemp = SST110_TEMPERATURE;
    break;
  case 2:
    minReceiverOutletTemp = SST120_TEMPERATURE;
    break;
  case 3:
    minReceiverOutletTemp = SST300_TEMPERATURE;
    break;
  case 4:
    minReceiverOutletTemp = SST400_TEMPERATURE;
    break;
  case 5:
    minReceiverOutletTemp = SST600_TEMPERATURE;
    break;
  case 6:
    minReceiverOutletTemp = SST700_TEMPERATURE;
    break;
  case 7:
    minReceiverOutletTemp = SST800_TEMPERATURE;
    break;
  case 8:
    minReceiverOutletTemp = SST900_TEMPERATURE;
    break;
  }
  
  if ( _centralReceiverOutletTemperature < minReceiverOutletTemp || _centralReceiverOutletTemperature > 995 )
    throw inputOutOfRange;
  
  if (_receiverNbOfTubes < 1)
    throw inputOutOfRange;
  
  if (_receiverNbOfTubes * _receiverTubesOutsideDiam > M_PI*_receiverApertureWidth / 2.)
    throw inputOutOfRange;
  
  if (_receiverInsulThickness < 0.01 || _receiverInsulThickness > 5)
    throw inputOutOfRange;
  
  if (_receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1)
    throw inputOutOfRange;
  
  if (_receiverTubesOutsideDiam < 0.0055 || _receiverTubesOutsideDiam > 0.1)
    throw inputOutOfRange;
  
  if (_receiverTubesOutsideDiam < _receiverTubesInsideDiam + 0.0005)
    throw inputOutOfRange;

   return true;
}

/*-------------------------------------------------*/
/*        validate problem maxHF_minCost (#8)      */
/*-------------------------------------------------*/
bool Scenario::validate_maxHF_minCost ( void ) const {

  invalid_argument inputOutOfRange ( "detected unacceptable value for one or more input parameters" );

  if (_heliostatHeight < 1 || _heliostatHeight > 40)
    throw inputOutOfRange;

  if (_heliostatWidth < 1 || _heliostatWidth > 40)
    throw inputOutOfRange;
  
  if (_towerHeight < 20 || _towerHeight > 250)
    throw inputOutOfRange;
  
  if (_towerHeight < 2 * _heliostatHeight)
    throw inputOutOfRange;
  
  if (_receiverApertureHeight < 1 || _receiverApertureHeight > 30)
    throw inputOutOfRange;
  
  if (_receiverApertureWidth < 1 || _receiverApertureWidth > 30)
    throw inputOutOfRange;
  
  if (_maxNumberOfHeliostats < 1)
    throw inputOutOfRange;
  
  if (_fieldMaxAngularSpan < 1 || _fieldMaxAngularSpan > 89)
    throw inputOutOfRange;
  
  if (_minimumDistanceToTower < 0 || _minimumDistanceToTower > 20)
    throw inputOutOfRange;
  
  if (_maximumDistanceToTower < 1 || _maximumDistanceToTower > 20)
    throw inputOutOfRange;
  
  if (_maximumDistanceToTower <= _minimumDistanceToTower)
    throw inputOutOfRange;
  
  if ( M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
       *(2 * _fieldMaxAngularSpan / 360.) > _cFieldSurface )
    throw inputOutOfRange;
  
  if (_receiverNbOfTubes < 1)
    throw inputOutOfRange;
  
  if (_receiverInsulThickness < 0.01 || _receiverInsulThickness > 5)
    throw inputOutOfRange;
  
  if (_receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1)
    throw inputOutOfRange;
  
  if (_receiverTubesOutsideDiam < 0.006 || _receiverTubesOutsideDiam > 0.1)
    throw inputOutOfRange;
  
  if (_receiverTubesOutsideDiam < _receiverTubesInsideDiam + 0.0005)
    throw inputOutOfRange;

   return true;
}

/*-------------------------------------------------*/
/*        validate problem maxNrg_minPar (#9)      */
/*-------------------------------------------------*/
bool Scenario::validate_maxNrg_minPar ( void ) const {

  invalid_argument inputOutOfRange ( "detected unacceptable value for one or more input parameters" );

  if ( _heliostatHeight < 1 || _heliostatHeight > 40 )
    throw inputOutOfRange;

  if ( _heliostatWidth  < 1 || _heliostatWidth > 40 )
    throw inputOutOfRange;

  if ( _towerHeight < 20 || _towerHeight > 250 )
    throw inputOutOfRange;

  if ( _towerHeight < 2 * _heliostatHeight )
    throw inputOutOfRange;
  
  if ( _receiverApertureHeight < 1 || _receiverApertureHeight > 30 )
    throw inputOutOfRange;
  
  if ( _receiverApertureWidth  < 1 || _receiverApertureWidth  > 30 )
    throw inputOutOfRange;
  
  if ( _fieldMaxAngularSpan < 1 || _fieldMaxAngularSpan > 89 )
    throw inputOutOfRange;

  if ( _minimumDistanceToTower < 0 || _minimumDistanceToTower > 20 )
    throw inputOutOfRange;

  if ( _maximumDistanceToTower < 1 || _maximumDistanceToTower > 20 )
    throw inputOutOfRange; 

  if ( _maximumDistanceToTower <= _minimumDistanceToTower )
    throw inputOutOfRange;

  if ( _maxNumberOfHeliostats < 1 )
    throw inputOutOfRange;
  
  if ( M_PI*(pow(_maximumDistanceToTower*_towerHeight, 2.) - pow(_minimumDistanceToTower*_towerHeight, 2.))
       *(2 * _fieldMaxAngularSpan / 360.) > _cFieldSurface )
    throw inputOutOfRange;

  double minReceiverOutletTemp = 0.;
  switch (_typeOfTurbine) {
  case 1:
    minReceiverOutletTemp = SST110_TEMPERATURE;
    break;
  case 2:
    minReceiverOutletTemp = SST120_TEMPERATURE;
    break;
  case 3:
    minReceiverOutletTemp = SST300_TEMPERATURE;
    break;
  case 4:
    minReceiverOutletTemp = SST400_TEMPERATURE;
    break;
  case 5:
    minReceiverOutletTemp = SST600_TEMPERATURE;
    break;
  case 6:
    minReceiverOutletTemp = SST700_TEMPERATURE;
    break;
  case 7:
    minReceiverOutletTemp = SST800_TEMPERATURE;
    break;
  case 8:
    minReceiverOutletTemp = SST900_TEMPERATURE;
    break;
  }
  
  if ( _centralReceiverOutletTemperature < minReceiverOutletTemp || _centralReceiverOutletTemperature > 995 )
    throw inputOutOfRange;

  if (_receiverNbOfTubes * _receiverTubesOutsideDiam > M_PI*_receiverApertureWidth / 2.)
    throw inputOutOfRange;

  if (_hotStorageHeight < 1 || _hotStorageHeight > 50)
    throw inputOutOfRange;

  if (_hotStorageDiameter < 1 || _hotStorageDiameter > 30)
    throw inputOutOfRange;

  if (_hotStorageInsulThickness < 0.01 || _hotStorageInsulThickness > 5)
    throw inputOutOfRange;

  if (_coldStorageInsulThickness < 0.01 || _coldStorageInsulThickness > 5)
    throw inputOutOfRange;

  if (_coldMoltenSaltMinTemperature < MELTING_POINT || _coldMoltenSaltMinTemperature > 650)
    throw inputOutOfRange;

  if (_receiverNbOfTubes < 1)
    throw inputOutOfRange;

  if (_receiverInsulThickness < 0.01 || _receiverInsulThickness > 5)
    throw inputOutOfRange;

  if (_receiverTubesInsideDiam < 0.005 || _receiverTubesInsideDiam > 0.1)
    throw inputOutOfRange;

  if (_receiverTubesOutsideDiam < 0.006 || _receiverTubesOutsideDiam > 0.1)
    throw inputOutOfRange;
  
  if (_receiverTubesOutsideDiam < _receiverTubesInsideDiam + 0.0001)
    throw inputOutOfRange;
  
  if (_exchangerTubesSpacing <= _exchangerTubesDout || _exchangerTubesSpacing < 0.007 || _exchangerTubesSpacing > 0.2)
    throw inputOutOfRange;

  if (_exchangerTubesLength < 0.5 || _receiverTubesOutsideDiam > 10)
    throw inputOutOfRange;

  if (_exchangerTubesDin < 0.005 || _exchangerTubesDin > 0.1)
    throw inputOutOfRange;
  
  if (_exchangerTubesDout < 0.006 || _exchangerTubesDout > 0.1)
    throw inputOutOfRange;

  if (_exchangerTubesDout < _exchangerTubesDin + 0.0005)
    throw inputOutOfRange;
  
  if (_exchangerNbOfBaffles < 2)
    throw inputOutOfRange;
  
  if (_exchangerBaffleCut < 0.15 || _exchangerBaffleCut > 0.4)
    throw inputOutOfRange;

  if (_exchangerNbOfTubes < 1)
    throw inputOutOfRange;

  if (_exchangerNbOfShells < 1 || _exchangerNbOfShells > 10)
    throw inputOutOfRange;
  
  if (_exchangerNbOfPassesPerShell < 1 || _exchangerNbOfPassesPerShell > 9)
    throw inputOutOfRange;

  if (_typeOfTurbine < 1 || _typeOfTurbine >  8)
    throw inputOutOfRange;

   return true;
}

/*****************************************************************************/
/*****************************************************************************/
/************************constructing model components************************/
/*****************************************************************************/
/*****************************************************************************/
void Scenario::construct_maxNrg_H1 ( void ) {

  _time = new Clock ( _numberOfTimeIncrements , 0 , _minutesPerTimeIncrement );
  _sun  = new Sun   ( _latitude , *_time , _day , _raysPerSquareMeters       );

  HeliostatField * field;
  Economics      * economics;

  field = new HeliostatField ( _maxNumberOfHeliostats  ,
			       _heliostatHeight        ,
			       _heliostatWidth         ,
			       _towerHeight            ,
			       _receiverApertureHeight ,
			       _receiverApertureWidth  ,
			       _minimumDistanceToTower ,
			       _maximumDistanceToTower ,
			       _fieldMaxAngularSpan    ,
			       *_sun                     );

  economics = new Economics();
  economics->set_heightOfTower(_towerHeight);
  economics->set_heightOfReceiverAperture(_receiverApertureHeight);
  economics->set_widthOfReceiverAperture(_receiverApertureWidth);
  economics->set_lengthOfHeliostats(_heliostatHeight);
  economics->set_widthOfHeliostats(_heliostatWidth);
  economics->set_exchangerModel ( _exchangerModel);

  _powerplant = new Powerplant ( *_time      ,
				 *_sun       ,
				 _model_type ,
				 field       ,
				 NULL        ,
				 NULL        ,
				 economics     );

  _powerplant->set_heliostatModel ( _heliostatsFieldModel );

}

void Scenario::construct_minCost_C1()
{
	_model_type = 2; // whole plant

	_time = new Clock(_numberOfTimeIncrements, 0, _minutesPerTimeIncrement);
	_sun = new Sun(_latitude, *_time, _day, _raysPerSquareMeters);
	fFillDemandVector();

	HeliostatField* field;
	HtfCycle* htfCycle;
	Powerblock* powerblock;
	Economics* economics;

	//Powerblock model
	powerblock = new Powerblock(_typeOfTurbine);

	//Constructing heliostats field
	field = new HeliostatField(_maxNumberOfHeliostats,
		_heliostatHeight,
		_heliostatWidth,
		_towerHeight,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_minimumDistanceToTower,
		_maximumDistanceToTower,
		_fieldMaxAngularSpan,
		*_sun);

	//Constructing Htf Cycle with desired steam generator model
	htfCycle = new HtfCycle(
		_centralReceiverOutletTemperature,
		_hotStorageHeight,
		_hotStorageDiameter,
		_hotStorageInsulThickness,
		_coldMoltenSaltMinTemperature,
		powerblock,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_receiverTubesInsideDiam,
		_receiverTubesOutsideDiam,
		_receiverNbOfTubes,
		_fixedPointsPrecision,
		_minutesPerTimeIncrement
		);
	htfCycle->initiateColdStorage();
	/*htfCycle->setStorage(_storageStartupCondition,
		0.98*_centralReceiverOutletTemperature,
		_coldMoltenSaltMinTemperature);*/

	//Investment cost model
	economics = new Economics();
	economics->set_heightOfTower(_towerHeight);
	economics->set_heightOfReceiverAperture(_receiverApertureHeight);
	economics->set_widthOfReceiverAperture(_receiverApertureWidth);
	economics->set_receiverNumberOfTubes(_receiverNbOfTubes);
	economics->set_receiverTubesDout(_receiverTubesOutsideDiam);
	economics->set_lengthOfHeliostats(_heliostatHeight);
	economics->set_widthOfHeliostats(_heliostatWidth);
	economics->set_hotStorageHeight(_hotStorageHeight);
	economics->set_storageDiameter(_hotStorageDiameter);
	economics->set_hotStorageInsulationThickness(_hotStorageInsulThickness);
	economics->set_coldStorageInsulationThickness(_coldStorageInsulThickness);
	economics->set_receiverInsulationThickness(_receiverInsulThickness);
	economics->set_turbineNominalPowerOutput(powerblock->get_powerOfTurbine());
	economics->set_exchangerModel ( _exchangerModel);

	_powerplant = new Powerplant ( *_time      ,
				       *_sun       ,
				       _model_type ,
				       field       ,
				       htfCycle    ,
				       powerblock  ,
				       economics     );
	
	//Setting demand vector
	_powerplant->set_demand              ( _demand               );
	_powerplant->set_heliostatModel      ( _heliostatsFieldModel );
}

void Scenario::construct_minCost_C2()
{
	_model_type = 2; // whole plant
	
	_time = new Clock(_numberOfTimeIncrements, 0, _minutesPerTimeIncrement);
	_sun = new Sun(_latitude, *_time, _day, _raysPerSquareMeters);
	fFillDemandVector();
	
	HeliostatField* field;
	HtfCycle* htfCycle;
	Powerblock* powerblock;
	Economics* economics;

	//Powerblock model
	powerblock = new Powerblock(_typeOfTurbine);

	//Constructing heliostats field
	field = new HeliostatField(_maxNumberOfHeliostats,
		_heliostatHeight,
		_heliostatWidth,
		_towerHeight,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_minimumDistanceToTower,
		_maximumDistanceToTower,
		_fieldMaxAngularSpan,
		*_sun);

	// TOTO: des new sont faits ici: sont-ils bien deletes quelque part
	
	
	//Constructing Htf Cycle with desired steam generator model
	htfCycle = new HtfCycle(
		_centralReceiverOutletTemperature,
		_hotStorageHeight,
		_hotStorageDiameter,
		_hotStorageInsulThickness,
		_coldMoltenSaltMinTemperature,
		powerblock,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_receiverTubesInsideDiam,
		_receiverTubesOutsideDiam,
		_receiverNbOfTubes,
		_fixedPointsPrecision,
		_minutesPerTimeIncrement,
		_exchangerTubesDin, 
		_exchangerTubesDout, 
		_exchangerTubesLength,
		_exchangerTubesSpacing,
		_exchangerBaffleCut,
		_exchangerNbOfBaffles,
		_exchangerNbOfTubes,
		_exchangerNbOfPassesPerShell,
		_exchangerNbOfShells);

	htfCycle->initiateColdStorage();
    
	//Investment cost model
	economics = new Economics();
	economics->set_heightOfTower(_towerHeight);
	economics->set_heightOfReceiverAperture(_receiverApertureHeight);
	economics->set_widthOfReceiverAperture(_receiverApertureWidth);
	economics->set_receiverNumberOfTubes(_receiverNbOfTubes);
	economics->set_receiverTubesDout(_receiverTubesOutsideDiam);
	economics->set_lengthOfHeliostats(_heliostatHeight);
	economics->set_widthOfHeliostats(_heliostatWidth);
	economics->set_hotStorageHeight(_hotStorageHeight);
	economics->set_storageDiameter(_hotStorageDiameter);
	economics->set_hotStorageInsulationThickness(_hotStorageInsulThickness);
	economics->set_coldStorageInsulationThickness(_coldStorageInsulThickness);
	economics->set_receiverInsulationThickness(_receiverInsulThickness);
	economics->set_turbineNominalPowerOutput(powerblock->get_powerOfTurbine());
	economics->set_exchangerModel ( _exchangerModel );



	
	_powerplant = new Powerplant ( *_time      ,
				       *_sun       ,
				       _model_type ,
				       field       ,
				       htfCycle    ,
				       powerblock  ,
				       economics     );

	//Setting demand vector
	_powerplant->set_demand(_demand);
	_powerplant->set_heliostatModel(_heliostatsFieldModel);
}

void Scenario::construct_minSurf_H1()
{
	_model_type = 2; // whole plant

	_time = new Clock(_numberOfTimeIncrements, 0, _minutesPerTimeIncrement);
	_sun = new Sun(_latitude, *_time, _day, _raysPerSquareMeters);
	fFillDemandVector();

	HeliostatField* field;
	HtfCycle* htfCycle;
	Powerblock* powerblock;
	Economics* economics;

	//Powerblock model
	powerblock = new Powerblock(_typeOfTurbine);

	//Constructing heliostats field
	field = new HeliostatField(_maxNumberOfHeliostats,
		_heliostatHeight,
		_heliostatWidth,
		_towerHeight,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_minimumDistanceToTower,
		_maximumDistanceToTower,
		_fieldMaxAngularSpan,
		*_sun);

	//Constructing Htf Cycle with desired steam generator model
	htfCycle = new HtfCycle(
		_centralReceiverOutletTemperature,
		_hotStorageHeight,
		_hotStorageDiameter,
		_hotStorageInsulThickness,
		_coldMoltenSaltMinTemperature,
		powerblock,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_receiverTubesInsideDiam,
		_receiverTubesOutsideDiam,
		_receiverNbOfTubes,
		_fixedPointsPrecision,
		_minutesPerTimeIncrement
		);

	htfCycle->setStorage(_storageStartupCondition,
		0.98*_centralReceiverOutletTemperature,
		_coldMoltenSaltMinTemperature);

	//Investment cost model
	economics = new Economics();
	economics->set_heightOfTower(_towerHeight);
	economics->set_heightOfReceiverAperture(_receiverApertureHeight);
	economics->set_widthOfReceiverAperture(_receiverApertureWidth);
	economics->set_receiverNumberOfTubes(_receiverNbOfTubes);
	economics->set_receiverTubesDout(_receiverTubesOutsideDiam);
	economics->set_lengthOfHeliostats(_heliostatHeight);
	economics->set_widthOfHeliostats(_heliostatWidth);
	economics->set_hotStorageHeight(_hotStorageHeight);
	economics->set_storageDiameter(_hotStorageDiameter);
	economics->set_hotStorageInsulationThickness(_hotStorageInsulThickness);
	economics->set_coldStorageInsulationThickness(_coldStorageInsulThickness);
	economics->set_receiverInsulationThickness(_receiverInsulThickness);
	economics->set_turbineNominalPowerOutput(powerblock->get_powerOfTurbine());
	economics->set_exchangerModel ( _exchangerModel);

	_powerplant = new Powerplant ( *_time      ,
				       *_sun       ,
				       _model_type ,
				       field       ,
				       htfCycle    ,
				       powerblock  ,
				       economics     );

	//Setting demand vector
	_powerplant->set_demand(_demand);
	_powerplant->set_heliostatModel(_heliostatsFieldModel);
}

void Scenario::construct_maxComp_HTF1 ( void ) {
  
	_model_type = 2; // whole plant

	_time = new Clock ( _numberOfTimeIncrements, 0, _minutesPerTimeIncrement );
	_sun  = new Sun   ( _latitude, *_time, _day, _raysPerSquareMeters        );
	fFillDemandVector();

	HeliostatField* field;
	HtfCycle      * htfCycle;
	Powerblock    * powerblock;
	Economics     * economics;

	//Powerblock model
	powerblock = new Powerblock(_typeOfTurbine);

	//Constructing heliostats field
	field = new HeliostatField(_maxNumberOfHeliostats,
	_heliostatHeight,
	_heliostatWidth,
	_towerHeight,
	_receiverApertureHeight,
	_receiverApertureWidth,
	_minimumDistanceToTower,
	_maximumDistanceToTower,
	_fieldMaxAngularSpan,
	*_sun);

	//Constructing Htf Cycle with desired steam generator model
	htfCycle = new HtfCycle(
		_centralReceiverOutletTemperature,
		_hotStorageHeight,
		_hotStorageDiameter,
		_hotStorageInsulThickness,
		_coldMoltenSaltMinTemperature,
		powerblock,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_receiverTubesInsideDiam,
		_receiverTubesOutsideDiam,
		_receiverNbOfTubes,
		_fixedPointsPrecision,
		_minutesPerTimeIncrement,
		_exchangerTubesDin,
		_exchangerTubesDout,
		_exchangerTubesLength,
		_exchangerTubesSpacing,
		_exchangerBaffleCut,
		_exchangerNbOfBaffles,
		_exchangerNbOfTubes,
		_exchangerNbOfPassesPerShell,
		_exchangerNbOfShells
		);

	htfCycle->initiateColdStorage();

	//Investment cost model
	economics = new Economics();
	economics->set_heightOfTower(_towerHeight);
	economics->set_heightOfReceiverAperture(_receiverApertureHeight);
	economics->set_widthOfReceiverAperture(_receiverApertureWidth);
	economics->set_receiverNumberOfTubes(_receiverNbOfTubes);
	economics->set_receiverTubesDout(_receiverTubesOutsideDiam);
	economics->set_lengthOfHeliostats(_heliostatHeight);
	economics->set_widthOfHeliostats(_heliostatWidth);
	economics->set_hotStorageHeight(_hotStorageHeight);
	economics->set_storageDiameter(_hotStorageDiameter);
	economics->set_hotStorageInsulationThickness(_hotStorageInsulThickness);
	economics->set_coldStorageInsulationThickness(_coldStorageInsulThickness);
	economics->set_receiverInsulationThickness(_receiverInsulThickness);
	economics->set_turbineNominalPowerOutput(powerblock->get_powerOfTurbine());
	economics->set_exchangerModel ( _exchangerModel );

	_powerplant = new Powerplant ( *_time      ,
				       *_sun       ,
				       _model_type ,
				       field       ,
				       htfCycle    ,
				       powerblock  ,
				       economics     );

	//Setting demand vector
	_powerplant->set_demand(_demand);
	_powerplant->set_heliostatModel(_heliostatsFieldModel);
		
	if  ( !_powerplant->set_heliostatFieldPowerOutput_MAXCOMP_HTF1() )
	  throw logic_error ( "error in the construction of the problem" );
}

void Scenario::construct_minCost_TS()
{
	_model_type = 2; // whole plant

	_time = new Clock(_numberOfTimeIncrements, 0, _minutesPerTimeIncrement);
	_sun = new Sun(_latitude, *_time, _day, _raysPerSquareMeters);
	fFillDemandVector();

	HeliostatField* field;
	HtfCycle* htfCycle;
	Powerblock* powerblock;
	Economics* economics;

	//Powerblock model
	powerblock = new Powerblock(_typeOfTurbine);

	//Constructing heliostats field
	field = NULL;

	//Constructing Htf Cycle with desired steam generator model
	htfCycle = new HtfCycle(
		_centralReceiverOutletTemperature,
		_hotStorageHeight,
		_hotStorageDiameter,
		_hotStorageInsulThickness,
		_coldMoltenSaltMinTemperature,
		powerblock,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_receiverTubesInsideDiam,
		_receiverTubesOutsideDiam,
		_receiverNbOfTubes,
		_fixedPointsPrecision,
		_minutesPerTimeIncrement
		);
	/*htfCycle->initiateColdStorage();*/
	htfCycle->setStorage(_storageStartupCondition,
		0.98*_centralReceiverOutletTemperature,
		_coldMoltenSaltMinTemperature);

	//Investment cost model
	economics = new Economics();
	economics->set_heightOfTower(_towerHeight);
	economics->set_heightOfReceiverAperture(_receiverApertureHeight);
	economics->set_widthOfReceiverAperture(_receiverApertureWidth);
	economics->set_receiverNumberOfTubes(_receiverNbOfTubes);
	economics->set_receiverTubesDout(_receiverTubesOutsideDiam);
	economics->set_lengthOfHeliostats(_heliostatHeight);
	economics->set_widthOfHeliostats(_heliostatWidth);
	economics->set_hotStorageHeight(_hotStorageHeight);
	economics->set_storageDiameter(_hotStorageDiameter);
	economics->set_hotStorageInsulationThickness(_hotStorageInsulThickness);
	economics->set_coldStorageInsulationThickness(_coldStorageInsulThickness);
	economics->set_receiverInsulationThickness(_receiverInsulThickness);
	economics->set_turbineNominalPowerOutput(powerblock->get_powerOfTurbine());
	economics->set_exchangerModel ( _exchangerModel );

	_powerplant = new Powerplant ( *_time      ,
				       *_sun       ,
				       _model_type ,
				       field       ,
				       htfCycle    ,
				       powerblock  ,
				       economics     );

	//Setting demand vector
	_powerplant->set_demand(_demand);
	_powerplant->set_heliostatModel(_heliostatsFieldModel);

	if  ( !_powerplant->set_heliostatFieldPowerOutput_MINCOST_TS() )
	  throw logic_error ( "error in the construction of the problem" );
}

void Scenario::construct_maxEff_RE()
{
	_model_type = 2; // whole plant

	_time = new Clock(_numberOfTimeIncrements, 0, _minutesPerTimeIncrement);
	_sun = new Sun(_latitude, *_time, _day, _raysPerSquareMeters);
	fFillDemandVector();

	HeliostatField* field;
	HtfCycle* htfCycle;
	Powerblock* powerblock;
	Economics* economics;

	//Powerblock model
	powerblock = new Powerblock(_typeOfTurbine);

	//Constructing heliostats field
	field = new HeliostatField(_maxNumberOfHeliostats,
		_heliostatHeight,
		_heliostatWidth,
		_towerHeight,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_minimumDistanceToTower,
		_maximumDistanceToTower,
		_fieldMaxAngularSpan,
		*_sun);

	//Constructing Htf Cycle with desired steam generator model
	htfCycle = new HtfCycle(
		_centralReceiverOutletTemperature,
		_hotStorageHeight,
		_hotStorageDiameter,
		_hotStorageInsulThickness,
		_coldMoltenSaltMinTemperature,
		powerblock,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_receiverTubesInsideDiam,
		_receiverTubesOutsideDiam,
		_receiverNbOfTubes,
		_fixedPointsPrecision,
		_minutesPerTimeIncrement
		);
	htfCycle->initiateColdStorage();

	//Investment cost model
	economics = new Economics();
	economics->set_heightOfTower(_towerHeight);
	economics->set_heightOfReceiverAperture(_receiverApertureHeight);
	economics->set_widthOfReceiverAperture(_receiverApertureWidth);
	economics->set_receiverNumberOfTubes(_receiverNbOfTubes);
	economics->set_receiverTubesDout(_receiverTubesOutsideDiam);
	economics->set_lengthOfHeliostats(_heliostatHeight);
	economics->set_widthOfHeliostats(_heliostatWidth);
	economics->set_hotStorageHeight(_hotStorageHeight);
	economics->set_storageDiameter(_hotStorageDiameter);
	economics->set_hotStorageInsulationThickness(_hotStorageInsulThickness);
	economics->set_coldStorageInsulationThickness(_coldStorageInsulThickness);
	economics->set_receiverInsulationThickness(_receiverInsulThickness);
	economics->set_turbineNominalPowerOutput(powerblock->get_powerOfTurbine());
	economics->set_exchangerModel ( _exchangerModel );

	_powerplant = new Powerplant ( *_time      ,
				       *_sun       ,
				       _model_type ,
				       field       ,
				       htfCycle    ,
				       powerblock  ,
				       economics     );

	//Setting demand vector
	_powerplant->set_demand(_demand);
	_powerplant->set_heliostatModel(_heliostatsFieldModel);
}

void Scenario::construct_maxHF_minCost()
{
	_model_type = 2; // whole plant

	_time = new Clock(_numberOfTimeIncrements, 0, _minutesPerTimeIncrement);
	_sun = new Sun(_latitude, *_time, _day, _raysPerSquareMeters);
	fFillDemandVector();

	HeliostatField* field;
	HtfCycle* htfCycle;
	Powerblock* powerblock;
	Economics* economics;

	//Powerblock model
	powerblock = new Powerblock(_typeOfTurbine);

	//Constructing heliostats field
	field = new HeliostatField(_maxNumberOfHeliostats,
		_heliostatHeight,
		_heliostatWidth,
		_towerHeight,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_minimumDistanceToTower,
		_maximumDistanceToTower,
		_fieldMaxAngularSpan,
		*_sun);

	//Constructing Htf Cycle with desired steam generator model
	htfCycle = new HtfCycle(
		_centralReceiverOutletTemperature,
		_hotStorageHeight,
		_hotStorageDiameter,
		_hotStorageInsulThickness,
		_coldMoltenSaltMinTemperature,
		powerblock,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_receiverTubesInsideDiam,
		_receiverTubesOutsideDiam,
		_receiverNbOfTubes,
		_fixedPointsPrecision,
		_minutesPerTimeIncrement
		);
	htfCycle->initiateColdStorage();

	//Investment cost model
	economics = new Economics();
	economics->set_heightOfTower(_towerHeight);
	economics->set_heightOfReceiverAperture(_receiverApertureHeight);
	economics->set_widthOfReceiverAperture(_receiverApertureWidth);
	economics->set_receiverNumberOfTubes(_receiverNbOfTubes);
	economics->set_receiverTubesDout(_receiverTubesOutsideDiam);
	economics->set_lengthOfHeliostats(_heliostatHeight);
	economics->set_widthOfHeliostats(_heliostatWidth);
	economics->set_hotStorageHeight(_hotStorageHeight);
	economics->set_storageDiameter(_hotStorageDiameter);
	economics->set_hotStorageInsulationThickness(_hotStorageInsulThickness);
	economics->set_coldStorageInsulationThickness(_coldStorageInsulThickness);
	economics->set_receiverInsulationThickness(_receiverInsulThickness);
	economics->set_turbineNominalPowerOutput(powerblock->get_powerOfTurbine());
	economics->set_exchangerModel ( _exchangerModel );

	_powerplant = new Powerplant ( *_time      ,
				       *_sun       ,
				       _model_type ,
				       field       ,
				       htfCycle    ,
				       powerblock  ,
				       economics     );

	//Setting demand vector
	_powerplant->set_demand(_demand);
	_powerplant->set_heliostatModel(_heliostatsFieldModel);
}

void Scenario::construct_maxNrg_minPar()
{
	_model_type = 2; // whole plant

	_time = new Clock(_numberOfTimeIncrements, 0, _minutesPerTimeIncrement);
	_sun = new Sun(_latitude, *_time, _day, _raysPerSquareMeters);
	fFillDemandVector();

	HeliostatField* field;
	HtfCycle* htfCycle;
	Powerblock* powerblock;
	Economics* economics;

	//Powerblock model
	powerblock = new Powerblock(_typeOfTurbine);

	//Constructing heliostats field
	field = new HeliostatField(_maxNumberOfHeliostats,
		_heliostatHeight,
		_heliostatWidth,
		_towerHeight,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_minimumDistanceToTower,
		_maximumDistanceToTower,
		_fieldMaxAngularSpan,
		*_sun);

	//Constructing Htf Cycle with desired steam generator model
	htfCycle = new HtfCycle(
		_centralReceiverOutletTemperature,
		_hotStorageHeight,
		_hotStorageDiameter,
		_hotStorageInsulThickness,
		_coldMoltenSaltMinTemperature,
		powerblock,
		_receiverApertureHeight,
		_receiverApertureWidth,
		_receiverTubesInsideDiam,
		_receiverTubesOutsideDiam,
		_receiverNbOfTubes,
		_fixedPointsPrecision,
		_minutesPerTimeIncrement,
		_exchangerTubesDin,
		_exchangerTubesDout,
		_exchangerTubesLength,
		_exchangerTubesSpacing,
		_exchangerBaffleCut,
		_exchangerNbOfBaffles,
		_exchangerNbOfTubes,
		_exchangerNbOfPassesPerShell,
		_exchangerNbOfShells);

	htfCycle->initiateColdStorage();

	//Investment cost model
	economics = new Economics();
	economics->set_heightOfTower(_towerHeight);
	economics->set_heightOfReceiverAperture(_receiverApertureHeight);
	economics->set_widthOfReceiverAperture(_receiverApertureWidth);
	economics->set_receiverNumberOfTubes(_receiverNbOfTubes);
	economics->set_receiverTubesDout(_receiverTubesOutsideDiam);
	economics->set_lengthOfHeliostats(_heliostatHeight);
	economics->set_widthOfHeliostats(_heliostatWidth);
	economics->set_hotStorageHeight(_hotStorageHeight);
	economics->set_storageDiameter(_hotStorageDiameter);
	economics->set_hotStorageInsulationThickness(_hotStorageInsulThickness);
	economics->set_coldStorageInsulationThickness(_coldStorageInsulThickness);
	economics->set_receiverInsulationThickness(_receiverInsulThickness);
	economics->set_turbineNominalPowerOutput(powerblock->get_powerOfTurbine());
	economics->set_exchangerModel ( _exchangerModel );

	_powerplant = new Powerplant ( *_time      ,
				       *_sun       ,
				       _model_type ,
				       field       ,
				       htfCycle    ,
				       powerblock  ,
				       economics     );

	//Setting demand vector
	_powerplant->set_demand(_demand);
	_powerplant->set_heliostatModel(_heliostatsFieldModel);
}

/*PRINT functions */
/*These functions will generate data files to analyze the solution's behavior*/

// TOTO: mettre a jour les fonctions suivantes pour un affichage final en mode verbose

//------------------------------------------------------------------------
//Function name : printColdStorage
//File name : DATA_ColdStorage.txt
//Will print in this format:
//Level (meters)   Temperature (K)
//------------------------------------------------------------------------
// void Scenario::printColdStorage ( void )
// {
// 	std::string output = "DATA_ColdStorage.txt";
// 	std::vector<double> level = _powerplant->get_coldStorageLevel();
// 	std::vector<double> temp = _powerplant->get_coldStorageTemp();

// 	ofstream coldStorage;
// 	coldStorage.open(output.c_str());
// 	coldStorage << "Level(m)" << "Temperature (K)" << endl;

// 	for (unsigned int i = 0; i < level.size(); i+=30)
// 	{
// 		coldStorage << level[i] << " " << temp[i] << endl;
// 	}
	
// 	coldStorage.close();
// }

//------------------------------------------------------------------------
//Function name : printHotStorage
//File name : DATA_HotStorage.txt
//Will print in this format:
//Level (meters)   Temperature (K)
//------------------------------------------------------------------------
// void Scenario::printHotStorage ( void )
// {
// 	std::string output = "DATA_HotStorage.txt";
// 	std::vector<double> level = _powerplant->get_hotStorageLevel();
// 	std::vector<double> temp = _powerplant->get_hotStorageTemp();

// 	ofstream hotStorage;
// 	hotStorage.open(output.c_str());
// 	hotStorage << "Level(m)" << "Temperature (K)" << endl;

// 	for (unsigned int i = 0; i < level.size(); i += 30)
// 	{
// 		hotStorage << level[i] << " " << temp[i] << endl;
// 	}

// 	hotStorage.close();
// }

//------------------------------------------------------------------------
//Function Name : printSteamGeneratorOutlet
//File Name : DATA_SteamGeneratorOutlet.txt
//Format:
//Mass flow (kg/s) Temperature(K) Heat Transfer Rate(W)
//------------------------------------------------------------------------
// void Scenario::printSteamGeneratorOutlet ( void )
// {
// 	std::string output = "DATA_SteamGeneratorOutlet.txt";
// 	std::vector<double> msRate = _powerplant->get_steamGeneratorInletFlow();
// 	std::vector<double> temp = _powerplant->get_steamGeneratorOutletTemperature();
// 	std::vector<double> HTrate = _powerplant->get_energyToPowerBlockWatts();

// 	ofstream steamGenOutlet;
// 	steamGenOutlet.open(output.c_str());
// 	steamGenOutlet << "Ms Rate (kg/s)" << " -Temperature (K)" << " -Heat Transfer Rate (W)" <<endl;

// 	for (unsigned int i = 0; i < msRate.size(); i += 30)
// 	{
// 		steamGenOutlet << msRate[i] << " " << temp[i] << " " << HTrate[i] << endl;
// 	}

// 	steamGenOutlet.close();
// }

//------------------------------------------------------------------------
//Function Name : printReceiver
//File Name : DATA_Receiver.txt
//Format:
//Mass flow (kg/s) SurfaceTemperature (K) HeatLosses (kW) Efficiency (0-1)
//------------------------------------------------------------------------
// void Scenario::printReceiver ( void )
// {
// 	std::string output = "DATA_Receiver.txt";
// 	std::vector<double> msRate = _powerplant->get_receiverMsRate();
// 	std::vector<double> temp = _powerplant->get_receiverSurfaceTemperature();
// 	std::vector<double> losses = _powerplant->get_receiverLosses();
// 	std::vector<double> eff = _powerplant->get_receiverEfficiency();

// 	ofstream Receiver;
// 	Receiver.open(output.c_str());

// 	Receiver << "Mass Flow(kg / s) " << "Surface Temperature (K) " << "Heat Losses (kW) " << "Efficiency" << std::endl;

// 	for (unsigned int i = 0; i < msRate.size(); i += 30)
// 	{
// 		Receiver << msRate[i] << " " << temp[i] << " " << losses[i] << " " << eff[i] << endl;
// 	}

// 	Receiver.close();
// }

//------------------------------------------------------------------------
//Function Name : printHeliostatsFieldLayout();
//File Name : DATA_HeliostatsLayout.txt
//Format:
//X(m) Y(m) Z(m)
//------------------------------------------------------------------------
// void Scenario::printHeliostatsFieldLayout ( void )
// {
// 	std::string output = "DATA_HeliostatsLayout.txt";

// 	_powerplant->get_heliostatField()->fOutputFieldLayout(output);

// 	output = "DATA_gridLayout.txt";
// 	_powerplant->get_heliostatField()->fOutputGridLayout(output);
// }

//------------------------------------------------------------------------
//Function Name : printHeliostatsFieldPerformance();
//File Name : DATA_HeliostatsPerformance.txt
//Format:
//Power output (W)
//------------------------------------------------------------------------
// void Scenario::printHeliostatsFieldPerformance ( void )
// {
// 	std::string output = "DATA_HeliostatsPerformance.txt";

// 	std::vector<double> powerOutput = _powerplant->get_heliostatFieldPowerOutput();

// 	ofstream powerFromField;
// 	foncteurPrintPower printPower(&powerFromField);

// 	powerFromField.open(output.c_str());

//     powerFromField << "Power To Receiver Aperture (W)" << std::endl;

// 	std::for_each(powerOutput.begin(), powerOutput.end(), printPower);

// 	powerFromField.close();
// }

//------------------------------------------------------------------------
//Function Name : printPowerOutput();
//File Name : DATA_PowerOutput.txt
//Format:
//Power output (W)
//------------------------------------------------------------------------
// void Scenario::printPowerOutput ( void )
// {
// 	std::string output = "DATA_PowerOutput.txt";

// 	std::vector<double> powerOutput = _powerplant->get_powerplantPowerOutput();

// 	ofstream ePower;
// 	ePower.open(output.c_str());

// 	ePower << "Electrical Power (W)" << endl;

// 	for (unsigned int i = 0; i < powerOutput.size(); i += 30)
// 	{
// 		ePower << powerOutput[i] << std::endl;
// 	}

// 	ePower.close();
// }

//------------------------------------------------------------------------
//Function Name : printDemand();
//File Name : DATA_Demand.txt
//Format:
//Power output (W)
//------------------------------------------------------------------------
// void Scenario::printDemand ( void )
// {
// 	std::string output = "DATA_Demand.txt";

// 	ofstream ePower;
// 	ePower.open(output.c_str());

// 	ePower << "Electrical Power (W)" << endl;

// 	for (unsigned int i = 0; i < _demand.size(); ++i)
// 	{
// 		ePower << _demand[i] << std::endl;
// 	}

// 	ePower.close();
// }

//------------------------------------------------------------------------
//Function Name : printRequiredThermalPower();
//File Name : DATA_Demand.txt
//Format:
//Power output (W)
//------------------------------------------------------------------------
// void Scenario::printRequiredThermalPower ( void )
// {
// 	std::string output = "DATA_RequiredThermal.txt";

// 	ofstream ePower;

// 	std::vector<double> requiredThermal = _powerplant->get_thermalPowerNeededFromBlock();
// 	ePower.open(output.c_str());

// 	ePower << "Thermal Power (W)" << std::endl;

// 	for (unsigned int i = 0; i < requiredThermal.size(); i += 30)
// 	{
// 		ePower << requiredThermal[i] << std::endl;
// 	}

// 	ePower.close();
// }

//------------------------------------------------------------------------
//Function Name : print();
//Calls all print functions
//------------------------------------------------------------------------
// TOTO: Mettre a jour cette fonction pour affichage final en mode verbose
void Scenario::print ( void ) const {
	// if (_problem == "MAXNRG_H1" || _problem == "MAXNRG_H1S")
	// {
  	//         printHeliostatsFieldLayout();
	// 	printHeliostatsFieldPerformance();
	// }
	// else if (_problem == "MAXCOMP_HTF1" || _problem == "MAXCOMP_HTF1S"
	// 	||_problem == "MINCOST_TS" || _problem == "MINCOST_TSS")
	// {
	// 	printColdStorage();
	// 	printHotStorage();
	// 	printSteamGeneratorOutlet();
	// 	printRequiredThermalPower();
	// 	printReceiver();
	// 	printHeliostatsFieldPerformance();
	// 	printPowerOutput();
	// 	printDemand();
	// }
	// else{
	// 	printColdStorage();
	// 	printHotStorage();
	// 	printSteamGeneratorOutlet();
	// 	printRequiredThermalPower();
	// 	printReceiver();
	// 	printHeliostatsFieldLayout();
	// 	printHeliostatsFieldPerformance();
	// 	printPowerOutput();
	// 	printDemand();
	// }
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
    for (int i = 0; i < _numberOfTimeIncrements; ++i) {
      time = timeVector[i] % 1440;
      if (time < _tStart || time > _tEnd)
	_demand.push_back(0.);
      else
	_demand.push_back(_maximumPowerDemand);
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

/*-----------------------------------------*/
/*              read inputs                */
/*-----------------------------------------*/
bool Scenario::read_x ( const std::string & x_file_name ) {

  std::ifstream in ( x_file_name );
  if ( in.fail() ) {
    in.close();
    return false;
  }
 
  for ( int i = 0 ; i < _n ; ++i ) {
    in >> _x[i];
    if ( in.fail() ) {
      in.close();
      return false;
    }
    in >> std::ws;
  }
  
  in.close();
  return true;
}

/*----------------------------------------------*/
/*        diplay vector of variables x          */
/*----------------------------------------------*/
void Scenario::display_x ( std::ostream & out ) const {
  for ( int i = 0 ; i < _n ; ++i )
    out << _x[i] << " ";
}

/*-------------------------------------*/
/*  delete vector of variables: _x     */
/*-------------------------------------*/
void Scenario::delete_x ( void ) {
  _n = 0;
  if ( _x ) {
    delete [] _x;
    _x = NULL;
  }
}
