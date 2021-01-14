#include "helpFunctions.hpp"
#include <iomanip>  // for setw

/*-----------------------------------------------------------*/
/*                   display list of problems                */
/*-----------------------------------------------------------*/
void display_problems ( std::ostream & out , const std::vector<Problem> & problems ) {

  out << "\t#\t" << std::setw(15) << "pb_id"
      << "\t"    << std::setw(37) << "obj.(f)"
      << "\t"    << std::setw(15) << "# of var.(n)"
      << "\t"    << std::setw(15) << "# of constr.(m)\n\n";

  for ( size_t i = 0 ; i < problems.size() ; ++i )
    out << "\t" << i+1
	<< "\t" << std::setw(15) << problems[i].get_pb_id()
	<< "\t" << std::setw(37) << problems[i].get_f_description()
	<< "\t" << std::setw(15) << problems[i].get_n()
      	<< "\t" << std::setw(15) << problems[i].get_m()
	<< std::endl;
}

/*-----------------------------------------------------------*/
/*                        display usage                      */
/*-----------------------------------------------------------*/
void display_usage ( std::ostream & out ) {
  out << "Run SOLAR : solar pb_id x.txt (add -v for verbose mode)" << std::endl
      << "Help(1)   : solar -h"               << std::endl
      << "Help(2)   : solar -h pb_id"         << std::endl
      << "Info      : solar -i"               << std::endl;
}

/*-----------------------------------------------------------*/
/*                        display info                       */
/*-----------------------------------------------------------*/
void display_info ( std::ostream & out , const std::string & version ) {
  out << std::endl << "SOLAR, the solar thermal power plant simulator, version "
      << version << std::endl << std::endl
      << "Contributors: M. Diago, S. Le Digabel, M. Lemyre-Garneau, B. Talgorn;"
      << " GERAD and Polytechnique Montreal"                      << std::endl << std::endl
      << "This code is distributed under the LGPL license"       << std::endl // TOTO: peut-etre autre licence
      << "https://github.com/bbopt/solar" << std::endl << std::endl
      << "Please report bugs to sebastien.le-digabel@polymtl.ca" << std::endl << std::endl;
}

/*-----------------------------------------------------------*/
/*              display (short) help for all problems        */
/*-----------------------------------------------------------*/
void display_help ( std::ostream & out , const std::vector<Problem> & problems ) {
  out << std::endl
      << "Run simulation    : solar pb_id x.txt"             << std::endl
      << "Help for a problem: solar pb_id or solar -h pb_id" << std::endl << std::endl
      << "\tpb_id: Problem ID"      << std::endl
      << "\tx.txt: Input vector (values separated with spaces, tabs, or line breaks)" << std::endl
      << std::endl
      << "List of problems:" << std::endl << std::endl;
  display_problems ( out , problems );
  out << std::endl;
}

/*-----------------------------------------------------------*/
/*           display (complete) help for one problem         */
/*-----------------------------------------------------------*/
void display_help ( std::ostream               & out      ,
		    const std::vector<Problem> & problems ,
		    const std::string          & pb_id      ) {
 
  out << "DISPLAY HELP FOR PROBLEM \"" << pb_id << "\"" << std::endl << std::endl;
  
  const Problem * pb = find_problem ( problems, pb_id );

  if ( pb ) {
   
    out << "PROBLEM: " << pb->get_pb_id() << " (solar" << pb->get_index() << ")"
	<< "\t"        << pb->get_f_description()
	<< "\tn="      << pb->get_n()
	<< "\tm="      << pb->get_m()
	<< std::endl;

    // #1:
    if ( pb->get_pb_id() == "MAXNRG_H1" )
      print_maxNrg_H1 ( out );

    // #2:
    else if ( pb->get_pb_id() == "MINSURF_H1" )
      print_minSurf_H1 ( out );

    // #3:
    else if ( pb->get_pb_id() == "MINCOST_C1" )
      print_minCost_C1 ( out);

    // #4:
    else if ( pb->get_pb_id() == "MINCOST_C2" )
      print_minCost_C2 ( out );

    // #5:
    else if ( pb->get_pb_id() == "MAXCOMP_HTF1" )
      print_maxComp_HTF1 ( out );

    // #6:
    else if ( pb->get_pb_id() == "MINCOST_TS" )
      print_minCost_TS ( out );

    // #7:
    else if ( pb->get_pb_id() == "MAXEFF_RE" )
      print_maxEff_RE ( out );

    // #8:
    else if ( pb->get_pb_id() == "MAXHF_MINCOST" )
      print_maxHF_minCost ( out );

    // #9:
    else if ( pb->get_pb_id() == "MAXNRG_MINPAR" )
      print_maxNrg_minPar ( out );

    else {
      out << "Cannot find detailed help for this problem" << std::endl;
      return;
    }
  }
  else {
    out << "This problem id is invalid" << std::endl;
    return;
  }
  out << "\n-----------------------------------------------------------------\n" << std::endl;
}

/*--------------  #1  ---------------------------*/
void print_maxNrg_H1 ( std::ostream & out ) {
/*-----------------------------------------------*/
  
  out << "\n-----------------------------------------------------------------\n"
      << "Parameters:\n"
      << "\tHeliostats field only\n"
      << "\tLatitude: 44.95 deg N\n"
      << "\tDay: April 10th\n"  // https://www.epochconverter.com/days/2019
      << "\tDuration: 24 hours\n"
      << "\tMaximum field surface: 195 hectares\n"
      << "\tBudget: $50M\n"
      << "\tMust provide 100% of the demand requirement\n"
      << std::endl
      << "Objective (first output)\n"
      << "\tMaximize the total solar energy concentrated on the receiver aperture through one day (kWh)\n"
      << std::endl
      << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\tx1: Heliostat height (m)         : Real in [ 1; 40]\n"
      << "\t\tx2: Heliostat width (m)          : Real in [ 1; 40]\n"
      << "\t\tx3: Tower height (m)             : Real in [20;250]\n"
      << "\t\tx4: Receiver aperture height (m) : Real in [ 1; 30]\n"
      << "\t\tx5: Receiver aperture width (m)  : Real in [ 1; 30]\n"
      << "\t\tx6: Maximum number of heliostats : Integer >= 1\n"
      << "\t\tx7: Field max. angular span (deg): Real in [1;89]\n"
      << "\t\tx8: Minimum distance to tower (%): Real in [0;20]\n"
      << "\t\tx9: Maximum distance to tower (%): Real in [0;20]\n"
      << std::endl
      << "Constraints (outputs 2 to 6 with format ci <= 0):\n"
      << "\tc1: Cost of field <= budget\n"
      << "\tc2: Field surface area\n"
      << "\tc3: Tower is at least twice as high as heliostats\n"
      << "\tc4: Min. distance to tower <= Max. distance to tower\n"
      << "\tc5: Max. number of heliostats\n"
      << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 9 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $1"     << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(    R    R     R    R    R   I    R    R    R )" << std::endl
      << "\tLOWER_BOUND      " << "(  1.0  1.0  20.0  1.0  1.0   1  1.0  0.0  0.0 )" << std::endl
      << "\tX0               " << "(  8.0  8.0 150.0  7.0  7.0 250 45.0  0.5  5.0 )" << std::endl
      << "\tUPPER_BOUND      " << "( 40.0 40.0 250.0 30.0 30.0   - 89.0 20.0 20.0 )" << std::endl;
}
  
/*--------------  #2  ---------------------------*/
void print_minSurf_H1 ( std::ostream & out ) {
/*-----------------------------------------------*/

  out << "\n-----------------------------------------------------------------\n"
      << "Parameters:\n"
      << "\tLatitude: 37.55 deg N\n"
      << "\tDay: June 29th\n" // https://www.epochconverter.com/days/2019
      << "\tDuration: 72 hours\n"
      << "\tMaximum field surface: 400 hectares\n"
      << "\tBudget: $300M\n"
      << "\tMust provide 100% of the demand requirement\n"
      << std::endl
      << "Objective (first output)\n"
      << "\tMminimize total heliostats field surface (m^2) with constraints to run a pre-determined powerplant\n"
      << std::endl
      << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\t x1: Heliostat height (m)         : Real in [ 1; 40]\n"
      << "\t\t x2: Heliostat width (m)          : Real in [ 1; 40]\n"
      << "\t\t x3: Tower height (m)             : Real in [20;250]\n"
      << "\t\t x4: Receiver aperture height (m) : Real in [ 1; 30]\n"
      << "\t\t x5: Receiver aperture width (m)  : Real in [ 1; 30]\n"
      << "\t\t x6: Maximum number of heliostats : Integer >= 1\n"
      << "\t\t x7: Field max. angular span (deg): Real in [1;89]\n"
      << "\t\t x8: Minimum distance to tower (%): Real in [0;20]\n"
      << "\t\t x9: Maximum distance to tower (%): Real in [1;20]\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx10: Central receiver outlet temperature (K): Real in [793;995]\n"
      << "\t\tx11: Receiver number of tubes               : Integer >= 1\n"
      << "\t\tx12: Receiver insulation thickness (m)      : Real in [0.01 ;5  ]\n"
      << "\t\tx13: Receiver tubes: inside diameter (m)    : Real in [0.005;0.1]\n"
      << "\t\tx14: Receiver tubes: outside diameter (m)   : Real in [0.005;0.1]\n"   
      << std::endl
      << "Constraints (outputs 2 to 14 with format ci <= 0):\n"
      << "\t c1: Field surface\n"
      << "\t c2: Demand compliance\n"
      << "\t c3: Budget\n"
      << "\t c4: Tower is at least twice as high as heliostats\n"
      << "\t c5: Min. distance to tower <= Max. distance to tower\n"
      << "\t c6: Max. number of heliostats\n"
      << "\t c7: Pressure in receiver tubes <= yield pressure\n"
      << "\t c8: Molten salt melting point <= hot storage lowest temperature\n"
      << "\t c9: Molten salt melting point <= cold storage lowest temperature\n"
      << "\tc10: Molten salt melting point <= steam generatour outlet temperature\n"
      << "\tc11: Receiver inside/outside diameter\n"
      << "\tc12: Tubes fit in receiver\n" 
      << "\tc13: Receiver outlet temperature > steam turbine inlet temperature\n";

  out << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 14 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $2" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(    R    R     R    R    R    I    R    R     R     R  I    R      R      R )" << std::endl
      << "\tLOWER_BOUND      " << "(  1.0  1.0  20.0  1.0  1.0    1  1.0  0.0   1.0 793.0  1 0.01  0.005 0.0050 )" << std::endl
      << "\tX0               " << "( 11.0 11.0 140.0 10.0 10.0 2650 89.0  0.5   5.0 838.0 36 0.30  0.020 0.0216 )" << std::endl
      << "\tUPPER_BOUND      " << "( 40.0 40.0 250.0 30.0 30.0    - 89.0 20.0  20.0 995.0  - 5.00  0.100 0.1000 )" << std::endl;
}

/*--------------  #3  ---------------------------*/
void print_minCost_C1 ( std::ostream & out ) {
/*-----------------------------------------------*/  
  out << "\n-----------------------------------------------------------------\n"
      << "Parameters:\n"
      << "\tWhole plant\n"
      << "\tLatitude: 35 deg N\n"
      << "\tDay: September 27th\n"
      << "\tDuration: 24 hours\n"
      << "\tDemand profile: 10 MW, starting at 3PM and ending at 9PM, 3 consecutive days\n"
      << "\tMaximum field surface: 800 000 m^2\n"
      << "\tMust provide 100% of the demand requirement\n"
      << std::endl
      << "Objective (first output)\n"
      << "\tMinimize Total investment cost ($)\n"
      << std::endl
      << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\t x1: Heliostat height            : Real in [ 1; 40]\n"
      << "\t\t x2: Heliostat width             : Real in [ 1; 40]\n"
      << "\t\t x3: Tower height                : Real in [20;250]\n"
      << "\t\t x4: Receiver aperture height    : Real in [ 1; 30]\n"
      << "\t\t x5: Receiver aperture width     : Real in [ 1; 30]\n"
      << "\t\t x6: Maximum number of heliostats: Integer >= 1\n"
      << "\t\t x7: Field max. angular span     : Real in [1;89]\n"
      << "\t\t x8: Minimum distance to tower   : Real in [0;20]\n"
      << "\t\t x9: Maximum distance to tower   : Real in [1;20]\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx10: Central receiver outlet temperature (K): Real in [793;995]\n"
      << "\t\tx11: Hot storage height (m)                 : Real in [1;50]\n"
      << "\t\tx12: Hot storage diameter (m)               : Real in [1;30]\n"
      << "\t\tx13: Hot storage insulation thickness (m)   : Real in [0.01;5]\n"
      << "\t\tx14: Cold storage insulation thickness (m)  : Real in [0.01;5]\n"
      << "\t\tx15: Cold molten salt min. temperature (K)  : Real in [495;650]\n"
      << "\t\tx16: Receiver number of tubes               : Integer >= 1\n"
      << "\t\tx17: Receiver insulation thickness (m)      : Real in [0.01 ;5  ]\n"
      << "\t\tx18: Receiver tubes: inside diameter (m)    : Real in [0.005;0.1]\n"
      << "\t\tx19: Receiver tubes: outside diameter (m)   : Real in [0.005;0.1]\n"
      << "\tPowerblock:\n"
      << "\t\tx20: Type of turbine: Integer in {1, 2, ..., 8}\n"
      << std::endl
      << "Constraints (outputs 2 to 14 with format ci <= 0):\n"
      << "\tc1: Surface area\n"
      << "\tc2: Compliance to demand\n"
      << "\tc3: Tower is at least twice as high as heliostats\n"
      << "\tc4: Min. distance to tower <= max. distance to tower\n"
      << "\tc5: Number of heliostats <= number of positions in the field\n"
      << "\tc6: Pressure in tubes does not exceed yield pressure\n"
      << "\tc7, c8, c9: Molten salt temperature does not fall below the melting point in storages and\n"
      << "\t            at the steam generator outlet\n"
      << "\tc10: Tubes inside diameter is smaller than outter diameter\n"
      << "\tc11: Tubes can fit inside the receiver\n"
      << "\tc12: Central receiver outlet temperature is higher than that required by the turbine\n"
      << "\tc13: Check if storage is back to initial conditions\n"
      << "\n-----------------------------------------------------------------\n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 20       << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $3" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(    R    R     R    R    R   I    R    R    R     R    R    R    R    R     R  I    R     R     R I )" << std::endl
      << "\tLOWER_BOUND      " << "(  1.0  1.0  20.0  1.0  1.0   1  1.0  0.0  1.0 793.0  1.0  1.0 0.01 0.01 495.0  1 0.01 0.005 0.005 1 )" << std::endl
      << "\tX0               " << "(  8.0  8.0 150.0  7.0  7.0 250 45.0  0.5  5.0 900.0  9.0  9.0 0.30 0.20 560.0 40 0.30 0.015 0.017 3 )" << std::endl
      << "\tUPPER_BOUND      " << "( 40.0 40.0 250.0 30.0 30.0   - 89.0 20.0 20.0 995.0 50.0 30.0 5.00 5.00 650.0  - 5.00 0.100 0.100 8 )" << std::endl;
}

/*--------------  #4  ---------------------------*/
void print_minCost_C2 ( std::ostream & out ) {
/*-----------------------------------------------*/
  out << "\n-----------------------------------------------------------------\n"
      << "Parameters:\n"
      << "\tWhole plant\n"
      << "\tLatitude: 35 deg N\n"
      << "\tDay: January 1st\n"
      << "\tDuration: 24 hours\n"
      << "\tDemand profile: 25 MW, starting at 3PM and ending at 9PM, 3 consecutive days\n"
      << "\tMaximum field surface: 2 000 000 m^2 (200 hectares)\n"
      << "\tMust provide 100% of the demand requirement\n"
      << std::endl

      << "Objective (first output)\n"
      << "\tMinimize the cost of powerplant to respect a given demand with a limited size of field ($)\n"
      << std::endl  

      << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\t x1: Heliostat height            : Real in [ 1; 40]\n"
      << "\t\t x2: Heliostat width             : Real in [ 1; 40]\n"
      << "\t\t x3: Tower height                : Real in [20;250]\n"
      << "\t\t x4: Receiver aperture height    : Real in [ 1; 30]\n"
      << "\t\t x5: Receiver aperture width     : Real in [ 1; 30]\n"
      << "\t\t x6: Maximum number of heliostats: Integer >= 1\n"
      << "\t\t x7: Field max. angular span     : Real in [1;89]\n"
      << "\t\t x8: Minimum distance to tower   : Real in [0;20]\n"
      << "\t\t x9: Maximum distance to tower   : Real in [1;20]\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx10: Central receiver outlet temperature (K): Real in [793;995]\n"
      << "\t\tx11: Hot storage height (m)                 : Real in [1;50]\n"
      << "\t\tx12: Hot storage diameter (m)               : Real in [1;30]\n"
      << "\t\tx13: Hot storage insulation thickness (m)   : Real in [0.01;5]\n"
      << "\t\tx14: Cold storage insulation thickness (m)  : Real in [0.01;5]\n"
      << "\t\tx15: Cold molten salt min. temperature (K)  : Real in [495;650]\n"
      << "\t\tx16: Receiver number of tubes               : Integer >= 1\n"
      << "\t\tx17: Receiver insulation thickness (m)      : Real in [0.01 ;5  ]\n"
      << "\t\tx18: Receiver tubes: inside diameter (m)    : Real in [0.005;0.1]\n"
      << "\t\tx19: Receiver tubes: outside diameter (m)   : Real in [0.006;0.1]\n"
      << "\tSteam generator:\n"
      << "\t\tx20: Tubes spacing (m)          : Real in [0.007;0.2]\n"
      << "\t\tx21: Tubes length (m)           : Real in [0.5;10]\n"
      << "\t\tx22: Tubes: inside diameter (m) : Real in [0.005;0.1]\n"
      << "\t\tx23: Tubes: outside diameter (m): Real in [0.006;0.1]\n"
      << "\t\tx24: Baffle cut                 : Real in [0.15;0.4]\n"
      << "\t\tx25: Number of baffles          : Integer >= 2\n"   
      << "\t\tx26: Number of tubes            : Integer >= 1\n"
      << "\t\tx25: Number of shells           : Integer in {1, 2, ..., 10}\n"
      << "\t\tx28: Number of passes per shell : Integer in {1, 2, ..., 9}\n"   
      << "\tPowerblock:\n"
      << "\t\tx29: Type of turbine: Integer in {1, 2, ..., 8}\n"
      << std::endl

      << "Constraints (outputs 2 to 17 with format ci <= 0):\n"
      << "\t c1: Maximum heliostats field surface\n"
      << "\t c2: Compliance to the demand\n"
      << "\t c3: Central tower higher than heliostats\n"
      << "\t c4: Min. distance to tower <= max. distance to tower\n"   
      << "\t c5: Number of heliostats <= number of positions in the grid\n"
      << "\t c6: Pressure in receiver tubes does not exceed yield pressure\n"
      << "\t c7: Hot storage temperature >= molten salt melting point\n"
      << "\t c8: Cold storage temperature >= molten salt melting point\n"
      << "\t c9: Steam generator outlet temperature >= molten salt melting point\n"
      << "\tc10: Receiver inside diameter <= receiver outside diameter\n"
      << "\tc11: Number of tubes in receiver fit inside receiver\n"
      << "\tc12: Receiver outlet temperature must exceed steam turbine inlet temperature\n"
      << "\tc13: Parasitics do not exceed 20% of energy production\n"
      << "\tc14: Steam generator outer tubes diameter <= tubes spacing\n"
      << "\tc15: Steam generator inside diameter <=  steam generator outside diameter\n"
      << "\tc16: Pressure in steam generator tubes <= yield pressure\n"

      << "\n-----------------------------------------------------------------\n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 29 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $4" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(    R    R     R    R    R    I    R    R    R     R    R    R    R    R     R   I    R      R     R     R    R      R     R    R I     I  I I I )" << std::endl
      << "\tLOWER_BOUND      " << "(  1.0  1.0  20.0  1.0  1.0    1  1.0  0.0  1.0 793.0  1.0  1.0 0.01 0.01 495.0   1 0.01 0.0050 0.006 0.007  0.5 0.0050 0.006 0.15 2     1  1 1 1 )" << std::endl
      << "\tX0               " << "(  9.0  9.0 150.0  6.0  8.0 1000 45.0  0.5  5.0 900.0  9.0  9.0 0.30 0.20 560.0 500 0.30 0.0165 0.018 0.017 10.0 0.0155 0.016 0.20 3 12000  1 2 2 )" << std::endl  
      << "\tUPPER_BOUND      " << "( 40.0 40.0 250.0 30.0 30.0    - 89.0 20.0 20.0 995.0 50.0 30.0 5.00 5.00 650.0   - 5.00 0.1000 0.100 0.200 10.0 0.1000 0.100 0.40 -     - 10 9 8 )" << std::endl;
}

/*--------------  #5  ---------------------------*/
void print_maxComp_HTF1 ( std::ostream & out ) {
/*-----------------------------------------------*/
  out << "\n-----------------------------------------------------------------\n"
      << "Parameters:\n"
      << "\tWhole plant\n"
      << "\tLatitude: 37.5581 deg N\n"
      << "\tDay: January 30th\n"  // https://www.epochconverter.com/days/2019
      << "\tDuration: 30 days\n"
      << "\tDemand: 12 MW\n"
      << "\tBudget: $300M\n"
      << "\tPre-determined sunlight input for a period of 1 month\n"
      << std::endl

      << "Objective (first output)\n"
      << "\tMaximize compliance to a demand profile\n"
      << std::endl  
      << "Variables:\n"
      << "\tHeat transfer loop:\n"
      << "\t\t x1: Central receiver outlet temperature (K): Real in [793;995]\n"
      << "\t\t x2: Hot storage height (m)                 : Real in [1;30]\n"
      << "\t\t x3: Hot storage diameter (m)               : Real in [1;30]\n"
      << "\t\t x4: Hot storage insulation thickness (m)   : Real in [0.01;2]\n"
      << "\t\t x5: Cold storage insulation thickness (m)  : Real in [0.01;2]\n"
      << "\t\t x6: Cold molten salt min. temperature (K)  : Real in [495;650]\n"
      << "\t\t x7: Receiver number of tubes               : Integer >= 1\n"
      << "\t\t x8: Receiver insulation thickness (m)      : Real in [0.100;2.0]\n"
      << "\t\t x9: Receiver tubes: inside diameter (m)    : Real in [0.005;0.1]\n"
      << "\t\tx10: Receiver tubes: outside diameter (m)   : Real in [0.005;0.1]\n"
      << "\tSteam generator:\n"
      << "\t\tx11: Tubes spacing (m)              : Real in [0.006; 0.2]\n"
      << "\t\tx12: Tubes length  (m)              : Real in [0.500;10.0]\n"
      << "\t\tx13: Tubes inner diameter (m)       : Real in [0.005; 0.1]\n"
      << "\t\tx14: Tubes outer diameter (m)       : Real in [0.006; 0.1]\n"
      << "\t\tx15: Baffles cut                    : Real in [0.150; 0.4]\n"
      << "\t\tx16: Number of baffles              : Integer >= 2\n"
      << "\t\tx17: Number of tubes                : Integer >= 1\n"
      << "\t\tx18: Number of shell passes         : Integer in {1, 2, ...,10}\n"
      << "\t\tx19: Number of tube passes per shell: Integer in {1, 2, ..., 9}\n"
      << "\tPowerblock:\n"
      << "\t\tx20: Type of turbine: Integer in {1, 2, ..., 8}\n"
      << std::endl;

  out << "Constraints (outputs 2 to 13 with format ci <= 0):\n"
      << "\t c1: Budget\n"
      << "\t c2: Pressure in receiver tubes <= yield pressure\n"
      << "\t c3: Molten salt melting point <= hot storage lowest temperature\n"
      << "\t c4: Molten salt melting point <= cold storage lowest temperature\n"
      << "\t c5: Molten salt melting point <= steam generatour outlet temperature\n"
      << "\t c6: Receiver inside diameter < receiver outside diameter\n"
      << "\t c7: Number of tubes in receiver fit inside receiver\n"
      << "\t c8: Receiver outlet temperature >= steam turbine inlet temperature\n"
      << "\t c9: Parasitic losses <= 18% of the generated output\n"
      << "\tc10: Steam generator tubes outer diameter <= tubes spacing\n"
      << "\tc11: Steam generator tubes inner diameter <= tubes outer diameter\n"
      << "\tc12: Pressure in steam generator tubes does not exceed yield pressure\n";

  out << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 20 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $5" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(     R    R    R    R    R   R  I    R     R     R      R    R       R     R    R I    I  I I I )" << std::endl
      << "\tLOWER_BOUND      " << "( 793.0  1.0  1.0 0.01 0.01 495  1 0.10 0.005 0.005  0.006  0.5   0.005 0.006 0.15 2    1  1 1 1 )" << std::endl
      << "\tX0               " << "( 900.0 10.0 12.0 0.15 0.10 560 24 0.35 0.020 0.023  0.050  8.0   0.020 0.023 0.20 2 5000  5 5 1 )" << std::endl
      << "\tUPPER_BOUND      " << "( 995.0 30.0 30.0 2.00 2.00 650  - 2.00 0.100 0.100  0.200 10.0   0.100 0.100 0.4  -    - 10 9 8 )" << std::endl;
}

/*--------------  #6  ---------------------------*/
void print_minCost_TS ( std::ostream & out ) {
/*-----------------------------------------------*/
    
  out << "\n-----------------------------------------------------------------\n"
      << "Parameters:\n"
      << "\tWhole plant\n"
      << "\tLatitude: 30.05 deg N\n"
      << "\tDay: January 1st\n"  // https://www.epochconverter.com/days/2019
      << "\tDuration: 24 hours\n"
      << "\tDemand: 120 MW\n"
      << "\tMust provide 100% of the demand requirement\n"
      << "\tThis instance runs a predetermined power plant using the molten salt cycle and power block models.\n"
      << "\tThe objective is to minimize the cost of the thermal storage units so that the power plant is able\n"
      << "\tto sustain a 120 MW electrical power output during 24 hours. Since the heliostat field is not being\n"
      << "\toptimized, its hourly power output is taken from prerecorded data instead of being simulated.\n"
      << std::endl;

  out << "Objective (first output)\n"
      << "\tMinimize the cost of storage\n"
      << std::endl;

  out << "Variables:\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx1: Central receiver outlet temperature (K): Real in [793;995]\n"
      << "\t\tx2: Hot storage height (m)                 : Real in [2;50]\n"
      << "\t\tx3: Hot storage diameter (m)               : Real in [2;30]\n"
      << "\t\tx4: Hot storage insulation thickness (m)   : Real in [0.01;5]\n"
      << "\t\tx5: Cold storage insulation thickness (m)  : Real in [0.01;5]\n"

      << "Constraints (outputs 2 to 7 with format ci <= 0):\n"
      << "\tc1: Compliance to the demand\n"
      << "\tc2: Pressure in receiver tubes does not exceed yield pressure\n"
      << "\tc3: Molten salt melting point <= hot storage lowest temperature\n"
      << "\tc4: Molten salt melting point <= cold storage lowest temperature\n"
      << "\tc5: Receiver outlet temperature >= steam turbine inlet temperature\n"
      << "\tc6: At midnight, storage must be at least at its original conditions\n"
    
      << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 5 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $6" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(     R    R    R    R    R )" << std::endl
      << "\tLOWER_BOUND      " << "( 793.0  2.0  2.0 0.01 0.01 )" << std::endl
      << "\tX0               " << "( 900.0 10.0 12.0 0.20 0.20 )" << std::endl
      << "\tUPPER_BOUND      " << "( 995.0 50.0 30.0 5.00 5.00 )" << std::endl;
}

/*--------------  #7  ---------------------------*/
void print_maxEff_RE ( std::ostream & out ) {
/*-----------------------------------------------*/

  out << "\n-----------------------------------------------------------------\n"
      << "Parameters:\n"
      << "\tWhole plant\n"
      << "\tLatitude: 30.05 deg N\n"
      << "\tDay: January 1st\n"  // https://www.epochconverter.com/days/2019
      << "\tDuration: 24 hours\n"
      << std::endl;
  
  out << "Objective (first output)\n"
      << "\tMaximize receiver efficiency, i.e the energy transfered to the molten salt\n"
      << std::endl;

  out << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\tx1: Receiver aperture height    : Real in [ 1; 30]\n"
      << "\t\tx2: Receiver aperture width     : Real in [ 1; 30]\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx3: Central receiver outlet temperature (K): Real in [793;995]\n"
      << "\t\tx4: Receiver number of tubes               : Integer >= 1\n"
      << "\t\tx5: Receiver insulation thickness (m)      : Real in [0.01  ;5.0]\n"
      << "\t\tx6: Receiver tubes: inside diameter (m)    : Real in [0.005 ;0.1]\n"
      << "\t\tx7: Receiver tubes: outside diameter (m)   : Real in [0.0055;0.1]\n"
     
      << "Constraints (outputs 2 to 7 with format ci <= 0):\n"
      << "\tc1: Budget\n"
      << "\tc2: Pressure in tubes does not exceed yield pressure\n"
      << "\tc3: Receiver inside diameter <= receiver outside diameter\n"
      << "\tc4: Receiver outlet temperature must exceed steam turbine inlet temperature\n"
      << "\tc5: Tubes fit in receiver\n"
      << "\tc6: Parasitics do not exceed 5% of the absorbed power\n"
    
      << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 7 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $7" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(    R    R     R  I    R     R      R )" << std::endl
      << "\tLOWER_BOUND      " << "(  1.0  1.0 793.0  1 0.01 0.005 0.0055 )" << std::endl
      << "\tX0               " << "(  7.0  7.0 850.0 40 0.20 0.010 0.0110 )" << std::endl
      << "\tUPPER_BOUND      " << "( 30.0 30.0 995.0  - 5.00 0.100 0.1000 )" << std::endl;
}

/*--------------  #8  ---------------------------*/
void print_maxHF_minCost ( std::ostream & out ) {
/*-----------------------------------------------*/
  
  out << "\n-----------------------------------------------------------------\n"
      << "Parameters:\n"
      << "\tWhole plant\n"
      << "\tLatitude: 45 deg N\n"
      << "\tDay: January 1st\n"  // https://www.epochconverter.com/days/2019
      << "\tDuration: 24 hours\n"
      << "\tMaximum field surface: 4x10^6 m^2\n"
      << "\tMinimum energy production: 400MWh\n"
      << std::endl;
    
  out << "Objectives (first and second outputs)\n"
      << "\tMaximize heliostat field performance (absorbed energy) and minimize cost of field, tower and receiver\n"
      << std::endl;

  out << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\t x1: Heliostat height            : Real in [ 1; 40]\n"
      << "\t\t x2: Heliostat width             : Real in [ 1; 40]\n"
      << "\t\t x3: Tower height                : Real in [20;250]\n"
      << "\t\t x4: Receiver aperture height    : Real in [ 1; 30]\n"
      << "\t\t x5: Receiver aperture width     : Real in [ 1; 30]\n"
      << "\t\t x6: Maximum number of heliostats: Integer >= 1\n"
      << "\t\t x7: Field max. angular span     : Real in [1;89]\n"
      << "\t\t x8: Minimum distance to tower   : Real in [0;20]\n"
      << "\t\t x9: Maximum distance to tower   : Real in [1;20]\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx10: Receiver number of tubes               : Integer >= 1\n"
      << "\t\tx11: Receiver insulation thickness (m)      : Real in [0.01 ;5  ]\n"
      << "\t\tx12: Receiver tubes: inside diameter (m)    : Real in [0.005;0.1]\n"
      << "\t\tx13: Receiver tubes: outside diameter (m)   : Real in [0.006;0.1]\n"

      << "Constraints (outputs 3 to 11 with format ci <= 0):\n"
      << "\tc1: Maximum heliostats field surface\n"
      << "\tc2: Tower is at least twice as high as heliostats\n"
      << "\tc3: Min. distance to tower <= max. distance to tower\n"   
      << "\tc4: Number of heliostats <= number of positions in the grid\n"
      << "\tc5: Pressure in receiver tubes <= yield pressure\n"
      << "\tc6: Receiver inside diameter <= receiver outside diameter\n"
      << "\tc7: Tubes fit in receiver\n"
      << "\tc8: Minimal acceptable energy production\n"
      << "\tc9: Parasitics must not exceed 8% of absorbed energy\n"
  
      << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "DIMENSION        " << 13 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $8" << std::endl
      << "BB_OUTPUT_TYPE   " << "OBJ OBJ CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "BB_INPUT_TYPE    " << "(    R    R     R    R    R    I    R    R    R  I    R     R      R )" << std::endl
      << "LOWER_BOUND      " << "(  1.0  1.0  20.0  1.0  1.0    1  1.0  0.0  1.0  1 0.01 0.005 0.0060 )" << std::endl
      << "X0               " << "( 11.0 11.0 200.0 10.0 10.0 2650 89.0  0.5  8.0 36 0.30 0.020 0.0216 )" << std::endl
      << "UPPER_BOUND      " << "( 40.0 40.0 250.0 30.0 30.0    - 89.0 20.0 20.0  - 5.00 0.100 0.1000 )" << std::endl;
}

/*--------------  #9  ---------------------------*/
void print_maxNrg_minPar ( std::ostream & out ) {
/*-----------------------------------------------*/
  
  out << "\n-----------------------------------------------------------------\n"
      << "Parameters:\n"
      << "\tWhole plant\n"
      << "\tLatitude: 25 deg N\n"
      << "\tDay: June 29th\n"  // https://www.epochconverter.com/days/2019
      << "\tDuration: 24 hours\n"
      << "\tMaximum field surface:  5M m^2\n"
      << "\tBudget: $1.2B\n"
      << "\tMinimum energy production: 250MWh\n"
      << std::endl;
  
  out << "Objectives (first and second outputs)\n"
      << "\tMaximize power and minimize losses\n"
      << std::endl;
  
  out << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\t x1: Heliostat height            : Real in [ 1; 40]\n"
      << "\t\t x2: Heliostat width             : Real in [ 1; 40]\n"
      << "\t\t x3: Tower height                : Real in [20;250]\n"
      << "\t\t x4: Receiver aperture height    : Real in [ 1; 30]\n"
      << "\t\t x5: Receiver aperture width     : Real in [ 1; 30]\n"
      << "\t\t x6: Maximum number of heliostats: Integer >= 1\n"
      << "\t\t x7: Field max. angular span     : Real in [1;89]\n"
      << "\t\t x8: Minimum distance to tower   : Real in [0;20]\n"
      << "\t\t x9: Maximum distance to tower   : Real in [1;20]\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx10: Central receiver outlet temperature (K): Real in [793;995]\n"
      << "\t\tx11: Hot storage height (m)                 : Real in [1;50]\n"
      << "\t\tx12: Hot storage diameter (m)               : Real in [1;30]\n"
      << "\t\tx13: Hot storage insulation thickness (m)   : Real in [0.01;5]\n"
      << "\t\tx14: Cold storage insulation thickness (m)  : Real in [0.01;5]\n"
      << "\t\tx15: Cold molten salt min. temperature (K)  : Real in [495;650]\n"
      << "\t\tx16: Receiver number of tubes               : Integer >= 1\n"
      << "\t\tx17: Receiver insulation thickness (m)      : Real in [0.01 ;5  ]\n"
      << "\t\tx18: Receiver tubes: inside diameter (m)    : Real in [0.005;0.1]\n"
      << "\t\tx19: Receiver tubes: outside diameter (m)   : Real in [0.006;0.1]\n"
      << "\t\tx20: Exchanger tubes spacing (m)            : Real in [0.007;0.2]\n"
      << "\t\tx21: Exchanger tubes length (m)             : Real in [0.5;10]\n"
      << "\t\tx22: Exchanger tubes: inside diameter (m)   : Real in [0.005;0.1]\n"
      << "\t\tx23: Exchanger tubes: outside diameter (m)  : Real in [0.006;0.1]\n"
      << "\t\tx24: Exchanger baffle cut (m)               : Real in [0.15;0.4]\n"
      << "\t\tx25: Exchanger number of baffles            : Integer >= 2\n"   
      << "\t\tx26: Exchanger number of tubes              : Integer >= 1\n"
      << "\t\tx25: Exchanger number of shells             : Integer in {1, 2, ..., 10}\n"
      << "\t\tx28: Exchanger number of passes per shell   : Integer in {1, 2, ..., 9}\n"   
      << "\tPowerblock:\n"
      << "\t\tx29: Type of turbine: Integer in {1, 2, ..., 8}\n"
      << std::endl;

  out << "Constraints (outputs 3 to 19 with format ci <= 0):\n"
      << "\t c1: Cost\n"
      << "\t c2: Minimum energy production is reached\n"
      << "\t c3: Maximum heliostats field surface\n"
      << "\t c4: Tower is at least twice as high as heliostats\n"
      << "\t c5: Min. distance to tower <= max. distance to tower\n"  
      << "\t c6: Number of heliostats <= number of positions in the field\n"
      << "\t c7: Pressure in receiver tubes <= yield pressure\n"
      << "\t c8: Molten salt melting point <= hot storage lowest temperature\n"
      << "\t c9: Molten salt melting point <= cold storage lowest temperature\n"
      << "\tc10: Molten salt melting point <= steam generatour outlet temperature\n"
      << "\tc11: Receiver inside diameter <= receiver outside diameter\n"
      << "\tc12: Tubes fit in receiver\n"
      << "\tc13: Receiver outlet temperature must exceed steam turbine inlet temperature\n"
      << "\tc14: Ratio parasitics vs power output <= 20%\n"
      << "\tc15: Steam generator outer tubes diameter <= tubes spacing\n"
      << "\tc16: Steam generator inside diameter <= steam generator outside diameter\n"
      << "\tc17: Pressure in steam generator tubes <= yield pressure\n";
  
  out << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "DIMENSION        " << 29 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $9" << std::endl
      << "BB_OUTPUT_TYPE   " << "OBJ OBJ CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "BB_INPUT_TYPE    " << "(    R    R     R    R    R    I    R    R    R     R    R    R    R    R     R  I    R      R     R     R    R      R     R    R I     I  I I I )" << std::endl
      << "LOWER_BOUND      " << "(  1.0  1.0  20.0  1.0  1.0    1  1.0  0.0  1.0 793.0  1.0  1.0 0.01 0.01 495.0  1 0.01 0.0050 0.006 0.007  0.5 0.0050 0.006 0.15 2     1  1 1 1 )" << std::endl
      << "X0               " << "(  9.0  9.0 150.0  6.0  8.0 1000 45.0  0.5  5.0 900.0  9.0  9.0 0.30 0.20 560.0 50 0.30 0.0165 0.018 0.017 10.0 0.0155 0.016 0.20 2 12000  1 2 2 )" << std::endl    
      << "UPPER_BOUND      " << "( 40.0 40.0 250.0 30.0 30.0    - 89.0 20.0 20.0 995.0 50.0 30.0 5.00 5.00 650.0  - 5.00 0.1000 0.100 0.200 10.0 0.1000 0.100 0.40 -     - 10 9 8 )" << std::endl;
}
