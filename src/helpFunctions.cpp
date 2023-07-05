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
#include "helpFunctions.hpp"

/*-----------------------------------------------------------*/
/*                  display best known values                */
/*-----------------------------------------------------------*/
void display_best_solutions ( std::ostream & out ) {
  out << "\tSOLAR1 \t-902,503.692418" << std::endl
      << "\tSOLAR2 \t841,839.671915"  << std::endl
      << "\tSOLAR3 \t70,813,885.0684" << std::endl
      << "\tSOLAR4 \t108,197,236.146" << std::endl
      << "\tSOLAR5 \t-28.8817193932"  << std::endl
      << "\tSOLAR6 \t43,954,935.1836" << std::endl  
      << "\tSOLAR7 \t-4,972.88703862" << std::endl
      << "\tSOLAR10\t42.447789"       << std::endl;
}

/*-----------------------------------------------------------*/
/*                   display list of problems                */
/*-----------------------------------------------------------*/
void display_problems ( std::ostream & out , const std::vector<Problem> & problems ) {

  out << "\t#\t"<< std::setw(22) << "pb_id"
      << "\t"   << std::setw(40) << "obj.(f)"
      << "\t"   << std::setw(15) << "# of objectives(p)"
      << "\t"   << std::setw(15) << "# of var.(n)"
      << "\t"   << std::setw(15) << "# of constr.(m)\n\n";
 
  for ( size_t i = 0 ; i < problems.size() ; ++i ) {
    std::ostringstream oss;
    oss << "SOLAR" << i+1;
    out << "\t" << oss.str()
	<< "\t" << std::setw(22) << problems[i].get_pb_id()
	<< "\t" << std::setw(40) << problems[i].get_f_description()
      	<< "\t" << std::setw(15) << problems[i].get_p()
	<< "\t" << std::setw(15) << problems[i].get_n()
      	<< "\t" << std::setw(15) << problems[i].get_m()
	<< std::endl;
  }
}

/*-----------------------------------------------------------*/
/*                        display usage                      */
/*-----------------------------------------------------------*/
void display_usage ( std::ostream & out ) {
  out << std::endl
      << "Run SOLAR (basic)   : solar pb_id x.txt (add -v for verbose mode)" << std::endl
      << "Run SOLAR (advanced): solar pb_id x.txt -seed=S -fid=F -rep=R -v"  << std::endl
      << " pb_id: Problem instance: integer in {1, 2, ..., 10}"              << std::endl
      << "     S: Random seed     : integer >=0 or \"diff\"; Default=0"      << std::endl
      << "     F: Fidelity        : real in ]0;1]; Default=1.0 (truth)"      << std::endl
      << "     R: Replications    : integer >= 1 ; Default=1\n"              << std::endl
      << "Validation: solar -check (can take several minutes)"               << std::endl
      << "Help(1)   : solar -h"                                              << std::endl
      << "Help(2)   : solar -h pb_id"                                        << std::endl
      << "Info      : solar -i\n"                                            << std::endl;
}

/*-----------------------------------------------------------*/
/*                        display info                       */
/*-----------------------------------------------------------*/
void display_info ( std::ostream & out , const std::string & version ) {
  out << std::endl << "SOLAR, the solar thermal power plant simulator, version "
      << version << std::endl << std::endl
      << "Contributors: M. Diago, S. Le Digabel, M. Lemyre-Garneau, B. Talgorn;"
      << " GERAD and Polytechnique Montreal"               << std::endl << std::endl
      << "This code is distributed under the LGPL license" << std::endl
      << "https://github.com/bbopt/solar" << std::endl << std::endl
      << "Please report bugs to sebastien.le-digabel@polymtl.ca" << std::endl << std::endl;
}

/*-----------------------------------------------------------*/
/*              display (short) help for all problems        */
/*-----------------------------------------------------------*/
void display_help ( std::ostream & out , const std::vector<Problem> & problems ) {
  out << std::endl
      << "Run simulation: solar pb_id x.txt -seed=S -fid=F -rep=R -v (optional)\n\n"
      << " pb_id: Problem instance (see list of problems below)\n\n"
      << " x.txt: Input vector: Point at which the simulator is evaluated\n"
      << "        Values separated with spaces\n"
      << "        It is possible to specify several vectors: Use one line for each\n\n"
      << "    -v: Verbose option\n\n"
      << "     S: Random seed:\n"
      << "          Some SOLAR instances are stochastic. This parameter impacts the value of stochastic outputs\n"
      << "          The seed is a natural integer\n"
      << "          If SOLAR is run twice at the same point with the same seed, it will give the same outputs\n"
      << "          The default value is 0\n"
      << "          Use -seed=diff to let SOLAR use a different random seed each time\n"
      << "          The random number generator can be validated by running 'solar -check'\n\n"
      << "     F: Fidelity of the simulator\n"
      << "          Real value in ]0;1]\n"
      << "          Default: 1.0 (full fidelity), which corresponds to the \"true blackbox\", or the \"truth\"\n"
      << "          Any value in ]0;1[ corresponds to a \"static surrogate\" of the truth\n"
      << "          The execution time increases with the fidelity\n"
      << "          A good default static surrogate is -fid=0.5\n\n"
      << "     R: Number of replications\n"
      << "          Integer >= 1, default=1\n"
      << "          Number of times that the simulator is run at the same point\n"
      << "          Each replication uses a different random seed dependent on the -seed option\n"
      << "          The mean value of stochastic outputs is displayed\n"
      << "          It is not possible to use R>1 with deterministic instances\n"
      << "\nHelp for a problem: solar pb_id or solar -h pb_id" << std::endl << std::endl
      << "List of problems:" << std::endl << std::endl;
  display_problems       ( out , problems );

  out << std::endl
      << "Best known values for single-objective instances (one replication, full fidelity, default seed of zero):"
      << std::endl << std::endl;
  display_best_solutions ( out );
  out << std::endl;
}

/*-----------------------------------------------------------*/
/*           display (complete) help for one instance        */
/*-----------------------------------------------------------*/
void display_help ( std::ostream               & out      ,
		    const std::vector<Problem> & problems ,
		    const std::string          & pb_id      ) {
 
  out << "Display help for Instance \"" << pb_id << "\":" << std::endl << std::endl;
  
  const Problem * pb = find_problem ( problems, pb_id );

  if ( pb ) {
   
    out << "Instance: " << pb->get_pb_id() << " (solar" << pb->get_index() << ")"
	<< "\t"         << pb->get_f_description()
	<< "\tn="       << pb->get_n()
	<< "\tm="       << pb->get_m()
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

    // #10:
    else if ( pb->get_pb_id() == "MINCOST_UNCONSTRAINED" )
      print_minCost_unconstrained ( out );

    else {
      out << "Cannot find detailed help for this instance" << std::endl;
      return;
    }
  }
  else {
    out << "This problem instance does not exist" << std::endl;
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
      << "\tMaximum field surface area: 195 hectares\n"
      << "\tBudget: $50M\n"
      << "\tMust provide 100% of the demand requirement\n"
      << "\tFidelity cannot be changed (must be 100%)\n"
      << std::endl
      << "Objective (first output, stochastic)\n"
      << "\tMaximize the total solar energy concentrated on the receiver aperture through one day (kWh)\n"
      << std::endl
      << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\tx1: Heliostats length (m)                          : Real in [ 1; 40]\n"
      << "\t\tx2: Heliostats width  (m)                          : Real in [ 1; 40]\n"
      << "\t\tx3: Tower height      (m)                          : Real in [20;250]\n"
      << "\t\tx4: Receiver aperture height (m)                   : Real in [ 1; 30]\n"
      << "\t\tx5: Receiver aperture width  (m)                   : Real in [ 1; 30]\n"
      << "\t\tx6: Number of heliostats to fit in the field       : Integer >= 1\n"
      << "\t\tx7: Field angular width (deg)                      : Real in [1;89]\n"
      << "\t\tx8: Minimum distance from tower (% of tower height): Real in [0;20]\n"
      << "\t\tx9: Maximum distance from tower (% of tower height): Real in [1;20]\n"
      << std::endl
      << "Constraints (outputs 2 to 6 with format ci <= 0):\n"
      << "\tc1: Cost of plant <= budget=$50M\n"
      << "\tc2: Field surface area: A priori constraint: PI*x3*x3(x9*x9-x8*x8) * x7/180 <= 1.95e6\n"
      << "\tc3: Tower is at least twice as high as heliostats       : A priori, linear constraint: 2x1 <= x3\n"
      << "\tc4: Min. distance from tower <= Max. distance from tower: A priori, linear constraint:  x8 <= x9\n"
      << "\tc5: Check that x6 heliostats can fit in the field\n"
      << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 9 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $1"     << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(    R    R     R    R    R   I    R    R    R )" << std::endl
      << "\tLOWER_BOUND      " << "(  1.0  1.0  20.0  1.0  1.0   1  1.0  0.0  1.0 )" << std::endl
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
      << "\tMaximum field surface area: 400 hectares\n"
      << "\tBudget: $300M\n"
      << "\tMust provide 100% of the demand requirement\n"
      << std::endl
      << "Objective (first output, analytic)\n"
      << "\tMinimize total heliostats field surface to run a pre-determined powerplant (square meters)\n"
      << "\tObjective = PI*x3*x3(x9*x9-x8*x8) * x7/180\n"
      << std::endl
      << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\t x1: Heliostats length (m)                          : Real in [ 1; 40]\n"
      << "\t\t x2: Heliostats width  (m)                          : Real in [ 1; 40]\n"
      << "\t\t x3: Tower height      (m)                          : Real in [20;250]\n"
      << "\t\t x4: Receiver aperture height (m)                   : Real in [ 1; 30]\n"
      << "\t\t x5: Receiver aperture width  (m)                   : Real in [ 1; 30]\n"
      << "\t\t x6: Number of heliostats to fit in the field       : Integer >= 1\n"
      << "\t\t x7: Field angular width (deg)                      : Real in [1;89]\n"
      << "\t\t x8: Minimum distance from tower (% of tower height): Real in [0;20]\n"
      << "\t\t x9: Maximum distance from tower (% of tower height): Real in [1;20]\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx10: Receiver outlet temperature   (K): Real in [793;995]\n"
      << "\t\tx11: Receiver number of tubes         : Integer in {1,2,...,9424}\n"
      << "\t\tx12: Receiver insulation thickness (m): Real in [0.01 ;5  ]\n"
      << "\t\tx13: Receiver tubes inner diameter (m): Real in [0.005;0.1]\n"
      << "\t\tx14: Receiver tubes outer diameter (m): Real in [0.005;0.1]\n"   
      << std::endl
      << "Constraints (outputs 2 to 13 with format ci <= 0):\n"
      << "\t c1: Field surface area: A priori constraint: PI*x3*x3(x9*x9-x8*x8) * x7/180 <= 4e6 or: objective <= 4e6\n"
      << "\t c2: Compliance to demand (stochastic)\n"
      << "\t c3: Cost of plant <= budget=$300M\n"
      << "\t c4: Tower is at least twice as high as heliostats       : A priori, linear constraint: 2x1 <= x3\n"
      << "\t c5: Min. distance from tower <= Max. distance from tower: A priori, linear constraint:  x8 <= x9\n"
      << "\t c6: Check that x6 heliostats can fit in the field\n"
      << "\t c7: Pressure in receiver tubes <= yield pressure                  (stochastic)\n"
      << "\t c8: Molten salt melting point  <= hot storage lowest temperature  (stochastic)\n"
      << "\t c9: Molten salt melting point  <= cold storage lowest temperature (stochastic)\n"
    // c10 removed in version 0.5.4:
    // << "\tc10: Molten salt melting point  <= steam generator outlet temperature\n"
      << "\tc10: Receiver tubes inside diameter <= outside diameter: A priori, linear constraint:     x13 <= x14\n"
      << "\tc11: Number of tubes in receiver fit inside receiver   : A priori constraint        : x11*x14 <= x5*PI/2\n"    
      << "\tc12: Receiver outlet temperature >= steam turbine inlet temperature\n";
  
  out << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 14 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $2" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(    R    R     R    R    R    I    R    R     R     R    I    R      R      R )" << std::endl
      << "\tLOWER_BOUND      " << "(  1.0  1.0  20.0  1.0  1.0    1  1.0  0.0   1.0 793.0    1 0.01  0.005 0.0050 )" << std::endl
      << "\tX0               " << "( 11.0 11.0 140.0 10.0 10.0 2650 89.0  0.5   5.0 838.0   36 0.30  0.020 0.0216 )" << std::endl
      << "\tUPPER_BOUND      " << "( 40.0 40.0 250.0 30.0 30.0    - 89.0 20.0  20.0 995.0 9424 5.00  0.100 0.1000 )" << std::endl;
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
      << "\tDemand profile: 10MW, starting at 3PM and ending at 9PM, 3 consecutive days\n"
      << "\tMaximum field surface area: 80 hectares\n"
      << "\tMust provide 100% of the demand requirement\n"
      << std::endl
      << "Objective (first output)\n"
      << "\tMinimize total investment cost ($)\n"
      << std::endl
      << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\t x1: Heliostats length (m)                          : Real in [ 1; 40]\n"
      << "\t\t x2: Heliostats width  (m)                          : Real in [ 1; 40]\n"
      << "\t\t x3: Tower height      (m)                          : Real in [20;250]\n"
      << "\t\t x4: Receiver aperture height (m)                   : Real in [ 1; 30]\n"
      << "\t\t x5: Receiver aperture width  (m)                   : Real in [ 1; 30]\n"
      << "\t\t x6: Number of heliostats to fit in the field       : Integer >= 1\n"
      << "\t\t x7: Field angular width (deg)                      : Real in [1;89]\n"
      << "\t\t x8: Minimum distance from tower (% of tower height): Real in [0;20]\n"
      << "\t\t x9: Maximum distance from tower (% of tower height): Real in [1;20]\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx10: Receiver outlet temperature (K)      : Real in [793;995]\n"
      << "\t\tx11: Hot storage height   (m)             : Real in [1;50]\n"
      << "\t\tx12: Hot storage diameter (m)             : Real in [1;30]\n"
      << "\t\tx13: Hot storage insulation thickness  (m): Real in [0.01;5]\n"
      << "\t\tx14: Cold storage insulation thickness (m): Real in [0.01;5]\n"
      << "\t\tx15: Mininum cold storage temperature  (K): Real in [495;650]\n"
      << "\t\tx16: Receiver number of tubes             : Integer in {1,2,...,9424}\n"
      << "\t\tx17: Receiver insulation thickness (m)    : Real in [0.01 ;5  ]\n"
      << "\t\tx18: Receiver tubes inner diameter (m)    : Real in [0.005;0.1]\n"
      << "\t\tx19: Receiver tubes outer diameter (m)    : Real in [0.005;0.1]\n"
      << "\tPowerblock:\n"
      << "\t\tx20: Type of turbine: Categorical: Integer in {1, 2, ..., 8}\n"
      << std::endl
      << "Constraints (outputs 2 to 14 with format ci <= 0):\n"
      << "\t c1: Field surface area: A priori constraint: PI*x3*x3(x9*x9-x8*x8) * x7/180 <= 800000\n"
      << "\t c2: Compliance to demand (stochastic)\n"
      << "\t c3: Tower is at least twice as high as heliostats       : A priori, linear constraint: 2x1 <= x3\n"
      << "\t c4: Min. distance from tower <= max. distance from tower: A priori, linear constraint:  x8 <= x9\n"
      << "\t c5: Check that x6 heliostats can fit in the field\n"
      << "\t c6: Pressure in receiver tubes <= yield pressure (stochastic)\n"
      << "\t c7: Molten salt melting point  <= hot storage lowest temperature  (stochastic)\n"
      << "\t c8: Molten salt melting point  <= cold storage lowest temperature (stochastic)\n"
    // c9 should be stochastic but this behavior was never observed
      << "\t c9: Molten salt melting point  <= steam generator outlet temperature\n"   
      << "\tc10: Receiver tubes inside diameter <= outside diameter: A priori, linear constraint:     x18 <= x19\n"
      << "\tc11: Number of tubes in receiver fit inside receiver   : A priori constraint        : x16*x19 <= x5*PI/2\n"
      << "\tc12: Receiver outlet temperature >= steam turbine inlet temperature\n"
      << "\tc13: Storage is back at least at its original conditions (stochastic)\n"
      << "\n-----------------------------------------------------------------\n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 20       << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $3" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(    R    R     R    R    R   I    R    R    R     R    R    R    R    R     R    I    R     R     R I )" << std::endl
      << "\tLOWER_BOUND      " << "(  1.0  1.0  20.0  1.0  1.0   1  1.0  0.0  1.0 793.0  1.0  1.0 0.01 0.01 495.0    1 0.01 0.005 0.005 1 )" << std::endl
      << "\tX0               " << "(  8.0  8.0 150.0  7.0  7.0 250 45.0  0.5  5.0 900.0  9.0  9.0 0.30 0.20 560.0   40 0.30 0.015 0.017 3 )" << std::endl
      << "\tUPPER_BOUND      " << "( 40.0 40.0 250.0 30.0 30.0   - 89.0 20.0 20.0 995.0 50.0 30.0 5.00 5.00 650.0 9424 5.00 0.100 0.100 8 )" << std::endl;
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
      << "\tDemand profile: 25MW, starting at 3PM and ending at 9PM, 3 consecutive days\n"
      << "\tMaximum field surface area: 200 hectares\n"
      << "\tMust provide 100% of the demand requirement\n"
      << std::endl

      << "Objective (first output)\n"
      << "\tMinimize the cost of powerplant to respect a given demand with a limited size of field ($)\n"
      << std::endl  

      << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\t x1: Heliostats length (m)                          : Real in [ 1; 40]\n"
      << "\t\t x2: Heliostats width  (m)                          : Real in [ 1; 40]\n"
      << "\t\t x3: Tower height      (m)                          : Real in [20;250]\n"
      << "\t\t x4: Receiver aperture height (m)                   : Real in [ 1; 30]\n"
      << "\t\t x5: Receiver aperture width  (m)                   : Real in [ 1; 30]\n"
      << "\t\t x6: Number of heliostats to fit in the field       : Integer >= 1\n"
      << "\t\t x7: Field angular width (deg)                      : Real in [1;89]\n"
      << "\t\t x8: Minimum distance from tower (% of tower height): Real in [0;20]\n"
      << "\t\t x9: Maximum distance from tower (% of tower height): Real in [1;20]\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx10: Receiver outlet temperature (K)      : Real in [793;995]\n"
      << "\t\tx11: Hot storage height   (m)             : Real in [1;50]\n"
      << "\t\tx12: Hot storage diameter (m)             : Real in [1;30]\n"
      << "\t\tx13: Hot storage insulation thickness  (m): Real in [0.01;5]\n"
      << "\t\tx14: Cold storage insulation thickness (m): Real in [0.01;5]\n"
      << "\t\tx15: Mininum cold storage temperature  (K): Real in [495;650]\n"
      << "\t\tx16: Receiver number of tubes             : Integer in {1,2,...,7853}\n"
      << "\t\tx17: Receiver insulation thickness     (m): Real in [0.01 ;5  ]\n"
      << "\t\tx18: Receiver tubes inner diameter     (m): Real in [0.005;0.1]\n"
      << "\t\tx19: Receiver tubes outer diameter     (m): Real in [0.006;0.1]\n"
      << "\tSteam generator:\n"
      << "\t\tx20: Tubes spacing (m)       : Real in [0.007;0.2]\n"
      << "\t\tx21: Tubes length  (m)       : Real in [0.5;10]\n"
      << "\t\tx22: Tubes inner diameter (m): Real in [0.005;0.1]\n"
      << "\t\tx23: Tubes outer diameter (m): Real in [0.006;0.1]\n"
      << "\t\tx24: Baffles cut             : Real in [0.15;0.4]\n"
      << "\t\tx25: Number of baffles       : Integer >= 2\n"   
      << "\t\tx26: Number of tubes         : Integer >= 1\n"
      << "\t\tx27: Number of shell passes  : Integer in {1, 2, ..., 10}\n"
      << "\t\tx28: Number of tubes passes  : Integer in {1, 2, ..., 9}\n"   
      << "\tPowerblock:\n"
      << "\t\tx29: Type of turbine: Categorical: Integer in {1, 2, ..., 8}\n"
      << std::endl
      << "Constraints (outputs 2 to 17 with format ci <= 0):\n"
      << "\t c1: Field surface area: A priori constraint: PI*x3*x3(x9*x9-x8*x8) * x7/180 <= 2e6\n"
      << "\t c2: Compliance to demand (stochastic)\n"
      << "\t c3: Tower is at least twice as high as heliostats       : A priori, linear constraint: 2x1 <= x3\n"
      << "\t c4: Min. distance from tower <= max. distance from tower: A priori, linear constraint:  x8 <= x9\n"   
      << "\t c5: Check that x6 heliostats can fit in the field\n"
      << "\t c6: Pressure in receiver tubes <= yield pressure (stochastic)\n"
      << "\t c7: Molten salt melting point  <= hot storage lowest temperature     (stochastic)\n"  
      << "\t c8: Molten salt melting point  <= cold storage lowest temperature    (stochastic)\n"
      << "\t c9: Molten salt melting point  <= steam generator outlet temperature (stochastic)\n"   
      << "\tc10: Receiver tubes inside diameter <= outside diameter: A priori, linear constraint: x18 <= x19\n"
      << "\tc11: Number of tubes in receiver fit inside receiver: A priori constraint: x16*x19 <= x5*PI/2\n"
      << "\tc12: Receiver outlet temperature >= steam turbine inlet temperature\n"
      << "\tc13: Parasitic losses <= 18% of the generated output (stochastic)\n"
      << "\tc14: Steam generator tubes outer diameter  <= tubes spacing: A priori, linear constraint: x23 <= x20\n"
      << "\tc15: Steam generator tubes inside diameter <= Steam generator tubes outer diameter: A priori, linear constraint: x22 <= x23\n"
      << "\tc16: Pressure in steam generator tubes <= yield pressure\n"

      << "\n-----------------------------------------------------------------\n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 29 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $4" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(    R    R     R    R    R    I    R    R    R     R    R    R    R    R     R    I    R      R     R     R    R      R     R    R I     I  I I I )" << std::endl
      << "\tLOWER_BOUND      " << "(  1.0  1.0  20.0  1.0  1.0    1  1.0  0.0  1.0 793.0  1.0  1.0 0.01 0.01 495.0    1 0.01 0.0050 0.006 0.007  0.5 0.0050 0.006 0.15 2     1  1 1 1 )" << std::endl
      << "\tX0               " << "(  9.0  9.0 150.0  6.0  8.0 1000 45.0  0.5  5.0 900.0  9.0  9.0 0.30 0.20 560.0  500 0.30 0.0165 0.018 0.017 10.0 0.0155 0.016 0.20 3 12000  1 2 2 )" << std::endl  
      << "\tUPPER_BOUND      " << "( 40.0 40.0 250.0 30.0 30.0    - 89.0 20.0 20.0 995.0 50.0 30.0 5.00 5.00 650.0 7853 5.00 0.1000 0.100 0.200 10.0 0.1000 0.100 0.40 -     - 10 9 8 )" << std::endl;
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
      << "\tDemand: 12MW\n"
      << "\tBudget: $100M\n"
      << "\tPre-determined sunlight input for a period of 1 month\n"
      << "\tNumber of heliostats to fit in the field: 3,800\n"
      << "\tDeterministic instance\n"
      << "\tFidelity cannot be changed (must be 100%)\n"
      << std::endl
      << "Objective (first output)\n"
      << "\tMaximize compliance to a demand profile\n"
      << std::endl  
      << "Variables:\n"
      << "\tHeat transfer loop:\n"
      << "\t\t x1: Receiver outlet temperature (K)      : Real in [793;995]\n"
      << "\t\t x2: Hot storage height   (m)             : Real in [1;30]\n"
      << "\t\t x3: Hot storage diameter (m)             : Real in [1;30]\n"
      << "\t\t x4: Hot storage insulation thickness  (m): Real in [0.01;2]\n"
      << "\t\t x5: Cold storage insulation thickness (m): Real in [0.01;2]\n"
      << "\t\t x6: Mininum cold storage temperature  (K): Real in [495;650]\n"
      << "\t\t x7: Receiver number of tubes             : Integer in {1,2,...,1884}\n"
      << "\t\t x8: Receiver insulation thickness     (m): Real in [0.100;2.0]\n"
      << "\t\t x9: Receiver tubes inner diameter     (m): Real in [0.005;0.1]\n"
      << "\t\tx10: Receiver tubes outer diameter     (m): Real in [0.005;0.1]\n"
      << "\tSteam generator:\n"
      << "\t\tx11: Tubes spacing (m)       : Real in [0.006; 0.2]\n"
      << "\t\tx12: Tubes length  (m)       : Real in [0.500;10.0]\n"
      << "\t\tx13: Tubes inner diameter (m): Real in [0.005; 0.1]\n"
      << "\t\tx14: Tubes outer diameter (m): Real in [0.006; 0.1]\n"
      << "\t\tx15: Baffles cut             : Real in [0.150; 0.4]\n"
      << "\t\tx16: Number of baffles       : Integer >= 2\n"
      << "\t\tx17: Number of tubes         : Integer >= 1\n"
      << "\t\tx18: Number of shell passes  : Integer in {1, 2, ...,10}\n"
      << "\t\tx19: Number of tube passes   : Integer in {1, 2, ..., 9}\n"
      << "\tPowerblock:\n"
      << "\t\tx20: Type of turbine: Categorical: Integer in {1, 2, ..., 8}\n"
      << std::endl
      << "Constraints (outputs 2 to 13 with format ci <= 0):\n"
      << "\t c1: Cost of plant <= budget=$100M\n"
      << "\t c2: Pressure in receiver tubes <= yield pressure\n"
      << "\t c3: Molten salt melting point  <= hot storage lowest temperature\n"
      << "\t c4: Molten salt melting point  <= cold storage lowest temperature\n"
      << "\t c5: Molten salt melting point  <= steam generator outlet temperature\n"
      << "\t c6: Receiver tubes inside diameter <= outside diameter: A priori, linear constraint: x9 <= x10\n"
      << "\t c7: Number of tubes in receiver fit inside receiver: A priori constraint: x7*x10 <= 3*PI\n"
      << "\t c8: Receiver outlet temperature >= steam turbine inlet temperature\n"
      << "\t c9: Parasitic losses <= 18% of the generated output \n"
      << "\tc10: Steam generator tubes outer diameter  <= tubes spacing: A priori, linear constraint: x14 <= x11\n"
      << "\tc11: Steam generator tubes inside diameter <= Steam generator tubes outer diameter: A priori, linear constraint: x13 <= x14\n"
      << "\tc12: Pressure in steam generator tubes <= yield pressure\n";

  out << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 20 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $5" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(     R    R    R    R    R   R    I    R     R     R      R    R       R     R    R I    I  I I I )" << std::endl
      << "\tLOWER_BOUND      " << "( 793.0  1.0  1.0 0.01 0.01 495    1 0.10 0.005 0.005  0.006  0.5   0.005 0.006 0.15 2    1  1 1 1 )" << std::endl
      << "\tX0               " << "( 900.0 10.0 12.0 0.15 0.10 560   24 0.35 0.020 0.023  0.050  8.0   0.020 0.023 0.20 2 5000  5 5 1 )" << std::endl
      << "\tUPPER_BOUND      " << "( 995.0 30.0 30.0 2.00 2.00 650 1884 2.00 0.100 0.100  0.200 10.0   0.100 0.100 0.4  -    - 10 9 8 )" << std::endl;
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
      << "\tDemand: 120MW\n"
      << "\tMust provide 100% of the demand requirement\n"
      << "\tNumber of heliostats to fit in the field: 12,232\n"  
      << "\tThis instance runs a predetermined power plant using the molten salt cycle and power block models.\n"
      << "\tThe objective is to minimize the cost of the thermal storage units so that the power plant is able\n"
      << "\tto sustain a 120MW electrical power output during 24 hours. Since the heliostat field is not being\n"
      << "\toptimized, its hourly power output is taken from prerecorded data instead of being simulated\n"
      << "\tDeterministic instance\n"
      << "\tFidelity cannot be changed (must be 100%)\n"
      << std::endl;

  out << "Objective (first output)\n"
      << "\tMinimize the cost of storage ($)\n"
      << std::endl;

  out << "Variables:\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx1: Receiver outlet temperature (K)      : Real in [793;995]\n"
      << "\t\tx2: Hot storage height   (m)             : Real in [2;50]\n"
      << "\t\tx3: Hot storage diameter (m)             : Real in [2;30]\n"
      << "\t\tx4: Hot storage insulation thickness  (m): Real in [0.01;5]\n"
      << "\t\tx5: Cold storage insulation thickness (m): Real in [0.01;5]\n"
      << std::endl
      << "Constraints (outputs 2 to 7 with format ci <= 0):\n"
      << "\tc1: Compliance to demand\n"
      << "\tc2: Pressure in receiver tubes  <= yield pressure\n"
      << "\tc3: Molten salt melting point   <= hot storage lowest temperature\n"
      << "\tc4: Molten salt melting point   <= cold storage lowest temperature\n"
      << "\tc5: Receiver outlet temperature >= steam turbine inlet temperature\n"
      << "\tc6: Storage is back at least at its original conditions\n"
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
      << "\tBudget: $45M\n"
      << "\tNumber of heliostats to fit in the field: 5,000\n"  
      << std::endl;
  
  out << "Objective (first output, stochastic)\n"
      << "\tMaximize receiver efficiency, i.e the energy transferred to the molten salt (J)\n"
      << std::endl;

  out << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\tx1: Receiver aperture height (m): Real in [ 1; 30]\n"
      << "\t\tx2: Receiver aperture width  (m): Real in [ 1; 30]\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx3: Receiver outlet temperature   (K): Real in [793;995]\n"
      << "\t\tx4: Receiver number of tubes         : Integer in {1,2,...,8567}\n"
      << "\t\tx5: Receiver insulation thickness (m): Real in [0.01  ;5.0]\n"
      << "\t\tx6: Receiver tubes inner diameter (m): Real in [0.005 ;0.1]\n"
      << "\t\tx7: Receiver tubes outer diameter (m): Real in [0.0055;0.1]\n"
      << std::endl
      << "Constraints (outputs 2 to 7 with format ci <= 0):\n"
      << "\tc1: Cost of plant <= budget=$45M\n"
      << "\tc2: Pressure in receiver tubes     <= yield pressure (stochastic)\n"
      << "\tc3: Receiver tubes inside diameter <= outside diameter: A priori, linear constraint: x6 <= x7\n"
      << "\tc4: Receiver outlet temperature    >= steam turbine inlet temperature\n"
      << "\tc5: Number of tubes in receiver fit inside receiver: A priori constraint: x4*x7 <= x2*PI/2\n"
      << "\tc6: Parasitic losses <= 3% of the generated output  (stochastic)\n"   
      << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 7 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $7" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(    R    R     R    I    R     R      R )" << std::endl
      << "\tLOWER_BOUND      " << "(  1.0  1.0 793.0    1 0.01 0.005 0.0055 )" << std::endl
      << "\tX0               " << "(  7.0  7.0 850.0   40 0.20 0.010 0.0110 )" << std::endl
      << "\tUPPER_BOUND      " << "( 30.0 30.0 995.0 8567 5.00 0.100 0.1000 )" << std::endl;
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
      << "\tMaximum field surface area: 400 hectares\n"
      << "\tMinimum energy production: 400MWh\n"
      << std::endl;
  
  out << "Objectives (first and second outputs; first objective is stochastic)\n"
      << "\tMaximize heliostat field performance (absorbed energy in Joules) and minimize cost of field, tower and receiver ($)\n"
      << std::endl;

  out << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\t x1: Heliostats length (m)                          : Real in [ 1; 40]\n"
      << "\t\t x2: Heliostats width  (m)                          : Real in [ 1; 40]\n"
      << "\t\t x3: Tower height      (m)                          : Real in [20;250]\n"
      << "\t\t x4: Receiver aperture height (m)                   : Real in [ 1; 30]\n"
      << "\t\t x5: Receiver aperture width  (m)                   : Real in [ 1; 30]\n"
      << "\t\t x6: Number of heliostats to fit in the field       : Integer >= 1\n"
      << "\t\t x7: Field angular width (deg)                      : Real in [1;89]\n"
      << "\t\t x8: Minimum distance from tower (% of tower height): Real in [0;20]\n"
      << "\t\t x9: Maximum distance from tower (% of tower height): Real in [1;20]\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx10: Receiver number of tubes         : Integer in {1,2,...,7853}\n"
      << "\t\tx11: Receiver insulation thickness (m): Real in [0.01 ;5  ]\n"
      << "\t\tx12: Receiver tubes inner diameter (m): Real in [0.005;0.1]\n"
      << "\t\tx13: Receiver tubes outer diameter (m): Real in [0.006;0.1]\n"
      << std::endl
      << "Constraints (outputs 3 to 11 with format ci <= 0):\n"
      << "\tc1: Field surface area: A priori constraint: PI*x3*x3(x9*x9-x8*x8) * x7/180 <= 4e6\n"
      << "\tc2: Tower is at least twice as high as heliostats       : A priori, linear constraint: 2x1 <= x3\n"
      << "\tc3: Min. distance from tower <= max. distance from tower: A priori, linear constraint:  x8 <= x9\n"   
      << "\tc4: Check that x6 heliostats can fit in the field\n"
      << "\tc5: Pressure in receiver tubes <= yield pressure (stochastic)\n"
      << "\tc6: Receiver tubes inside diameter <= outside diameter: A priori, linear constraint: x12 <= x13\n"
      << "\tc7: Number of tubes in receiver fit inside receiver: A priori constraint: x10*x13 <= x5*PI/2\n"
      << "\tc8: Minimal acceptable energy production (lower bound on the first objective, stochastic)\n"
      << "\tc9: Parasitic losses <= 8% of the generated output  (stochastic)\n"
  
      << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 13 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $8" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ OBJ CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(    R    R     R    R    R    I    R    R    R    I    R     R      R )" << std::endl
      << "\tLOWER_BOUND      " << "(  1.0  1.0  20.0  1.0  1.0    1  1.0  0.0  1.0    1 0.01 0.005 0.0060 )" << std::endl
      << "\tX0               " << "( 11.0 11.0 200.0 10.0 10.0 2650 89.0  0.5  8.0   36 0.30 0.020 0.0216 )" << std::endl
      << "\tUPPER_BOUND      " << "( 40.0 40.0 250.0 30.0 30.0    - 89.0 20.0 20.0 7853 5.00 0.100 0.1000 )" << std::endl;
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
      << "\tMaximum field surface area: 500 hectares\n"
      << "\tBudget: $1.2B\n"
      << "\tMinimum energy production: 120MWh\n"
      << std::endl;
  
  out << "Objectives (first and second outputs; the two objectives are stochastic)\n"
      << "\tMaximize power (Wh) and minimize losses (Wh)\n"
      << std::endl;
  
  out << "Variables:\n"
      << "\tHeliostats Field:\n"
      << "\t\t x1: Heliostats length (m)                          : Real in [ 1; 40]\n"
      << "\t\t x2: Heliostats width  (m)                          : Real in [ 1; 40]\n"
      << "\t\t x3: Tower height      (m)                          : Real in [20;250]\n"
      << "\t\t x4: Receiver aperture height (m)                   : Real in [ 1; 30]\n"
      << "\t\t x5: Receiver aperture width  (m)                   : Real in [ 1; 30]\n"
      << "\t\t x6: Number of heliostats to fit in the field       : Integer >= 1\n"
      << "\t\t x7: Field angular width (deg)                      : Real in [1;89]\n"
      << "\t\t x8: Minimum distance from tower (% of tower height): Real in [0;20]\n"
      << "\t\t x9: Maximum distance from tower (% of tower height): Real in [1;20]\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx10: Receiver outlet temperature (K)      : Real in [793;995]\n"
      << "\t\tx11: Hot storage height   (m)             : Real in [1;50]\n"
      << "\t\tx12: Hot storage diameter (m)             : Real in [1;30]\n"
      << "\t\tx13: Hot storage insulation thickness  (m): Real in [0.01;5]\n"
      << "\t\tx14: Cold storage insulation thickness (m): Real in [0.01;5]\n"
      << "\t\tx15: Mininum cold storage temperature  (K): Real in [495;650]\n"
      << "\t\tx16: Receiver number of tubes             : Integer in {1,2,...,7853}\n"
      << "\t\tx17: Receiver insulation thickness (m)    : Real in [0.01 ;5  ]\n"
      << "\t\tx18: Receiver tubes inner diameter (m)    : Real in [0.005;0.1]\n"
      << "\t\tx19: Receiver tubes outer diameter (m)    : Real in [0.006;0.1]\n"
      << "\tSteam generator:\n"
      << "\t\tx20: Tubes spacing (m)       : Real in [0.007;0.2]\n"
      << "\t\tx21: Tubes length  (m)       : Real in [0.5;10]\n"
      << "\t\tx22: Tubes inner diameter (m): Real in [0.005;0.1]\n"
      << "\t\tx23: Tubes outer diameter (m): Real in [0.006;0.1]\n"
      << "\t\tx24: Baffles cut             : Real in [0.15;0.4]\n"
      << "\t\tx25: Number of baffles       : Integer >= 2\n"   
      << "\t\tx26: Number of tubes         : Integer >= 1\n"
      << "\t\tx27: Number of shell passes  : Integer in {1, 2, ..., 10}\n"
      << "\t\tx28: Number of tube passes   : Integer in {1, 2, ..., 9}\n"   
      << "\tPowerblock:\n"
      << "\t\tx29: Type of turbine: Categorical: Integer in {1, 2, ..., 8}\n"
      << std::endl
      << "Constraints (outputs 3 to 19 with format ci <= 0):\n"
      << "\t c1: Cost of plant <= budget=$1.2B\n"
      << "\t c2: Minimal acceptable energy production (stochastic)\n"
      << "\t c3: Field surface area: A priori constraint: PI*x3*x3(x9*x9-x8*x8) * x7/180 <= 5e6\n"
      << "\t c4: Tower is at least twice as high as heliostats       : A priori, linear constraint: 2x1 <= x3\n"
      << "\t c5: Min. distance from tower <= max. distance from tower: A priori, linear constraint:  x8 <= x9\n"  
      << "\t c6: Check that x6 heliostats can fit in the field\n"
      << "\t c7: Pressure in receiver tubes <= yield pressure                     (stochastic)\n"
      << "\t c8: Molten salt melting point  <= hot storage lowest temperature     (stochastic)\n"
      << "\t c9: Molten salt melting point  <= cold storage lowest temperature    (stochastic)\n"
      << "\tc10: Molten salt melting point  <= steam generator outlet temperature (stochastic)\n"
      << "\tc11: Receiver tubes inside diameter <= outside diameter: A priori, linear constraint: x18 <= x19\n"
      << "\tc12: Number of tubes in receiver fit inside receiver: A priori constraint:  x16*x19 <= x5*PI/2\n"
      << "\tc13: Receiver outlet temperature >= steam turbine inlet temperature\n"
      << "\tc14: Parasitic losses <= 20% of the generated output  (stochastic)\n"
      << "\tc15: Steam generator tubes outer diameter  <= tubes spacing: A priori, linear constraint: x23 <= x20\n"
      << "\tc16: Steam generator tubes inside diameter <= Steam generator tubes outer diameter: A priori, linear constraint: x22 <= x23\n"
      << "\tc17: Pressure in steam generator tubes <= yield pressure\n";

  out << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 29 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $9" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ OBJ CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR CSTR" << std::endl
      << "\tBB_INPUT_TYPE    " << "(    R    R     R    R    R    I    R    R    R     R    R    R    R    R     R    I    R      R     R     R    R      R     R    R I     I  I I I )" << std::endl
      << "\tLOWER_BOUND      " << "(  1.0  1.0  20.0  1.0  1.0    1  1.0  0.0  1.0 793.0  1.0  1.0 0.01 0.01 495.0    1 0.01 0.0050 0.006 0.007  0.5 0.0050 0.006 0.15 2     1  1 1 1 )" << std::endl
      << "\tX0               " << "(  9.0  9.0 150.0  6.0  8.0 1000 45.0  0.5  5.0 900.0  9.0  9.0 0.30 0.20 560.0  500 0.30 0.0165 0.018 0.017 10.0 0.0155 0.016 0.20 3 12000  1 2 2 )" << std::endl   
      << "\tUPPER_BOUND      " << "( 40.0 40.0 250.0 30.0 30.0    - 89.0 20.0 20.0 995.0 50.0 30.0 5.00 5.00 650.0 7853 5.00 0.1000 0.100 0.200 10.0 0.1000 0.100 0.40 -     - 10 9 8 )" << std::endl;
}

/*--------------  #10 -------------------------------------*/
void print_minCost_unconstrained ( std::ostream & out ) {
/*---------------------------------------------------------*/

  out << "\n-----------------------------------------------------------------\n"
      << "Parameters:\n"
      << "\tWhole plant\n"
      << "\tLatitude: 30.05 deg N\n"
      << "\tDay: January 1st\n"  // https://www.epochconverter.com/days/2019
      << "\tDuration: 24 hours\n"
      << "\tDemand: 120MW\n"
      << "\tMust provide 100% of the demand requirement\n"
      << "\tNumber of heliostats to fit in the field: 12,232\n"
      << "\tDeterministic instance\n"
      << "\tThis instance is the unconstrained version of instance #6\n"
      << std::endl;

  out << "Objective (first output)\n"
      << "\tMinimize the cost of storage + penalties on the 6 constraints of Instance #6:\n"
      << "\tUnconstrained objective = 1E-6 f + 0.5 ( (c1+)^2 + (2E-6 c2+)^2 + (c3+)^2 + (c4+)^2 + (c5+)^2 + (c6+)^2 )\n"
      << "\twith f, c1, c2, ..., c6 the outputs of Instance #6 and cj+ = max{0,cj}, j=1,2,...,6\n"
      << std::endl;

  out << "Variables:\n"
      << "\tHeat transfer loop:\n"
      << "\t\tx1: Receiver outlet temperature (K)      : Real in [793;995]\n"
      << "\t\tx2: Hot storage height   (m)             : Real in [2;50]\n"
      << "\t\tx3: Hot storage diameter (m)             : Real in [2;30]\n"
      << "\t\tx4: Hot storage insulation thickness  (m): Real in [0.01;5]\n"
      << "\t\tx5: Cold storage insulation thickness (m): Real in [0.01;5]\n"
      << std::endl
      << "\n----------------------------------------------------------------- \n"
      << "NOMAD parameters:\n\n"
      << "\tDIMENSION        " << 5 << std::endl
      << "\tBB_EXE           " << "$SOLAR_HOME/bin/solar $10" << std::endl
      << "\tBB_OUTPUT_TYPE   " << "OBJ" << std::endl
      << "\tBB_INPUT_TYPE    " << "(     R    R    R    R    R )" << std::endl
      << "\tLOWER_BOUND      " << "( 793.0  2.0  2.0 0.01 0.01 )" << std::endl
      << "\tX0               " << "( 900.0 10.0 12.0 0.20 0.20 )" << std::endl
      << "\tUPPER_BOUND      " << "( 995.0 50.0 30.0 5.00 5.00 )" << std::endl;
}
