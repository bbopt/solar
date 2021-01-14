/*-------------------------------------------------------------------------------*/
/*  SOLAR - The solar thermal power plant simulator - version 0.0                */
/*  https://github.com/bbopt/solar                                               */
/*                                                                               */
/*  2021-01-11                                                                   */
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

// TOTO: le test dans toto genere des fichiers: a virer



// TOTO: Checker la violation des contraintes EB. Ca devrait faire planter la simulation des que possible (meme de facon artificielle)


// TOTO: mettre license LGPL dans le bon repertoire; Ou autre licence;
// TOTO: hardcoder les fichiers de data
// TOTO: tester le truc de Bastien (TODO3.eml)
// TOTO: arranger l'interface; que les neuf instances soient mieux identifiees
// TOTO: enlever les commentaires en francais du main
// TOTO: mettre les surrogates

#include "helpFunctions.hpp"
#include "Global.hpp"
#include <iostream>
#include <ostream>
#include <fstream>
#include <algorithm>

// Version:
const std::string VERSION = "0.1, 2021-01-11";

/*-----------------------------------------------------------*/
/*                       main function                       */
/*-----------------------------------------------------------*/
int main ( int argc , char ** argv ) {
  
  // create problem descriptions:
  std::vector<Problem> problems;
  create_problems ( problems );

  
  /*---------------------------------------------------------------*/
  /*  check arguments and call appropriate display/help functions  */
  /*---------------------------------------------------------------*/
  {
    // ./solar
    if ( argc == 1 ) {
      display_usage ( std::cout );
      return 1;
    }

    std::string arg1 = toupper(argv[1]);

    // solar -i
    if ( arg1 == "-I" ) {
      display_info ( std::cout , VERSION );
      return 1;
    }
    
    if ( argc != 2 && argc != 3 && argc != 4 ) {
      display_usage ( std::cout );
      return 1;
    }
  
    if ( argc == 2 ) {

      if ( arg1 == "-H" ) {

	// display help for all problems: <solar -h> :
	display_help ( std::cout , problems );
	return 1;
      }
    
      // display help for one problem: <solar pb_id> :
      display_help ( std::cout , problems, argv[1] );
      return 1;
    }

    // display help for one problem: <solar -h pb_id> :
    if ( toupper(argv[1]) == "-H" ) {
      display_help ( std::cout, problems, argv[2] );
      return 1;
    }

    // display help for one problem: <solar pb_id -v> or <solar pb_id -h> :
    if ( toupper(argv[2]) == "-V" || toupper(argv[2]) == "-H" ) {
      display_help ( std::cout, problems, argv[1] );
      return 1;
    }
  }

  /*--------------------------------------------------------*/
  /*  run simulation with argv[1]==pb_id and argv[2]=x.txt  */
  /*--------------------------------------------------------*/

  // TOTO: gerer les surrogates statiques
  
  std::string pb_id   = toupper(argv[1]);
  std::string x_file  = argv[2];
  bool        verbose = ( argc == 4 && toupper(argv[3])=="-V" ) ? true : false;

  // check the problem id:
  {
    const Problem * pb = find_problem ( problems, pb_id );
    if ( pb )
      pb_id = pb->get_pb_id();
    else {
      std::cerr << "The problem id \""<< pb_id << "\" is invalid" << std::endl;
      return 1;
    }
  }

  if ( verbose )
    std::cout << "RUN SIMULATION FOR pb_id=" << pb_id << " and x_file=" << x_file << std::endl;

  if ( verbose )
    display_help ( std::cout , problems, pb_id );

  // launch the simulation:
  bool simulation_completed = false, cnt_eval = false;
  
  try{
    Scenario scenario ( pb_id, x_file.c_str() );
    simulation_completed = scenario.simulate ( std::cout , cnt_eval , verbose );
    if ( verbose ) {
      std::cout << std::endl
		<< "simulation completed: " << simulation_completed << std::endl
		<< "count evaluation    : " << cnt_eval             << std::endl
		<< std::endl;
      // scenario.print();  // TOTO permettre en mode verbose
    }
  }
  catch ( invalid_argument &e ) {
    // problem with simulation (e.g. bad input): Simulation cannot be used by solver;
    // error message must be printed since outputs are not displayed.
    if ( verbose )
      std::cout << std::endl
		<< "simulation completed: " << simulation_completed << std::endl
		<< "count evaluation    : " << cnt_eval             << std::endl
		<< std::endl;
    std::cout << e.what() << endl;
    return 1;
  }
  // (controlled) problem with the simulation (e.g. a priori constraint violated): Simulation should be usable by solver;
  // some of the outputs are displayed with 1e20 values.
  catch ( logic_error &e ) {
    if ( verbose ) {
      std::cout << std::endl
		<< "simulation completed: " << simulation_completed << std::endl
		<< "count evaluation    : " << cnt_eval             << std::endl
		<< std::endl;
      std::cout << e.what() << endl; 
    }

    // simulation_completed may be equal to 'false' but since some outputs are valid,
    // the simulation can still be used by a solver: 'return 0' is then appropriate.
    return 0;
  }
  catch (...) {
    if ( verbose ) {
      std::cout << std::endl
		<< "simulation completed: " << simulation_completed << std::endl
		<< "count evaluation    : " << cnt_eval             << std::endl
		<< std::endl;
    }
    std::cerr << "Error: simulation was interrupted for an unknown reason" << endl;
    return 1;
  }
  
  return 0;
}
