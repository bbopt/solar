/*-------------------------------------------------------------------------------*/
/*  SOLAR - The solar thermal power plant simulator - version 0.5.7              */
/*  https://github.com/bbopt/solar                                               */
/*                                                                               */
/*  2023-07-05                                                                   */
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
/*                                                                               */
/* Best known values for single-objective instances                              */
/* (one replication, full fidelity, default seed of zero):                       */
/*                                                                               */
/* SOLAR1    -902,503.692418                                                     */
/* SOLAR2     841,839.671915                                                     */
/* SOLAR3  70,813,885.0684                                                       */
/* SOLAR4 108,197,236.146                                                        */
/* SOLAR5         -28.8817193932                                                 */
/* SOLAR6  43,954,935.1836                                                       */
/* SOLAR7      -4,972.88703862                                                   */
/* SOLAR10         42.447789                                                     */
/*                                                                               */
/*-------------------------------------------------------------------------------*/
#include "Evaluator.hpp"
#include "sampling.hpp"

// version:
const std::string VERSION = "0.5.7, 2023-07-05";

// validation functions:
bool check ( bool fast );

bool check_eval ( const std::string & pb_id        ,
		  int                 seed         ,
		  double              fidelity     ,
		  int                 replications ,
		  const double      * x            ,
		  const std::string & expected     ,
		  std::string       & error          );

// to parse and get the options:
bool get_options ( int           argc         ,
		   char       ** argv         ,
		   std::string & x_file       ,
 		   bool        & verbose      ,
		   int         & seed         ,
		   double      & fidelity     ,
		   int         & replications   );

/*------------------------------------------------------------------------------*/
/*                                  main function                               */
/*------------------------------------------------------------------------------*/
/* + another version of the main function is available after this one           */
/* + it is called main_minimal and is a minimal example of a single evaluation  */
/*------------------------------------------------------------------------------*/
int main ( int argc , char ** argv ) {
  
  // create problem descriptions:
  std::vector<Problem> problems;
  create_problems ( problems );
  
  /*---------------------------------------------------------------*/
  /*  check arguments and call appropriate display/help functions  */
  /*---------------------------------------------------------------*/
  {
    // ./solar
    if ( argc == 1 || argc > 7 ) {
      display_usage ( std::cout );
      return 1;
    }

    std::string arg1 = toupper(argv[1]);

    // <solar -c> or <solar -check>: call validation function:
    if ( arg1 == "-C" || arg1 == "-CHECK" || arg1 == "--C" || arg1 == "--CHECK" )
      return check(false) ? 0 : 1;
    
    // <solar -i> or <solar -v> or <solar -info>:
    if ( arg1 ==  "-I" || arg1 ==  "-V" || arg1 ==  "-INFO" ||
	 arg1 == "--I" || arg1 == "--V" || arg1 == "--INFO"    ) {
      display_info ( std::cout , VERSION );
      return 1;
    }
   
    if ( argc == 2 ) {

      if ( arg1 == "-H" || arg1 == "-HELP" || arg1 == "--H" || arg1 == "--HELP" ) {

	// display help for all problems: <solar -h> :
	display_help ( std::cout , problems );
	return 1;
      }
    
      // display help for one problem: <solar pb_id> :
      display_help ( std::cout , problems, argv[1] );
      return 1;
    }

    // display help for one problem: <solar -h pb_id> :
    if ( arg1 == "-H" || arg1 == "-HELP" || arg1 == "--H" || arg1 == "--HELP") {
      display_help ( std::cout, problems, argv[2] );
      return 1;
    }

    // display help for one problem: <solar pb_id -v> or <solar pb_id -h> :
    if ( toupper(argv[2]) ==  "-V" || toupper(argv[2]) ==  "-H" ||  toupper(argv[2]) ==  "-HELP" ||
	 toupper(argv[2]) == "--V" || toupper(argv[2]) == "--H" ||  toupper(argv[2]) == "--HELP"    ) {
      display_help ( std::cout, problems, argv[1] );
      return 1;
    }
  }

  /*--------------------------------------------------------*/
  /*  run simulation with argv[1]==pb_id and argv[2]=x.txt  */
  /*--------------------------------------------------------*/
 
  std::string pb_id  = toupper(argv[1]);
  
  // get the options:
  std::string x_file;
  bool        verbose      = false;
  int         seed         = 0;
  double      fidelity     = 1.0;
  int         replications = 1;

  if ( !get_options ( argc, argv, x_file, verbose, seed, fidelity, replications ) ) {
    display_usage ( std::cout );
    return 1;
  }
 
  // random seed:
  bool change_seed = false;
  if ( seed < 0 ) {
    change_seed = true;
    seed = RNG::get_pid();
  }
   
  if ( verbose ) {
    std::cout << std::endl
	      << "Seed        : " << seed;
    if ( change_seed )
      std::cout << " (diff)";
    std::cout << std::endl
	      << "Fidelity    : " << fidelity     << std::endl
	      << "Replications: " << replications << std::endl
	      << std::endl;
  }
   
  // check the problem id:
  const Problem * pb = find_problem ( problems, pb_id );
  if ( pb )
    pb_id = pb->get_pb_id();
  else {
    std::cerr << "The problem id \""<< pb_id << "\" is invalid" << std::endl;
    return 1;
  }

  // simulation ready to launch:
  if ( verbose ) {
    std::cout << "Run simulation for Problem " << pb_id << " and input file " << x_file << "\n\n";
    display_help ( std::cout , problems, pb_id );
  }

  bool simulation_completed = false, cnt_eval = false;

  Clock clock1;
  
  // creation of the evaluator:
  Evaluator evaluator ( *pb, std::cout );
 
  // read inputs:
  if ( !evaluator.read_x ( x_file ) ) {
    std::cerr << "Problem with input file \""<< x_file << "\"" << std::endl;
    return 1;
  }

  // For each point in input file:
  bool error = false;
  
  for ( int i = 0; i < evaluator.get_nb_input_points(); ++i ) {

    Clock       clock2;
    std::string err_msg;
   
    // random seed:
    if ( change_seed && i > 0 ) {
      seed += RNG::rand()%100;
      if ( verbose )
	std::cout << std::endl << "New seed: " << seed << std::endl;
    }
    
    // evaluation(s):
    if ( !evaluator.eval_x ( i, seed, fidelity, replications , simulation_completed, cnt_eval, err_msg, verbose ) ) {

      error = true;
      
      // displays:
      if ( verbose && !err_msg.empty()) {
	std::cout << err_msg << std::endl;
	std::cout << std::endl
		  << "Simulation completed: " << simulation_completed   << std::endl
		  << "Count evaluation    : " << cnt_eval               << std::endl
		  << "CPU time            : " << clock2.get_CPU_time()  << std::endl
		  << "real time           : " << clock2.get_real_time() << std::endl;
	if ( evaluator.has_intermediate_outputs() ) {
	  std::cout << "intermediate outputs: ";
	  evaluator.display_intermediate_outputs();
	  std::cout << std::endl;
	}
	std::cout << "Outputs             : ";
      }

      // should display only infinite values:
      evaluator.display_outputs();
      std::cout << std::endl;
      
      continue;
    }
  
    // displays:
    if ( verbose && !err_msg.empty() )
      std::cout << "\n" << err_msg << std::endl;
    
    if ( verbose ) {
      std::cout << std::endl
		<< "Simulation completed: " << simulation_completed   << std::endl
		<< "Count evaluation    : " << cnt_eval               << std::endl
		<< "CPU time            : " << clock2.get_CPU_time()  << std::endl
		<< "real time           : " << clock2.get_real_time() << std::endl;
      if ( evaluator.has_intermediate_outputs() ) {
	std::cout << "intermediate outputs: ";
	evaluator.display_intermediate_outputs();
	std::cout << std::endl;
      }
      std::cout << "Outputs             : ";
    }
    evaluator.display_outputs();
    std::cout << std::endl;
  }

  if ( verbose && evaluator.get_nb_input_points() > 1 )
    std::cout << std::endl
	      << "Total CPU time : " << clock1.get_CPU_time()  << std::endl
	      << "Total real time: " << clock1.get_real_time() << std::endl;
 
  return error ? 1 : 0;
}

/*------------------------------------------------------------------------*/
/*                   minimal example of a single evaluation               */
/*------------------------------------------------------------------------*/
/*  + this simple example is equivalent to <solar 10 x.txt>               */
/*  + it allows to retrieve the simulation_completed and cnt_eval flags   */
/*  + to use it, rename main_minimal by main and rename the "true" solar  */
/*    main function                                                       */
/*  + see main_minimal2() for a variation of this example                 */
/*------------------------------------------------------------------------*/
int main_minimal ( void ) {

  // set the problem SOLAR10:
  std::vector<Problem> problems;
  create_problems ( problems );
  const Problem * pb = find_problem ( problems, "10" );

  // will contain the outputs: they are not available as real values but only
  // through an output stream (to ensure the right number of decimals):
  std::ostringstream output_stream;

  // creation of the evaluator:
  Evaluator evaluator ( *pb, output_stream );
  
  // set the input point:
  double x[] = { 900, 10, 12, 0.2, 0.2 };
  evaluator.set_x(x);

  // parameters of the evaluation:
  int         x_index              = 0;     // only one point
  int         seed                 = 0;     // random seed
  double      fidelity             = 1.0;   // full fidelity
  int         replications         = 1;     // one replication
  bool        simulation_completed = false; // output flag
  bool        cnt_eval             = false; // output flag
  bool        verbose              = false; // no display
  std::string error_msg;                    // error message
  
  // the evaluation:
  evaluator.eval_x ( x_index, seed, fidelity, replications, simulation_completed, cnt_eval, error_msg, verbose );
  
  // get the outputs: this puts the outputs into output_stream:
  // (SOLAR10 has only one output, the objective)
  evaluator.display_outputs();

  // to get the output as a real:
  // double f = std::atof(output_stream.str().c_str());
  
  // display the outputs: objective and flags:
  std::cout << "output              : " << output_stream.str()  << std::endl
	    << "simulation completed: " << simulation_completed << std::endl
	    << "count evaluation    : " << cnt_eval             << std::endl;
  
  return 0;
}

/*---------------------------------------------------------------------------------*/
/*  another minimal example of a single evaluation (with inputs read from a file)  */
/*---------------------------------------------------------------------------------*/
int main_minimal2 ( int argc, char ** argv ) {

  if ( argc != 2 )
    return 1;

  // set the problem: SOLAR5
  std::vector<Problem> problems;
  create_problems ( problems );
  const Problem * pb = find_problem ( problems, "5" );
  int n = 20;
  
  std::ostringstream output_stream;
  Evaluator evaluator ( *pb, output_stream );
  
  // set the input point:
  double * x = new double[n];
  std::ifstream in ( argv[1] );

  if ( in.fail() ) {
    delete [] x;
    in.close();
    return 1;
  }
   
  for ( int i = 0 ; i < n ; ++i )
    in >> x[i];
  in.close();
  evaluator.set_x(x);
  delete [] x;
  
  // parameters of the evaluation:
  int         x_index              = 0;     // only one point
  int         seed                 = 0;     // random seed
  double      fidelity             = 1.0;   // full fidelity
  int         replications         = 1;     // one replication
  bool        simulation_completed = false; // output flag
  bool        cnt_eval             = false; // output flag
  bool        verbose              = false; // no display
  std::string error_msg;                    // error message
  
  // the evaluation:
  evaluator.eval_x ( x_index, seed, fidelity, replications, simulation_completed, cnt_eval, error_msg, verbose );
  
  // get the outputs: this puts the outputs into output_stream:
  evaluator.display_outputs();
  
  // display the outputs: objective and flags:
  std::cout << output_stream.str()  << " " << cnt_eval << std::endl;
  
  return 0;
}

/*----------------------------------------------------------*/
/*                parse and get the options                 */
/*                                                          */
/*  Example: solar ID  X.txt -v  -fid=0.5 -seed=0 -rep=10  */
/*----------------------------------------------------------*/
bool get_options ( int           argc      ,
		   char       ** argv      ,
		   std::string & x_file    ,
 		   bool        & verbose   ,
		   int         & seed      ,
		   double      & fidelity  ,
		   int         & replications ) {
  
  verbose      = false;
  seed         = 0;
  fidelity     = 1.0;
  replications = 1;

  x_file.clear();
  
  bool bseed = false;
  bool bfid  = false;
  bool brep  = false;
  bool chk   = false;
  
  std::string arg, sarg, sarg2;
  
  for ( int i = 2 ; i < argc ; ++i ) {

    arg = toupper(argv[i]);
    chk = false;
   
    // verbose:
    // --------
    if ( arg == "-V" || arg == "--V" ) {
      if ( verbose )
	return false;
      verbose = true;
      chk     = true;
    }
        
    else {
    
      // random seed:
      // ------------
      sarg = arg.substr(0,6);
      if ( sarg == "-SEED=" ) {
	if ( bseed )
	  return false;
	bseed = true;
	sarg2 = arg.substr(6,arg.size()-6);
	if ( sarg2 == "DIFF" || sarg2 == "-1" )
	  seed = -1;
	else if ( !string_to_int ( sarg2 , seed ) || seed < 0 )
	  return false;
	chk = true;
      }

      // fidelity:
      // ---------
      sarg = arg.substr(0,5);
      if ( sarg == "-FID=" ) {
	if ( bfid )
	  return false;
	bfid  = true;
	sarg2 = arg.substr(5,arg.size()-5);
	if ( !string_to_double(sarg2,fidelity) || fidelity <= 0.0 || fidelity > 1.0 )
	  return false;
	chk = true;
      }
      
      // number of replications:
      // -----------------------
      if ( sarg == "-REP=" ) {
	if ( brep )
	  return false;
	brep = true;
	sarg2 = arg.substr(5,arg.size()-5);
	if ( !string_to_int ( sarg2 , replications ) || replications <= 0 )
	  return false; 
	chk = true;
      }
    }
   
    if ( !chk ) {     
      if ( x_file.empty() )
	x_file = argv[i];
      else
	return false;
    }
  }
  return true;
}

/*----------------------------------------------------------------*/
/*  validation functions: called by "solar -check" or "solar -c"  */
/*----------------------------------------------------------------*/
bool check_eval ( const std::string & pb_id        ,
		  int                 seed         ,
		  double              fidelity     ,
		  int                 replications ,
		  const double      * x            ,
		  const std::string & expected     ,
		  std::string       & error          ) {

  std::vector<Problem> problems;
  create_problems ( problems );

  bool               simulation_completed, cnt_eval;
  std::string        err_msg;
  const Problem    * pb = find_problem ( problems, pb_id );    
  std::ostringstream oss;
  Evaluator          evaluator ( *pb , oss );
  
  evaluator.set_x(x);

  bool chk = false;
 
  if ( evaluator.eval_x ( 0, seed, fidelity, replications, simulation_completed, cnt_eval, err_msg, false ) ) {

    evaluator.display_outputs();
    if ( oss.str() == expected ) {
      chk = true;
    }
    else {
      chk    = false;
      error += "Error with " + pb_id + ":\n" +
	"\tOutput  : " + oss.str() + "\n\tExpected: " + expected + "\n";
    }
  }
  else {
    chk    = false;
    error += "Error with " + pb_id + ":\n\t" + err_msg + "\n";
  }
  
  return chk; 
}

/*----------------------------------------------------------------*/
bool check ( bool fast ) {
/*----------------------------------------------------------------*/
 
  Clock       clock1 , clock2;
  std::string error;
  bool        chk = true;
  
  std::cout << "\nValidation tests (can take several minutes):\n\n";
  
  // check random number generator (RNG):
  // ------------------------------------
  {
    std::ostringstream oss;
    std::string expected
      = "31.4297651712 72.7963156003 25.580539374 33.5965603203 50.053351454 55.2808740072 50.8702460562 ";
    RNG::set_seed(17);
    oss << std::setprecision(12);
    for ( int i = 0 ; i < 7 ; ++i )
      oss << RNG::rand(0,100) << " ";
    std::cout << "\tRNG test  ( 1/ 2) ..." << std::flush;
    if ( expected == oss.str() )
      std::cout << "... Ok";
    else {
      std::cout << "... Fail";
      error += "Error with RNG(1): Output=" + oss.str() + "\n";
      chk = false;
    }
    std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    clock2.reset();
    oss.str("");
    RNG::set_seed(17);
    for ( int i = 0 ; i < 7 ; ++i )
      oss << RNG::rand(0,100) << " ";
    std::cout << "\tRNG test  ( 2/ 2) ..." << std::flush;
    if ( expected == oss.str() )
      std::cout << "... Ok";
    else {
      std::cout << "... Fail";
      error += "Error with RNG(2): Output=" + oss.str() + "\n";
      chk = false;
    }
    std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
  }

  // check evaluations:
  // ------------------
  {

    int nb_eval_tests = 26;
    std::string expected_output;
      
    // tests for SOLAR1:
    {
      clock2.reset();
      std::cout << "\tEval test ( 1/" << nb_eval_tests << ") ..." << std::flush;
      double x[9] = { 8, 8, 150, 7, 7, 250, 45, 0.5, 5 };
      expected_output = "-122505.5978 -10881140.57 -1512631.39776 -134 -4.5 0";
      if ( check_eval ( "SOLAR1", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
    {
      clock2.reset();
      std::cout << "\tEval test ( 2/" << nb_eval_tests << ") ..." << std::flush;
      double x[9] = { 14.972466423, 13.2463292833, 109.7450505, 15.35568477,
    		      12.167043, 328, 81.23243276, 1.1838879569, 7.11915913 };
      expected_output = "-536680.901352 -425.8820987 -1108498.72275 -79.800117654 -5.9352711731 0";     
      if ( check_eval ( "SOLAR1", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;    
    }

    // tests for SOLAR2:
    if ( !fast ) {
      clock2.reset();
      std::cout << "\tEval test ( 3/" << nb_eval_tests << ") ..." << std::flush;
      double x[14] = { 11, 11, 140, 10, 10, 2650, 89, 0.5, 5, 838, 36, 0.3, 0.02, 0.0216 };
      expected_output =
    	"753526.705927 -3246473.29407 50.0941358025 -152038072.932 -118 -4.5 1527 178330543.218 -325.41767 -59.324704 -0.0016 -14.9303632679 -45";
      if ( check_eval ( "SOLAR2", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
    if ( !fast ) {
      clock2.reset();
      std::cout << "\tEval test ( 4/" << nb_eval_tests << ") ..." << std::flush;
      double x[14] = { 10.0849, 8.46, 138.941, 10.1341, 29.099, 2610,
    		       88.5164, 0.587, 6.9972, 989.693, 28, 0.2024, 0.033, 0.0562 };
      expected_output =
    	"1449917.55379 -2550082.44621 1.45717592593 -105818274.511 -118.7712 -6.4102 0 -108957877.298 -473.92953 -59.515714 -0.0232 -44.1350023134 -196.693";
      if ( check_eval ( "SOLAR2", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }

    if ( !fast ) {
      clock2.reset();
      std::cout << "\tEval test ( 5/" << nb_eval_tests << ") ..." << std::flush;
      double x[14] = { 12.9018637, 7.74041232, 206.6668478, 14.6224053, 14.737138, 2108, 65.82325952,
    		       1.559523919, 4.750141173, 991.8610965, 142, 3.43340143, 0.0193787, 0.02588297 };
      expected_output =
    	"987823.606284 -3012176.39372 0 -97967207.3177 -180.8631204 -3.190617254 0 -167813199.552 -476.05114 -59.529626 -0.00650427 -19.4736604979 -198.8610965";
      if ( check_eval ( "SOLAR2", 15478, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }  

    // tests for SOLAR3:
    {
      clock2.reset();
      std::cout << "\tEval test ( 6/" << nb_eval_tests << ") ..." << std::flush;
      double x[20] = { 8, 8, 150, 7, 7, 250, 45, 0.5, 5, 900, 9, 9, 0.3, 0.2, 560, 40, 0.3, 0.015, 0.017, 3 };
      expected_output =
    	"107541652.033 -362631.397758 73.0119047619 -134 -4.5 0 -138720188.816 -385.43902 -64.222892 -65 -0.002 -10.3155742876 -107 -0.05007853702";
      if ( check_eval ( "SOLAR3", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
    {
      clock2.reset();
      std::cout << "\tEval test ( 7/" << nb_eval_tests << ") ..." << std::flush;
      double x[20] = { 9.30925125, 9.39667813, 117.117247, 9.338565, 9.58816529,
    		       439, 26.19944197, 0.54375869, 7.0398380584, 856.7626312, 9.653758656, 10.91743348,
    		       0.59864555, 0.03119655, 529.0704, 72, 0.0222818, 0.02203051, 0.030181399, 2 };
      expected_output =
    	"70946739.3175 -491014.803208 2.94444444444 -98.4987445 -6.4960793684 0 -215134002.847 -357.94231 -33.703392 -34.0704 -0.008150889 -12.8879940902 -53.7626312 -0.05001958666";
      if ( check_eval ( "SOLAR3", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
    {
      clock2.reset();
      std::cout << "\tEval test ( 8/" << nb_eval_tests << ") ..." << std::flush;
      double x[20] = { 13, 9, 150, 7, 8, 249, 46, 0, 6, 900, 10, 8, 0.3, 0.2, 560, 42, 0.3, 0.015, 0.017, 2 };
      expected_output =
    	// "1e+20 -149690.320707 1e+20 -124 -6 1e+20 1e+20 1e+20 1e+20 1e+20 -0.002 -11.8523706144 -97 1e+20";
        "1e+20 -149690.320707 1e+20 -124 -6 1e+20 1e+20 1e+20 1e+20 1e+20 -0.002 -11.8523706144 1e+20 1e+20";
      if ( check_eval ( "SOLAR3", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
    {
      clock2.reset();
      std::cout << "\tEval test ( 9/" << nb_eval_tests << ") ..." << std::flush;
      double x[20] = { 9.31557819, 9.46006023, 116.559148749, 9.1581308098, 9.5177389864,
    		       438, 26.1305053825, 0.5575474642, 7.2069156184, 957.962200399, 9.5079960916,
    		       9.6607452026, 0.14149852, 0.01, 526.0248528820, 69, 0.8741599092, 0.01586339,
    		       0.0722243104, 2 };
      expected_output =
    	"65367065.2237 -480103.329801 4.17857142857 -97.927992369 -6.6493681542 0 -372150822.587 -441.45854 -29.493412 -31.024853 -0.0563609204 -9.96695202163 -154.962200399 -0.0500492024";
      if ( check_eval ( "SOLAR3", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }

    {
      clock2.reset();
      std::cout << "\tEval test (10/" << nb_eval_tests << ") ..." << std::flush;
      double x[20] = { 6.98582, 15.0696, 111.158, 6.87621, 15.2899, 366, 38.7892, 0.451441, 8.92151, 971.883,
    		       20.1932, 6.52063, 0.126425, 0.01, 526.466, 323, 0.633657, 0.00795822, 0.0371892, 2 };
      expected_output =
    	"63038104.0469 -135900.836066 0 -97.18636 -8.470069 0 -369541611.54 -472.9671 -29.453781 -31.466 -0.02923098 -12.0052071571 -168.883 -0.7244717952";
      if ( check_eval ( "SOLAR3", 17, 0.75, 10, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }

    {
      clock2.reset();
      std::cout << "\tEval test (11/" << nb_eval_tests << ") ..." << std::flush;
      double x[20] = { 14, 13, 150, 15, 14, 260, 64, 1.5, 4, 910, 29, 9, 0.01, 0.01, 650, 56, 3.97, 0.0143, 0.1, 2 };

      expected_output =
    	"93179363.113 -454424.808105 100 -122 -2.5 0 -451010886.47 -415 -146.77014 -155 -0.0857 -16.3911485751 -107 0";
      
      if ( check_eval ( "SOLAR3", 33, 0.27, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }  
    
    // tests for SOLAR4:
    {
      clock2.reset();
      std::cout << "\tEval test (12/" << nb_eval_tests << ") ..." << std::flush;
      double x[29] = { 9, 9, 150, 6, 8, 1000, 45, 0.5, 5, 900, 9, 9, 0.3, 0.2, 560, 500, 0.3, 0.0165, 0.018, 0.017, 10, 0.0155, 0.016, 0.2, 3, 12000, 1, 2, 2 };
      expected_output =
    	"78622308.4162 -1562631.39776 84.4142857143 -132 -4.5 3 -203757427.579 -386.09689 -63.656929 -65 -0.0015 -3.56637061436 -97 -0.1486950791 -0.001 -0.0005 -9206349.20635";
      if ( check_eval ( "SOLAR4", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
    {
      clock2.reset();
      std::cout << "\tEval test (13/" << nb_eval_tests << ") ..." << std::flush;
      double x[29] = { 4.796025798, 7.853658908, 125.7281067, 5.447879269, 7.141820154, 978, 44.821162, 1.497980897,
    		       7.039382352, 984.4806143, 8.587450041, 10.59030397, 0.1330657894, 0.0192747504, 649.9999406,
    		       295, 0.4625754212, 0.006425283473, 0.01258886031, 0.01946272908, 4.107100165, 0.00512445631,
    		       0.01667123691, 0.399935942, 8, 11998, 2, 1, 3, };
      expected_output =
    	"99048342.3651 -1414982.84004 55.5753968254 -116.136055104 -5.541401455 0 -159476828.703 -445.45938 -15.502131 4.6080534 -0.006163576837 -7.50463107308 -191.4806143 0.297471648 -0.00279149217 -0.0115467806 -307268627.816";
      if ( check_eval ( "SOLAR4", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
   
    // tests for SOLAR5:
    if ( !fast) {
      clock2.reset();
      std::cout << "\tEval test (14/" << nb_eval_tests << ") ..." << std::flush;
      double x[20] = { 847.015, 27.2039, 29.9981, 1.231, 0.0101, 650, 23, 0.2211, 0.026504, 0.037583, 0.045, 9.38, 0.010004, 0.0209, 0.204, 2, 54700, 3, 1, 2 };      
      expected_output =
     	"-28.8817193932 -26173.114929 -192920008.529 -344.9493 -28.349607 -0.08269015 -0.011079 -8.56036896077 -44.015 -0.0430727657 -0.0241 -0.010896 -204459343.324";
      if ( check_eval ( "SOLAR5", 0, 1.0, 1, x, expected_output, error ) )
     	std::cout << "... Ok";
      else {
     	std::cout << "... Fail";
     	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }  
    
    if ( !fast ) {
      clock2.reset();
      std::cout << "\tEval test (15/" << nb_eval_tests << ") ..." << std::flush;
      double x[20] = { 900, 10, 12, 0.15, 0.1, 560, 24, 0.35, 0.02, 0.023, 0.05, 8, 0.02, 0.023, 0.2, 2, 5000, 5, 5, 1 };
      expected_output =
    	"-30.0614936989 -68427260.159 3469782.67137 -343.36191 0 136.36512 -0.003 -8.87277796077 -97 1.383897976 -0.027 -0.003 -40342363.2245";
      if ( check_eval ( "SOLAR5", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
  
    // tests for SOLAR6:
    {
      clock2.reset();
      std::cout << "\tEval test (16/" << nb_eval_tests << ") ..." << std::flush;
      double x[5] = { 900, 10, 12, 0.2, 0.2 };
      expected_output =
    	"4136232.10319 41.5648148148 -90697638.1288 -368.75421 -34.73144 -87 44.99652957";
      if ( check_eval ( "SOLAR6", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
    if ( !fast ) {
      clock2.reset();
      std::cout << "\tEval test (17/" << nb_eval_tests << ") ..." << std::flush;
      double x[5] = { 994.947, 46.1, 21.8312, 0.01, 0.01 };
      expected_output =
    	"44298455.5682 0 -81139035.9655 -451.54499 -31.908255 -181.947 -13.8870143";
      if ( check_eval ( "SOLAR6", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
    if ( !fast ) {
      clock2.reset();
      std::cout << "\tEval test (18/" << nb_eval_tests << ") ..." << std::flush;
      double x[5] = { 994.84699999999998, 46.090000000000003, 21.841200000000004, 0.01, 0.11 };
      expected_output =
    	"46257817.4945 0 -81103180.0025 -451.4926 -31.908287 -181.847 -13.90461209";
      if ( check_eval ( "SOLAR6", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
    if ( !fast ) {
      clock2.reset();
      std::cout << "\tEval test (19/" << nb_eval_tests << ") ..." << std::flush;
      double x[5] = { 994.847, 46.09, 21.8412, 0.01, 0.11 };
      expected_output =
    	"46257817.4945 0 -81103180.0025 -451.4926 -31.908287 -181.847 -13.90461209";
      if ( check_eval ( "SOLAR6", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
    // tests for SOLAR7:
    {
      clock2.reset();
      std::cout << "\tEval test (20/" << nb_eval_tests << ") ..." << std::flush;
      double x[7] = { 7, 7, 850, 40, 0.2, 0.01, 0.011 };
      expected_output =
    	"-2768.29493993 -30016849.1891 5656779685 -0.001 -47 -10.5555742876 6.491331149";
      if ( check_eval ( "SOLAR7", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }

    {
      clock2.reset();
      std::cout << "\tEval test (21/" << nb_eval_tests << ") ..." << std::flush;
      double x[7] = { 11.118612, 14.0293, 922.245, 403, 2.704774, 0.02526764, 0.02948456 };
      expected_output =
    	"-4939.4070342 -22502463.2016 -212401980.756 -0.00421692 -119.245 -10.1548952275 -0.02886186514";
      if ( check_eval ( "SOLAR7", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
  
    // tests for SOLAR8:
    {
      clock2.reset();
      std::cout << "\tEval test (22/" << nb_eval_tests << ") ..." << std::flush;
      double x[13] = { 11, 11, 200, 10, 10, 2650, 89, 0.5, 8, 36, 0.3, 0.02, 0.0216 };
      expected_output =
    	"-4.63128915645e+12 124663605.86 -38975.2625989 -178 -7.5 0 460968232.465 -0.0016 -14.9303632679 -3.19128915645e+12 0.5478427658";
      if ( check_eval ( "SOLAR8", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
   
    // tests for SOLAR9:
    {
      clock2.reset();
      std::cout << "\tEval test (23/" << nb_eval_tests << ") ..." << std::flush;
      double x[29] = { 9, 9, 150, 6, 8, 1000, 45, 0.5, 5, 900, 9, 9, 0.3, 0.2, 560, 50, 0.3, 0.0165, 0.018, 0.017, 10, 0.0155, 0.016, 0.2, 2, 12000, 1, 2, 2 };
      expected_output =
    	"1e+20 1e+20 1e+20 1e+20 -4562631.39776 -132 -4.5 1e+20 1e+20 1e+20 1e+20 1e+20 -0.0015 -11.6663706144 1e+20 1e+20 -0.001 -0.0005 1e+20";
      if ( check_eval ( "SOLAR9", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }
    {
      clock2.reset();
      std::cout << "\tEval test (24/" << nb_eval_tests << ") ..." << std::flush;
      double x[29] = { 10.41264945, 9.891229103, 118.7647519, 10.66521456, 11.24343755, 788, 68.3347142, 0.5023836989, 9.091052143, 934.0300707,
    		       22.53631097, 10.24318735, 0.06012616095, 0.01003464624, 649.9997812, 332, 3.670187562, 0.01526149657, 0.02659821692,
    		       0.01829685361, 7.568243894, 0.006819263246, 0.01705667356, 0.399353877, 11, 3390, 2, 5, 3 };
      expected_output =
    	"-763500000001 4.645627836e+10 -1088133975.31 -331500000001 -3613902.13451 -97.939453 -8.5886684441 0 -263540971.585 -386.410223 -55.33238453 -1.1455984573 -0.01133672035 -8.83054238665 -141.0300707 -0.13915353195 -0.00124018005 -0.010237410314 -220956570.629";
      
      if ( check_eval ( "SOLAR9", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }

    // tests for SOLAR10:
    {     
      clock2.reset();
      std::cout << "\tEval test (25/" << nb_eval_tests << ") ..." << std::flush;
      double x[5] = { 883.41, 23.63, 29.26, 1.81, 0.01 };
      expected_output = "93.557882";
      if ( check_eval ( "SOLAR10", 0, 1.0, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }

    {     
      clock2.reset();
      std::cout << "\tEval test (26/" << nb_eval_tests << ") ..." << std::flush;
      double x[5] = { 995, 50, 20.48, 0.02, 0.01 };
      expected_output = "42.488721";
      if ( check_eval ( "SOLAR10", 0, 0.65, 1, x, expected_output, error ) )
    	std::cout << "... Ok";
      else {
    	std::cout << "... Fail";
    	chk = false;
      }
      std::cout << " \tTime: CPU=" << clock2.get_CPU_time() << " \treal=" << clock2.get_real_time() << std::endl;
    }   
  
  } // end of validation tests

  if ( chk )
    std::cout << "\nThis version of SOLAR is valid" << std::endl;
  else
    std::cout << "\nValidation failed: This version of SOLAR is not valid\n\n"
	      << error << std::endl;
  
  std::cout << std::endl
	    << "CPU time : " << clock1.get_CPU_time () << "s" << std::endl
	    << "Real time: " << clock1.get_real_time() << "s" << std::endl;
  
  return chk;
}
