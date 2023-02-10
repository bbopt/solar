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
#include "Evaluator.hpp"

/*-----------------------------------------*/
/*                destructor               */
/*-----------------------------------------*/
Evaluator::~Evaluator ( void ) {
  delete_x();
  delete [] _outputs;
  if ( _intermediate_outputs )
    delete [] _intermediate_outputs;
}

/*-----------------------------------------*/
/*          display one input point        */
/*-----------------------------------------*/
void Evaluator::display_x ( int x_index ) const {
  for ( int i = 0 ; i < _problem.get_n() ; ++i )
    _out << _x[x_index][i] << " ";
}

/*-----------------------------------------*/
/*             display all inputs          */
/*-----------------------------------------*/ 
void Evaluator::display_x ( void ) const {
  _out << std::setprecision(12);
  for ( size_t k = 0 ; k < _x.size() ; ++k ) {
    _out << "input point #" << k+1 << ": ( ";
    display_x(static_cast<int>(k));
    _out << ")" << std::endl;
  }
}

/*-----------------------------------------*/
/*              clear inputs               */
/*-----------------------------------------*/
void Evaluator::delete_x ( void ) {
  for ( size_t i = 0 ; i < _x.size() ; ++i )
    delete [] _x[i];
  _x.clear();
}

/*-----------------------------------------*/
/*               read inputs               */
/*-----------------------------------------*/
bool Evaluator::read_x ( const std::string & x_file_name ) {

  delete_x();
  
  std::ifstream in ( x_file_name.c_str() );
  if ( in.fail() ) {
    in.close();
    return false;
  }

  int n = _problem.get_n();
    
  while ( !in.eof() ) {

    double * x = new double[n];

    for ( int i = 0 ; i < n ; ++i ) {
      in >> x[i];
      if ( in.fail() ) {
	in.close();
	delete [] x;
	delete_x();
	return false;
      }
      in >> std::ws;
    }

    _x.push_back(x);
  }
  in.close();
  return true;
}

/*-----------------------------------------*/
/*          define one input point         */
/*-----------------------------------------*/
bool Evaluator::set_x ( const double * x ) {
  delete_x();
  int n = _problem.get_n();
  double * xx = new double[n];
  for ( int i = 0 ; i < n ; ++i )
    xx[i] = x[i];
  _x.push_back(xx);
  return true;
}

/*------------------------------------------------------------*/
/*   Compute the mean and the coefficient of variation        */
/*   for one vector values of the same output                 */
/*   The best 80% of the values are considered                */
/*------------------------------------------------------------*/
bool Evaluator::compute_mean_var ( std::vector<double> & output       ,
				   int                   output_index ,
				   double              & mean         ,
				   double              & var            ) const {

  mean = var = 0.0;
 
  if ( output.empty() ) {
    mean = var = 1e20;
    return false;
  }
    
  if ( output.size() == 1 ) {
    mean = output[0];
    var  = 0.0;
    return true;
  }

  size_t i;  
  size_t n = output.size();

  mean = 0.0;

  int nb_fails = 0;
  
  for ( i = 0 ; i < n ; ++i ) {
    if ( output[i] == 1e20 ) {
      ++nb_fails;
      
      // we do not accept more than 20% of failures:
      if ( nb_fails > n*0.2 ) {
	mean = var = 1e20;
	return true;
      }
      
    }
    else
      mean += output[i];
  }
  mean = mean / (1.0*(n-nb_fails));
    
  var = 0.0;
  for ( i = 0 ; i < n ; ++i ) {
    if ( output[i] != 1e20 )
      var += pow(output[i]-mean,2.0);
  }
  
  var = sqrt(var/(n-nb_fails-1.0));
  if ( mean != 0.0 )
    var /= fabs(mean);
  
  if ( !_problem.is_stochastic(output_index) && fabs(var) > 1e-10 ) {
    // std::cout << " i=" << output_index << " var=" << var << std::endl; // DEBUG
    return false;
  }
      
  return true;
}

/*------------------------------------------------------------*/
/*  multiple evaluations (replications with different seeds)  */
/*------------------------------------------------------------*/
bool Evaluator::eval_x ( int           x_index              ,
			 int           seed                 ,
			 double        fidelity             ,
			 int           replications         ,
			 bool        & simulation_completed ,
			 bool        & cnt_eval             ,
			 std::string & err_msg              ,
			 bool          verbose                ) {
  
  if ( !_problem.is_stochastic() ) {
    if ( replications > 1 ) {
      err_msg = "Number or replications > 1 for a deterministic instance";
      return false;
    }
    if ( seed != 0 ) {
      err_msg = "Seed != 0 for a deterministic instance";
      return false;
    }
  }
  
  if ( replications == 1 )
    return eval_x ( x_index, seed, fidelity, simulation_completed, cnt_eval, err_msg, verbose );
  
  reset_outputs(1e20);
  cnt_eval = simulation_completed = false;
  
  if ( replications <= 1 ) {
    err_msg = "Number or replications must be > 1";
    return false;
  }

  if ( x_index < 0 || x_index >= static_cast<int>(_x.size()) ) {
    err_msg = "Evaluator::eval_x(): Bad input index";
    return false;
  }
  
  if ( verbose ) {
    _out << std::endl
	 << "Begin simulation of " << _problem.get_pb_id() << " for ";
    if ( _x.size() > 1 )
      _out << "point #" << x_index+1 << "/" << _x.size() << " = (";
    else
      _out << "x=( ";
    display_x ( x_index );
    _out <<  ") and " << replications << " replications:" << std::endl;
  }

  simulation_completed = true;
  bool simc, cnt;
  int nbo = _problem.get_nb_outputs();
  
  std::vector<double> * output_matrix = new std::vector<double>[nbo];
  double              * means         = new double[nbo];
  double              * vars          = new double[nbo];
 
  for ( int o = 0; o < nbo; ++o )
    means[o] = vars[o] = 1e20;
 
  for ( int r = 0 ; r < replications ; ++r ) {

    if ( verbose )
      _out << "\tRun #" << std::setw(3) << r+1 << " (seed=" << seed << ") ..." << std::flush;

    reset_outputs(1e20);

    // blackbox evaluation is here:
    if ( !eval_x ( x_index, seed, fidelity, simc, cnt, err_msg, false ) ) {
      reset_outputs(1e20);
      delete [] output_matrix;
      delete [] means;
      delete [] vars;
      if ( verbose )
	_out << " Fail\n";
      if ( cnt )
	cnt_eval = true;
      if ( !simc )
	simulation_completed = false;
      return false;
    }

    if ( verbose ) {
      if ( !err_msg.empty() )
	_out << " [warning: " << err_msg << "] ";
      _out << "... outputs: ";
      display_outputs();
    }
    
    // compute mean and coeff. of variation:
    for ( int o = 0; o < nbo; ++o ) {
      
      output_matrix[o].push_back(_outputs[o]);
      
      if ( verbose ) {
	
	if  ( !compute_mean_var ( output_matrix[o], o, means[o], vars[o] ) ) {
	  reset_outputs(1e20);
	  delete [] output_matrix;
	  delete [] means;
	  delete [] vars;
	  _out << " Fail\n";

	  std::ostringstream err;
	  err << "problem with the computation of mean and coeff. of variation for output #" << o;
	  err_msg = err.str();
	  
	  return false;
	}
	
	if ( _problem.is_stochastic(o) )
	  _out << "\t m" << o << "=" << means[o]; //  << " v" << o << "=" << vars[o];
      }
    }
    
    if ( verbose )
      _out << std::endl;
  					  
    err_msg.clear();

    if ( cnt )
      cnt_eval = true;
    if ( !simc )
      simulation_completed = false;

    // Seed update:
    {
      double seed_tmp = seed + 11.0 + RNG::rand()%298712;
      if ( seed_tmp > 1.0*INT_MAX )
	seed_tmp = RNG::rand()%179841;
      seed = static_cast<int>(floor(seed_tmp));
    }
    
  }

  bool error = false;
  
  if ( !verbose ) {
    for ( int o = 0; o < nbo; ++o ) {
      if  ( !compute_mean_var ( output_matrix[o], o, means[o], vars[o] ) ) {
	error = true;
	std::ostringstream err;
	err << "problem with the computation of mean and coeff. of variation for output #" << o;
	err_msg = err.str();
      }
    }
  }
	  
  for ( int o = 0; o < nbo; ++o )
    _outputs[o] = means[o];
  
  delete [] output_matrix;
  delete [] means;
  delete [] vars;

  return !error;
}

/*---------------------------------------------*/
/*   one evaluation (interface to simulator)   */
/*---------------------------------------------*/
bool Evaluator::eval_x ( int           x_index              ,
			 int           seed                 ,
			 double        fidelity             ,
			 bool        & simulation_completed ,
			 bool        & cnt_eval             ,
			 std::string & err_msg              ,
			 bool          verbose                ) {

  reset_outputs(1e20);
  
  RNG::set_seed(seed);
  
  simulation_completed = cnt_eval = false;

  if ( x_index < 0 || x_index >= static_cast<int>(_x.size()) ) {
    err_msg = "Evaluator::eval_x(): Bad input index";
    return false;
  }
  
  bool status = true;
  err_msg.clear();
  
  try { 
     
    Scenario scenario ( _problem.get_pb_id(), fidelity );
   
    cnt_eval = scenario.set_x ( _x[x_index] );

    if ( verbose ) {
      _out << std::endl
	   << "Begin simulation of " << _problem.get_pb_id() << " for ";
      if ( _x.size() > 1 )
	_out << "point #" << x_index+1 << "/" << _x.size() << " = (";
      else
	_out << "x=( ";
      display_x ( x_index );
      _out <<  ") ..." << std::endl;
    }

    simulation_completed = scenario.simulate ( _outputs , _intermediate_outputs, fidelity, cnt_eval );

    if ( verbose )
      _out << "... done" << std::endl;
  }

  catch ( const std::invalid_argument & e ) {
    // problem with simulation (e.g. bad input): Simulation cannot be used by solver;
    // error message must be printed since outputs are not displayed.     
    err_msg = e.what();   
    status  = false;
  }
  // (controlled) problem with the simulation (e.g. a priori constraint violated):
  // Simulation should be usable by solver;
  // Some of the outputs are displayed with 1e20 values.
  // Simulation_completed may be equal to 'false' but since some outputs are valid,
  // the simulation can still be used by a solver: 'return 0' is then appropriate.
  catch ( const Simulation_Interruption & e ) {
    err_msg = e.what();
    status  = true;
  }
  catch (...) {
    err_msg = "Error: simulation was interrupted for an unknown reason";
    status  = false;
  }

  return status;
}

/*-------------------------------------------------------*/
/*           reset all outputs to a specific value       */
/*-------------------------------------------------------*/
void Evaluator::reset_outputs ( double v ) {
  int i;
  for ( i = 0; i < _problem.get_nb_outputs(); ++i )
    _outputs[i] = v;
  if ( _intermediate_outputs ) {
    delete [] _intermediate_outputs;
    _intermediate_outputs = NULL;
  }
  if ( _problem.get_pb_id() == "MINCOST_UNCONSTRAINED" ) {
    _intermediate_outputs = new double[7];
    for ( i = 0; i < 7; ++i )
      _intermediate_outputs[i] = v;
  }
}
  
/*--------------------------------------------------------*/
/*                      display outputs                   */
/*  (each output is displayed with a specific precision)  */
/*--------------------------------------------------------*/
void Evaluator::display_outputs ( void ) const {
 
  const std::string pb_id = _problem.get_pb_id();
 
  // solar 1:
  if ( pb_id == "MAXNRG_H1" )
    _out << std::setprecision(12) << _outputs[0] << " "
	 << std::setprecision(10) << _outputs[1] << " " // New in Version 0.4.2
	 << std::setprecision(12) << _outputs[2] << " "
	 << _outputs[3] << " "
	 << _outputs[4] << " "
	 << _outputs[5] << std::setprecision(12);
  
  // solar 2:
  else if ( pb_id == "MINSURF_H1" )
    _out << std::setprecision(12) << _outputs[0] << " "
	 << _outputs[1] << " "
	 << _outputs[2] << " "
	 << _outputs[3] << " "
	 << _outputs[4] << " "
	 << _outputs[5] << " "
	 << _outputs[6] << " "
	 << _outputs[7] << " "
	 << std::setprecision( 8) << _outputs[ 8] << " "
	 << std::setprecision( 8) << _outputs[ 9] << " "
         // << std::setprecision( 8) << _outputs[10] << " " c10 removed in version 0.5.4
	 << std::setprecision(12) << _outputs[10] << " "
	 << _outputs[11] << " "
	 << _outputs[12] << std::setprecision(12);

  // solar 3:
  else if ( pb_id == "MINCOST_C1" )
    _out << std::setprecision(12) << _outputs[0] << " "
	 << _outputs[1] << " "
	 << _outputs[2] << " "
	 << _outputs[3] << " "
	 << _outputs[4] << " "
	 << _outputs[5] << " "
	 << _outputs[6] << " "
	 << std::setprecision( 8) << _outputs[7] << " "
	 << std::setprecision( 8) << _outputs[8] << " "
	 << std::setprecision( 8) << _outputs[9] << " "
	 << std::setprecision(12) << _outputs[10] << " "
	 << _outputs[11] << " "
	 << _outputs[12] << " "
	 << std::setprecision(10) << _outputs[13]
	 << std::setprecision(12);
  
  // solar 4:
  else if ( pb_id == "MINCOST_C2" )
    _out << std::setprecision(12) << _outputs[0] << " "
	 << _outputs[1] << " "
	 << _outputs[2] << " "
	 << _outputs[3] << " "
	 << _outputs[4] << " "
	 << _outputs[5] << " "
	 << _outputs[6] << " "
	 << std::setprecision( 8) << _outputs[7] << " "
	 << std::setprecision( 8) << _outputs[8] << " "
	 << std::setprecision( 8) << _outputs[9] << " "
	 << std::setprecision(12) << _outputs[10] << " "
	 << _outputs[11] << " "
	 << _outputs[12] << " "
	 << std::setprecision(10) << _outputs[13] << " "
	 << std::setprecision(12) << _outputs[14] << " "
	 << _outputs[15] << " "
	 << _outputs[16] << std::setprecision(12);

  // solar 5:
  else if ( pb_id == "MAXCOMP_HTF1" )
    _out << std::setprecision(12) << _outputs[0] << " "
	 << std::setprecision(11) << _outputs[1] << " " // New in Version 0.4.2
	 << std::setprecision(12) << _outputs[2] << " "
	 << std::setprecision( 8) << _outputs[3] << " "
	 << std::setprecision( 8) << _outputs[4] << " "
	 << std::setprecision( 8) << _outputs[5] << " "
	 << std::setprecision(12) << _outputs[6] << " "
	 << _outputs[7] << " "
	 << _outputs[8] << " "
	 << std::setprecision(10) << _outputs[ 9] << " "
	 << std::setprecision(12) << _outputs[10] << " "
	 << _outputs[11] << " "
	 << _outputs[12] << std::setprecision(12);
  
  // solar 6:
  else if ( pb_id == "MINCOST_TS" )
    _out << std::setprecision(12) << _outputs[0] << " "
	 << _outputs[1] << " "
	 << _outputs[2] << " "
	 << std::setprecision( 8) << _outputs[3] << " "
	 << std::setprecision( 8) << _outputs[4] << " "
	 << std::setprecision(12) << _outputs[5] << " "
	 << std::setprecision(10) << _outputs[6]
	 << std::setprecision(12);
  
  // solar 7:
  else if ( pb_id == "MAXEFF_RE" )
    _out << std::setprecision(12) << _outputs[0] << " "
	 << _outputs[1] << " "
	 << _outputs[2] << " "
	 << _outputs[3] << " "
	 << _outputs[4] << " "
	 << _outputs[5] << " "
	 << std::setprecision(10) << _outputs[6]
	 << std::setprecision(12);
  
  // solar 8:
  else if ( pb_id == "MAXHF_MINCOST" )
    _out << std::setprecision(12) << _outputs[0] << " "
	 << _outputs[1] << " "
	 << _outputs[2] << " "
	 << _outputs[3] << " "
	 << _outputs[4] << " "
	 << _outputs[5] << " "
	 << _outputs[6] << " "
	 << _outputs[7] << " " 
	 << _outputs[8] << " "
	 << _outputs[9] << " "
	 << std::setprecision(10) << _outputs[10]
	 << std::setprecision(12);
  
  // solar 9:
  else if ( pb_id == "MAXNRG_MINPAR" )
    _out << std::setprecision(12) << _outputs[ 0] << " "
	 << std::setprecision(10) << _outputs[ 1] << " " // New in Version 0.4.2
	 << std::setprecision(12) << _outputs[ 2] << " "
	 << _outputs[ 3] << " "
	 << _outputs[ 4] << " "
	 << _outputs[ 5] << " "
	 << _outputs[ 6] << " "
	 << _outputs[ 7] << " "
	 << _outputs[ 8] << " "
	 << std::setprecision( 9) << _outputs[ 9] << " " // New in Version 0.4.2
	 << std::setprecision(10) << _outputs[10] << " " // New in Version 0.4.2
	 << std::setprecision(12) << _outputs[11] << " "
	 << _outputs[12] << " "
	 << _outputs[13] << " "
	 << _outputs[14] << " "
	 << std::setprecision(11) << _outputs[15] << " " // New in Version 0.4.2
	 << std::setprecision(12) << _outputs[16] << " "
	 << _outputs[17] << " "
	 << _outputs[18] << std::setprecision(12);

  // solar 10:
  else if ( pb_id == "MINCOST_UNCONSTRAINED" )
    _out << std::setprecision(8) << _outputs[0];
}

void Evaluator::display_intermediate_outputs ( void ) const {
 
  // solar 10:
  if ( _problem.get_pb_id() == "MINCOST_UNCONSTRAINED" ) {
    _out << std::setprecision(12) << _intermediate_outputs[0] << " "
	 << _intermediate_outputs[1] << " "
	 << _intermediate_outputs[2] << " "
	 << std::setprecision( 8) << _intermediate_outputs[3] << " "
	 << std::setprecision( 8) << _intermediate_outputs[4] << " "
	 << std::setprecision(12) << _intermediate_outputs[5] << " "
	 << std::setprecision(10) << _intermediate_outputs[6]
	 << std::setprecision(12);
  }
}
