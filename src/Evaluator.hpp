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
#ifndef __EVALUATOR__
#define __EVALUATOR__

#include "helpFunctions.hpp"
#include "Global.hpp"
#include "Clock.hpp"
#include "Problem.hpp"
#include <iostream>
#include <ostream>
#include <fstream>
#include <algorithm>

class Evaluator {
        
private:

  const Problem         & _problem;
  std::vector<double *>   _x;
  double                * _outputs;
  double                * _intermediate_outputs;
  std::ostream          & _out;
  
  void display_x     ( int x_index ) const;
  void delete_x      ( void );
  void reset_outputs ( double v );

  bool compute_mean_var ( std::vector<double> & output       ,
			  int                   output_index ,
			  double              & mean         ,
			  double              & var            ) const;

  // one evaluation:
  bool eval_x ( int           x_index              ,
		int           seed                 ,
		double        fidelity             ,
		bool        & simulation_completed ,
		bool        & cnt_eval             ,
		std::string & err_msg              ,
		bool          verbose                );
  
public:
        
  Evaluator ( const Problem & pb, std::ostream & out ) :
    _problem              ( pb                              ) ,
    _outputs              ( new double[pb.get_nb_outputs()] ) ,
    _intermediate_outputs ( NULL                            ) ,
    _out                  ( out                             )   { reset_outputs(1e20); }

  ~Evaluator ( void );

  int get_nb_input_points ( void ) const { return static_cast<int>(_x.size()); }
  
  bool read_x ( const std::string & x_file_name );

  bool set_x ( const double * x );
  
  void display_x ( void ) const;

  bool has_intermediate_outputs ( void ) const { return _intermediate_outputs; }
  
  // multiple evaluations:
  bool eval_x ( int           x_index              ,
		int           seed                 ,
		double        fidelity             ,
		int           replications         ,
		bool        & simulation_completed ,
		bool        & cnt_eval             ,
		std::string & err_msg              ,
		bool          verbose                );

  void display_outputs              ( void ) const;
  void display_intermediate_outputs ( void ) const;
};

#endif
