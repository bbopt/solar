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
#ifndef __PROBLEM_H__
#define __PROBLEM_H__

#include <vector>
#include "Global.hpp"

class Problem {

private:

  std::string _pb_id;
  std::string _f_description;
  int         _index , _p , _n , _m;

public:

  Problem ( const std::string & pb_id, const std::string & f, int index , int p , int n, int m )
    : _pb_id(pb_id), _f_description(f), _index(index) , _p(p) , _n(n), _m(m) {}
  
  Problem ( const Problem & pb ) : _pb_id         ( pb._pb_id         ) ,
				   _f_description ( pb._f_description ) ,
				   _index         ( pb._index         ) ,
				   _p             ( pb._p             ) ,
				   _n             ( pb._n             ) ,
				   _m             ( pb._m             )   {}
  ~Problem ( void ) {}

  const std::string & get_pb_id         ( void ) const { return _pb_id; }
  const std::string & get_f_description ( void ) const { return _f_description; }

  bool is_stochastic ( int output_index = -1 ) const;
  
  int get_index      ( void ) const { return _index; }
  int get_p          ( void ) const { return _p; }  // number of objectives
  int get_n          ( void ) const { return _n; }  // number of variables
  int get_m          ( void ) const { return _m; }  // number of constraints
  int get_nb_outputs ( void ) const { return _p+_m; }
};


// create the list of problems:
void create_problems ( std::vector<Problem> & problems );

// find a problem from its pb_id:
const Problem * find_problem ( const std::vector<Problem> & problems , const std::string & pb_id );

#endif
