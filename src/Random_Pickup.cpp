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
#include "Random_Pickup.hpp"

Random_Pickup::Random_Pickup ( int n ) : _n0   ( n          ) ,
					 _n    ( n          ) ,
					 _elts ( new int[n] )   {
  for ( int i = 0 ; i < n ; ++i )
    _elts[i] = i;
}

void Random_Pickup::reset ( void ) {
  _n = _n0;
  for ( int i = 0 ; i < _n ; ++i )
    _elts[i] = i;
}

int Random_Pickup::pickup ( void ) {
  if ( _n == 0 )
    return 0;
  int ind = RNG::rand()%_n;
  int tmp = _elts[ind];
  if ( ind < _n - 1 ) {
    _elts[ind ] = _elts[_n-1];
    _elts[_n-1] = tmp;
  }
  --_n;
  return tmp;
}

void Random_Pickup::cancel_last_pickup ( void ) {
  if ( _n < _n0 )
    ++_n;
}

