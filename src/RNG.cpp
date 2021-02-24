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
#include "RNG.hpp"


// Default values for the provided number seed:
int RNG::_s = 0;


uint32_t RNG::x_def = 123456789;
uint32_t RNG::y_def = 362436069;
uint32_t RNG::z_def = 521288629;
uint32_t RNG::_x = x_def;
uint32_t RNG::_y = y_def;
uint32_t RNG::_z = z_def;

int RNG::get_pid ( void ) {
#ifdef _MSC_VER
  return _getpid();
#else
  return getpid();
#endif
}

void RNG::set_seed ( int s ) {
  
  if( s<=INT_MAX && s>=0 )
    _s=s;
  else
    throw std::invalid_argument ( "RNG::set_seed(): invalid seed. Seed should be in [0,INT_MAX]" );
  RNG::reset_private_seed_to_default();
  for ( int i = 0 ; i < _s ; i++)
    RNG::rand();
}

// http://madrabbit.org/~ray/code/xorshf96.c //period 2^96-1
uint32_t RNG::rand ( void ) {
  
  uint32_t t;
  _x ^= _x << 16;
  _x ^= _x >> 5;
  _x ^= _x << 1;
  t = _x;
  _x = _y;
  _y = _z;
  _z = t ^ _x ^ _y;
  
  return _z;
}
