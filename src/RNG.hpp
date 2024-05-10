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

/*-----------------------------------------------*/
/*  RNG class imported from NOMAD                */
/*-----------------------------------------------*/
#ifndef __RNG__
#define __RNG__

#include <cmath>
#include <stdexcept>
#include <climits>

// use of 'access' or '_access', and getpid() or _getpid():
#ifdef _MSC_VER
#include <io.h>
#include <process.h>
#else
#include <unistd.h>
#endif

#if !defined(UINT32_MAX)
typedef unsigned int uint32_t;
#define UINT32_MAX    0xffffffff
#endif

class RNG {

private:
  
  static uint32_t x_def,y_def,z_def,_x,_y,_z;
  static int _s;

  // reset seed to its default value:
  static void reset_private_seed_to_default ( void ) {
    _x=x_def;
    _y=y_def;
    _z=z_def;
  }
  
public:  
  
  static int  get_seed ( void ) { return static_cast<int>(_s); }
  static void set_seed ( int s );
  static int  get_pid  ( void );

  static uint32_t rand ( void );
  static double   rand ( double a, double b) { return a+((b-a)*RNG::rand())/UINT32_MAX; }

};

#endif
