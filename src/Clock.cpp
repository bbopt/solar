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
#include "Clock.hpp"

/*-----------------------------------*/
/*   static members initialization   */
/*-----------------------------------*/
const double Clock::_D_CLOCKS_PER_SEC = static_cast<double>(CLOCKS_PER_SEC);

/*---------------------------------------------------------*/
/*  compute the wall-clock time (real time) elapsed since  */
/*  the construction of the Clock object                   */
/*---------------------------------------------------------*/
int Clock::get_real_time ( void ) const
{
    time_t t2;
    time  (&t2);
    int dti= static_cast<int>( difftime ( t2 , _real_t0 ) );
    return dti;
}
