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
#ifndef __MOLTEN_SALT_H__
#define __MOLTEN_SALT_H__

#include "Constants.hpp"
#include <cmath>

//Simple class object containing basic thermodynamical properties of molten salt

class MoltenSalt {

private:
  
  double _temperature; // in K
  double _enthalpy;    // in J/kg*K
  double _pressure;    // in kPa
  double _massFlow;    // in kg/s
  double _viscosity;   // in kg/ms

public:
  
  MoltenSalt ( double, double, double );
  MoltenSalt ( double, double );
  MoltenSalt ( MoltenSalt& );

  double get_temperature ( void ) const { return _temperature; }
  double get_enthalpy    ( void ) const{ return _enthalpy; }
  double get_pressure    ( void ) const{ return _pressure; }
  double get_massFlow    ( void ) const{ return _massFlow; }
  double get_viscosity   ( void ) const { return _viscosity; }
  
  void set_pressure      ( double pres) { _pressure = pres; }
  void set_massFlow      ( double masf) { _massFlow = masf; }
  void set_temperature   ( double );
  void set_enthalpy      ( double );
  void fComputeViscosity ( void   );
  void fModifyEnergy     ( double );

  static double fComputeViscosity ( double );
};

#endif
