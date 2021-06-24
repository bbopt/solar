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
#include "Powerblock.hpp"

/*--------------------------------------------------------------------------*/
/*                                 constructor                              */
/*--------------------------------------------------------------------------*/
Powerblock::Powerblock ( int typeOfTurbine ) : _turbine(typeOfTurbine) {

  int n = _turbine._nbOfStages;

  // tuning partial load correction coefficients:
  a = A1 + B1*n + C1*pow(n,2) + D1*pow(n,3);
  b = A2 + B2*n + C2*pow(n,2) + D2*pow(n,3);
  c = A3 + B3*n + C3*pow(n,2) + D3*pow(n,3);
  d = A4 + B4*n + C4*pow(n,2) + D4*pow(n,3);
  
  _powerOutput.reserve         (86400); // 24 * 60 * 60
  _requiredThermalPower.reserve(86400); // 24 * 60 * 60
}

/*--------------------------------------------------------------------------*/
double Powerblock::fComputeRequiredThermalEnergy ( double & Pout ) {
/*--------------------------------------------------------------------------*/

  _Pout = 0.;

  double requiredThermalEnergy;
  
  if ( Pout > 0 ) {

    double f, correctionFactor, efficiency;
      
    if (Pout >= _turbine._maxPower) {
      _Pout = _turbine._maxPower;
      f     = 100;
    }
    else if (Pout < _turbine._minPower) {
      _Pout = _turbine._minPower;
      f     = _Pout / _turbine._maxPower;
    }
    else {
      _Pout = Pout;
      f     = 100 * (_Pout / _turbine._maxPower);
    }

    correctionFactor = exp(a + f*b + f*f*c + f*f*f*d);

    efficiency = _turbine._basicEfficiency * correctionFactor;

    requiredThermalEnergy = _Pout / efficiency;
  }
  else {
    _Pout = 0.;
    requiredThermalEnergy = 0.;
  }

  return requiredThermalEnergy;
}

/*-------------------------------------------------------------------------------------*/
void Powerblock::adjustPowerData ( double thermalTransferred , double thermalNeeded ) {
/*-------------------------------------------------------------------------------------*/
  if ( fabs(thermalTransferred - thermalNeeded) > 1000.0 ) {
    _powerOutput.pop_back();
    _powerOutput.push_back(0.0);
  }
  else
    _powerOutput.push_back(thermalTransferred);
}
