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
#ifndef __POWERBLOCK_H__
#define __POWERBLOCK_H__

#include "Global.hpp"
#include "Turbine.hpp"
#include <vector>

class Powerblock {
  
private:

  Turbine             _turbine;
  std::vector<double> _requiredThermalPower;
  std::vector<double> _powerOutput;
  double              _Pout;
  double              _steamRate;

public:

  Powerblock ( int );

  int    get_typeOfTurbine         ( void  ) const { return _turbine._ID;               }
  double get_hotEnthalpy           ( void  ) const { return _turbine._inletEnthalpy;    }
  double get_coldEnthalpy          ( void  ) const { return _turbine._outletEnthalpy;   }
  double get_temperature           ( void  ) const { return _turbine._inletTemperature; }
  double get_pressure              ( void  ) const { return _turbine._inletPressure;    }
  double get_powerOfTurbine        ( void  ) const { return _turbine._maxPower;         }
  double get_powerOutput           ( int i ) const { return _powerOutput[i];            }
  double get_Pout                  ( void  ) const { return _Pout;                      }
  double get_turbineOutletEnthalpy ( void  ) const { return _turbine._outletEnthalpy;   }
  double get_steamRate             ( void  ) const { return _steamRate;                 }

  void adjustPowerData ( double, double );
  void set_steamRate   ( double x ) { _steamRate = x; }
  
  double fComputeRequiredThermalEnergy(double& );

  void add_requiredThermalPower ( double x ) { _requiredThermalPower.push_back(x); }
  void reserve ( size_t n ) { _powerOutput.reserve(n); }
};

#endif
