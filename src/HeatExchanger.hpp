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
#ifndef __HEAT_EXCHANGER_H__
#define __HEAT_EXCHANGER_H__

#include "MoltenSalt.hpp"
#include "Global.hpp"
#include "Constants.hpp"
#include "Powerblock.hpp"
#include "Turbine.hpp"
#include <stdexcept>
#include <vector>
#include <algorithm>

class HeatExchanger {

private:
  MoltenSalt * _input;
  MoltenSalt * _output;
  Powerblock * _powerblock;

  double _inletWaterTemperature;
  double _inletWaterPressure;
  double _outletSteamTemperature;
  double _outletSteamPressure;
  int    _exchangerModel;
		
  double _tubesLength;
  double _tubesDin;
  double _tubesDout;
  double _tubesSpacing;
  double _baffleCut;
  int    _nbOfBaffles;
  int    _nbOfTubes;
  int    _nbOfPassesPerShell;
  int    _nbOfShells;

  // corollary attributes
  double _baffleSpacing;
  double _shellWidth;
  double _shellCrossSection;
  int    _crossFlowRows;
  double _L_E;
  double _windowAngle;
  double _windowArea_FG;
  double _windowArea_FR;
  double _windowArea_F;
  double _window_Nrows;
  int    _totalRows;
  double _window_Ntubes;
  double _window_EqDiameter;
  double _window_EqPerimeter;
  double _longitudinalPitch;
  double _bypassArea;
  double _bundleArea;
  double _bundle_EqDiameter;
  double _nozzlesDiameter;
  double _nozzlesArea;
  double _a, _b, _c, _e;

  // Data gathering
  std::vector<double> _heatTransferred;
  
  //geometric parameters to add
  double fComputeEpsilon(void);
  double fComputeC1(double, double) const;
  double fComputeC2(void) const;
  double fComputeM(double, double) const;

  double fComputeWaterEnthalpy ( double T_o ) const {
    return 0.146*pow(T_o, 2.0) + 2182.4*T_o + 2.0*pow(10.0,6);
  }

public:

  HeatExchanger ( MoltenSalt * input  ,
		  MoltenSalt * output ,
		  Powerblock * powerblock );

  HeatExchanger ( MoltenSalt * input          ,
		  MoltenSalt * output         ,
		  Powerblock * powerblock     ,
		  double       tubesLength    ,
		  double       tubesDin       ,
		  double       tubesDout      ,
		  double       tubesSpacing   ,
		  double       baffleCut      ,
		  int          nbOfBaffles    ,
		  int          nbOfTubes      ,
		  int          passesPerShell ,
		  int          nbOfShells       );
  
  ~HeatExchanger ( void ) {}

  double fComputeRequiredMoltenSaltMassFlow ( double, double) const;
  void   fCalculateEnergyTransferred        ( void );
  double fComputeRequiredMoltenSaltMassFlow ( double, double, double);
  double fEnergyToPowerBlock                ( int );

  double computeYieldPressure ( void ) const {
    return (_tubesDout - _tubesDin) * SS316_YIELD_PRESSURE / (0.5*(_tubesDin + _tubesDout));
  }

  double computePressureInTubes  ( double ) const;
  double computePressureInShells ( void   ) const;
  
  const std::vector<double> & get_heatTransferred ( void ) const { return _heatTransferred; }

  int    get_exchangerModel     ( void ) const { return _exchangerModel;     }
  double get_tubesSpacing       ( void ) const { return _tubesSpacing;       }
  double get_tubesDin           ( void ) const { return _tubesDin;           }
  double get_tubesDout          ( void ) const { return _tubesDout;          }
  double get_tubesLength        ( void ) const { return _tubesLength;        }
  int    get_nbOfTubes          ( void ) const { return _nbOfTubes;          }
  int    get_nbOfPassesPerShell ( void ) const { return _nbOfPassesPerShell; }
  int    get_nbOfShells         ( void ) const { return _nbOfShells;         }

};

#endif
