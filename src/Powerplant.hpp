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
#ifndef __POWERPLANT_H__
#define __POWERPLANT_H__

#include "HtfCycle.hpp"
#include "MoltenSalt.hpp"
#include "HeliostatField.hpp"
#include "Economics.hpp"
#include "Powerblock.hpp"
#include "Constants.hpp"
#include "Global.hpp"
#include "Time_Manager.hpp"

#include <cmath>
#include <numeric>
#include <algorithm>
#include <sstream>
#include <fstream>

class Powerplant
{
private:

  // Simulation components
  Time_Manager _time;
  int  _model_type; // 1:heliostats field; 2: whole plant
  int  _heliostatsFieldModel;

  // Powerplant components
  HtfCycle       * _moltenSaltLoop;
  HeliostatField * _heliostatsField;
  Powerblock     * _powerblock;
  Economics      * _investmentCost;

  std::vector<double> _heliostatFieldPowerOutput;
  
  // Include day, starting hour, demand info, storage starting conditions, etc.
	
  std::vector<double> _receiverOutletFlow;
  std::vector<double> _receiverPumpHead;
  std::vector<double> _pressureShellSide;
  std::vector<double> _pressureTubesSide;
  std::vector<double> _steamRate;
  std::vector<double> _msRateSteamGen;
  std::vector<double> _heliostatsFieldPar;
  std::vector<double> _hotStoragePowerLosses;
  std::vector<double> _coldStoragePowerLosses;
  std::vector<double> _powerplantPowerOutput;
  std::vector<double> _demand;

  double _reflectiveSurface;
  double _costOfHeliostatsField;
  double _totalEnergyConcentrated;
  double _maximumPressureInReceiver;
  double _maximumPressureInExchanger;
  double _yieldPressureReceiver;
  double _yieldPressureExchanger;
  double _overallComplianceToDemand;

  void clean ( void );
  
public:

  Powerplant ( const Time_Manager&, int, HeliostatField*, HtfCycle*, Powerblock*, Economics* );
  ~Powerplant ( void ) { clean(); }

  const HtfCycle * get_moltenSaltLoop ( void ) const { return _moltenSaltLoop;  }
  HeliostatField * get_heliostatField ( void ) const { return _heliostatsField; }
  Economics      * get_investmentCost ( void ) const { return _investmentCost;  }

  void   fSimulatePowerplant         ( bool low_fid );
  void   fSimulateHeliostatField     ( void   );
  double fComputeSteamRate           ( double );
  double fComputePressureInExchanger ( void   );
  double fComputeParasiticLosses     ( void   ) const;

  double fComputeParasiticsForPb7 ( void ) const {
    return std::inner_product(_receiverPumpHead.begin(), _receiverPumpHead.end(),
			      _receiverOutletFlow.begin(), 0.0) / MS_DENSITY;
  }
  
  double fComputeParasiticsForPb3 ( void ) const;
  double fComputeParasiticsForPb9 ( void ) const;

  const std::vector<double> & get_powerplantPowerOutput ( void ) const { return _powerplantPowerOutput; }

  double fComputeThermalEnergy ( double steamRate ) const {
    return steamRate * (_powerblock->get_hotEnthalpy() - WATER_300K_1ATM_ENTHALPY);
  }
  
  double get_reflectiveSurface            ( void ) const { return _reflectiveSurface; }
  double get_costOfHeliostatsField        ( void ) const { return _costOfHeliostatsField; }
  double get_totalEnergyConcentrated      ( void ) const { return _totalEnergyConcentrated; }
  double get_costOfHeliostatField         ( void ) const { return _investmentCost->evaluateCostOfField();          }
  double get_costOfReceiver               ( void ) const { return _investmentCost->evaluateCostOfReceiver();       }
  double get_costOfTower                  ( void ) const { return _investmentCost->evaluateCostOfTower();          }
  double get_costOfSteamGenerator         ( void ) const { return _investmentCost->evaluateCostOfSteamGenerator(); }
  double get_costOfPowerblock             ( void ) const { return _investmentCost->evaluateCostOfPowerblock();     }
  double get_costOfStorage                ( void ) const { return _investmentCost->evaluateCostOfStorage();        }
  double get_steamTurbineInletTemperature ( void ) const { return _powerblock->get_temperature();                  }
  double get_overallComplianceToDemand    ( void ) const { return _overallComplianceToDemand; }
  double get_maximumPressureInReceiver    ( void ) const { return _maximumPressureInReceiver; }
  double get_yieldPressureInReceiver      ( void ) const { return _yieldPressureReceiver;                   }
  double get_yieldPressureInExchanger     ( void ) const { return _yieldPressureExchanger;                  }
  double get_maximumPressureInExchanger   ( void ) const { return _maximumPressureInExchanger;              }
  double get_minColdStorageTemp           ( void ) const { return _moltenSaltLoop->get_minColdStorageTemp();}
  double get_minHotStorageTemp            ( void ) const { return _moltenSaltLoop->get_minHotStorageTemp(); }
  double get_minSteamGenTemp              ( void ) const { return _moltenSaltLoop->get_minSteamGenTemp();   }

  void set_demand         ( const std::vector<double>& demandVector ) { _demand = demandVector;     }
  void set_heliostatModel ( int                        hm           ) { _heliostatsFieldModel = hm; }

  bool set_heliostatFieldPowerOutput_MINCOST_TS   ( void );
  bool set_heliostatFieldPowerOutput_MAXCOMP_HTF1 ( void );
  
};

#endif
