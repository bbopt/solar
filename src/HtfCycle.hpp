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
#ifndef __HTF_CYCLE_H__
#define __HTF_CYCLE_H__

#include "CentralReceiver.hpp"
#include "ThermalStorage.hpp"
#include "HeatExchanger.hpp"
#include "MoltenSalt.hpp"
#include "Constants.hpp"
#include "Powerblock.hpp"
#include "Global.hpp"

#include <cmath>

class HtfCycle {

private:
  
  CentralReceiver _centralReceiver;
  ThermalStorage  _hotStorage;
  ThermalStorage  _coldStorage;
  HeatExchanger   _steamGenerator;
  MoltenSalt      _centralReceiverInlet;
  MoltenSalt      _centralReceiverOutlet;
  MoltenSalt      _steamGeneratorInlet;
  MoltenSalt      _steamGeneratorOutlet;
  int             _timeInterval; //minutes

  // Data gathering
  std::vector<double> _steamGenOutletMsRate;
  std::vector<double> _steamGenOutletTemp;
  double              _minColdStorageTemp;
  double              _minHotStorageTemp;
  double              _minSteamGenOutletTemp;
  std::vector<double> _storageHeat;

public:

  HtfCycle ( double       receiverTemp           ,
	     double       storageHeight          ,
	     double       storageDiameter        ,
	     double       insulationThickness    ,
	     double       exchangerExitTemp      ,
	     Powerblock * powerblock             ,
	     double       apertureHeight         ,
	     double       apertureWidth          ,
	     double       receiverTubesDin       ,
	     double       receiverTubesThickness ,
	     int          receiverNbTubes        ,
	     int          timeInterval             );

  HtfCycle ( double       receiverTemp                ,
	     double       storageHeight               ,
	     double       storageDiameter             ,
	     double       insulationThickness         ,
	     double       exchangerExitTemp           ,
	     Powerblock * powerblock                  ,
	     double       apertureHeight              ,
	     double       apertureWidth               ,
	     double       receiverTubesDin            ,
	     double       receiverTubesThickness      ,
	     int          receiverNbTubes             ,
	     int          timeInterval                ,
	     double       exchangerTubesDin           ,
	     double       exchangerTubesDout          ,
	     double       exchangerTubesLength        ,
	     double       exchangerTubesSpacing       ,
	     double       baffleCut                   ,
	     int          nbOfBaffles                 ,
	     int          exchangerNbOfTubes          ,
	     int          exchangerNbOfPassesPerShell ,
	     int          exchangerNbOfShells           );

  ~HtfCycle ( void ) {}

  double compute_CR_YieldPressure   ( void ) const { return _centralReceiver.computeYieldPressure();   }
  double compute_CR_PressureInTubes ( void ) const { return _centralReceiver.computePressureInTubes(); }
  
  double compute_SG_YieldPressure   ( void     ) const { return _steamGenerator.computeYieldPressure();    }
  double compute_SG_PressureInTubes ( double s ) const { return _steamGenerator.computePressureInTubes(s); }
  int    get_exchangerModel         ( void     ) const { return _steamGenerator.get_exchangerModel();      }
  double computePressureInShells    ( void     ) const { return _steamGenerator.computePressureInShells(); }
  
  const	ThermalStorage & get_hotStorage            ( void ) const { return _hotStorage;            }
  const	MoltenSalt     & get_centralReceiverOutlet ( void ) const { return _centralReceiverOutlet; }
  const	MoltenSalt     & get_steamGeneratorInlet   ( void ) const { return _steamGeneratorInlet;   }
  const	MoltenSalt     & get_steamGeneratorOutlet  ( void ) const { return _steamGeneratorOutlet;  }
  
  const	std::vector<double> & get_steamGenOutletMsRate ( void ) const { return _steamGenOutletMsRate; }
  const	std::vector<double> & get_steamGenOutletTemp   ( void ) const { return _steamGenOutletTemp;   }
  const	std::vector<double> & get_storageHeatV         ( void ) const { return _storageHeat;          }
  
  double get_storageHeat ( int i ) const { return _storageHeat[i]; }
  
  double get_minColdStorageTemp ( void ) const { return _minColdStorageTemp;    }
  double get_minHotStorageTemp  ( void ) const { return _minHotStorageTemp;     }
  double get_minSteamGenTemp    ( void ) const { return _minSteamGenOutletTemp; }
  
  void initiateColdStorage ( void );
  void setStorage          ( double, double, double );
  void fOperateCycle       ( int    timeInSeconds       ,
			     double energyFromField     ,
			     double requiredPowerOutput ,
			     bool   low_fid               );

};

#endif
