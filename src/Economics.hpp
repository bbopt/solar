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
#ifndef __ECONOMICS_H__
#define __ECONOMICS_H__

#include "Global.hpp"

class Economics {
  
private:

  int    _nbOfHeliostats;
  double _hotStorageInsulationThickness;
  double _coldStorageInsulationThickness;
  double _hotStorageHeight;
  double _coldStorageHeight;   // New in V2 (P.B., 2025-07)
  double _hotStorageDiameter;
  double _coldStorageDiameter; // New in V2 (SLD, 2025-07-24)
  double _heightOfTower;
  double _lengthOfHeliostats;
  double _widthOfHeliostats;
  double _turbineNominalPowerOutput;
  double _heightOfReceiverAperture;
  double _widthOfReceiverAperture;  
  int    _receiverNumberOfTubes;
  double _receiverTubesDout;
  double _reflectiveArea;
  int    _exchangerModel;
  double _exchangerTubesOutterDiameter;
  double _exchangerTubesLength;
  int    _exchangerNumberOfTubes;
  int    _exchangerTubePassesPerShell;
  int    _exchangerNumberOfShell;
  double _costOfField;
  double _costPerHeliostat;
  double _costOfTower;
  double _costOfStorage;
  double _costOfPowerblock;
  double _costOfSteamGenerator;
  double _costOfReceiver;
  double _totalCost;

public:
  
  // constructor:
  Economics ( int    numberOfHeliostats          ,
	      double hotStorageInsul             ,
	      double coldStorageInsul            ,
	      double hotStorageHeight            ,
	      double coldStorageHeight           , // new in V2 (P.B)
	      double hotStorageDiameter          ,
	      double coldStorageDiameter         , // new in V2 (SLD)
	      double heightOfTower               ,
	      double lengthOfHeliostats          ,
	      double widthOfHeliostats           ,
	      double turbinePower                ,
	      double heightOfReceiverAperture    ,
	      double widthOfAperture             ,
	      int    receiverNbOfTubes           , // new in V2 (SLD)
	      double receiverTubesOutsideDiam    , // new in V2 (SLD)
	      int    exchangerModel              , // new in V2 (SLD)
	      double exchangerTubesDiameter      ,
	      double exchangerTubesLength        , 
	      int    exchangerNumberOfTubes      ,
	      int    exchangerTubePassesPerShell ,
	      int    exchangerNumberOfShell        );

  // compute costs:
  double evaluateCostOfField          ( void );
  double evaluateCostOfHeliostat      ( void );
  double evaluateCostOfTower          ( void );
  double evaluateCostOfStorage        ( void );
  double evaluateCostOfReceiver       ( void );
  double evaluateCostOfPowerblock     ( void );
  double evaluateCostOfSteamGenerator ( void );
  double evaluateTotalInvestmentCost  ( void );
  
  // set methods:
  void set_nbOfHeliostats ( int    n ) { _nbOfHeliostats = n; }
  void set_reflectiveArea ( double r ) { _reflectiveArea = r; }

  // get methods:
  double get_costOfField          ( void ) const { return _costOfField;          }
  double get_costOfTower          ( void ) const { return _costOfTower;          }
  double get_costOfStorage        ( void ) const { return _costOfStorage;        }
  double get_costOfReceiver       ( void ) const { return _costOfReceiver;       }
  double get_costOfPowerblock     ( void ) const { return _costOfPowerblock;     }
  double get_costOfSteamGenerator ( void ) const { return _costOfSteamGenerator; }

};

#endif
