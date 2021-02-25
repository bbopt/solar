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

  double _hotStorageInsulationThickness;
  double _coldStorageInsulationThickness;
  double _hotStorageHeight;
  double _storageDiameter;
  double _receiverInsulationThickness;
  double _heightOfTower;
  double _heightOfReceiverAperture;
  double _widthOfReceiverAperture;
  int    _receiverNumberOfTubes;
  double _receiverTubesDout;
  double _lengthOfHeliostats;
  double _widthOfHeliostats;
  int    _nbOfHeliostats;
  double _reflectiveArea;
  double _totalMoltenSaltMass;
  double _turbineNominalPowerOutput;
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

  // heliostats Field only:
  Economics ( void );
  
  // full plant:
  Economics ( int    numberOfHeliostats          ,
	      double hotStorageInsul             ,
	      double coldStorageInsul            ,
	      double hotStorageHeight            ,
	      double storageDiameter             ,
	      double receiverInsul               ,
	      double heightOfTower               ,
	      double lengthOfHeliostats          ,
	      double widthOfHeliostats           ,
	      double turbinePower                ,
	      double heightOfReceiverAperture    ,
	      double widthOfAperture             ,
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
  void set_heightOfTower                 ( double x ) { _heightOfTower                  = x; }
  void set_heightOfReceiverAperture      ( double x ) { _heightOfReceiverAperture       = x; }
  void set_widthOfReceiverAperture       ( double x ) { _widthOfReceiverAperture        = x; }
  void set_lengthOfHeliostats            ( double x ) { _lengthOfHeliostats             = x; }
  void set_widthOfHeliostats             ( double x ) { _widthOfHeliostats              = x; }
  void set_exchangerModel                ( int    x ) { _exchangerModel                 = x; }
  void set_hotStorageInsulationThickness ( double x ) { _hotStorageInsulationThickness  = x; }
  void set_coldStorageInsulationThickness( double x ) { _coldStorageInsulationThickness = x; }
  void set_hotStorageHeight              ( double x ) { _hotStorageHeight               = x; }
  void set_storageDiameter               ( double x ) { _storageDiameter                = x; }
  void set_receiverInsulationThickness   ( double x ) { _receiverInsulationThickness    = x; }
  void set_receiverNumberOfTubes         ( int    x ) { _receiverNumberOfTubes          = x; }
  void set_receiverTubesDout             ( double x ) { _receiverTubesDout              = x; }
  void set_nbOfHeliostats                ( int    x ) { _nbOfHeliostats                 = x; }
  void set_reflectiveArea                ( double x ) { _reflectiveArea                 = x; }
  void set_totalMoltenSaltMass           ( double x ) { _totalMoltenSaltMass            = x; }
  void set_turbineNominalPowerOutput     ( double x ) { _turbineNominalPowerOutput      = x; }
  void set_exchangerTubesOutterDiameter  ( double x ) { _exchangerTubesOutterDiameter   = x; }
  void set_exchangerTubesLength          ( double x ) { _exchangerTubesLength           = x; }
  void set_exchangerNumberOfTubes        ( int    x ) { _exchangerNumberOfTubes         = x; }
  void set_exchangerTubePassesPerShell   ( int    x ) { _exchangerTubePassesPerShell    = x; }
  void set_exchangerNumberOfShell        ( int    x ) { _exchangerNumberOfShell         = x; }

  // get methods:
  double get_costOfField          ( void ) const { return _costOfField;          }
  double get_costOfTower          ( void ) const { return _costOfTower;          }
  double get_costOfStorage        ( void ) const { return _costOfStorage;        }
  double get_costOfReceiver       ( void ) const { return _costOfReceiver;       }
  double get_costOfPowerblock     ( void ) const { return _costOfPowerblock;     }
  double get_costOfSteamGenerator ( void ) const { return _costOfSteamGenerator; }

};

#endif
