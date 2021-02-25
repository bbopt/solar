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
#ifndef __THERMAL_STORAGE_H__
#define __THERMAL_STORAGE_H__

#include "MoltenSalt.hpp"
#include "Constants.hpp"
#include "Global.hpp"
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <sstream>

class ThermalStorage {
	
private:
  double            _storedMass;
  double            _storedTemperature;
  const MoltenSalt* _inputHTF;
  MoltenSalt      * _outputHTF;

  //Dimensions parameters in meters
  double _heightOfStorage;
  double _diameterOfStorage;
  double _thicknessOfInsulation;
	
  //height of volume stored
  double _heightOfVolumeStored;

public:

  ThermalStorage ( const MoltenSalt* input, MoltenSalt* output, double height, double diameter, double thickness ) :
    _storedMass            ( 0.0       ) ,
    _storedTemperature     ( 0.0       ) ,
    _inputHTF              ( input     ) ,
    _outputHTF             ( output    ) ,
    _heightOfStorage       ( height    ) ,
    _diameterOfStorage     ( diameter  ) ,
    _thicknessOfInsulation ( thickness ) ,
    _heightOfVolumeStored  ( 0.0       )   {}
  
  double fInitialStorageTemperature ( int ) const;
  double fComputeStorageTemperature ( int );

  double fComputeStorageMass ( int timeInterval ) const {
    return _storedMass + timeInterval * 60 * (_inputHTF->get_massFlow() - _outputHTF->get_massFlow());
  }

  double fComputeStorageLevel ( void );

  double fComputeStorageLevel ( double storedMass ) const {
    return (storedMass / MS_DENSITY) / (PI*pow(0.5*_diameterOfStorage, 2.0));
  }

  double fInitialStorageMass ( int timeInterval ) const {
    return _storedMass + timeInterval * 60 * _inputHTF->get_massFlow();
  }
  
  double fComputeEnergyLosses          ( double, double);
  double fComputeInsideRadiationLosses ( void );

  double fSolveForT   ( double, double, double, double, double, double ) const;
  double fSolveForT_i ( double, double, double, double, double, double ) const;

  void set_storage  ( double, double );
  void set_storage2 ( double, double );

  double get_heightOfVolumeStored  ( void ) const { return _heightOfVolumeStored; }
  double get_storedMass            ( void ) const { return _storedMass; }
  double get_storedTemperature     ( void ) const { return _storedTemperature; }
  double get_heightOfStorage       ( void ) const { return _heightOfStorage; }
  double get_diameterOfStorage     ( void ) const { return _diameterOfStorage; }
  double get_thicknessOfInsulation ( void ) const { return _thicknessOfInsulation; }
  const MoltenSalt * get_inputHTF  ( void ) const { return _inputHTF; }
  const MoltenSalt * get_outputHTF ( void ) const { return _outputHTF; }
};

#endif
