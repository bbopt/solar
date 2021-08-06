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
#ifndef __HELIOSTATFIELD_H__
#define __HELIOSTATFIELD_H__

#include "Heliostat.hpp"
#include "Global.hpp"
#include "Sun.hpp"

#include <list>
#include <vector>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

class HeliostatField {

private:

  size_t                             _nbOfHeliostats;
  double                             _heliostatLength;
  double                             _heliostatWidth;
  double                             _towerHeight;          // height of aim point
  double                             _apertureHeight;
  double                             _apertureWidth;
  double                             _minDistanceFromTower; // as a multiple of tower height
  double                             _maxDistanceFromTower; // as a multiple of tower height
  double                             _maxAngularDeviation;  // degrees
  std::vector<Heliostat*>            _listOfHeliostats;
  std::vector< std::vector<double> > _gridLayoutAngularCoordinates;
  std::vector< std::vector<double> > _gridLayoutCartesianCoordinates;
  std::vector< std::vector<double> > _fieldLayoutCartesian;
  std::vector<double>                _fieldEfficiency;
  std::vector<double>                _powerOutput;
  Sun                                _sun;
  
  void   fComputeStaggeredGridLayout    ( void ); 
  void   fComputeEfficiency             ( void );
  void   fComputeAtmosphericAttenuation ( void );
  void   fComputeCosineAndSpillage      ( void );		
  void   fConfigureField                ( void );
  double fComputeFieldEfficiency        ( void );
  double fEvaluateFieldSurface          ( void ) const;
  
  void   delete_heliostats              ( void );

public:
  
  HeliostatField  ( size_t, double, double, double, double, double, double, double, double, const Sun & );
  ~HeliostatField ( void ) { delete_heliostats(); }

  void   fGenerateField              ( void );
  double fCalculateTotalEnergyOutput ( void );
  void   fGenerateSunrays            ( void );
  void   fTimeIncrement              ( void ) { _sun.fTimeIncrement(); }
 
  double get_heliostatLength      ( void ) const { return _heliostatLength; }
  double get_heliostatWidth       ( void ) const { return _heliostatWidth; }
  double get_towerHeight          ( void ) const { return _towerHeight; }
  double get_minDistanceFromTower ( void ) const { return _minDistanceFromTower; }
  double get_maxDistanceFromTower ( void ) const { return _maxDistanceFromTower; }
  double get_maxAngularDeviation  ( void ) const { return _maxAngularDeviation; }
  size_t get_nb_heliostats        ( void ) const { return _listOfHeliostats.size(); }
};

bool compareDistanceToSun ( const Heliostat * h1, const Heliostat* h2 );
bool comparePositions     ( const std::vector<double> & firstPosition, const std::vector<double> & secondPosition );

#endif
