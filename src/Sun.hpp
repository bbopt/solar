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
#ifndef __SUN_H__
#define __SUN_H__

#include "Time_Manager.hpp"
#include "Sunray.hpp"
#include "RNG.hpp"
#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>

class Sun {

private:

  double               _latitude;
  long double          _elevation;
  long double          _azimuth;
  double               _raysPerSquareMeters;
  Time_Manager         _time;
  int                  _day;
  double               _sunDeclination;
  std::vector<Sunray*> _listOfSunrays;

public:

  Sun ( double latitude, const Time_Manager & time, int day, double raysPerSquareMeters ) :
    _latitude            ( latitude            ) ,
    _elevation           ( 0.0                 ) ,
    _azimuth             ( 0.0                 ) ,
    _raysPerSquareMeters ( raysPerSquareMeters ) ,
    _time                ( time                ) ,
    _day                 ( day                 ) ,
    _sunDeclination      ( DECLINATION * DEG_TO_RAD * sin((360 * (_day + 284) / 365)*DEG_TO_RAD) ) { fComputeSunPosition(); }

  Sun  ( const Sun & sun);
  ~Sun ( void );

  void fComputeSunPosition ( void );
  void fAddNewSunray       ( double Rmax, double thetaMax, double Z );

  double get_elevation           ( void ) const { return _elevation;           }
  double get_azimuth             ( void ) const { return _azimuth;             }
  double get_raysPerSquareMeters ( void ) const { return _raysPerSquareMeters; }
  
  void fTimeIncrement ( void ) { _time.fTimeIncrement(); }
  void fResetTime     ( void ) { _time.fResetTime();     }

  int get_incrementsCounter  ( void ) const { return _time.get_incrementsCounter (); }
  int get_numberOfIncrements ( void ) const { return _time.get_numberOfIncrements(); }

  void projectTarget     ( int i ) { _listOfSunrays[i]->projectTarget(); }
  
  void set_isIntercepted ( int i, bool b ) { _listOfSunrays[i]->set_isIntercepted(b); }
  
  bool computeCollision  ( int i, const Heliostat & h ) { return _listOfSunrays[i]->computeCollision(h); }
  
  size_t get_nb_sunrays ( void ) const { return _listOfSunrays.size(); }
  
};

#endif
