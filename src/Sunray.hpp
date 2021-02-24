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
#ifndef __SUNRAY_H__
#define __SUNRAY_H__

#include "Heliostat.hpp"
#include <vector>
#include <iostream>

class Sunray {

private:

  static double _azimuth;
  static double _elevation;
  static double _minDistance;
  
  double _xTarget;
  double _yTarget;
  double _zTarget;
  
  std::vector<double> _projectedTarget;
  
  bool _isIntercepted;
  int  _interceptedBy; 

public:

  Sunray ( double x, double y, double z );

  Sunray ( const Sunray & s ) :
    _xTarget         ( s._xTarget         ) ,
    _yTarget         ( s._yTarget         ) ,
    _zTarget         ( s._zTarget         ) ,
    _projectedTarget ( s._projectedTarget ) ,
    _isIntercepted   ( s._isIntercepted   ) ,
    _interceptedBy   ( s._interceptedBy   )   {}

  void projectTarget    ( void );
  bool computeCollision ( const Heliostat & heliostat );

  void set_xTarget       ( double x  ) { _xTarget = x; }
  void set_yTarget       ( double y  ) { _yTarget = y; }
  void set_zTarget       ( double z  ) { _zTarget = z; }
  void set_isIntercepted ( bool   is ) { _isIntercepted = is; }
  void set_interceptedBy ( int    id ) { _interceptedBy = id; }

  static void set_azimuth     ( double az) { Sunray::_azimuth   = az; }
  static void set_elevation   ( double el) { Sunray::_elevation = el; }
  static void set_minDistance ( void     ) {
    Sunray::_minDistance = sqrt(pow(Heliostat::get_width(), 2.0) + pow(Heliostat::get_length(), 2.0));
  }
  
  double get_xTarget ( void ) const { return _xTarget; }
  double get_yTarget ( void ) const { return _yTarget; }
  double get_zTarget ( void ) const { return _zTarget; }

  bool get_isIntercepted ( void ) const { return _isIntercepted; }
  int  get_interceptedBy ( void ) const { return _interceptedBy; }
};

#endif
