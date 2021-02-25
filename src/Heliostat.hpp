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
#ifndef __HELIOSTAT_H__
#define __HELIOSTAT_H__

#include "Constants.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

class Sun;
class Sunray;

class Heliostat {
  
private:
  
  static int _IDmax;
  int        _ID;

  static double _width;
  static double _length;

  // Position of the pilar:
  double _x;
  double _y;
  double _z;

  // Projected position of the pilar:
  double _xProj;
  double _yProj;
  double _zProj;

  // Projected position of heliostats corners
  // on a plane perpendicular to sun rays:
  std::vector<double> _cTopLeft;
  std::vector<double> _cTopRight;
  std::vector<double> _cBottomLeft;
  std::vector<double> _cBottomRight;

  std::vector<double> _cTopLeftProj;
  std::vector<double> _cTopRightProj;
  std::vector<double> _cBottomLeftProj;
  std::vector<double> _cBottomRightProj;

  // Projected vectors representing the top and left sides of the heliostat
  // as seen from the projection plan:
  std::vector<double> _cTL_to_TR;
  std::vector<double> _cTL_to_BL;

  // Optical efficiency factors:
  double _cosineEfficiency;
  double _atmosphericAttenuation;
  
  // Orientation (horizontal coordinates)
  // For optimal reflection, the flat mirror should always be pointing directly half way
  // between the sun and the receptor
  double _azimuth;
  double _elevation;

  double _azimuthToAimpoint;
  double _elevationToAimpoint;

  int _sunraysCount;

  Heliostat ( const Heliostat & );
  
public:

  Heliostat  ( double, double, double, double, double, double );
  ~Heliostat ( void ) {}

  void computePilarProjection    ( const Sun & );
  void computeCornersPositions   ( void );
  void computeCornersProjections ( const Sun & );
  void computeAngles             ( const Sun &, double );

  int    get_ID    ( void ) const { return _ID; }
  double get_x     ( void ) const { return _x;  }
  double get_y     ( void ) const { return _y;  }
  double get_z     ( void ) const { return _z;  }
  double get_xProj ( void ) const { return _xProj; }
  double get_yProj ( void ) const { return _yProj; }
  double get_zProj ( void ) const { return _zProj; }

  double get_cosineEfficiency       ( void ) const { return _cosineEfficiency; }
  double get_atmosphericAttenuation ( void ) const { return _atmosphericAttenuation; }
  double get_azimuth                ( void ) const { return _azimuth;   }
  double get_elevation              ( void ) const { return _elevation; }
  int    get_sunraysCount           ( void ) const { return _sunraysCount; }
  
  double get_cTopLeftProj    ( int i ) const { return _cTopLeftProj[i];    }
  double get_cTopRightProj   ( int i ) const { return _cTopRightProj[i];   }
  double get_cBottomLeftProj ( int i ) const { return _cBottomLeftProj[i]; }
  double get_cTL_to_TR       ( int i ) const { return _cTL_to_TR[i];       }
  double get_cTL_to_BL       ( int i ) const { return _cTL_to_BL[i];       }
 
  void increase_sunraysCount ( void ) { ++_sunraysCount;   }
  void clear_sunraysCount    ( void ) { _sunraysCount = 0; }
 
  double fComputeSpillage ( double, double ) const;

  static void   set_width  ( double w ) { _width  = w;    }
  static void   set_length ( double l ) { _length = l;    }
  static double get_width  ( void     ) { return _width;  }
  static double get_length ( void     ) { return _length; }
};

#endif
