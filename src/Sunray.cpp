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
#include "Sunray.hpp"

double Sunray::_azimuth     = 0.0;
double Sunray::_elevation   = 0.0;
double Sunray::_minDistance = 0.0;

/*------------------------------------------------------------------*/
/*                              constructor                         */
/*------------------------------------------------------------------*/
Sunray::Sunray ( double x, double y, double z ) :
  _xTarget(x), _yTarget(y), _zTarget(z), _isIntercepted(false), _interceptedBy(0) {

  _projectedTarget = std::vector<double>(3, 0.);

  long double x_p0 = cos(_azimuth*DEG_TO_RAD);
  long double x_p1 = sin(_azimuth*DEG_TO_RAD);

  long double y_p0 = -cos(_elevation*DEG_TO_RAD)*sin(_azimuth*DEG_TO_RAD);
  long double y_p1 =  cos(_elevation*DEG_TO_RAD)*cos(_azimuth*DEG_TO_RAD);
  long double y_p2 =  sin(_elevation*DEG_TO_RAD);

  long double z_p0 =  sin(_elevation*DEG_TO_RAD)*sin(_azimuth*DEG_TO_RAD);
  long double z_p1 = -sin(_elevation*DEG_TO_RAD)*cos(_azimuth*DEG_TO_RAD);
  long double z_p2 =  cos(_elevation*DEG_TO_RAD);

  _projectedTarget[0] = _xTarget * x_p0 + _yTarget * x_p1;
  _projectedTarget[1] = _xTarget * y_p0 + _yTarget * y_p1 + _zTarget * y_p2;
  _projectedTarget[2] = _xTarget * z_p0 + _yTarget * z_p1 + _zTarget * z_p2;
}

/*------------------------------------------------------------------*/
bool Sunray::computeCollision ( const Heliostat & heliostat ) {
/*------------------------------------------------------------------*/

  if ( !_isIntercepted ) {
     
    if (fabs(_projectedTarget[0] - heliostat.get_xProj()) < _minDistance
	&& fabs(_projectedTarget[2] - heliostat.get_zProj()) < _minDistance) {
      double target0 = 0.0 , target1 = 0.0;
      double u_1, u_2, detTopLeft;
	  
      target0 = _projectedTarget[0] - heliostat.get_cTopLeftProj(0);
      target1 = _projectedTarget[2] - heliostat.get_cTopLeftProj(2);
	  
      detTopLeft = heliostat.get_cTL_to_TR(0) * heliostat.get_cTL_to_BL(1)
	- heliostat.get_cTL_to_TR(1) * heliostat.get_cTL_to_BL(0);

      if (detTopLeft != 0.) {
	u_2 = (target1 * heliostat.get_cTL_to_TR(0) - target0 * heliostat.get_cTL_to_TR(1)) / detTopLeft;
	u_1 = (target0 - u_2*heliostat.get_cTL_to_BL(0)) / heliostat.get_cTL_to_TR(0);
	if ( u_1 >= 0. && u_1 <= 1. && u_2 >= 0. && u_2 <= 1.0 ) {
	  _isIntercepted = true;
	  _interceptedBy = heliostat.get_ID();
	  return true;
	}
      }
    }
  }
  return false;
}

/*------------------------------------------------------------------*/
void Sunray::projectTarget ( void ) {
/*------------------------------------------------------------------*/
  
  long double cosElev = cos(_elevation*DEG_TO_RAD);
  long double cosAzm  = cos(_azimuth*DEG_TO_RAD);
  long double sinElev = sin(_elevation*DEG_TO_RAD);
  long double sinAzm  = sin(_azimuth*DEG_TO_RAD);
  
  _projectedTarget[0] = _xTarget *  cosAzm + _yTarget * sinAzm;
  _projectedTarget[1] = _xTarget * (-cosElev*sinAzm) + _yTarget * cosElev*cosAzm + _zTarget * sinElev;
  _projectedTarget[2] = _xTarget * sinElev*sinAzm + _yTarget * (-sinElev*cosAzm) + _zTarget * cosElev;
}
