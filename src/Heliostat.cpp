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
#include "Heliostat.hpp"
#include "Sun.hpp"

double Heliostat::_width  = 0.0;
double Heliostat::_length = 0.0;
int    Heliostat::_IDmax  = 0;

/*-----------------------------------------------------*/
/*                     constructor                     */
/*-----------------------------------------------------*/
Heliostat::Heliostat ( double X                ,
		       double Y                ,
		       double Z                ,
		       double cosineEfficiency ,
		       double atmAttenuation   ,
		       double towerHeight        ) :
  _x                      ( X                ) ,
  _y                      ( Y                ) ,
  _z                      ( Z                ) ,
  _cosineEfficiency       ( cosineEfficiency ) ,
  _atmosphericAttenuation ( atmAttenuation   ) ,
  _sunraysCount           ( 0                )   {
  
  ++Heliostat::_IDmax;
  _ID = Heliostat::_IDmax;

  _cTopRight       = std::vector<double>(3, 0.0);
  _cTopLeft        = std::vector<double>(3, 0.0);
  _cBottomLeft     = std::vector<double>(3, 0.0);
  _cTopRightProj   = std::vector<double>(3, 0.0);
  _cTopLeftProj    = std::vector<double>(3, 0.0);
  _cBottomLeftProj = std::vector<double>(3, 0.0);
  _cTL_to_TR       = std::vector<double>(2, 0.0);
  _cTL_to_BL       = std::vector<double>(2, 0.0);
  _azimuthToAimpoint   = std::atan2(_x, -_y);
  _elevationToAimpoint = std::atan((towerHeight - _z) / std::sqrt(_x * _x + _y * _y));
}

/*--------------------------------------------------------------------------------*/
/*  Computes the projection of heliostat pilar position on a plane perpendicular  */
/*  to current radiations using the same axis convention as for sun position;     */
/*  Rotated frame is such that y' axis points directly towards the sun;           */
/*  Rotation is done following this order:                                        */
/*  1: rotate frame around z axis by an angle equivalent to the sun azimuth;      */
/*  2: rotate frame around x' axis by an angle equivalent to the sun elevation;   */
/*  Unit vectors for x', y' and z' are the following:                             */
/*    x' = ( cos(az), sin(az), 0 )                                                */
/*    y' = ( -cos(el)(sin(az), cos(el)cos(az), sin(el) )                          */
/*    z' = ( sin(el)sin(az), sin(el)cos(az), cos(el) )                            */
/*--------------------------------------------------------------------------------*/
void Heliostat::computePilarProjection ( const Sun & sun ) {

  long double el = sun.get_elevation() * DEG_TO_RAD;
  long double az = sun.get_azimuth() * DEG_TO_RAD;

  long double x_p0 = std::cos(az);
  long double x_p1 = std::sin(az);
	
  long double y_p0 = -std::cos(el)*std::sin(az);
  long double y_p1 =  std::cos(el)*std::cos(az);
  long double y_p2 =  std::sin(el);
  
  long double z_p0 =  std::sin(el)*std::sin(az);
  long double z_p1 = -std::sin(el)*std::cos(az);
  long double z_p2 =  std::cos(el);
  
  _xProj = _x * x_p0 + _y * x_p1;
  _yProj = _x * y_p0 + _y * y_p1 + _z * y_p2;
  _zProj = _x * z_p0 + _y * z_p1 + _z * z_p2;

  // DEBUG:
  // ------
  // if ( fabs(_x-14.1408)<1e-1  && fabs(_y+463.795)<1e-1 ) {
  //   std::cout << "HELIOSTAT #" << _ID << ":\n";
  //   std::cout << "\t(x,y,z)=" << _x << " " << _y << " " << _x << std::endl;
  //   std::cout << "\t(el,az)=" << el << " " << az << std::endl;
  //   std::cout << "\txp=" << x_p0 << " " << x_p1 << std::endl;
  //   std::cout << "\typ=" << y_p0 << " " << y_p1 << " " << y_p2 << std::endl;
  //   std::cout << "\tzp=" << z_p0 << " " << z_p1 << " " << z_p2 << std::endl;
  //   std::cout << "\tPROJ="<< _xProj << " " << _yProj << " " << _zProj << std::endl;
  // }
  
}

/*--------------------------------------------------------------------------------*/
/*  Calculates values of x, y and z for each of the 4 corners of the heliostat    */
/*  as a function of  the two angles of inclination (elevation and azimuth of     */
/*  the vector normal to its reflective surface)                                  */
/*--------------------------------------------------------------------------------*/
void Heliostat::computeCornersPositions ( void ) {

  _cTopRight[0] = _x + (_width / 2.0)*std::cos(_azimuth) + (_length / 2.0)*(std::sin(_elevation)*std::sin(_azimuth));
  _cTopRight[1] = _y + (_width / 2.0)*std::sin(_azimuth) - (_length / 2.0)*(std::sin(_elevation)*std::cos(_azimuth));
  _cTopRight[2] = _z + (_length / 2.0)*std::cos(_elevation);

  _cTopLeft[0] = _x - (_width / 2.0)*std::cos(_azimuth) + (_length / 2.0)*std::sin(_elevation)*std::sin(_azimuth);
  _cTopLeft[1] = _y - (_width / 2.0)*std::sin(_azimuth) - (_length / 2.0)*std::sin(_elevation)*std::cos(_azimuth);
  _cTopLeft[2] = _z + (_length / 2.0) * std::cos(_elevation);

  _cBottomLeft[0] = _x - (_width / 2.0)*std::cos(_azimuth) - (_length / 2.0) * std::sin(_elevation)*std::sin(_azimuth);
  _cBottomLeft[1] = _y - (_width / 2.0)*std::sin(_azimuth) + (_length / 2.0) * std::sin(_elevation)*std::cos(_azimuth);
  _cBottomLeft[2] = _z - (_length / 2.0)*std::cos(_elevation);	
}

/*--------------------------------------------------------------*/
void Heliostat::computeCornersProjections ( const Sun & sun ) {
/*--------------------------------------------------------------*/

  computeCornersPositions();
  double el = sun.get_elevation() * DEG_TO_RAD;
  double az = sun.get_azimuth() * DEG_TO_RAD;

  std::vector<double > x_p(3, 0.0);
  std::vector<double > y_p(3, 0.0);
  std::vector<double > z_p(3, 0.0);

  x_p[0] = std::cos(az);
  x_p[1] = std::sin(az);

  y_p[0] = -std::cos(el)*std::sin(az);
  y_p[1] =  std::cos(el)*std::cos(az);
  y_p[2] =  std::sin(el);

  z_p[0] =  std::sin(el)*std::sin(az);
  z_p[1] = -std::sin(el)*std::cos(az);
  z_p[2] =  std::cos(el);

  _cTopRightProj[0] = _cTopRight[0] * x_p[0] + _cTopRight[1] * x_p[1];
  _cTopRightProj[1] = _cTopRight[0] * y_p[0] + _cTopRight[1] * y_p[1] + _cTopRight[2] * y_p[2];
  _cTopRightProj[2] = _cTopRight[0] * z_p[0] + _cTopRight[1] * z_p[1] + _cTopRight[2] * z_p[2];

  _cTopLeftProj[0] = _cTopLeft[0] * x_p[0] + _cTopLeft[1] * x_p[1];
  _cTopLeftProj[1] = _cTopLeft[0] * y_p[0] + _cTopLeft[1] * y_p[1] + _cTopLeft[2] * y_p[2];
  _cTopLeftProj[2] = _cTopLeft[0] * z_p[0] + _cTopLeft[1] * z_p[1] + _cTopLeft[2] * z_p[2];
  
  _cBottomLeftProj[0] = _cBottomLeft[0] * x_p[0] + _cBottomLeft[1] * x_p[1];
  _cBottomLeftProj[1] = _cBottomLeft[0] * y_p[0] + _cBottomLeft[1] * y_p[1] + _cBottomLeft[2] * y_p[2];
  _cBottomLeftProj[2] = _cBottomLeft[0] * z_p[0] + _cBottomLeft[1] * z_p[1] + _cBottomLeft[2] * z_p[2];

  // update collision vectors:
  _cTL_to_TR[0] = _cTopRightProj[0] - _cTopLeftProj[0];
  _cTL_to_TR[1] = _cTopRightProj[2] - _cTopLeftProj[2];

  _cTL_to_BL[0] = _cBottomLeftProj[0] - _cTopLeftProj[0];
  _cTL_to_BL[1] = _cBottomLeftProj[2] - _cTopLeftProj[2];
}

/*----------------------------------------------------------------------*/
void Heliostat::computeAngles ( const Sun & sun, double towerHeight ) {
/*----------------------------------------------------------------------*/
  double az_to_receiver = -std::atan2(-_x, -_y);
  double el_to_receiver =  std::atan2(towerHeight - (_length / 2.0), std::sqrt(_x*_x + _y*_y));
  _azimuth   = (az_to_receiver + sun.get_azimuth()*DEG_TO_RAD) / 2.0;
  _elevation = (el_to_receiver + sun.get_elevation()*DEG_TO_RAD) / 2.0;
}

/*------------------------------------------------------------------------------------------*/
double Heliostat::fComputeSpillage ( double apertureHeight, double apertureWidth ) const {
/*------------------------------------------------------------------------------------------*/
  double noSpillageRatio = 0.0;
  double projectedWidth  = std::abs(_cTopLeftProj[0] - _cTopRightProj[0]);
  double projectedHeight = std::abs(_cTopLeftProj[2] - _cBottomLeftProj[2]);

  if ( projectedWidth > apertureWidth * std::cos(_azimuthToAimpoint) )
    noSpillageRatio = (apertureWidth*std::cos(_azimuthToAimpoint))/projectedWidth;
  else
    noSpillageRatio = 1.0;
  if ( projectedHeight > apertureHeight*std::cos(_elevationToAimpoint) )
    noSpillageRatio *= (apertureHeight*std::cos(_elevationToAimpoint))/projectedHeight;
  else
    noSpillageRatio *= 1.0;
  return noSpillageRatio;
}
