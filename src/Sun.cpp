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
#include "Sun.hpp"

/*------------------------------------------------------*/
/*                    copy constructor                  */
/*------------------------------------------------------*/
Sun::Sun ( const Sun & sun) :
    _latitude            ( sun._latitude            ) ,
    _elevation           ( sun._elevation           ) ,
    _azimuth             ( sun._azimuth             ) ,
    _raysPerSquareMeters ( sun._raysPerSquareMeters ) ,
    _time                ( sun._time                ) ,
    _day                 ( sun._day                 ) ,
    _sunDeclination      ( sun._sunDeclination      )   {

  for ( size_t k = 0 ; k < sun._listOfSunrays.size(); ++k )
    _listOfSunrays.push_back ( new Sunray( *sun._listOfSunrays[k] ) );
    
  fComputeSunPosition();
}

/*------------------------------------------------------*/
/*                       destructor                     */
/*------------------------------------------------------*/
Sun::~Sun ( void ) {
  for ( size_t k = 0 ; k < _listOfSunrays.size(); ++k )
    delete _listOfSunrays[k];
}

/*------------------------------------------------------*/
void Sun::fComputeSunPosition ( void ) {
/*------------------------------------------------------*/

  //Assuming that the sun azimuth is due south at 12:00
  //Using  ESO's FITS convention where azimuth is measured from the south increasing towards the west Thus it is 0 at 12:00
  //and at this point elevation is 90 degrees - latitude
  //Neglecting the earth's precession axis inclination (thus supposing the earth is perfectly straight on its orbit

  //radians - put source. Calculating at summer solstice
  double HRA = (1.0*(_time.get_currentTime() - 720)*EARTH_OMEGA) * DEG_TO_RAD; // Local Solar Time
  double latitude_in_rad = _latitude * DEG_TO_RAD;
	
  _elevation = asin(sin(latitude_in_rad)*sin(_sunDeclination) + cos(latitude_in_rad)*cos(_sunDeclination)*cos(HRA)) * RAD_TO_DEG;
  long double azimuth_argument = (sin(_elevation * DEG_TO_RAD)*sin(latitude_in_rad) - sin(_sunDeclination )) /
    (cos(_elevation * DEG_TO_RAD)*cos(latitude_in_rad));
  if ( azimuth_argument > 1.0 )
    azimuth_argument = 1.0;
  else if (azimuth_argument < -1.0)
    azimuth_argument = -1.0;
  
  if ( HRA >= 0.0 )
    _azimuth = acos(azimuth_argument) * RAD_TO_DEG;
  else
    _azimuth = (- acos(azimuth_argument)) * RAD_TO_DEG;

  Sunray::set_azimuth   ( _azimuth   );
  Sunray::set_elevation ( _elevation );
}

/*------------------------------------------------------------------*/
void Sun::fAddNewSunray ( double Rmax, double thetaMax, double Z ) {
/*------------------------------------------------------------------*/
  
  long double fRandomR     = Rmax*sqrt(RNG::rand(0,1));
  long double fRandomTheta = (2.0 * thetaMax*RNG::rand(0,1) - thetaMax)*DEG_TO_RAD;
  long double fRandomZ     = Z * RNG::rand(0,1);

  _listOfSunrays.push_back ( new Sunray ( fRandomR*sin(fRandomTheta), -fRandomR*cos(fRandomTheta), fRandomZ ) );
}
