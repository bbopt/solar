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
#include "HeliostatField.hpp"

/*------------------------------------------------------------------------*/
/*                              constructor                               */
/*------------------------------------------------------------------------*/
HeliostatField::HeliostatField ( size_t      nbOfHeliostats       ,
				 double      heliostatLength      ,
				 double      heliostatWidth       ,
				 double      towerHeight          ,
				 double      apertureHeight       ,
				 double      apertureWidth        ,
				 double      minDistanceFromTower ,
				 double      maxDistanceFromTower ,
				 double      maxAngle             ,
				 const Sun & sun                    ) :
  _nbOfHeliostats       ( nbOfHeliostats       ) ,
  _heliostatLength      ( heliostatLength      ) ,
  _heliostatWidth       ( heliostatWidth       ) ,
  _towerHeight          ( towerHeight          ) ,
  _apertureHeight       ( apertureHeight       ) ,
  _apertureWidth        ( apertureWidth        ) ,
  _minDistanceFromTower ( minDistanceFromTower ) ,
  _maxDistanceFromTower ( maxDistanceFromTower ) ,
  _maxAngularDeviation  ( maxAngle             ) ,
  _sun                  ( sun                  )   {

  fComputeStaggeredGridLayout();
  fComputeEfficiency();
  fConfigureField();

  Heliostat::set_width(heliostatWidth);
  Heliostat::set_length(heliostatLength);
  Sunray::set_minDistance();

  _powerOutput.reserve     ( 24 );
  _fieldEfficiency.reserve ( 24 );
}

void HeliostatField::delete_heliostats ( void ) {
  for ( size_t k = 0 ; k < _listOfHeliostats.size() ; ++k )
    delete _listOfHeliostats[k];
  _listOfHeliostats.clear();
}

/*---------------------------------------------------------*/
/*  Determine the coordinates of all potential heliostats  */
/*  within the boundaries of the field                     */
/*---------------------------------------------------------*/
void HeliostatField::fComputeStaggeredGridLayout ( void ) {

  // initial check: R_min cannot be equal to zero:
  double R_min = _minDistanceFromTower * _towerHeight;  // Hidden constraint: SOLAR1: x3*x8 != 0.0

  if ( R_min == 0.0 )
    throw Simulation_Interruption ( "R_min equal to zero" );
  
  double z_0     = _heliostatLength / 2.0; // heigth of heliostat center from the ground
  double l_m     = _heliostatLength;
  double w_m     = _heliostatWidth;
  double H_t     = _towerHeight;
  double l_r     = l_m;
  double Phi_max = _maxAngularDeviation * DEG_TO_RAD;
  double R_max   = _maxDistanceFromTower * _towerHeight;
  double Beta_L  = 0.0; // terrain slope will be assumed to be zero

  double f      = w_m / l_m; //ratio of heliostat width over length
  double z_1    = H_t - 0.5*l_r;
  double r_m    = 0.5*l_m;
  double dS_min = 2.0*f - sqrt(1.0 + pow(f, 2.0)); //eq 2
  double dS     = dS_min;
  
  // Now that all variables have been identified, do the following sequence:
  // 1- Determine the position of all heliostats on the first ring
  // 2- Determine the minimum radial spacing
  // 3- Determine the radius of the first staggered ring.
  //    If the difference between the 2nd radius and the first radius is smaller
  //    than the minimal radial spacing, define the new radius as R_0 + min_delta_R
  // 4- Determine the position of all heliostats on the 2nd ring
  // 5- Determine the next two radiuses using the no blocking method.
  //    Check that the difference between each radius is always bigger than the minimal radial spacing
  // 6- Check if filling the 3rd radius using the current angular spacing will amount to
  //    a better heliostat/field area than skipping the 3rd radius and filling the 4th completely 
  //   as a new essential ring.

  double DM = ( l_m * (sqrt(1.0 + pow(f, 2.0)) + dS) >= 2.0*w_m) ?
    l_m * (sqrt(1. + pow(f, 2.0)) + dS) :
    2.0*w_m;

  long double tmp1 = 30.0*PI / 180.0;
  long double tmp2 = Beta_L* (PI / 180.0);
	
  double Delta_r_min = DM * cos(tmp1) * cos(tmp2);

  int i = 0, j = 0, n = 0, im1 = 0, np1 = 0;
	
  std::vector< std::vector<double> > R;
  std::vector<double> gamma;

  // The radius for the first ring is defined already by R_min
  R.push_back(std::vector<double>(1, R_min));
  gamma.push_back(DM / (2.0*R[0][0]));
	
  //Determine the number of heliostats on each ring of this group
  int N_essential, N_staggered, N_new_group, N_new_circle;

  N_essential = 1 + 2 * myround(floor(Phi_max / (2.0 * gamma[j])));
  N_staggered = 2 + 2 * myround(floor((Phi_max - gamma[j]) / (2.0 * gamma[j])));

  double z_m, a, b, c, y_r, A, B, C, y_m1, y_m2;
	
  while ( R[j][i] <= R_max ) {
	  
    //Now that the first ring has been filled we will generate the next one
    if ( i == 0 ) {

      long double tmp3 = gamma[j];
      long double tmp4 = Beta_L * (PI / 180.0);
		  
      R[j].push_back ( R[j][i] * cos(tmp3)
		       + sqrt(pow(DM*cos(tmp4), 2.0) + pow(R[j][i] * sin(tmp3), 2.0)) );	
      ++i;
    }
    else {

      long double tmp5 = Beta_L*(PI / 180.0);

      z_m = z_0 + R[j][i] * tan(tmp5);
			
      a = pow(z_1, 2.0)*(pow(r_m, 2.0) - pow(R[j][i], 2.0));
      b = 2 * R[j][i] * z_1*(z_1 - z_m);
      c = pow(r_m, 2.) - pow(z_m - z_1, 2.0);

      y_r = (-b + sqrt(b*b - 4.0 * a*c)) / (2.0 * a);
			
      A = -((2 * z_1*y_r + tan(tmp5))*tan(tmp5) + pow(z_1*y_r, 2.0));
      B = 2 * (z_1 - z_0)*(z_1*y_r + tan(tmp5));
			
      C = r_m*r_m*(1 + pow(z_1*y_r, 2.0)) - pow(z_1 - z_0, 2.0);

      // Check if the new radius will place heliostats far enough from the last defined circle
      R[j][i] + Delta_r_min >= (-B - sqrt(B*B - 4.0 * A*C)) / (2.0 * A) ?
	(y_m2 = R[j][i] + Delta_r_min) : (y_m2 = (-B - sqrt(B*B - 4.0 * A*C)) / (2.0 * A));
			
      N_new_group = 1 + 2 * myround(floor(Phi_max / (2.0 * (DM / (2.0 * y_m2)))));

      im1 = i - 1;
      z_m = z_0 + R[j][im1] * tan(tmp5);
      a = pow(z_1, 2.0)*(pow(r_m, 2.0) - pow(R[j][im1], 2.0));
      b = 2 * R[j][im1] * z_1*(z_1 - z_m);
      c = pow(r_m, 2.0) - pow(z_m - z_1, 2.0);

      y_r = (-b + sqrt(b*b - 4.0 * a*c)) / (2.0 * a);

      A = -((2 * z_1*y_r + tan(tmp5))*tan(tmp5) + pow(z_1*y_r, 2.0));
			
      B = 2 * (z_1 - z_0)*(z_1*y_r + tan(tmp5));
      C = r_m*r_m*(1 + pow(z_1*y_r, 2.0)) - pow(z_1 - z_0, 2.0);
			
      // Check if the new radius will place heliostats far enough from the last defined circle
      R[j][i] + Delta_r_min >= (-B - sqrt(B*B - 4.0 * A*C)) / (2.0 * A) ?
	(y_m1 = R[j][i] + Delta_r_min) : (y_m1 = (-B - sqrt(B*B - 4.0 * A*C)) / (2.0 * A));
			
      i + 1 % 2 == 0 ? N_new_circle = N_essential : N_new_circle = N_staggered;

      if ( N_new_group / (PI*(y_m2*y_m2 - R[j][i] * R[j][i])) >= N_new_circle / (PI*(y_m1*y_m1 - R[j][i] * R[j][i]))
	   && y_m2 + (DM / 2.0) <= R_max ) {

	//Create a new group
	++j;
	R.push_back(std::vector<double>(1, y_m2));
	gamma.push_back(DM / (2.0 * y_m2));
	i = 0;
	//Nombre d'heliostat pour les cercles du groupe j
	N_essential = 1 + 2 * myround(floor(Phi_max / (2.0 * gamma[j])));
	N_staggered = 2 + 2 * myround(floor((Phi_max - gamma[j]) / (2.0 * gamma[j])));
      }
      else {
	y_m2 = y_m1;
	if (y_m2 + (DM / 2.0) <= R_max) {
	  R[j].push_back(y_m2);
	  ++i;
	}
	else
	  break;
      }
    }
  }

  std::vector<double> positionPolaire     ( 7, 0.0 );
  std::vector<double> positionCartesienne ( 7, 0.0 );

  //[2] -> z
  //[3] -> average cosine efficiency
  //[4] -> atmospheric transmittivity
  //[5] -> average spillage

  for (j = 0; j < static_cast<int>(R.size()); ++j) {
    for (i = 0; i < static_cast<int>(R[j].size()); i++) {
      positionPolaire[0] = R[j][i];
      long double tmp6 = Beta_L*(PI / 180.0);
      positionPolaire[2] = z_0 + R[j][i] * tan(tmp6);
      			
      if (i % 2 == 0) {
	N_essential = 1 + 2 * static_cast<int>(floor(Phi_max / (2 * gamma[j])));
	for ( n = 0; n < N_essential; n += 2 )
	  {
	    if (n == 0) {
	      positionPolaire[1] = 0.;
	      _gridLayoutAngularCoordinates.push_back(positionPolaire);
	    }
	    else {
	      positionPolaire[1] = n*gamma[j];
	      _gridLayoutAngularCoordinates.push_back(positionPolaire);
	      positionPolaire[1] = -positionPolaire[1];
	      _gridLayoutAngularCoordinates.push_back(positionPolaire);
	    }
	  }
      }
      if (i % 2 == 1) {
	N_staggered = 2 + 2 * (int)floor((Phi_max - gamma[j]) / (2 * gamma[j]));
	for (n = 0; n < N_staggered; n += 2) {
	  np1 = n + 1;
	  positionPolaire[1] = np1*gamma[j];
	  _gridLayoutAngularCoordinates.push_back(positionPolaire);
	  positionPolaire[1] = -positionPolaire[1];
	  _gridLayoutAngularCoordinates.push_back(positionPolaire);
	}
      }
    }
  }

  for ( size_t i = 0; i < _gridLayoutAngularCoordinates.size(); ++i) {

    long double tmp7 = _gridLayoutAngularCoordinates[i][1];
	  
    positionCartesienne[0] =  _gridLayoutAngularCoordinates[i][0] * sin(tmp7); // x
    positionCartesienne[1] = -_gridLayoutAngularCoordinates[i][0] * cos(tmp7); // y
    positionCartesienne[2] =  _gridLayoutAngularCoordinates[i][2];             // z
    _gridLayoutCartesianCoordinates.push_back(positionCartesienne);
  }
}

/*-----------------------------------------------------------------*/
/*  Determines overall efficiency of a given position in the grid  */
/*-----------------------------------------------------------------*/
void HeliostatField::fComputeEfficiency ( void ) {

  //1- Compute maximum azimuth for sun position as a function of latitude
  //2- Compute integral for average cosine efficiency through the day for the whole grid
  //3- Compute total efficiency of given heliostat by factoring in atmospheric attenuation
  //4- dump data in _gridLayoutAngularCoordinates[n][3]
  fComputeAtmosphericAttenuation();
  fComputeCosineAndSpillage();
}

/*--------------------------------------------------------------------------------------------*/
/*  Determines atmospheric attenuation between a given position in the grid and the receiver  */
/*--------------------------------------------------------------------------------------------*/
/* Calculates the atmospheric attenuation losses for every position in the grid layout        */
/* inserts the value of transmitivity (1 - losses) in _gridLayoutAngularCoordinates[n][4]     */
/* Equations taken according to MODTRAN for a clear day/environment                           */
/*--------------------------------------------------------------------------------------------*/
void HeliostatField::fComputeAtmosphericAttenuation ( void ) {
  double slantRange;
  for ( size_t i = 0; i < _gridLayoutAngularCoordinates.size(); ++i ) {
    slantRange = sqrt(pow(_towerHeight, 2.0) + pow(_gridLayoutAngularCoordinates[i][0], 2.0)) / 1000.0; //km
    // Transmissivite; clear day
    _gridLayoutAngularCoordinates[i][4]
      = 1 - (0.29544 + 15.22128*slantRange - 1.8598 * pow(slantRange, 2.0) + 0.15182 * pow(slantRange, 3.0)) / 100.0;
  }
}

/*-----------------------------------------------------------------------------*/
/*  Determines losses due to the angle of the receiver to the heliostat cross  */
/*  section of receiver aperture seen from the heliostat may be smaller than   */
/*  the heliostat surface, causing some of the parallel light reflected from   */
/*  the heliostat to miss the aperture.                                        */
/*-----------------------------------------------------------------------------*/
/*  calculates the average losses due to spillage through the day for every    */
/*  position in the grid layout. Uses the Weiler-Atherton polygon clipping     */
/*  algorithm.                                                                 */
/*-----------------------------------------------------------------------------*/
void HeliostatField::fComputeCosineAndSpillage ( void ) {

  //angles en rad
  long double alpha_s, gamma_s, theta;
  long double summation;
  long double slantRange;
  std::vector<long double> vector_I ( 3, 0.0 );
  std::vector<long double> vector_R ( 3, 0.0 );
  int n = 0;

  Sun tempSun ( _sun );

  // Assuming that the sun azimuth is due south at 12:00
  // Using  ESO's FITS convention where azimuth is measured from the south increasing towards the west Thus it is 0 at 12:00
  // and at this point elevation is 90 degrees - latitude
  // Neglecting the earth's precession axis inclination (thus supposing the earth is perfectly straight on its orbit

  std::vector<long double> alpha;

  long double theta_argument;
  long double azimuthToAimpoint;
  long double elevationToAimpoint;
  long double heliostatElevation;
  long double heliostatAzimuth;
  std::vector<long double> vector_A = std::vector<long double>(3, 0.);
  long double vector_A_norm     = 0.0;
  long double summationSpillage = 0.0;
  long double nonSpilledRatio   = 0.0;

  for ( size_t i = 0; i < _gridLayoutCartesianCoordinates.size(); ++i ) {

    tempSun.fResetTime();
    n = 0;
    summation = 0.0;
    summationSpillage = 0.0;
    slantRange = sqrt(_gridLayoutCartesianCoordinates[i][0] * _gridLayoutCartesianCoordinates[i][0] +
		      _gridLayoutCartesianCoordinates[i][1] * _gridLayoutCartesianCoordinates[i][1] +
		      (_towerHeight - _gridLayoutCartesianCoordinates[i][2]) *
		      (_towerHeight - _gridLayoutCartesianCoordinates[i][2]));

    long double tmp1 = _gridLayoutCartesianCoordinates[i][0];
    long double tmp2 = _gridLayoutCartesianCoordinates[i][1];
		
    azimuthToAimpoint = atan2(tmp1, -tmp2);

    long double tmp3 = (_towerHeight - _gridLayoutCartesianCoordinates[i][2]) /
      sqrt(_gridLayoutCartesianCoordinates[i][0] * _gridLayoutCartesianCoordinates[i][0] +
	   _gridLayoutCartesianCoordinates[i][1] * _gridLayoutCartesianCoordinates[i][1]);
    
    elevationToAimpoint = atan( tmp3 );

    vector_R[0] = -_gridLayoutCartesianCoordinates[i][0] / slantRange;
    vector_R[1] = -_gridLayoutCartesianCoordinates[i][1] / slantRange;
    vector_R[2] = (_towerHeight - _gridLayoutCartesianCoordinates[i][2]) / slantRange;

    while ( tempSun.get_incrementsCounter() < tempSun.get_numberOfIncrements() ) {

      tempSun.fComputeSunPosition();

      alpha_s = tempSun.get_elevation()*DEG_TO_RAD;

      if ( alpha_s > 0.0 ) {

	gamma_s = tempSun.get_azimuth()*DEG_TO_RAD;
	
	long double tmp10 = alpha_s;
	long double tmp11 = gamma_s;
				
	vector_I[0] = -cos(tmp10)*sin(tmp11);
	vector_I[1] =  cos(tmp10)*cos(tmp11);
	vector_I[2] =  sin(tmp10);
	
	theta_argument = (vector_I[0] * vector_R[0]
			  + vector_I[1] * vector_R[1]
			  + vector_I[2] * vector_R[2]);
	
	if ( theta_argument >  1.0 )
	  theta_argument = 1.0;
	if ( theta_argument < -1.0 )
	  theta_argument = -1.0;

	long double tmp12 = theta_argument;
	
	theta = 0.5*acos(tmp12);
	long double tmp13 = theta;
				
	summation += cos(tmp13);

	// Computing spillage (approximation):
	vector_A[0] = 0.5*(vector_I[0] + vector_R[0]);
	vector_A[1] = 0.5*(vector_I[1] + vector_R[1]);
	vector_A[2] = 0.5*(vector_I[2] + vector_R[2]);
	vector_A_norm = sqrt(pow(vector_A[0], 2.0) + pow(vector_A[1], 2.0) + pow(vector_A[2], 2.0));

	long double tmp3 = vector_A[2] / vector_A_norm;
	heliostatElevation = asin(tmp3);
		
	long double tmp4 = -vector_A[0], tmp5 = vector_A[1];
	heliostatAzimuth = atan2 ( tmp4, tmp5 );

	tmp3 = azimuthToAimpoint - heliostatAzimuth;
	tmp4 = azimuthToAimpoint;
				
	if ( _heliostatWidth*cos(tmp3) > _apertureWidth*cos(tmp4) ) {
	  nonSpilledRatio = _apertureWidth*cos(tmp4) / _heliostatWidth*cos(tmp3);
	}
	else
	  nonSpilledRatio = 1.0;

	long double tmp17 = elevationToAimpoint - heliostatElevation;
	long double tmp18 = elevationToAimpoint;
				
	if ( _heliostatLength*cos(tmp17) > _apertureHeight*cos(tmp18) )
	  nonSpilledRatio *= _apertureHeight*cos(tmp18) / _heliostatLength*cos(tmp17);
	
	alpha.push_back(nonSpilledRatio);
	summationSpillage += nonSpilledRatio;
	
	++n;
      }
      tempSun.fTimeIncrement();
    }

    // Rounding for different machine precisions:
    long double tmp21 = roundl ( summationSpillage * 1e16 ) / 1e16;
    summationSpillage = tmp21;

    _gridLayoutCartesianCoordinates[i][3] = summation / n;
    _gridLayoutCartesianCoordinates[i][5] = summationSpillage / n;
  }
}

/*--------------------------------------------------------*/
/*  Uses the grid layout to determine the best potential  */
/*  heliostats to be placed on the field                  */
/*--------------------------------------------------------*/
bool comparePositions ( const std::vector<double> & firstPosition, const std::vector<double> & secondPosition ) {
  // Used to break equalities and have consistency accross platforms:
  if ( fabs ( firstPosition[6] - secondPosition[6] ) < 1e-13 )
    return firstPosition[0] > secondPosition[0];
  return firstPosition[6] > secondPosition[6]; 
}

/*------------------------------------------------------------------------*/
bool compareDistanceToSun ( const Heliostat * h1, const Heliostat* h2 ) {
/*------------------------------------------------------------------------*/
  return ( h1->get_yProj() > h2->get_yProj() );
}

/*-------------------------------------------------------------*/
/*  Determines the position of the _nbOfheliostats heliostats  */
/*-------------------------------------------------------------*/
void HeliostatField::fConfigureField ( void ) {

  size_t i = 0;
  std::vector<double> fieldPosition(6, 0.0);
	
  for ( i = 0; i < _gridLayoutAngularCoordinates.size(); ++i ) {
    _gridLayoutCartesianCoordinates[i][4] = _gridLayoutAngularCoordinates[i][4];
    _gridLayoutCartesianCoordinates[i][6] = _gridLayoutCartesianCoordinates[i][3] * //cosine efficiency
      _gridLayoutCartesianCoordinates[i][4] * //atmospheric transmissivity
      _gridLayoutCartesianCoordinates[i][5]; // non-spilled radiation
  }

  if ( _nbOfHeliostats >= (_gridLayoutCartesianCoordinates.size()) ) {
    _fieldLayoutCartesian = _gridLayoutCartesianCoordinates;
    _nbOfHeliostats = static_cast<int>(_fieldLayoutCartesian.size());
  }
  else {

    // DEBUG:
    // {
    //   std::cout.precision(20);
    //   size_t kk = _gridLayoutCartesianCoordinates.size();
    //   for ( size_t t = 0 ; t < kk ; ++t  ) {
    //  	std::cout << "BEFORE[" << t << "]={ ";
    //  	for ( size_t k = 0 ; k < _gridLayoutCartesianCoordinates[t].size() ; ++k )
    // 	  std::cout << _gridLayoutCartesianCoordinates[t][k] << ", ";
    // 	std::cout << "}\n";
    //   }
    // }
	     
    std::sort(_gridLayoutCartesianCoordinates.begin(), _gridLayoutCartesianCoordinates.end(), comparePositions);

    // DEBUG:
    // {
    //   size_t kk = _gridLayoutCartesianCoordinates.size();
    //   for (size_t t = 0 ; t < kk ; ++t  ) {
    // 	std::cout << "AFTER[" << t << "]={ ";
    // 	for ( size_t k = 0 ; k < _gridLayoutCartesianCoordinates[t].size() ; ++k )
    // 	  std::cout << _gridLayoutCartesianCoordinates[t][k] << ", ";
    // 	std::cout << "}\n";
    //   }
    // }
											   
    for ( size_t i = 0; i < _gridLayoutCartesianCoordinates.size(); ++i ) {
      fieldPosition[0] = _gridLayoutCartesianCoordinates[i][0];
      fieldPosition[1] = _gridLayoutCartesianCoordinates[i][1];
      fieldPosition[2] = _gridLayoutCartesianCoordinates[i][2];
      fieldPosition[3] = _gridLayoutCartesianCoordinates[i][3];
      fieldPosition[4] = _gridLayoutCartesianCoordinates[i][4];
      fieldPosition[5] = _gridLayoutCartesianCoordinates[i][5];	
      _fieldLayoutCartesian.push_back(fieldPosition);
    }
  }
}

/*----------------------------------------------------------*/
/*  Constructs the _nbOfHeliostats heliostats of the field  */
/*----------------------------------------------------------*/
void HeliostatField::fGenerateField ( void ) {
  
  delete_heliostats();
  
  size_t n = ( _nbOfHeliostats < _fieldLayoutCartesian.size() ) ?
    _nbOfHeliostats : _fieldLayoutCartesian.size();
  
  for ( size_t i = 0; i < n; ++i )
    _listOfHeliostats.push_back ( new Heliostat ( _fieldLayoutCartesian[i][0] ,
						  _fieldLayoutCartesian[i][1] ,
						  _fieldLayoutCartesian[i][2] ,
						  _fieldLayoutCartesian[i][3] ,
						  _fieldLayoutCartesian[i][4] ,
						  _towerHeight                  ) );
}

/*----------------------------------------------------------------------------*/
/* Evaluate the total visible surface of the field as projected onto a plane  */
/* perpendicular to solar radiation                                           */
// 0- top surface 1- west surface 2- east surface
/*----------------------------------------------------------------------------*/
double HeliostatField::fEvaluateFieldSurface ( void ) const {

  std::vector<double> vTopArea  (3, 0.0);
  std::vector<double> vWestSide (3, 0.0);
  std::vector<double> vEastSide (3, 0.0);
  std::vector<double> vWBackSide(3, 0.0);
  std::vector<double> vEBackSide(3, 0.0);
  std::vector<double> vSunrays  (3, 0.0);

  double tmp, maxRadius = 0.0;

  if ( _listOfHeliostats.size() < _gridLayoutAngularCoordinates.size() ) {
    for ( size_t i = 0; i < _listOfHeliostats.size(); ++i ) {
      tmp = sqrt(pow(_listOfHeliostats[i]->get_x(), 2.0) + pow(_listOfHeliostats[i]->get_y(), 2.0));
      if ( tmp > maxRadius )
	maxRadius = tmp;
    }
  }
  else
    maxRadius = _maxDistanceFromTower * _towerHeight;

  long double tmp2 = _maxAngularDeviation*PI / 180.0;
	
  double fTopArea  = pow(maxRadius, 2.0)*_maxAngularDeviation* DEG_TO_RAD;
  double fSideArea = _heliostatLength * maxRadius;
  double fBackArea = sqrt( pow(maxRadius * std::sin(tmp2), 2.0) +
			   pow(maxRadius * (1.0 - std::cos(tmp2)), 2.0)) *_heliostatLength;

  double fProjectedArea      = 0.0;
  double fScalarProductTop   = 0.0;
  double fScalarProductWest  = 0.0;
  double fScalarProductEast  = 0.0;
  double fScalarProductWBack = 0.0;
  double fScalarProductEBack = 0.0;
  
  vTopArea[0] = 0.0;
  vTopArea[1] = 0.0;
  vTopArea[2] = 1.0;

  tmp2 = _maxAngularDeviation*PI / 180.0;
	
  vWestSide[0] = -std::cos(tmp2);
  vWestSide[1] =  std::sin(tmp2);
  vWestSide[2] = 0.0;

  tmp2 /= 2.0;
  
  vWBackSide[0] = -std::cos(tmp2);
  vWBackSide[1] = -std::sin(tmp2);
  vWBackSide[2] = 0.0;
  vEastSide [0] = -vWestSide[0];
  vEastSide [1] =  vWestSide[1];
  vEastSide [2] =  vWestSide[2];
  vEBackSide[0] = -vWBackSide[0];
  vEBackSide[1] =  vWBackSide[1];
  vEBackSide[2] =  vWBackSide[2];

  long double tmp3 = _sun.get_elevation()*DEG_TO_RAD;
  long double tmp4 = _sun.get_azimuth()  *DEG_TO_RAD;
	
  vSunrays[0] = -std::cos(tmp3)*std::sin(tmp4);
  vSunrays[1] =  std::cos(tmp3)*std::cos(tmp4);
  vSunrays[2] =  std::sin(tmp3);
  
  for ( size_t i = 0; i < 3; ++i ) {
    vTopArea  [i] = vTopArea  [i] * fTopArea;
    vWestSide [i] = vWestSide [i] * fSideArea;
    vEastSide [i] = vEastSide [i] * fSideArea;
    vWBackSide[i] = vWBackSide[i] * fBackArea;
    vEBackSide[i] = vEBackSide[i] * fBackArea;
    fScalarProductTop   += vSunrays[i] * vTopArea[i];
    fScalarProductWest  += vSunrays[i] * vWestSide[i];
    fScalarProductEast  += vSunrays[i] * vEastSide[i];
    fScalarProductWBack += vSunrays[i] * vWBackSide[i];
    fScalarProductEBack += vSunrays[i] * vEBackSide[i];
  }
  
  // Total projected area is equal to the sum of the projection of the scalar product between radiations and
  // normal vectors, for positive scalar product. If the scalar product is negative it means that this face should
  // be hit from the inside of the field, and so it is not visible.
  if ( fScalarProductTop   > 0.0 ) { fProjectedArea += fScalarProductTop;   }
  if ( fScalarProductWest  > 0.0 ) { fProjectedArea += fScalarProductWest;  }
  if ( fScalarProductEast  > 0.0 ) { fProjectedArea += fScalarProductEast;  }
  if ( fScalarProductWBack > 0.0 ) { fProjectedArea += fScalarProductWBack; }
  if ( fScalarProductEBack > 0.0 ) { fProjectedArea += fScalarProductEBack; }

  return fProjectedArea;
}

/*------------------------------------------------*/
void HeliostatField::fGenerateSunrays ( void ) {
/*------------------------------------------------*/

  double rMax     = 0.0;
  double distance = 0.0;
  
  if ( _nbOfHeliostats < _gridLayoutAngularCoordinates.size() ) {

    for ( size_t i = 0; i < _listOfHeliostats.size(); ++i ) {
      distance = sqrt (   _listOfHeliostats[i]->get_x()*_listOfHeliostats[i]->get_x()
			+ _listOfHeliostats[i]->get_y()*_listOfHeliostats[i]->get_y()  );
      if ( distance > rMax )
	rMax = distance;
    }
  }
  else
    rMax = _maxDistanceFromTower * _towerHeight;

  double fSurfaceArea = (_maxAngularDeviation * PI / 90.0) * pow(rMax, 2.0);

  int nb_sunrays = static_cast<int>(ceil(fSurfaceArea*_sun.get_raysPerSquareMeters()));

  for ( int i = 0; i < nb_sunrays ; ++i )
    _sun.fAddNewSunray ( rMax, _maxAngularDeviation, _heliostatLength );
}

/*----------------------------------------------------------------------------------*/
/*  Compute overall field efficiency:                                               */
/*    total energy reflected to receiver / total energy flux accross field surface  */
/*----------------------------------------------------------------------------------*/
double HeliostatField::fComputeFieldEfficiency ( void ) {

  _sun.fComputeSunPosition();

  int nb_sunrays = static_cast<int>(_sun.get_nb_sunrays());
  
  if (_sun.get_elevation() >= 0.) {
	  
    for ( int i = 0; i < nb_sunrays ; ++i )
      _sun.projectTarget(i);
	  
    for ( size_t i = 0; i < _listOfHeliostats.size(); ++i) {
      _listOfHeliostats[i]->computeAngles(_sun, _towerHeight);
      _listOfHeliostats[i]->computePilarProjection(_sun);
      _listOfHeliostats[i]->computeCornersProjections(_sun);
    }

    // DEBUG
    // size_t kk = _listOfHeliostats.size();
    // if ( kk > 15 ) kk = 15;
    // for ( size_t i = 0; i < kk ; ++i) {
    //   std::cout << "BEFORE: HELIOSTAT #" << _listOfHeliostats[i]->get_ID()
    // 		<< ": x=" << _listOfHeliostats[i]->get_x()
    // 		<< "  y=" << _listOfHeliostats[i]->get_y()
    // 		<< " PROJ=" << _listOfHeliostats[i]->get_xProj() << " " << _listOfHeliostats[i]->get_yProj()
    // 		<< " " <<_listOfHeliostats[i]->get_zProj()
    // 		<< std::endl;
    // }
    
    // Sort heliostats in increasing order of projected Y coordinate
    std::sort(_listOfHeliostats.begin(), _listOfHeliostats.end(), compareDistanceToSun);

    // DEBUG
    // std::cout << std::endl;
    // for ( size_t i = 0; i < kk ; ++i) {
    //   std::cout << "AFTER: HELIOSTAT #" << _listOfHeliostats[i]->get_ID()
    // 		<< ": x=" << _listOfHeliostats[i]->get_x()
    // 		<< "  y=" << _listOfHeliostats[i]->get_y()
    // 		<< " PROJ=" << _listOfHeliostats[i]->get_xProj() << " " << _listOfHeliostats[i]->get_yProj()
    // 		<< " " <<_listOfHeliostats[i]->get_zProj()
    // 		<< std::endl;
    // }
    
    // compute collisions:
    int listSize      = static_cast<int>(_sun.get_nb_sunrays());
    int nb_heliostats = static_cast<int>(_listOfHeliostats.size());
    
    for ( int i = 0; i < nb_heliostats; ++i )
      for ( int j = 0; j < listSize; ++j )
	if ( _sun.computeCollision(j,*_listOfHeliostats[i]) )
	  _listOfHeliostats[i]->increase_sunraysCount();

    // compute resulting efficiency:
    double sumOfEfficiencies = 0.0;
   
    for ( size_t i = 0; i < _listOfHeliostats.size(); ++i ) {
    
      sumOfEfficiencies += _listOfHeliostats[i]->get_sunraysCount()
	*_listOfHeliostats[i]->get_atmosphericAttenuation()
	*_listOfHeliostats[i]->fComputeSpillage(_apertureHeight, _apertureWidth);

      _listOfHeliostats[i]->clear_sunraysCount();      
    }
   
    int nb_sunrays = static_cast<int>(_sun.get_nb_sunrays());
    for ( int i = 0; i < nb_sunrays; ++i )
      _sun.set_isIntercepted ( i, false );
	  
    _fieldEfficiency.push_back(sumOfEfficiencies / listSize*1.0);
   
    return sumOfEfficiencies / (listSize*1.0);
  }
  
  _fieldEfficiency.push_back ( 0.0 );
  return 0.0;
}

/*---------------------------------------------------------------------------*/
/*  Assigns a shadowing efficiency coefficlong doubleient to each heliostat  */
/*---------------------------------------------------------------------------*/
double HeliostatField::fCalculateTotalEnergyOutput ( void ) {
  double fieldEfficiency           = fComputeFieldEfficiency();
  double fieldPerpendicularSurface = fEvaluateFieldSurface();
  double energyPerSquareMeter      = EXTRATERRESTRIAL_INSOLATION*(1.0 - ATMOSPHERE_ATTENUATION / 100.0);
  double totalOutput               = energyPerSquareMeter*fieldPerpendicularSurface*fieldEfficiency;
  _powerOutput.push_back(totalOutput);
  return totalOutput;
}
