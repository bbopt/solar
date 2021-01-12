#include "HeliostatField.hpp"

HeliostatField::HeliostatField(unsigned int nbOfHeliostats,
	double heliostatLength, double heliostatWidth,
	double towerHeight, double apertureHeight, double apertureWidth,
	double minDistanceFromTower, double maxDistanceFromTower, double maxAngle,
	Sun& sun)
	:_nbOfHeliostats(nbOfHeliostats),
	_heliostatLength(heliostatLength), _heliostatWidth(heliostatWidth),
	_towerHeight(towerHeight), _apertureHeight(apertureHeight), _apertureWidth(apertureWidth),
	_minDistanceFromTower(minDistanceFromTower),
	_maxDistanceFromTower(maxDistanceFromTower),
	_maxAngularDeviation(maxAngle),
	_sun(&sun)
{
	fComputeStaggeredGridLayout();
	fComputeEfficiency();
	fConfigureField();

	Heliostat::set_width(heliostatWidth);
	Heliostat::set_length(heliostatLength);
	Sunray::set_minDistance();

	_powerOutput.reserve(24);
	_fieldEfficiency.reserve(24);
}

HeliostatField::~HeliostatField()
{
}

void HeliostatField::fComputeStaggeredGridLayout() {

  // initial check: R_min cannot be equal to zero:
  double R_min = _minDistanceFromTower * _towerHeight;
  if ( R_min == 0.0 )
    throw logic_error ( "Simulation could not go through: R_min equal to zero" );
  
	//Initial values --------------------------------------------------------
	double z_0 = _heliostatLength / 2; //Zo = heigth of heliostat center from the ground.
	double l_m = _heliostatLength;
	double w_m = _heliostatWidth;
	double H_t = _towerHeight;
	double l_r = l_m; //receiver height assumed to be the same height as the heliostats since they are flat
	//^^ incorrect assumption. See spillage.
	double Phi_max = _maxAngularDeviation * DEG_TO_RAD;
	double R_max = _maxDistanceFromTower * _towerHeight;
	// double f_a = 1.0; // ratio of net reflecting surface area
	double Beta_L = 0.; //terrain slope will be assumed to be zero -- Make this a parameter

	//Step 1
	double f = w_m / l_m; //ratio of heliostat width over length
	// double A_m = f*f_a*pow(l_m, 2.0); // net reflective area
	double z_1 = H_t - 0.5*l_r;
	double r_m = 0.5*l_m;

	double dS_min = 2.*f - sqrt(1. + pow(f, 2.0)); //eq 2
	double dS = dS_min;

	//Now that all variables have been identified, do the following sequence:
	//1- Determine the position of all heliostats on the first ring
	//2- Determine the minimum radial spacing
	//3- Determine the radius of the first staggered ring. If the difference between the 2nd radius and the first radius is smaller than the minimal radial spacing,
	//   define the new radius as R_0 + min_delta_R
	//4- Determine the position of all heliostats on the 2nd ring
	//5- Determine the next two radiuses using the no blocking method. Verify that the difference between each radius is always bigger than the minimal radial spacing
	//6- Verify if filling the 3rd radius using the current angular spacing will amount to a better heliostat/field area than skipping the 3rd radius and filling the 4th completely 
	//   as a new essential ring.


	//Step 2
	double DM;
	if (l_m * (sqrt(1. + pow(f, 2.0)) + dS) >= 2.*w_m) { DM = l_m * (sqrt(1. + pow(f, 2.0)) + dS); }
	else { DM = 2.*w_m; }
	double Delta_r_min = DM * cos(30.*M_PI / 180.) * cos(Beta_L* (M_PI / 180));

	//Step 3
	/*Note: on ne gèrera pas de "groupes" mais plutôt chaque cercle individuellement. À chaque cercle on vérifiera si le nouveau cercle devrait
	être créé en conservant le même espacement angulaire ou en l'éloignant et en redéfinissant l'espacement angulaire.*/
	unsigned int i, j, n;
	i = 0; j = 0; n = 0;

	vector < vector<double> > R;			 //Rayons des cercles.
	vector < double> gamma;				 //Espacements angulaires

	//The radius for the first ring is defined already by R_min
	R.push_back(vector<double>(1, R_min));			//Radius for the initial circle
	gamma.push_back(DM / (2.*R[0][0]));				//Et espacement angulaire gamma[0]fDM
	
	//Determine the number of heliostats on each ring of this group ------
	unsigned int N_essential, N_staggered, N_new_group, N_new_circle;
	// bool isNewGroup = true;

	//Nombre d'heliostat pour les cercles du groupe j
	N_essential = 1 + 2 * (int)floor(Phi_max / (2 * gamma[j]));
	N_staggered = 2 + 2 * (int)floor((Phi_max - gamma[j]) / (2 * gamma[j]));

	//
	double z_m, a, b, c, y_r, A, B, C, y_m1, y_m2;
	
	while (R[j][i] <= R_max)
	{
	  
		//Now that the first ring has been filled we will generate the next one
		if (i == 0)
		{
			R[j].push_back ( R[j][i] * cos(gamma[j])
					 + sqrt(pow(DM*cos(Beta_L * (M_PI / 180.)), 2.) + pow(R[j][i] * sin(gamma[j]), 2.)) );	
			++i;
		}
		else
		{
			//a partir du 3e cercle du groupe, on doit verifier s'il est avantageux
			//de poursuivre avec ce groupe ou de creer un nouveau groupe.
			//Puisque les miroirs sont identiques d'un groupe a l'autre, on peut simplement
			//compter la densite de miroirs totale.
			//Si on cree un nouveau groupe, on devra l'espacer a partir du dernier anneau.
			//Si on conserve le meme groupe, on devra l'espacer a partir du dernier anneau
			//de meme type (essential or staggered
			//Dans les deux cas, on doit mesurer le rayon a partir du dernier cercle pour calculer l'aire

			z_m = z_0 + R[j][i] * tan(Beta_L*(M_PI / 180.));
			a = pow(z_1, 2.)*(pow(r_m, 2.) - pow(R[j][i], 2.));
			b = 2 * R[j][i] * z_1*(z_1 - z_m);
			c = pow(r_m, 2.) - pow(z_m - z_1, 2.);

			y_r = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);

			A = -((2 * z_1*y_r + tan(Beta_L *(M_PI / 180.)))*tan(Beta_L*(M_PI / 180.)) + pow(z_1*y_r, 2.));
			B = 2 * (z_1 - z_0)*(z_1*y_r + tan(Beta_L*(M_PI / 180.)));
			C = r_m*r_m*(1 + pow(z_1*y_r, 2.)) - pow(z_1 - z_0, 2.);

			//Verify if the new radius will place heliostats far enough from the last defined circle
			R[j][i] + Delta_r_min >= (-B - sqrt(B*B - 4 * A*C)) / (2 * A) ?
				(y_m2 = R[j][i] + Delta_r_min) : (y_m2 = (-B - sqrt(B*B - 4 * A*C)) / (2 * A));
			
			N_new_group = 1 + 2 * (int)floor(Phi_max / (2 * (DM / (2 * y_m2))));

			//---------------------------------------------------------------------------------------------
			z_m = z_0 + R[j][i - 1] * tan(Beta_L*(M_PI / 180.));
			a = pow(z_1, 2.)*(pow(r_m, 2.) - pow(R[j][i - 1], 2.));
			b = 2 * R[j][i - 1] * z_1*(z_1 - z_m);
			c = pow(r_m, 2.) - pow(z_m - z_1, 2.);

			y_r = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);

			A = -((2 * z_1*y_r + tan(Beta_L *(M_PI / 180.)))*tan(Beta_L*(M_PI / 180.)) + pow(z_1*y_r, 2.));
			B = 2 * (z_1 - z_0)*(z_1*y_r + tan(Beta_L*(M_PI / 180.)));
			C = r_m*r_m*(1 + pow(z_1*y_r, 2.)) - pow(z_1 - z_0, 2.);
			
			//Verify if the new radius will place heliostats far enough from the last defined circle
			R[j][i] + Delta_r_min >= (-B - sqrt(B*B - 4 * A*C)) / (2 * A) ?
				(y_m1 = R[j][i] + Delta_r_min) : (y_m1 = (-B - sqrt(B*B - 4 * A*C)) / (2 * A));
			
			i + 1 % 2 == 0 ? N_new_circle = N_essential : N_new_circle = N_staggered;

			if (N_new_group / (M_PI*(y_m2*y_m2 - R[j][i] * R[j][i])) >= N_new_circle / (M_PI*(y_m1*y_m1 - R[j][i] * R[j][i]))
				&& y_m2 + (DM / 2) <= R_max)
			{
				//Create a new group
	  
				++j;
				R.push_back(vector<double>(1, y_m2));
				gamma.push_back(DM / (2 * y_m2));
				i = 0;
				//Nombre d'heliostat pour les cercles du groupe j
				N_essential = 1 + 2 * (int)floor(Phi_max / (2 * gamma[j]));
				N_staggered = 2 + 2 * (int)floor((Phi_max - gamma[j]) / (2 * gamma[j]));
			}
			else
			{
				y_m2 = y_m1;
				if (y_m2 + (DM / 2.) <= R_max)
				{
					R[j].push_back(y_m2);
					++i;
				}
				else{ break; }
			}
		}
	}

	
	//Maintenant que toute les valeurs de rayons et de gamma sont determinees
	//On trouve toutes les positions de la grille.
	vector<double> positionPolaire(7, 0.);
	vector<double> positionCartesienne(7, 0.);
	//[2] -> z
	//[3] -> average cosine efficiency
	//[4] -> atmospheric transmittivity
	//[5] -> average spillage


	for (j = 0; j < R.size(); ++j)
	{
		for (i = 0; i < R[j].size(); i++)
		{
			positionPolaire[0] = R[j][i];		
			positionPolaire[2] = z_0 + R[j][i] * tan(Beta_L*(M_PI / 180.));
			if (i % 2 == 0)
			{
				N_essential = 1 + 2 * (int)floor(Phi_max / (2 * gamma[j]));
				for (n = 0; n < N_essential; n += 2)
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
			if (i % 2 == 1)
			{
				N_staggered = 2 + 2 * (int)floor((Phi_max - gamma[j]) / (2 * gamma[j]));
				for (n = 0; n < N_staggered; n += 2)
				{
					positionPolaire[1] = (n + 1)*gamma[j];
					_gridLayoutAngularCoordinates.push_back(positionPolaire);
					positionPolaire[1] = -positionPolaire[1];
					_gridLayoutAngularCoordinates.push_back(positionPolaire);
				}
			}
		}
	}

	// Puis la grille cartesienne pour tracage

	for (unsigned int i = 0; i < _gridLayoutAngularCoordinates.size(); ++i)
	{
		positionCartesienne[0] = _gridLayoutAngularCoordinates[i][0] * sin(_gridLayoutAngularCoordinates[i][1]); // x
		positionCartesienne[1] = -_gridLayoutAngularCoordinates[i][0] * cos(_gridLayoutAngularCoordinates[i][1]); // y
		positionCartesienne[2] = _gridLayoutAngularCoordinates[i][2];                         // z
		_gridLayoutCartesianCoordinates.push_back(positionCartesienne);
	}

	//for_each(_gridLayoutAngularCoordinates.begin(), _gridLayoutAngularCoordinates.end(),
	//	[&](const vector<double> coordPolaires){
	//	positionCartesienne[0] = coordPolaires[0] * sin(coordPolaires[1]); // x
	//	positionCartesienne[1] = -coordPolaires[0] * cos(coordPolaires[1]); // y
	//	positionCartesienne[2] = coordPolaires[2];                         // z
	//	_gridLayoutCartesianCoordinates.push_back(positionCartesienne);
	//});
}

void HeliostatField::fComputeEfficiency()
{
	//1- Compute maximum azimuth for sun position as a function of latitude
	//2- Compute integral for average cosine efficiency through the day for the whole grid
	//3- Compute total efficiency of given heliostat by factoring in atmospheric attenuation
	//4- dump data in _gridLayoutAngularCoordinates[n][3]
	fComputeAtmosphericAttenuation();
	fComputeCosineAndSpillage();
}

//-----------------------------------------------------------------------
//Calculates the atmospheric attenuation losses for every position in the
// grid layout.
// inserts the value of transmitivity (1 - losses) in 
// _gridLayoutAngularCoordinates[n][4]
//-----------------------------------------------------------------------
void HeliostatField::fComputeAtmosphericAttenuation()
{
	double slantRange;
	//Equations taken according to MODTRAN for a clear day/environment

	for (size_t i = 0; i < _gridLayoutAngularCoordinates.size(); ++i)
	{
		slantRange = sqrt(pow(_towerHeight, 2.) + pow(_gridLayoutAngularCoordinates[i][0], 2.)) / 1000; //km
		//Transmissivite --- clear day
		_gridLayoutAngularCoordinates[i][4] = 1 - (0.29544
			+ 15.22128*slantRange
			- 1.8598 * pow(slantRange, 2.)
			+ 0.15182 * pow(slantRange, 3.)) / 100;
	}

	//for_each(_gridLayoutAngularCoordinates.begin(), _gridLayoutAngularCoordinates.end(),
	//	[&](vector<double>& gridPosition)
	//{
	//	slantRange = sqrt(pow(_towerHeight, 2.) + pow(gridPosition[0], 2.)) / 1000; //km
	//	//Transmissivite --- clear day
	//	gridPosition[4] = 1 - (0.29544
	//		+ 15.22128*slantRange
	//		- 1.8598 * pow(slantRange, 2.)
	//		+ 0.15182 * pow(slantRange, 3.)) / 100;

	//	////Transmissivite --- hazy day
	//	//slantRange = 1 - (0.77941
	//	//	+ 55.49083 * slantRange
	//	//	- 14.78875 * pow(slantRange, 2.)
	//	//	+ 1.53718  * pow(slantRange, 3.));
	//});
}

//-----------------------------------------------------------------------
//calculates the average losses due to spillage through the day for every
//position in the grid layout.
//Uses the Weiler-Atherton polygon clipping algorithm
//-----------------------------------------------------------------------
void HeliostatField::fComputeCosineAndSpillage()
{
	//angles en rad
	double alpha_s, gamma_s, theta;
	double latitude, declin;
	double timeStep;
	double summation;
	double slantRange;
	vector<double> vector_I(3, 0.);
	vector<double> vector_R(3, 0.);
	int time = 0; int n = 0;

	latitude = _sun->get_latitude()*DEG_TO_RAD;
	declin = _sun->get_sunDeclination();
	timeStep = _sun->get_time().get_sizeOfIncrements();
	Sun tempSun(*_sun);
	//Assuming that the sun azimuth is due south at 12:00
	//Using  ESO's FITS convention where azimuth is measured from the south increasing towards the west Thus it is 0 at 12:00
	//and at this point elevation is 90 degrees - latitude
	//Neglecting the earth's precession axis inclination (thus supposing the earth is perfectly straight on its orbit
	//On calcule "cosine efficiency" sur toute la journee

	//***
	int K = 0;
	vector<double> alpha;
	ofstream alphaDebug;
	//***
	double theta_argument;
	double azimuthToAimpoint;
	double elevationToAimpoint;
	double heliostatElevation;
	double heliostatAzimuth;
	double H_apparent;
	double W_apparent;
	vector<double> vector_A = vector<double>(3, 0.);
	double vector_A_norm = 0.;
	double summationSpillage;
	double nonSpilledRatio = 0.;

	for (unsigned int i = 0; i < _gridLayoutCartesianCoordinates.size(); ++i)
	{
		tempSun.get_time().fResetTime();
		n = 0;
		time = tempSun.get_time().get_currentTime();
		summation = 0.;
		summationSpillage = 0.;
		slantRange = sqrt(_gridLayoutCartesianCoordinates[i][0] * _gridLayoutCartesianCoordinates[i][0] +
			_gridLayoutCartesianCoordinates[i][1] * _gridLayoutCartesianCoordinates[i][1] +
			(_towerHeight - _gridLayoutCartesianCoordinates[i][2])* (_towerHeight - _gridLayoutCartesianCoordinates[i][2]));

		azimuthToAimpoint = atan2(_gridLayoutCartesianCoordinates[i][0], -_gridLayoutCartesianCoordinates[i][1]);
		elevationToAimpoint = atan((_towerHeight - _gridLayoutCartesianCoordinates[i][2]) /
			sqrt(_gridLayoutCartesianCoordinates[i][0] * _gridLayoutCartesianCoordinates[i][0] + _gridLayoutCartesianCoordinates[i][1] * _gridLayoutCartesianCoordinates[i][1]));

		H_apparent = _apertureHeight*cos(elevationToAimpoint);
		W_apparent = _apertureWidth*cos(azimuthToAimpoint);

		//Verifier car ici les x et y sont calcules en considerant y positif
		// vers le nord mais l'article prend y positif vers le sud
		vector_R[0] = -_gridLayoutCartesianCoordinates[i][0] / slantRange;
		vector_R[1] = -_gridLayoutCartesianCoordinates[i][1] / slantRange; //article met -heliostat[1]
		vector_R[2] = (_towerHeight - _gridLayoutCartesianCoordinates[i][2]) / slantRange;

		while (tempSun.get_time().get_incrementsCounter() < tempSun.get_time().get_numberOfIncrements())
		{
			tempSun.fComputeSunPosition();

			/*HRA = (1.*(time-720)*EARTH_OMEGA) * DEG_TO_RAD;*/
			alpha_s = tempSun.get_elevation()*DEG_TO_RAD;

			if (alpha_s > 0.) // on ne calcule la moyenne que pendant les periodes d'ensoleillement
			{
				gamma_s = tempSun.get_azimuth()*DEG_TO_RAD;

				vector_I[0] = -cos(alpha_s)*sin(gamma_s);
				vector_I[1] = cos(alpha_s)*cos(gamma_s);
				vector_I[2] = sin(alpha_s);

				theta_argument = (vector_I[0] * vector_R[0]
					+ vector_I[1] * vector_R[1]
					+ vector_I[2] * vector_R[2]);

				if (theta_argument > 1.0) { theta_argument = 1.0; }
				if (theta_argument < -1.0){ theta_argument = -1.0; }

				theta = 0.5*acos(theta_argument);

				summation += cos(theta);

				//Computing spillage (approximation) -----------------------------------
				vector_A[0] = 0.5*(vector_I[0] + vector_R[0]);
				vector_A[1] = 0.5*(vector_I[1] + vector_R[1]);
				vector_A[2] = 0.5*(vector_I[2] + vector_R[2]);
				vector_A_norm = sqrt(pow(vector_A[0], 2.) + pow(vector_A[1], 2.) + pow(vector_A[2], 2.));
				heliostatElevation = asin(vector_A[2] / vector_A_norm);
				heliostatAzimuth = atan2(-vector_A[0], vector_A[1]);

				if (_heliostatWidth*cos(azimuthToAimpoint - heliostatAzimuth) > _apertureWidth*cos(azimuthToAimpoint))
				{
					nonSpilledRatio = _apertureWidth*cos(azimuthToAimpoint) /
						_heliostatWidth*cos(azimuthToAimpoint - heliostatAzimuth);
				}
				else{ nonSpilledRatio = 1.0; }

				if (_heliostatLength*cos(elevationToAimpoint - heliostatElevation) > _apertureHeight*cos(elevationToAimpoint))
				{
					nonSpilledRatio *= _apertureHeight*cos(elevationToAimpoint) /
						_heliostatLength*cos(elevationToAimpoint - heliostatElevation);
				}
				else{ nonSpilledRatio *= 1.0; }
				alpha.push_back(nonSpilledRatio);
				summationSpillage += nonSpilledRatio;
				//-----------------------------------------------------------------------------------

				++n;
			}
			tempSun.get_time().fTimeIncrement();
			time = tempSun.get_time().get_currentTime();
		}
		//***
		K = 1;
		//***
		_gridLayoutCartesianCoordinates[i][3] = summation / n;
		_gridLayoutCartesianCoordinates[i][5] = summationSpillage / n;
	}
}

//Uses the grid layout to determine the best potential heliostats to be placed on the field
//Predicate
bool comparePositions(vector<double> firstPosition, vector<double> secondPosition)
{
	return firstPosition[6] > secondPosition[6];
}

void HeliostatField::fConfigureField()
{
	unsigned int i = 0;
	vector<double> fieldPosition(6, 0.);

	for (i = 0; i < _gridLayoutAngularCoordinates.size(); ++i)
	{
		_gridLayoutCartesianCoordinates[i][4] = _gridLayoutAngularCoordinates[i][4];
		_gridLayoutCartesianCoordinates[i][6] = _gridLayoutCartesianCoordinates[i][3] * //cosine efficiency
			_gridLayoutCartesianCoordinates[i][4] * //atmospheric transmissivity
			_gridLayoutCartesianCoordinates[i][5]; // non-spilled radiation		
	}

	if ( _nbOfHeliostats >= (_gridLayoutCartesianCoordinates.size()) )	{
		_fieldLayoutCartesian = _gridLayoutCartesianCoordinates;
		_nbOfHeliostats = _fieldLayoutCartesian.size();


		
	}
	else{
		std::sort(_gridLayoutCartesianCoordinates.begin(), _gridLayoutCartesianCoordinates.end(),
			comparePositions);

		for (unsigned int i = 0; i < _gridLayoutCartesianCoordinates.size(); ++i)
		{
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

void HeliostatField::fGenerateField()
{
	int n;
	_nbOfHeliostats < _fieldLayoutCartesian.size() ?
		n = _nbOfHeliostats : n = _fieldLayoutCartesian.size();

	for (int i = 0; i < n; ++i)
	{
		_listOfHeliostats.push_back(new Heliostat(_fieldLayoutCartesian[i][0], _fieldLayoutCartesian[i][1], _fieldLayoutCartesian[i][2], _fieldLayoutCartesian[i][3], _fieldLayoutCartesian[i][4], _towerHeight));
	}
}

//Function to evaluate the total visible surface of the field as projected onto a plane
//perpendicular to solar radiation.
//NOTE ------ 
//How should we consider the area. Should the whole available area (grid layout) be used or
//should only the field layout area be used. 
double HeliostatField::fEvaluateFieldSurface() const
{
	//Vecteur normaux aux trois surfaces, pointant vers l'exterieur du champs
	//Le vecteur de al surface arriere est neglige parce que le soleil ne devrait jamais se trouver derriere le champs
	vector<double> vTopArea(3, 0.);
	vector<double> vWestSide(3, 0.);
	vector<double> vEastSide(3, 0.);
	vector<double> vWBackSide(3, 0.);
	vector<double> vEBackSide(3, 0.);
	vector<double> vSunrays(3, 0.);
	double maxRadius = 0.;

	if (_listOfHeliostats.size() < _gridLayoutAngularCoordinates.size())
	{
		for (unsigned int i = 0; i < _listOfHeliostats.size(); ++i)
		{
			if (sqrt(pow(_listOfHeliostats[i]->get_x(), 2.) + pow(_listOfHeliostats[i]->get_y(), 2.)) > maxRadius)
			{
				maxRadius = sqrt(pow(_listOfHeliostats[i]->get_x(), 2.) + pow(_listOfHeliostats[i]->get_y(), 2.));
			}
		}
	}
	else { maxRadius = _maxDistanceFromTower * _towerHeight; }
	double fTopArea = pow(maxRadius, 2.)*_maxAngularDeviation* DEG_TO_RAD;
	double fSideArea = _heliostatLength * maxRadius;
	double fBackArea = sqrt(
		pow(maxRadius * sin(_maxAngularDeviation*M_PI / 180.), 2.) +
		pow(maxRadius * (1. - cos(_maxAngularDeviation*M_PI / 180.)), 2.))
		*_heliostatLength;

	double fProjectedArea = 0.;
	double fScalarProductTop = 0.;
	double fScalarProductWest = 0.;
	double fScalarProductEast = 0.;
	double fScalarProductWBack = 0.;
	double fScalarProductEBack = 0.;



	vTopArea[0] = 0.;
	vTopArea[1] = 0.;
	vTopArea[2] = 1.;
	vWestSide[0] = -cos(_maxAngularDeviation*M_PI / 180.);
	vWestSide[1] = sin(_maxAngularDeviation*M_PI / 180.);
	vWestSide[2] = 0.;
	vWBackSide[0] = -cos(_maxAngularDeviation*0.5*M_PI / 180.);
	vWBackSide[1] = -sin(_maxAngularDeviation*0.5*M_PI / 180.);
	vWBackSide[2] = 0.;
	vEastSide[0] = -vWestSide[0];
	vEastSide[1] = vWestSide[1];
	vEastSide[2] = vWestSide[2];
	vEBackSide[0] = -vWBackSide[0];
	vEBackSide[1] = vWBackSide[1];
	vEBackSide[2] = vWBackSide[2];

	//Same definition as Vector_I in fComputeCosineEfficiency
	vSunrays[0] = -cos(_sun->get_elevation()*DEG_TO_RAD)*sin(_sun->get_azimuth()*DEG_TO_RAD);
	vSunrays[1] = cos(_sun->get_elevation()*DEG_TO_RAD)*cos(_sun->get_azimuth()*DEG_TO_RAD);
	vSunrays[2] = sin(_sun->get_elevation()*DEG_TO_RAD);

	for (unsigned int i = 0; i < 3; ++i)
	{
		vTopArea[i] = vTopArea[i] * fTopArea;
		vWestSide[i] = vWestSide[i] * fSideArea;
		vEastSide[i] = vEastSide[i] * fSideArea;
		vWBackSide[i] = vWBackSide[i] * fBackArea;
		vEBackSide[i] = vEBackSide[i] * fBackArea;
		fScalarProductTop += vSunrays[i] * vTopArea[i];
		fScalarProductWest += vSunrays[i] * vWestSide[i];
		fScalarProductEast += vSunrays[i] * vEastSide[i];
		fScalarProductWBack += vSunrays[i] * vWBackSide[i];
		fScalarProductEBack += vSunrays[i] * vEBackSide[i];
	}

	//Total projected area is equal to the sum of the projection of the scalar product between radiations and
	//normal vectors, for positive scalar product. If the scalar product is negative it means that this face should
	//be hit from the inside of the field, and so it is not visible.

	if (fScalarProductTop > 0) { fProjectedArea += fScalarProductTop; }
	if (fScalarProductWest > 0) { fProjectedArea += fScalarProductWest; }
	if (fScalarProductEast > 0) { fProjectedArea += fScalarProductEast; }
	if (fScalarProductWBack > 0) { fProjectedArea += fScalarProductWBack; }
	if (fScalarProductEBack > 0) { fProjectedArea += fScalarProductEBack; }

	return fProjectedArea;
}

void HeliostatField::fGenerateSunrays()
{
	double rMax = 0.;
	double distance = 0.;

	//Si le champs n'est pas rempli par defaut alors on reduira la repartition des rayons
	//a la zone occupee par les heliostats pour une meilleure precision
	if (_nbOfHeliostats < (_gridLayoutAngularCoordinates.size()))
	{
		for (unsigned int i = 0; i < _listOfHeliostats.size(); ++i)
		{
			distance = sqrt(_listOfHeliostats[i]->get_x()*_listOfHeliostats[i]->get_x() + _listOfHeliostats[i]->get_y()*_listOfHeliostats[i]->get_y());
			if (distance > rMax) { rMax = distance; }
		}
	}
	else{ rMax = _maxDistanceFromTower * _towerHeight; }

	double fSurfaceArea = (_maxAngularDeviation*2. * M_PI / 180.)
		*(pow(rMax, 2.)/* - pow(_towerHeight*_minDistanceFromTower, 2.)*/);
	_sun->get_listOfSunrays().reserve((int)ceil(fSurfaceArea*_sun->get_raysPerSquareMeters()));

	for (unsigned int i = 0; i < ceil(fSurfaceArea*_sun->get_raysPerSquareMeters()); ++i)
	{
		_sun->fAddNewSunray(rMax, _maxAngularDeviation, _heliostatLength);
	}
}

//------------------------------------------------------------------------------------------------
//Determine si un heliostat est plus pres du soleil (plan perpendiculaires aux radiations)
//en comparant la valeur de leur coordonnee _yProj, qui est par definition l'axe pointant
//directement vers le soleil. On fait une compariason > plutot que < parce que par definition
//l'axe est positif vers le soleil, les valeurs les plus faibles de Y sont donc les plus eloignees
//-------------------------------------------------------------------------------------------------
bool compareDistanceToSun(Heliostat* heliostat1, Heliostat* heliostat2)
{
	return heliostat1->get_yProj() > heliostat2->get_yProj();
}

double HeliostatField::fComputeFieldEfficiency() {

  _sun->fComputeSunPosition();

  if (_sun->get_elevation() >= 0.) {
	  
    for (unsigned int i = 0; i < _sun->get_listOfSunrays().size(); ++i)
      _sun->get_sunray(i)->projectTarget();
	  
    for (unsigned int i = 0; i < _listOfHeliostats.size(); ++i) {
	_listOfHeliostats[i]->computeAngles(*_sun, _towerHeight);
	_listOfHeliostats[i]->computePilarProjection(*_sun);
	_listOfHeliostats[i]->computeCornersProjections(*_sun);
      }

    //Sort heliostats in increasing order of projected Y coordinate
    std::sort(_listOfHeliostats.begin(), _listOfHeliostats.end(), compareDistanceToSun);

    //3 - compute collisions;
    unsigned int i = 0;
    unsigned int j = 0;
    unsigned int listSize = _sun->get_listOfSunrays().size();
    unsigned int remainingSize = listSize;
	  
    for (i = 0; i < _listOfHeliostats.size(); ++i) {
      remainingSize = _sun->get_listOfSunrays().size();
      for (j = 0; j < remainingSize; ++j)
	_sun->get_sunray(j)->computeCollision(_listOfHeliostats[i]);
    }

    //4 - compute resulting efficiency;
    double sumOfEfficiencies = 0.;
    for (unsigned int i = 0; i < _listOfHeliostats.size(); ++i) {
      sumOfEfficiencies += _listOfHeliostats[i]->get_sunraysCount()
	*_listOfHeliostats[i]->get_atmosphericAttenuation()
	*_listOfHeliostats[i]->fComputeSpillage(_apertureHeight, _apertureWidth);    
      _listOfHeliostats[i]->clear_sunraysCount();
    }

    for (unsigned int i = 0; i < _sun->get_listOfSunrays().size(); ++i)
      _sun->get_sunray(i)->set_isIntercepted(false);
	  
    _fieldEfficiency.push_back(sumOfEfficiencies / listSize*1.);
   
    return sumOfEfficiencies / (listSize*1.);
  }
  else {
    _fieldEfficiency.push_back(0.);
    return 0.;
  }
}

double HeliostatField::fCalculateTotalEnergyOutput()
{
	double fieldEfficiency = fComputeFieldEfficiency();
	double fieldPerpendicularSurface = fEvaluateFieldSurface();
	double energyPerSquareMeter = EXTRATERRESTRIAL_INSOLATION*(1. - ATMOSPHERE_ATTENUATION / 100.);
	double totalOutput = energyPerSquareMeter*fieldPerpendicularSurface*fieldEfficiency;
	_powerOutput.push_back(totalOutput);
	return totalOutput;
}

//Notes:
//Avant de determiner le "spillage", compléter une fonction d'affichage qui insert les deux
//facteurs d'efficacité déjà calculés dans le vecteur cartésien et qui imprime le vecteur
//dans un fichier. Ensuite, utiliser MatLab pour tracer trois graphiques:
//1- la grille en vue superposée
//2-3- la grille avec une valeur de fonction égale à l'efficacité pour chaque point. Le graphique
//     devrait montrer la répartition de l'efficacité sur la grille.

struct foncteurPrint
{
	ofstream* _stream;
	foncteurPrint(ofstream* stream)
	{
		_stream = stream;
	}

	void operator() (vector<double> position)
	{
		for (size_t i = 0; i < position.size(); ++i)
		{
			*_stream << position[i] << " ";
		}
		*_stream << "\n";
	}
};

void HeliostatField::fOutputGridLayout(string outputFile) const
{
	ofstream gridLayout;
	foncteurPrint printPosition(&gridLayout);
	gridLayout.open(outputFile.c_str());
	gridLayout << "X     " << "Y     " << "Z     " << "efficiency" << endl;
	for_each(_gridLayoutCartesianCoordinates.begin(), _gridLayoutCartesianCoordinates.end(),
		printPosition);

	gridLayout.close();
}

struct foncteurPrintXYZ
{
	ofstream* _stream;
	foncteurPrintXYZ(ofstream* stream)
	{
		_stream = stream;
	}

	void operator()(vector<double> vecteur)
	{
		*_stream << vecteur[0] << " " << vecteur[1] << " " << vecteur[2] << endl;
	}
};

void HeliostatField::fOutputFieldLayout(string outputFile) const
{
	ofstream fieldLayout;
	foncteurPrintXYZ printPosition(&fieldLayout);
	fieldLayout.open(outputFile.c_str());
	fieldLayout << "X     " << "Y     " << "Z     " << endl;
	for_each(_fieldLayoutCartesian.begin(), _fieldLayoutCartesian.end(),
		printPosition);

	fieldLayout.close();
}

struct foncteurPrintAtmAtn
{
	ofstream* _stream;

	foncteurPrintAtmAtn(ofstream* stream)
	{
		_stream = stream;
	}

	void operator() (Heliostat* h)
	{
		*_stream << h->get_atmosphericAttenuation() << " " << h->get_sunraysCount() << "\n";
	}
};

void HeliostatField::fOutputAtmosphericAttenuation(string outputFile) const
{
	ofstream atmAttenuation;
	foncteurPrintAtmAtn printAttenuation(&atmAttenuation);
	atmAttenuation.open(outputFile.c_str());
	atmAttenuation << "Atmospheric attenuation " << "sunray count" << endl;
	for_each(_listOfHeliostats.begin(), _listOfHeliostats.end(),
		printAttenuation);

	atmAttenuation.close();
}



//void HeliostatField::fPrintProjections(string outputFile) const
//{
//	ofstream projections;
//	projections.open(outputFile);
//	projections << "_xProj _yProj _zProj" << endl;
//	for_each(_listOfHeliostats.begin(), _listOfHeliostats.end(),
//		[&](Heliostat* heliostat)
//	{
//		projections << heliostat->get_xProj() << " " << heliostat->get_yProj() << " " << heliostat->get_zProj() << "\n";
//	});
//
//	projections << "projections de sunray target" << endl;
//	vector<Sunray*> list = _sun->get_listOfSunrays();
//	for_each(list.begin(), list.end(),
//		[&](Sunray* sunray)
//	{
//		projections << sunray->get_projectedTarget()[0] << " " << sunray->get_projectedTarget()[1] << " " << sunray->get_projectedTarget()[2] << endl;
//	});
//
//	projections.close();
//}

//void HeliostatField::fPrintSunraysTargetProjections(string outputFile) const
//{
//	ofstream projections;
//	projections.open(outputFile);
//	projections << "_xProj _yProj _zProj" << endl;
//
//	vector<Sunray*> list = _sun->get_listOfSunrays();
//	for_each(list.begin(), list.end(),
//		[&](Sunray* sunray)
//	{
//		projections << sunray->get_projectedTarget()[0] << " " << sunray->get_projectedTarget()[1] << " " << sunray->get_projectedTarget()[2] << endl;
//	});
//	projections.close();
//}

//void HeliostatField::fOutputRayInterception(string outputFile) const
//{
//	ofstream intercept;
//	intercept.open(outputFile);
//	intercept << "_x _y _nbOfSunrays" << endl;
//
//	for_each(_listOfHeliostats.begin(), _listOfHeliostats.end(),
//		[&](Heliostat* heliostat)
//	{
//		intercept << heliostat->get_x() << " " << heliostat->get_y() << " " << heliostat->get_sunraysCount() << endl;
//	});
//	intercept.close();
//}
