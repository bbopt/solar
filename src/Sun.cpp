#include "Sun.hpp"

Sun::Sun(double latitude, Clock time, int day, double raysPerSquareMeters)
  : _latitude(latitude),
    _raysPerSquareMeters(raysPerSquareMeters),
    _time(time),
    _day(day)
{
	_sunDeclination = DECLINATION * DEG_TO_RAD * sin((360 * (_day + 284) / 365)*DEG_TO_RAD);
	fComputeSunPosition();
}

Sun::Sun(const Sun& sun)
  : _latitude(sun._latitude),
    _raysPerSquareMeters(sun._raysPerSquareMeters),
    _time(sun._time),
    _day(sun._day)
{
	_sunDeclination = DECLINATION * DEG_TO_RAD * sin((360 * (_day + 284) / 365)*DEG_TO_RAD);
	fComputeSunPosition();
}

Sun::~Sun()
{
  clearListOfSunrays();
  // for_each(_listOfSunrays.begin(), _listOfSunrays.end(),
  // 	[](Sunray* sunray)
  // {
  // 	delete sunray;
  // });
}

void Sun::fComputeSunPosition()
{
	//Assuming that the sun azimuth is due south at 12:00
	//Using  ESO's FITS convention where azimuth is measured from the south increasing towards the west Thus it is 0 at 12:00
	//and at this point elevation is 90 degrees - latitude
	//Neglecting the earth's precession axis inclination (thus supposing the earth is perfectly straight on its orbit

		//radians - put source. Calculating at summer solstice
	double HRA = (1.0*(_time.get_currentTime() - 720)*EARTH_OMEGA) * DEG_TO_RAD;			// Local Solar Time
	double latitude_in_rad = _latitude * DEG_TO_RAD;
	
	_elevation = asin(sin(latitude_in_rad)*sin(_sunDeclination) + cos(latitude_in_rad)*cos(_sunDeclination)*cos(HRA)) * RAD_TO_DEG;
	double azimuth_argument;
	azimuth_argument = (sin(_elevation * DEG_TO_RAD)*sin(latitude_in_rad) - sin(_sunDeclination )) / (cos(_elevation * DEG_TO_RAD)*cos(latitude_in_rad));
	//cout << "azimuth argument " << azimuth_argument << endl;
	if (azimuth_argument > 1.0) { azimuth_argument = 1.0; }
	if (azimuth_argument < -1.0) { azimuth_argument = -1.0; }
	if (HRA >= 0.) {
		_azimuth = acos(azimuth_argument) * RAD_TO_DEG;
	}
	else{
		_azimuth = (- acos(azimuth_argument)) * RAD_TO_DEG;
	}


	fUpdateSunrays();
}

void Sun::fAddNewSunray(double Rmax, double thetaMax, double Z)
{
	double fRandomR = Rmax*sqrt(1.*rand()/RAND_MAX);
	double fRandomTheta = (2. * thetaMax*(1.*rand() / RAND_MAX) - thetaMax)*DEG_TO_RAD;
	double fRandomZ = Z * (1.*rand() / RAND_MAX);
	_listOfSunrays.push_back(new Sunray(fRandomR*sin(fRandomTheta), -fRandomR*cos(fRandomTheta), fRandomZ));
}

void Sun::fUpdateSunrays()
{
	Sunray::set_azimuth(_azimuth);
	Sunray::set_elevation(_elevation);
}

void Sun::removeSunray(unsigned int j)
{
	_listOfSunrays.erase(_listOfSunrays.begin() + j);
}

bool predInterceptedSunrays(Sunray* sunray)
{
	return sunray->get_isIntercepted();
}

void Sun::clearInterceptedSunrays()
{
	int interceptedSunrays = count_if(_listOfSunrays.begin(), _listOfSunrays.end(), predInterceptedSunrays);
	int j = 0;
	vector<Sunray*>::iterator i;
	while (j < interceptedSunrays)
	{
		i = find_if(_listOfSunrays.begin(), _listOfSunrays.end(), predInterceptedSunrays);
		delete _listOfSunrays[i - _listOfSunrays.begin()];
		_listOfSunrays.erase(i);
		j++;
	}
}

void deleteSunrays(Sunray* sunray)
{
	delete sunray;
}

void Sun::clearListOfSunrays()
{
	for_each(_listOfSunrays.begin(), _listOfSunrays.end(),
		deleteSunrays);
	_listOfSunrays.clear();
}

//-----------------------------------------------------------------------
//Affichage
struct foncteurPrintSunrays
{
	ofstream* _stream;

	foncteurPrintSunrays(ofstream* stream)
	{
		_stream = stream;
	}

	void operator()(Sunray* sunray)
	{
		*_stream << sunray->get_xTarget() << " " << sunray->get_yTarget() << "\n";
	}
};

void Sun::fOutputListOfSunrays(string outputFile) const
{
	ofstream listOfSunrays;
	foncteurPrintSunrays printSunray(&listOfSunrays);
	listOfSunrays.open(outputFile.c_str());
	listOfSunrays << "X     " << "Y     " <<  endl;
	for_each(_listOfSunrays.begin(), _listOfSunrays.end(), printSunray);

	listOfSunrays.close();
}
