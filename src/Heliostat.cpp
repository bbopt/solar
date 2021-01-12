#include "Heliostat.hpp"
#include "Sun.hpp"

double Heliostat::_width = 0.;
double Heliostat::_length = 0.;
int Heliostat::_nbOfHeliostats = 0;
//Constructor
Heliostat::Heliostat(double X, double Y, double Z, double cosineEfficiency, double atmAttenuation, double& towerHeight)
:_x(X), _y(Y), _z(Z),
_cosineEfficiency(cosineEfficiency),
_atmosphericAttenuation(atmAttenuation),
_sunraysCount(0)
{
	++_nbOfHeliostats;
	_ID = _nbOfHeliostats; //L'indice 0 n'est associe a aucun heliostat

	_cTopRight = vector<double>(3, 0.);
	_cTopLeft = vector<double>(3, 0.);
	//_cBottomRight = vector<double>(3, 0.);
	_cBottomLeft = vector<double>(3, 0.);
	_cTopRightProj = vector<double>(3, 0.);
	_cTopLeftProj = vector<double>(3, 0.);
	//_cBottomRightProj = vector<double>(3, 0.);
	_cBottomLeftProj = vector<double>(3, 0.);
	_cTL_to_TR = vector < double>(2, 0.);
	_cTL_to_BL = vector <double>(2, 0.);

	_azimuthToAimpoint = atan2(_x, -_y);
	_elevationToAimpoint = atan((towerHeight - _z) /
		sqrt(_x * _x + _y * _y));
}

//Heliostat::Heliostat(double X, double Y, double Z)
//:_x(X), _y(Y), _z(Z)
//{
//	_cTopRight = vector<double>(3, 0.);
//	_cTopLeft = vector<double>(3, 0.);
//	_cBottomRight = vector<double>(3, 0.);
//	_cBottomLeft = vector<double>(3, 0.);
//	_cTopRightProj = vector<double>(3, 0.);
//	_cTopLeftProj = vector<double>(3, 0.);
//	_cBottomRightProj = vector<double>(3, 0.);
//	_cBottomLeftProj = vector<double>(3, 0.);
//	_cTL_to_TR = vector < double>(2, 0.);
//	_cTL_to_BL = vector <double>(2, 0.);
//}

//Heliostat& Heliostat::operator=(const Heliostat& heliostat)
//{
//	*this = Heliostat(heliostat._x, heliostat._y, heliostat._z);
//	return *this;
//}

//Computes the projection of heliostat pilar position on a plane perpendicular to current radiations.
//Using the same axis convention as for sun position
//Rotated frame is such that y' axis points directly towards the sun.
//Rotation is done following this order :
//1- rotate frame around z axis by an angle equivalent to the sun azimuth.
//2- rotate frame around x' axis by an angle equivalent to the sun elevation.
//unit vectors for x', y' and z' are the following:
//x' = ( cos(az), sin(az), 0 )
//y' = ( -cos(el)(sin(az), cos(el)cos(az), sin(el) )
//z' = ( sin(el)sin(az), sin(el)cos(az), cos(el) )

void Heliostat::computePilarProjection(Sun& sun)
{
	double el = sun.get_elevation() * DEG_TO_RAD;
	double az = sun.get_azimuth() * DEG_TO_RAD;

	vector<double > x_p(3, 0.);
	vector<double > y_p(3, 0.);
	vector<double > z_p(3, 0.);

	x_p[0] = cos(az);
	x_p[1] = sin(az);
	
	y_p[0] = -cos(el)*sin(az);
	y_p[1] = cos(el)*cos(az);
	y_p[2] = sin(el);

	z_p[0] = sin(el)*sin(az);
	z_p[1] = -sin(el)*cos(az);
	z_p[2] = cos(el);

	_xProj = _x *  x_p[0] + _y * x_p[1];
	_yProj = _x * y_p[0] + _y * y_p[1] + _z * y_p[2];
	_zProj = _x * z_p[0] + _y * z_p[1] + _z * z_p[2];
}

//Calculates de values of x, y and z for each of the 4 corners of the heliostat as a function of 
//its two angles of inclination (elevation and azimuth of the vector normal to its reflective surface)

void Heliostat::computeCornersPositions()
{

	_cTopRight[0] = _x + (_width / 2.)*cos(_azimuth) + (_length / 2.)*(sin(_elevation)*sin(_azimuth));  //x
	_cTopRight[1] = _y + (_width / 2.)*sin(_azimuth) - (_length / 2.)*(sin(_elevation)*cos(_azimuth));  //y
	_cTopRight[2] = _z + (_length / 2.)*cos(_elevation);												//z

	//_cBottomRight[0] = _x + (_width / 2.)*cos(_azimuth) - (_length / 2.)*sin(_elevation)*sin(_azimuth);
	//_cBottomRight[1] = _y + (_width / 2.)*sin(_azimuth) + (_length / 2.)*sin(_elevation)*cos(_azimuth);
	//_cBottomRight[2] = _z - (_length / 2.)*cos(_elevation);

	_cTopLeft[0] = _x - (_width / 2.)*cos(_azimuth) + (_length / 2.)*sin(_elevation)*sin(_azimuth);
	_cTopLeft[1] = _y - (_width / 2.)*sin(_azimuth) - (_length / 2.)*sin(_elevation)*cos(_azimuth);
	_cTopLeft[2] = _z + (_length / 2.) * cos(_elevation);

	_cBottomLeft[0] = _x - (_width / 2.)*cos(_azimuth) - (_length / 2.) * sin(_elevation)*sin(_azimuth);
	_cBottomLeft[1] = _y - (_width / 2.)*sin(_azimuth) + (_length / 2.) * sin(_elevation)*cos(_azimuth);
	_cBottomLeft[2] = _z - (_length / 2.)*cos(_elevation);	
}

void Heliostat::computeCornersProjections(Sun& sun)
{
	computeCornersPositions();
	double el = sun.get_elevation() * DEG_TO_RAD;
	double az = sun.get_azimuth() * DEG_TO_RAD;

	vector<double > x_p(3, 0.);
	vector<double > y_p(3, 0.);
	vector<double > z_p(3, 0.);

	x_p[0] = cos(az);
	x_p[1] = sin(az);

	y_p[0] = -cos(el)*sin(az);
	y_p[1] = cos(el)*cos(az);
	y_p[2] = sin(el);

	z_p[0] = sin(el)*sin(az);
	z_p[1] = -sin(el)*cos(az);
	z_p[2] = cos(el);

	_cTopRightProj[0] = _cTopRight[0] * x_p[0] + _cTopRight[1] * x_p[1];
	_cTopRightProj[1] = _cTopRight[0] * y_p[0] + _cTopRight[1] * y_p[1] + _cTopRight[2] * y_p[2];
	_cTopRightProj[2] = _cTopRight[0] * z_p[0] + _cTopRight[1] * z_p[1] + _cTopRight[2] * z_p[2];

	_cTopLeftProj[0] = _cTopLeft[0] * x_p[0] + _cTopLeft[1] * x_p[1];
	_cTopLeftProj[1] = _cTopLeft[0] * y_p[0] + _cTopLeft[1] * y_p[1] + _cTopLeft[2] * y_p[2];
	_cTopLeftProj[2] = _cTopLeft[0] * z_p[0] + _cTopLeft[1] * z_p[1] + _cTopLeft[2] * z_p[2];

	//_cBottomRightProj[0] = _cBottomRight[0] * x_p[0] + _cBottomRight[1] * x_p[1];
	//_cBottomRightProj[1] = _cBottomRight[0] * y_p[0] + _cBottomRight[1] * y_p[1] + _cBottomRight[2] * y_p[2];
	//_cBottomRightProj[2] = _cBottomRight[0] * z_p[0] + _cBottomRight[1] * z_p[1] + _cBottomRight[2] * z_p[2];

	_cBottomLeftProj[0] = _cBottomLeft[0] * x_p[0] + _cBottomLeft[1] * x_p[1];
	_cBottomLeftProj[1] = _cBottomLeft[0] * y_p[0] + _cBottomLeft[1] * y_p[1] + _cBottomLeft[2] * y_p[2];
	_cBottomLeftProj[2] = _cBottomLeft[0] * z_p[0] + _cBottomLeft[1] * z_p[1] + _cBottomLeft[2] * z_p[2];

	//update collision vectors
	_cTL_to_TR[0] = _cTopRightProj[0] - _cTopLeftProj[0];
	_cTL_to_TR[1] = _cTopRightProj[2] - _cTopLeftProj[2];

	_cTL_to_BL[0] = _cBottomLeftProj[0] - _cTopLeftProj[0];
	_cTL_to_BL[1] = _cBottomLeftProj[2] - _cTopLeftProj[2];
}

void Heliostat::computeAngles(Sun& sun, double& towerHeight)
{
	double az_to_receiver = -atan2(-_x, -_y);
	double el_to_receiver = atan2(towerHeight - (_length / 2.), sqrt(_x*_x + _y*_y));

	_azimuth = (az_to_receiver + sun.get_azimuth()*DEG_TO_RAD) / 2.;
	_elevation = (el_to_receiver + sun.get_elevation()*DEG_TO_RAD) / 2.;
}

void Heliostat::increase_sunraysCount()
{
	++_sunraysCount;
}

void Heliostat::clear_sunraysCount()
{
	_sunraysCount = 0;
}

double Heliostat::fComputeSpillage(double& apertureHeight, double& apertureWidth)
{
	double noSpillageRatio = 0.;
	double projectedWidth = abs(_cTopLeftProj[0] - _cTopRightProj[0]);
	double projectedHeight = abs(_cTopLeftProj[2] - _cBottomLeftProj[2]);

	if ( projectedWidth > apertureWidth * cos(_azimuthToAimpoint))
	{
		noSpillageRatio = (apertureWidth*cos(_azimuthToAimpoint))/projectedWidth;
	}
	else{ noSpillageRatio = 1.0; }

	if (projectedHeight > apertureHeight*cos(_elevationToAimpoint))
	{
		noSpillageRatio *= (apertureHeight*cos(_elevationToAimpoint))/projectedHeight;
	}
	else{ noSpillageRatio *= 1.0; }

	return noSpillageRatio;
}

void Heliostat::fOutputCorners(string outputFile) const
{
	ofstream projections;
	projections.open(outputFile.c_str());
	projections << "Originaux " << endl;
	projections << _cTopRight[0] << " " << _cTopRight[1] << " " << _cTopRight[2] << endl;
	projections << _cTopLeft[0] << " " << _cTopLeft[1] << " " << _cTopLeft[2] << endl;
	//projections << _cBottomRight[0] << " " << _cBottomRight[1] << " " << _cBottomRight[2] << endl;
	projections << _cBottomLeft[0] << " " << _cBottomLeft[1] << " " << _cBottomLeft[2] << endl;

	projections << "Projections " << endl;
	projections << _cTopRightProj[0] << " " << _cTopRightProj[1] << " " << _cTopRightProj[2] << endl;
	projections << _cTopLeftProj[0] << " " << _cTopLeftProj[1] << " " << _cTopLeftProj[2] << endl;
	//projections << _cBottomRightProj[0] << " " << _cBottomRightProj[1] << " " << _cBottomRightProj[2] << endl;
	projections << _cBottomLeftProj[0] << " " << _cBottomLeftProj[1] << " " << _cBottomLeftProj[2] << endl;

	projections.close();
}
