#include "Sunray.hpp"

double Sunray::_azimuth = 0.;
double Sunray::_elevation = 0.;
double Sunray::_minDistance = 0.;

Sunray::Sunray(double x, double y, double z)
: _xTarget(x),
_yTarget(y),
_zTarget(z),
_isIntercepted(false),
_interceptedBy(0)
{
	_projectedTarget = std::vector<double>(3, 0.);
	std::vector<double > x_p(3, 0.);
	std::vector<double > y_p(3, 0.);
	std::vector<double > z_p(3, 0.);

	x_p[0] = cos(_azimuth*DEG_TO_RAD);
	x_p[1] = sin(_azimuth*DEG_TO_RAD);

	y_p[0] = -cos(_elevation*DEG_TO_RAD)*sin(_azimuth*DEG_TO_RAD);
	y_p[1] = cos(_elevation*DEG_TO_RAD)*cos(_azimuth*DEG_TO_RAD);
	y_p[2] = sin(_elevation*DEG_TO_RAD);

	z_p[0] = sin(_elevation*DEG_TO_RAD)*sin(_azimuth*DEG_TO_RAD);
	z_p[1] = -sin(_elevation*DEG_TO_RAD)*cos(_azimuth*DEG_TO_RAD);
	z_p[2] = cos(_elevation*DEG_TO_RAD);

	_projectedTarget[0] = _xTarget *  x_p[0] + _yTarget * x_p[1];
	_projectedTarget[1] = _xTarget * y_p[0] + _yTarget * y_p[1] + _zTarget * y_p[2];
	_projectedTarget[2] = _xTarget * z_p[0] + _yTarget * z_p[1] + _zTarget * z_p[2];
}

//-------------------------------------------------------------------------------------------------
//On verifie si la cible du point arrive a l'interieur de la projection de l'heliostat sur XZ
//La projection de l'heliostat est un parallelograme, on utilise donc la condition suivante :
//Soit Top et Left deux vecteurs.
//Target est le vecteur reliant TopLeft au rayon
//Top est le vecteur reliant TopLeft a TopRight
//Left est le vecteur reliant TopLeft a BottomLeft
//alors le rayon impacte heliostat si :
//0 < produit scalaire de Target et Top < Top*Top 
//&&
//0 < produit scalaire de Target et Left < Left*Left
//-------------------------------------------------------------------------------------------------

bool Sunray::computeCollision(Heliostat* heliostat)
{
	if (_isIntercepted == false)
	{
		if (fabs(_projectedTarget[0] - heliostat->get_xProj()) < _minDistance
			&& fabs(_projectedTarget[2] - heliostat->get_zProj()) < _minDistance)
		{
			std::vector<double> Target(2, 0.);
			double u_1, u_2, detTopLeft;
			
			Target[0] = _projectedTarget[0] - heliostat->get_cTopLeftProj()[0];
			Target[1] = _projectedTarget[2] - heliostat->get_cTopLeftProj()[2];
			
			//cout << top_target << " " << left_target << " " << left_left << " " << top_top << endl;
			detTopLeft = heliostat->get_cTL_to_TR()[0] * heliostat->get_cTL_to_BL()[1]
					- heliostat->get_cTL_to_TR()[1] * heliostat->get_cTL_to_BL()[0];
			if (detTopLeft != 0.)
			{
				u_2 = (Target[1] * heliostat->get_cTL_to_TR()[0] - Target[0] * heliostat->get_cTL_to_TR()[1])
					/ detTopLeft;
				u_1 = (Target[0] - u_2*heliostat->get_cTL_to_BL()[0]) / heliostat->get_cTL_to_TR()[0];
				if (u_1 >= 0. && u_1 <= 1. && u_2 >= 0. && u_2 <= 1.)
				{
					_isIntercepted = true;
					_interceptedBy = heliostat->get_ID();
					heliostat->increase_sunraysCount();
					return true;
				}
			}
		}
	}
	return false;
}

void Sunray::projectTarget()
{
	_cosElev = cos(_elevation*DEG_TO_RAD);
	_cosAzm = cos(_azimuth*DEG_TO_RAD);
	_sinElev = sin(_elevation*DEG_TO_RAD);
	_sinAzm = sin(_azimuth*DEG_TO_RAD);

	_projectedTarget[0] = _xTarget *  _cosAzm 
		+ _yTarget * _sinAzm;
	_projectedTarget[1] = _xTarget * (-_cosElev*_sinAzm)
		+ _yTarget * _cosElev*_cosAzm
		+ _zTarget * _sinElev;
	_projectedTarget[2] = _xTarget * _sinElev*_sinAzm
		+ _yTarget * (-_sinElev*_cosAzm)
		+ _zTarget * _cosElev;
}
