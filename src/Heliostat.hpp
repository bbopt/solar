#ifndef _HELIOSTAT_H_
#define _HELIOSTAT_H_

#include "Constants.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

class Sun;
class Sunray;

using namespace std;

class Heliostat
{
private:

	//ID
	static int _nbOfHeliostats;
	int _ID;

	//Heliostats dimentions
	static double _width;
	static double _length;

	//Position of the pilar
	double _x;
	double _y;
	double _z;

	//Projected position of the pilar
	double _xProj;
	double _yProj;
	double _zProj;

	//Projected position of heliostats corners
	//On a plane perpendicular to sun rays
	vector<double> _cTopLeft;
	vector<double> _cTopRight;
	vector<double> _cBottomLeft;
	vector<double> _cBottomRight;

	vector<double> _cTopLeftProj;
	vector<double> _cTopRightProj;
	vector<double> _cBottomLeftProj;
	vector<double> _cBottomRightProj;

	//Projected vectors representing the top and left sides of the heliostat
	//as seen from the projection plan
	vector<double> _cTL_to_TR;
	vector<double> _cTL_to_BL;

	//Optical efficiency factors
	double _cosineEfficiency;
	double _atmosphericAttenuation;
	double _spillageEfficiency;

	//Orientation (horizontal coordinates)
	//For optimal reflection, the flat mirror should always be pointing directly half way
	//between the sun and the receptor
	double _azimuth;
	double _elevation;

	double _azimuthToAimpoint;
	double _elevationToAimpoint;

	int _sunraysCount;

public:
	Heliostat(double, double, double, double, double, double&);
	/*Heliostat(double, double, double);*/
	~Heliostat();

	//-----------------------------------------------------------------------
	void computePilarProjection(Sun&);
	void computeCornersPositions();
	void computeCornersProjections(Sun&);
	void computeAngles(Sun&, double&);

	int get_ID()						 { return _ID; }
	double& get_x()						 { return _x; }
	double& get_y()						 { return _y; }
	double& get_z()						 { return _z; }
	double& get_xProj()					 { return _xProj; }
	double& get_yProj()					 { return _yProj; }
	double& get_zProj()					 { return _zProj; }
	double& get_cosineEfficiency()		 { return _cosineEfficiency; }
	double& get_atmosphericAttenuation() { return _atmosphericAttenuation; }
	double& get_azimuth()				 { return _azimuth; }
	double& get_elevation()				 { return _elevation; }
	int& get_sunraysCount()				 { return _sunraysCount; }
	vector<double>& get_cTopLeft() { return _cTopLeft; }
	vector<double>& get_cTopRight() { return _cTopRight; }
	vector<double>& get_cBottomLeft() { return _cBottomLeft; }
	vector<double>& get_cTopLeftProj() { return _cTopLeftProj; }
	vector<double>& get_cTopRightProj() { return _cTopRightProj; }
	vector<double>& get_cBottomLeftProj() { return _cBottomLeftProj; }
	vector<double>& get_cTL_to_TR() { return _cTL_to_TR; }
	vector<double>& get_cTL_to_BL() { return _cTL_to_BL; }
	//Heliostat& operator=(const Heliostat&);

	void increase_sunraysCount();
	void clear_sunraysCount();
	double fComputeSpillage(double&, double&);

	static void set_width(double width) { _width = width; }
	static void set_length(double length) { _length = length; }
	static double& get_width() { return _width; }
	static double& get_length() { return _length; }

	void fOutputCorners(string) const;
};

#endif
