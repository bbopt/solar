#ifndef _SUNRAY_H_
#define _SUNRAY_H_

#include "Heliostat.hpp"

#include <vector>
#include <iostream>

class Sunray
{
private:
	double _xTarget;
	double _yTarget;
	double _zTarget;

	std::vector<double> _projectedTarget;

	static double _azimuth;
	static double _elevation;
	static double _minDistance;

	bool _isIntercepted;
	int _interceptedBy;

	//for calculations
	double _cosElev, _cosAzm, _sinElev, _sinAzm;

public:

	Sunray(double, double, double);

	//-----------------------------------------------------------------------
	void projectTarget();
	bool computeCollision(Heliostat*);

	void set_xTarget(double x) { _xTarget = x; }
	void set_yTarget(double y) { _yTarget = y; }
	void set_zTarget(double z) { _zTarget = z; }
	void set_isIntercepted(bool is) { _isIntercepted = is; }
	void set_interceptedBy(int heliostatID) { _interceptedBy = heliostatID; }

	static void set_azimuth(double az) { _azimuth = az; }
	static void set_elevation(double el) { _elevation = el; }
	static void set_minDistance(){ _minDistance = sqrt(pow(Heliostat::get_width(), 2.) + pow(Heliostat::get_length(), 2.)); }

	double get_xTarget() const { return _xTarget; }
	double get_yTarget() const { return _yTarget; }
	double get_zTarget() const { return _zTarget; }
	vector<double> get_projectedTarget() const { return _projectedTarget; }
	bool get_isIntercepted() const { return _isIntercepted; }
	int get_interceptedBy() const { return _interceptedBy; }
};

#endif
