#ifndef _HELIOSTATFIELD_H_
#define _HELIOSTATFIELD_H_

#include "Heliostat.hpp"
#include "Sun.hpp"

#include <list>
#include <vector>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;

class HeliostatField
{

private:
	unsigned int _nbOfHeliostats;
	double _heliostatLength;
	double _heliostatWidth;
	double _towerHeight; //height of aim point
	double _apertureHeight;
	double _apertureWidth;
	double _minDistanceFromTower; //as a multiple of tower height
	double _maxDistanceFromTower; //as a multiple of tower height
	double _maxAngularDeviation;  //degrees
	vector < Heliostat* > _listOfHeliostats;
	vector < vector<double> > _gridLayoutAngularCoordinates;
	vector <vector <double> > _gridLayoutCartesianCoordinates;
	vector <vector <double> > _fieldLayoutCartesian;
	Sun* _sun;

	vector<double> _fieldEfficiency;
	vector<double> _powerOutput;


	void fComputeStaggeredGridLayout(); //Determines the coordinates of all potential heliostats within the boundaries of the field

	void fComputeEfficiency();				//Determines overall efficiency of a given position in the grid
	void fComputeCosineEfficiency();		//Determines Cosine efficiency for each position in the staggered grid layout
	void fComputeAtmosphericAttenuation();	//Determines atmospheric attenuation between a given position in the grid and the receiver
	void fComputeCosineAndSpillage();		//Determines losses due to the angle of the receiver to the heliostat cross section of receiver aperture
	// seen from the heliostat may be smaller than the heliostat surface, causing some of the parallel light 
	// reflected from the heliostat to miss the aperture.

	double fEvaluateFieldSurface() const;   // 0- top surface 1- west surface 2- east surface


public:
	//Construct / destruct
	HeliostatField(unsigned int, double, double, double, double, double, double, double, double, Sun&);
	~HeliostatField();

	//-----------------------------------------------------------------------
	void fConfigureField();         //Determines the position of the _nbOfheliostats heliostats
	void fGenerateField();			//Constructs the _nbOfHeliostats heliostats of the field

	void fSortByDistanceToSun();	//For a given sun azimuth, sorts heliostats in increasing order of distance to the sun

	void fComputeFluxToHeliostat(); //Calculates collision of incoming rays from the sun to heliostats
	//Assigns a shadowing efficiency coefficient to each heliostat
	//This requires to find a way to manage sunrays. SunRay object? How to generate quickly?
	double fCalculateTotalEnergyOutput();

	void fGenerateSunrays();
	void fOutputGridLayout(string) const;
	void fOutputFieldLayout(string) const;
	void fOutputAtmosphericAttenuation(string) const;
	void fPrintProjections(string) const;
	void fPrintSunraysTargetProjections(string) const;
	void fOutputRayInterception(string) const;

	double fComputeFieldEfficiency(); //Compute overall field efficiency (total energy reflected to receiver / total energy flux accross field surface)


	//-----------------------------------------------------------------------
	double get_heliostatLength() const { return _heliostatLength; }
	double get_heliostatWidth() const { return _heliostatWidth; }
	double get_towerHeight() const { return _towerHeight; }
	double get_minDistanceFromTower() const { return _minDistanceFromTower; }
	double get_maxDistanceFromTower() const { return _maxDistanceFromTower; }
	double get_maxAngularDeviation() const { return _maxAngularDeviation; }
	vector<Heliostat* > get_listOfHeliostats() const { return _listOfHeliostats; }
	vector<vector<double> > get_gridLayoutCartesianCoordinates() const { return _gridLayoutCartesianCoordinates; }
	vector<vector<double> > get_fieldLayoutCartesian() const { return _fieldLayoutCartesian; }
	vector<double>& get_fieldEfficiency(){ return _fieldEfficiency; }
	vector<double>& get_powerOutput(){ return _powerOutput; }
	Sun* get_sun() const { return _sun; }

	friend class heliostat;
};


#endif
