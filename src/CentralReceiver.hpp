#ifndef _CENTRAL_RECEIVER_H_
#define _CENTRAL_RECEIVER_H_

#include "MoltenSalt.hpp"
#include "Constants.hpp"
#include <vector>
#include <stdexcept>
#include <exception>
#include <iostream>
#include <cmath>

class CentralReceiver
{
private:
	//fluid conditions at inlet and outlet
	MoltenSalt* _input;
	MoltenSalt* _output;

	//input design parameters
	double _apertureHeight;
	double _apertureWidth;
	double _insulationThickness;
	double _tubesInsideDiameter;
	double _tubesOutsideDiameter;
	int _numberOfTubes;

	//consequential attributes
	int _numberOfPasses;
	double _receiverSurfaceArea;
	double _receiverEfficiency;

	double computeReflectionLosses(double);
	double computeEmissionLosses(double );
	double computeConvectionLosses(double );
	double computeConductionLosses(double );
	double fSolveForT(double, double, double, double, double);

	//Simulation data
	std::vector<double> _losses;
	std::vector<double> _efficiency;
	std::vector<double> _surfaceTemperature;
	std::vector<double> _msRate;
	

public:
	CentralReceiver(MoltenSalt*, MoltenSalt*, double, double, double, double, double, int);
	~CentralReceiver();

	double& get_receiverEfficiency(){ return _receiverEfficiency; }

	//Methods
	
	double computeEnergyToFluid(double ); //function shall provide an option to calculate if the output temperature
							 //is imposed or if the mass flow is imposed.
								 
								 //for imposed temperature the mass flow will vary but varying the mass flow
								 //effectively changes the convection transfer... Or not, because the
								 // temperature being fixed means that the temperature distribution should be
								 //the same wether the flow is fast or slow.

								 //For imposted flow the output temperature will vary and the higher it is,
								 //the less heat transfer will go on. The temperature distribution must be
								 //determined along the tubes so that the steady-state transfer may be determined.
								 //that is, the amount of energy absorbed to the fluid isn't the same in the transient
								 //regime as in the steady-state regime. Though the difference shouldn't be much...

	double computeYieldPressure();
	double computePressureInTubes();

	std::vector<double>& get_losses(){ return _losses; }
	std::vector<double>& get_efficiency(){ return _efficiency; }
	std::vector<double>& get_surfaceTemperature(){ return _surfaceTemperature; }
	std::vector<double>& get_msRate(){ return _msRate; }

};

#endif
