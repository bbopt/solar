#include "CentralReceiver.hpp"

CentralReceiver::CentralReceiver(MoltenSalt* input, MoltenSalt* output,
	double apertureHeight, double apertureWidth, double insulationThickness, double tubeDi, double tubeThickness, int nbTubes)
	:_input(input),
	_output(output),
	_apertureHeight(apertureHeight),
	_apertureWidth(apertureWidth),
	_insulationThickness(insulationThickness),
	_tubesInsideDiameter(tubeDi),
	_tubesOutsideDiameter(tubeDi + tubeThickness),
	_numberOfTubes(nbTubes)
{

	if (_numberOfTubes * _tubesOutsideDiameter > M_PI*_apertureWidth/2.)
	{
		_numberOfTubes = (int)floor(M_PI*_apertureWidth / _tubesOutsideDiameter);
		_numberOfPasses = 1;
	}
	else
	{
		_numberOfPasses = (int)floor((M_PI*_apertureWidth / 2.) / 
			(_numberOfTubes*_tubesOutsideDiameter));
	}
	_receiverSurfaceArea = M_PI*_apertureWidth*_apertureHeight/2;
	_receiverEfficiency = 1.;
	_losses.reserve(24);
	_efficiency.reserve(24);
	_surfaceTemperature.reserve(24);
}

//-----------------------------------------------------------------------
//Returns the amount of molten salt that can be heated to the design point
//conditions with the provided amount of energy.
//-----------------------------------------------------------------------
double CentralReceiver::computeEnergyToFluid(double Q_in)
{
	//The procedure will go through if energy is inputed from the field
	if (Q_in <= 0.){ 
		_losses.push_back(0.);
		_efficiency.push_back(0.);
		_surfaceTemperature.push_back(0.);
		_msRate.push_back(0.);
		return 0.; }

	double To = _output->get_temperature();
	double Ti  = _input->get_temperature();
	double eff = 0.7; //default efficiency assumed at 70%
	double m_dot_1 = 0.;
	double T_ms = (To + Ti) / 2.;
	double Del_T;
	double mu = MoltenSalt::fComputeViscosity(T_ms);
	double m_dot_2, h_ms, T_re_sur,Nus, Re, Pr, V, Q_loss, Q_no_reflect, f;
	double q_no_reflect, Q_cond, Q_em, Q_ref, Q_tot1, Q_tot2, Q_abs;
	// double q_in = Q_in / _receiverSurfaceArea; //Q_in per tubes unit length
	int count = 0;
	int count2 = 0;
	// Use fluid property at average T
	// double tubesPerMeter = 1. / _tubesOutsideDiameter;

	Q_no_reflect = Q_in - computeReflectionLosses(Q_in);
	q_no_reflect = Q_no_reflect / _receiverSurfaceArea;
	m_dot_2 = eff * Q_in / (HEAT_CAPACITY * (To - Ti));

	try{
		while (fabs(m_dot_2 - m_dot_1) > fabs(0.001*m_dot_1) && count < 150)
		{
			Q_tot1 = 0;
			m_dot_1 = m_dot_2;
			/*  remove */
			V = m_dot_1 / (1.* _numberOfTubes * MS_DENSITY * M_PI * pow(_tubesInsideDiameter, 2.) / 4.);
			Re = MS_DENSITY * V * _tubesInsideDiameter / mu;
			Pr = HEAT_CAPACITY * mu / MS_CONDUCTIVITY;

			if (Re <= 3000.)
			{
				//for low Re the flow will be laminar and we can't use this
				//equation. For laminar flow in circular tubes with constant
				//heat flux on the surface Nus =  4.36
				Nus = 4.36;
			}
			else
			{
				//for high Re the flow is assumed to be turbulent
				f = pow(0.790*log(Re) - 1.64, -2.);
				Nus = ((f / 8.)*(Re - 1000)*Pr)
					/ (1. + 12.7*sqrt(f / 8.)*(pow(Pr, 2. / 3.) - 1.));
			}


			h_ms = Nus * MS_CONDUCTIVITY / _tubesInsideDiameter;
			T_re_sur = (eff * Q_in / (_numberOfTubes*_numberOfPasses * _apertureHeight) /*/ tubesPerMeter*/)*(
				(_tubesOutsideDiameter / (h_ms * _tubesInsideDiameter))
				+ ((_tubesOutsideDiameter / (2.*SS_COND))
				* log(_tubesOutsideDiameter / _tubesInsideDiameter))
				)
				+ T_ms;

			Del_T = T_re_sur/2.;
			count2 = 0;
			while (fabs(Q_tot1 - Q_in) > 100)
			{
				if (count2 > 0)
				{
					if (Q_tot1 < Q_in && Q_tot2 < Q_in)
					{
						T_re_sur += Del_T;
					}
					else if (Q_tot1 > Q_in && Q_tot2 < Q_in)
					{
						Del_T /= 5.;
						T_re_sur -= Del_T;
					}
					else if (Q_tot1 < Q_in && Q_tot2 > Q_in)
					{
						Del_T /= 5.;
						T_re_sur += Del_T;
					}
					else if (Q_tot1 > Q_in && Q_tot2 > Q_in)
					{
						T_re_sur -= (Del_T*(2./3.));
					}
				}

				Q_tot2 = Q_tot1;
				Q_abs = (T_re_sur - T_ms)*(_numberOfTubes*_numberOfPasses * _apertureHeight)
					/ (
					(_tubesOutsideDiameter / (h_ms * _tubesInsideDiameter))
					+ ((_tubesOutsideDiameter / (2.*SS_COND))
					* log(_tubesOutsideDiameter / _tubesInsideDiameter))
					);

				//compute all losses type
				Q_ref = computeReflectionLosses(Q_in);
				Q_em = computeEmissionLosses(T_re_sur);
				Q_cond = computeConductionLosses(T_re_sur);
				Q_loss = Q_ref + Q_em + Q_cond;

				Q_tot1 = Q_loss + Q_abs;

				++count2;
			}

			
			eff = 1 - (Q_loss / Q_in);
			m_dot_2 = eff * Q_in / (HEAT_CAPACITY * (To - Ti));

			++count;
		}

		if (count >= 150){
			std::runtime_error noConvergence("Could not find converging value for central receiver absorbed energy");
			throw noConvergence;
		}

		if (eff <= 0.){
			eff = 0.;
			m_dot_2 = 0.;
		}
	}
	catch (std::runtime_error& e)
	{
		m_dot_2 = 0.;
		eff = 0.;
		Q_loss = 0.;
		T_re_sur = 0.;
		//std::cerr << e.what() << std::endl;
	}
	
	_input->set_massFlow(m_dot_2);
	_output->set_massFlow(m_dot_2);
	_receiverEfficiency = eff;

	//Gathering data
	_losses.push_back(Q_loss);
	_efficiency.push_back(eff);
	_surfaceTemperature.push_back(T_re_sur);
	_msRate.push_back(m_dot_2);

	return m_dot_2;
}

CentralReceiver::~CentralReceiver()
{}

double CentralReceiver::computeReflectionLosses(double Q_in)
{
	return Q_in * RECEIVER_SURF_REFLECTIVITY * ((_apertureHeight*_apertureWidth)/_receiverSurfaceArea);
}



double CentralReceiver::computeEmissionLosses(double T_sur)
{
	double losses, Fr, eps_avg;
	Fr = (_apertureHeight * _apertureWidth) / (_receiverSurfaceArea + M_PI*pow(_apertureWidth/2.,2.));
	eps_avg = EPSILON_RECEIVER_SURF / (EPSILON_RECEIVER_SURF + (1 - EPSILON_RECEIVER_SURF)*Fr);

	losses = eps_avg * BOLTZMANN * (pow(T_sur, 4.) - pow(T_ATM, 4.)) 
		* Fr * _receiverSurfaceArea;

	return losses;
}

double CentralReceiver::computeConductionLosses(double T)
{
	double k_0 = 0.043; //W/mK
	double k_1 = 1.3*pow(10., -4); //W/mK^2
	double k_insul = k_0 + k_1*(((T + T_ATM) / 2.) - 273); //W/mK
	double A_out = _apertureHeight*M_PI*_apertureWidth/2.;
	double t_insul = _insulationThickness;
	double Re = WIND_VELOCITY * (_apertureWidth + 2.*t_insul) / AIR_VISCOSITY;
	double Pr = 0.707; //1 atm 30 deg C
	double Nu;
	int count = 0;

	//Hilpert's correlation for Nu, finding h
	double C, m;
	if (Re >= 4. && Re < 40.){ C = 0.911; m = 0.385; }
	if (Re >= 40. && Re < 4000.) { C = 0.683; m = 0.466; }
	if (Re >= 4000. && Re < 40000.) { C = 0.193; m = 0.618; }
	if (Re >= 40000. /*  should be upperbound  */){ C = 0.027; m = 0.805; }
	Nu = C*pow(Re, m)*pow(Pr, 1. / 3.);
	double h_out = AIR_CONDUCTIVITY*Nu / (_apertureWidth  + 2.*t_insul);

	//The overall heat transfer coefficient including outside surface convection
	double UA = (M_PI*0.5*_apertureWidth*_apertureHeight) / (
		+ (0.5*_apertureWidth / k_insul)* log((0.5*_apertureWidth + t_insul) / (0.5*_apertureWidth))
		+ ((0.5*_apertureWidth) / ((0.5*_apertureWidth + t_insul)*h_out))
		);

	//The overall heat transfer coefficient for conduction through the tank wall
	double UA_cond = (M_PI*0.5*_apertureWidth*_apertureHeight) / 
		((0.5*_apertureWidth / k_insul)*log((0.5*_apertureWidth + t_insul) / (0.5*_apertureWidth))
		);

	//The heat transfer resistance for convection 
	double R_conv = 1 / (h_out*A_out);
	double k_conv = 1 / R_conv;
	double k_rad_wet = BOLTZMANN*EPSILON_OUT*A_out;

	//First approximation of heat transfer to outer air excluding radiation losses
	double q_out = UA*(T - T_ATM);
	double T_surf_out = T - q_out / UA_cond;

	double q_1, q_2;
	q_1 = 0.;
	q_2 = q_out;

	try{
		while (fabs(q_2 - q_1) / q_2 > 0.001  && count < 150)
		{
			q_1 = q_2;

			T_surf_out = fSolveForT(k_rad_wet, k_conv, T, q_1, 0.01);

			k_insul = k_0 + k_1*(((T_surf_out + T) / 2.) - 273); //W/mK

			UA_cond = (M_PI*0.5*_apertureWidth*_apertureHeight) / (
				+(0.5*_apertureWidth / k_insul)*log((0.5*_apertureWidth + t_insul) / (0.5*_apertureWidth))
				);

			q_2 = (T - T_surf_out)*UA_cond;
			++count;
		}

		if (count >= 150){ std::runtime_error noConvergence("Could not find converging value for receiver conduction losses rate");
		throw noConvergence;
		}

	}
	catch (std::runtime_error& e){
		q_2 = q_out;
		//std::cerr << e.what() << std::endl;
	}

	q_out = q_2;
	return q_out;
}

//This function solves the typical k1*T^4 + k2*T - q = 0 equation using Newton<s method.
double CentralReceiver::fSolveForT(double coef_T4, double coef_T, double T_max, double q, double eps)
{
	double T_1, T_2;
	double g_k, Dg_k;
	int count = 0;

	T_1 = 0.;
	T_2 = T_max;

	try{
		while (fabs(T_2 - T_1) > eps && count < 150)
		{
			T_1 = T_2;
			g_k = coef_T4*pow(T_1, 4.) + coef_T*T_1 - (q + coef_T4*pow(T_ATM, 4.) + coef_T*T_ATM);
			Dg_k = 4.*coef_T4*pow(T_1, 3.) + coef_T;
			//**** insert throw if Dg_k = 0.
			T_2 = T_1 - (g_k / Dg_k);

			++count;

			if (T_2 < T_ATM && T_1 < T_ATM) 
			{ 
				std::range_error invalidTemperature("Newton's method gives impossible result for external surface temperature (Receiver) Setting to T_max");
				throw invalidTemperature;
			}
		}

		if (count >= 150){ std::runtime_error noConvergence("Newton's method could not converge to an external surface temperature (Receiver)"); 
		throw noConvergence;
		}
	}
	catch (std::runtime_error& e){
		std::cout << e.what() << std::endl;
		T_2 = T_max;
	}
	return T_2;
}

double CentralReceiver::computeYieldPressure()
{
	double YT = SS316_YIELD_PRESSURE;
	double yieldPressure = (_tubesOutsideDiameter - _tubesInsideDiameter)* YT / (0.5*(_tubesInsideDiameter + _tubesOutsideDiameter));

	return yieldPressure;
}

double CentralReceiver::computePressureInTubes()
{
	double A_tubes = M_PI*pow(_tubesInsideDiameter / 2., 2.); //m^2
	double V = _input->get_massFlow() / (MS_DENSITY * A_tubes * _numberOfTubes); //m/s
	double mu = MoltenSalt::fComputeViscosity(_input->get_temperature());
	double Re = V * _tubesInsideDiameter / mu;
	double Lambda, S, Pm, Dh, P, L;
	std::out_of_range description("Reynolds number out of range");

	if (V > 0){
		if (Re < 2300){
			Lambda = 64. / Re;
		}
		else if (Re < 4000){
			Lambda = 0.5*(64. / Re + 0.3164*pow(Re, -0.25));
		}
		else if(Re < 100000){
			Lambda = 0.3164*pow(Re, -0.25);
		}
		else{
			throw description;
		}
		S = A_tubes;
		Pm = M_PI*_tubesInsideDiameter;
		Dh = 4 * S / Pm;
		L = _apertureHeight * _numberOfPasses;
		P = Lambda * L * MS_DENSITY * V*V / (2 * Dh);
	}
	else { P = 0; }

	

	return P;
}
