#include "HeatExchanger.hpp"

HeatExchanger::HeatExchanger(MoltenSalt* input, MoltenSalt* output, Powerblock* powerblock)
	:
  _input(input),
  _output(output),
  _powerblock(powerblock)
{
	//exchangerModel ought to be 1 or 2
	//1- simple energy balance model with specified outlet conditions
	//2- efficiency model with floating conditions for the molten salt outlet
	_exchangerModel = 1;
	_inletWaterTemperature = T_ATM;
	_inletWaterPressure = P_ATM;
	_outletSteamTemperature = powerblock->get_temperature();
	_outletSteamPressure = powerblock->get_pressure();

	_heatTransfered.reserve(24 * 60 * 60);
}

HeatExchanger::HeatExchanger(MoltenSalt* input, MoltenSalt* output, Powerblock* powerblock,
	double tubesLength, double tubesDin, double tubesDout, double tubesSpacing,
	double baffleCut, int nbOfBaffles,
	int nbOfTubes, int passesPerShell, int nbOfShells)
:
  _input(input), 
  _output(output),
  _powerblock(powerblock),
  _tubesLength(tubesLength),
  _tubesDin(tubesDin),
  _tubesDout(tubesDout),
  _tubesSpacing(tubesSpacing),
  _baffleCut(baffleCut),
  _nbOfBaffles(nbOfBaffles),
  _nbOfTubes(nbOfTubes),
  _nbOfPassesPerShell(passesPerShell),
  _nbOfShells(nbOfShells)
{
	_exchangerModel = 2;

	_inletWaterTemperature = T_ATM;
	_inletWaterPressure = P_ATM;
	_outletSteamTemperature = powerblock->get_temperature();
	_outletSteamPressure = powerblock->get_pressure();

	//Determining Baffles spacing
	_baffleSpacing = _tubesLength / (_nbOfBaffles + 1);

	//Determining Shells Diameter
	//To compute shell diameter, we suppose an hexagonal arrangement. The tubes are stacked
	//in hexagonal arrangement and each tiny hexagon representing a tube has a side length
	//equal to r_hex. We then find the "radius" of the big hexagon that contains all of them
	//such that the amount of small hexagons inside the big hexagon is > than the total nb of tubes.
	//We find a first radius for the shell by supposing that all tubes are organized in purely hexagonal
	//arrangements, leaving empty spaces on the side. We then also find a radius r_h such that
	//the section would have the exact same area as the sum of all small hexagons.
	//Then we do a weighted average between the two configurations, with the weight being a function
	//of the number of small hexagons (Tubes). The larger the number of tubes, the more it is possible to
	//approximate a circle with them.
	double r_h, A_h, R_H, r_H, u;
	int N_tot, N_H;

	//r_h is the hexagons side length
	//A_h is the small hexagons area
	//A_H is the big hexagon area
	//R_H is the shell's largest radius estimate
	//r_H is the shell's smallest radius estimate
	//u   is the average's weight factor
	//N_tot is the total pieces of tubes (tubes * passes)
	//N_H is the number of hexagon layers

	N_tot = _nbOfPassesPerShell * _nbOfTubes;
	//Find the side length and area of small hexagons
	r_h = (_tubesSpacing/2) * sqrt(4. / 3.);
	A_h = r_h *r_h * 3 * sqrt(3. / 4.);

	//find the number of hexagonal layers to get enough tubes
	bool flag = false;
	N_H = 1;
	while (!flag)
	{
		if (1 + 6 * N_H*(N_H - 1) / 2 >= N_tot)
		{
			flag = true;
		}
		else{ ++N_H; }
	}
	R_H = _tubesSpacing * (2 * N_H - 1);
	r_H = sqrt(N_tot * A_h / M_PI);

	u = 1 / sqrt(N_tot);
	_shellWidth = 2 * (u*R_H + (1 - u)*r_H);

	if (_baffleSpacing / _shellWidth < 0.2 || _baffleSpacing / _shellWidth > 1.0)
	  {
	    std::invalid_argument error("_baffleSpacing inappropriate value for this shell diameter.");
	    throw error;
	  }

	_shellCrossSection = M_PI*pow(_shellWidth / 2., 2.);
	_longitudinalPitch = _tubesSpacing * sin(M_PI / 3.);
	_totalRows = floor(_shellWidth / _longitudinalPitch);

	_bundleArea = N_tot*A_h;
	_bypassArea = M_PI*pow(_shellWidth / 2., 2.) - _bundleArea;
	if (_bypassArea < 0.){ _bypassArea = 0.; }
	_bundle_EqDiameter = sqrt(_shellWidth*_shellWidth - 4. * _bypassArea / M_PI);

	//cross flow section properties
	_crossFlowRows = ceil((1 - 2 * _baffleCut)*_totalRows);
	_L_E = (_tubesSpacing - _tubesDout) * (_crossFlowRows + 1);
	//computing window section properties
	_windowAngle = 2 * acos(1 - 2 * _baffleCut);
	_windowArea_FG = (M_PI*pow(_shellWidth, 2.) / 4.)*(_windowAngle / (2 * M_PI))
		- (_shellWidth - 2 * _shellWidth*_baffleCut)*_shellWidth * sin(_windowAngle / 2.) / 4;
	_window_Nrows = ceil(_shellWidth*_baffleCut / _longitudinalPitch);
	_window_Ntubes = ceil(_windowArea_FG *N_tot / _shellCrossSection);
	_windowArea_FR = _window_Ntubes * M_PI*pow(_tubesDout / 2., 2.);
	_windowArea_F = _windowArea_FG - _windowArea_FR;
	_window_EqPerimeter = _shellWidth * (_windowAngle / 2.) + M_PI*_tubesDout*_window_Ntubes;
	_window_EqDiameter = 4. * _windowArea_F / _window_EqPerimeter;
	//Other geometric properties
	
	_a = _tubesSpacing / _tubesDout;
	_b = _longitudinalPitch / _tubesDout;
	_c = 0;
	_e = (_a - 1) * _tubesDout;
	if (_b < 0.5*sqrt(2 * _a + 1))
	{
		_c = sqrt(pow(a / 2., 2.) + _b*_b);
		_e = (_c - 1)*_tubesDout;
	}
	
	//We assume nozzle diameters to be so that the outlet/inlet area is the same as the
	//smallest between a window section and the crossflow sections.
	_baffleSpacing * _L_E < _windowArea_F ?
		_nozzlesArea = _baffleSpacing*_L_E/5. : _nozzlesArea = _windowArea_F/5.;
	
	_nozzlesDiameter = sqrt(4 * _nozzlesArea / M_PI);
}

HeatExchanger::~HeatExchanger()
{}

double HeatExchanger::fComputeRequiredMoltenSaltMassFlow(double energyOutputRequired, double inputMSTemp)
{
	double msInletInternalEnergy, msOutletInternalEnergy, msMassFlow;
	if (inputMSTemp > MELTING_POINT && energyOutputRequired > 0)
	{
		msInletInternalEnergy = inputMSTemp * HEAT_CAPACITY;
		msOutletInternalEnergy = _output->get_temperature() * HEAT_CAPACITY;
		msMassFlow = -energyOutputRequired / (msOutletInternalEnergy - msInletInternalEnergy);
	}
	else{
		msMassFlow = 0.;
	}
	
	return msMassFlow;
}

double HeatExchanger::fEnergyToPowerBlock(int timeInSeconds)
{

	double inputMass = timeInSeconds * _input->get_massFlow();
	double inputTemperature = _input->get_temperature();
	double outputTemperature = _output->get_temperature();

	double Q_transfered = HEAT_CAPACITY * inputMass * (inputTemperature - outputTemperature);

	if (_exchangerModel == 1){ _heatTransfered.push_back(Q_transfered); }

	return Q_transfered;
}


/**********************************************/
/**********************************************/
/**********************************************/
/**********************************************/
/**********************************************/
/**********************************************/
/**********************************************/

double HeatExchanger::fComputeRequiredMoltenSaltMassFlow(double energyOutputRequired, double inputMSTemp, double maximumFlow)
{
	if (energyOutputRequired == 0.){ 
		_heatTransfered.push_back(0.);
		return 0.; }
	//Subscripts ms and w go for "molten salt" and "water" respectively
	double m_dot_ms1, m_dot_ms2, m_dot_w, Del_m; //kg/s
	double V_dot_ms, V_dot_w; //m^3/s
	double vel_w, S_D, S_T, V, V_max;
	// double vel_ms
	double shellCrossArea, C1, C2, m;
	// int N_L;
	double Re_w, Re_ms;
	double Pr_w, Pr_ms;
	double h_ms_i, h_w_i, h_w_o; //specific enthalpy
	// double h_ms_o;
	double h_ms, h_w; //convection coefficients
	double T_in_ms, T_in_w, T_out_ms, T_out_w;
	double C_min, C_max, C_r; //m_dot * c   C_ms  C_w
	double c_ms, c_w;
	double eps, eps1, NTU1, eps_requis; // NTU
	double shellTransferArea, UA_shell;  // transferArea
	// double q_trans;
	double Nus_w, Nus_ms, f;
	double visc_ms; // visc_w;
	double Q_to_water, Q1;
	int count = 0;

	T_in_ms = inputMSTemp;
	T_out_ms = _output->get_temperature();
	T_in_w = _inletWaterTemperature;
	T_out_w = _powerblock->get_temperature();
	//ADD enthalpy calculation. total energy must comprise phase change energy
	//Using enthalpy value, fine m_dot_w
	//values for enthalpy must be provided by turbine model.
	h_w_i = WATER_300K_1ATM_ENTHALPY;
	h_w_o = _powerblock->get_hotEnthalpy();
	c_ms = HEAT_CAPACITY;

	h_ms_i = HEAT_CAPACITY*T_in_ms;
	m_dot_w = _powerblock->get_steamRate();
	Q_to_water = m_dot_w * (h_w_o - h_w_i);
	//effective c_w
	c_w = (h_w_o - h_w_i) / (T_out_w - T_in_w);

	//Finding first draft for m_dot
	m_dot_ms1 = Q_to_water / (c_ms*(T_in_ms - T_out_ms));

	/*****************************************************/
	if (m_dot_ms1 > maximumFlow) { 
		_heatTransfered.push_back(0.);
		return 0.; }
	m_dot_ms2 = 0.;

	//Determining convection coefficient in the tube
	V_dot_w = m_dot_w / WATER_DENSITY;
	vel_w = V_dot_w / (_nbOfTubes*M_PI*pow(_tubesDin / 2., 2.));
	Re_w = WATER_DENSITY*vel_w*_tubesDin / WATER_300K_1ATM_VISCOSITY;
	Pr_w = WATER_HEAT_CAPACITY*WATER_300K_1ATM_VISCOSITY / WATER_300K_1ATM_CONDUCTIVITY; /**/
	if (Re_w <= 3000.)
	{
		//for low Re the flow will be laminar and we can't use this
		//equation. For laminar flow in circular tubes with constant
		//heat flux on the surface Nus =  4.36
		Nus_w = 4.36;
	}
	else
	{
		//for high Re the flow is assumed to be turbulent
		f = pow(0.790*log(Re_w) - 1.64, -2.);
		Nus_w = ((f / 8.)*(Re_w - 1000)*Pr_w)
			/ (1. + 12.7*sqrt(f / 8.)*(pow(Pr_w, 2. / 3.) - 1.));
	}
	h_w = Nus_w * WATER_300K_1ATM_CONDUCTIVITY / _tubesDin;

	//Findins m_dot_ms required---------------------------------------------------------------------

	//Verify if it is possible to reach the desired heat transfer with maximum flow
	m_dot_w * c_w < maximumFlow*c_ms ? (C_min = m_dot_w*c_w, C_max = maximumFlow*c_ms) :
		(C_min = maximumFlow*c_ms, C_max = m_dot_w*c_w);

	//epsilon that is required
	eps_requis = Q_to_water / (C_min * (T_in_ms - T_in_w));

	if (eps_requis > 0.95)
	{
		_output->set_massFlow(0);
		_input->set_massFlow(0);
		_heatTransfered.push_back(0);
		return 0;
	}

	C_r = C_min / C_max;
	//Determining coefficient outside the tubes
	V_dot_ms = maximumFlow / MS_DENSITY;
	S_T = _tubesSpacing;
	S_D = S_T*sqrt(5) / 2.;
	shellCrossArea = _tubesLength*((_totalRows + 1) * S_T)/(_nbOfBaffles + 1.);

	V = V_dot_ms / shellCrossArea;

	//Incropera P.439
	if (S_T*(sqrt(5) - 1.) < _tubesDout)
	{
		V_max = V*S_T / (2.*(S_D - _tubesDout));
	}
	else{
		V_max = V*S_T / (S_T / _tubesDout);
	}

	visc_ms = MoltenSalt::fComputeViscosity(0.5*(T_in_ms + T_out_ms));

	Re_ms = MS_DENSITY*V_max*_tubesDout / visc_ms;
	Pr_ms = HEAT_CAPACITY * visc_ms / MS_CONDUCTIVITY;
	C1 = fComputeC1(S_T, _tubesDout);
	m = fComputeM(S_T, _tubesDout);

	Nus_ms = 1.13 * C1 * pow(Re_ms, m) * pow(Pr_ms, 1. / 3.);

	if (_totalRows < 10)
	{
		C2 = fComputeC2();
		Nus_ms *= C2;
	}
	else{ C2 = 1.; }
	h_ms = Nus_ms*MS_CONDUCTIVITY / _tubesDout;

	//According to Incropera P.688, for a shell-tubes exchanger, NUT1
	//is assumed identical for every shell with NTU = n(NTU1)
	//Note that the tubes thickness and thermal resistance is neglected.
	shellTransferArea = _nbOfTubes *_nbOfPassesPerShell*
		2. * M_PI*(_tubesDin / 2.)*
		_tubesLength;

	UA_shell = shellTransferArea*pow(
		(_tubesDin / (_tubesDout * h_ms)) + (1. / h_w) +
		((_tubesDin / (2.*SS_COND))* log(_tubesDout / _tubesDin)),
		-1.);

	NTU1 = UA_shell / C_min;
	eps1 = 2. * pow(
		1. + C_r + sqrt(1. + C_r*C_r)
		*(1 + exp(-NTU1*sqrt(1. + C_r*C_r)))
		/ (1 - exp(-NTU1*sqrt(1 + C_r*C_r)))
		, -1.);

	_nbOfShells == 1 ? eps = eps1 :
		eps = (pow((1. - eps1*C_r) / (1. - eps1), _nbOfShells) - 1)
		* pow(pow((1. - eps1*C_r) / (1. - eps1), _nbOfShells) - C_r, -1.);

	//if (eps_requis > eps){ 
	//	return 0.; }

	//**********************************************
	//**********************************************
	//**********************************************
	Q_to_water = 0.;
	Q1 = 0.;
	Del_m = m_dot_ms1;
	try{
		while ((fabs(Q_to_water - energyOutputRequired) > 1e2
			|| fabs(m_dot_ms1 - m_dot_ms2)/m_dot_ms1 > 0.001)
			&& count < 500)
		{
			////Initializing
			//if (count > 0)
			//{
			//	m_dot_ms2 = m_dot_ms1;
			//	if (Q_to_water < energyOutputRequired
			//		&& Q1 < energyOutputRequired)
			//	{
			//		m_dot_ms1 += Del_m;
			//	}
			//	else if (Q_to_water > energyOutputRequired
			//		&& Q1 < energyOutputRequired)
			//	{
			//		Del_m /= 5.;
			//		m_dot_ms1 -= Del_m;
			//	}
			//	else if (Q_to_water < energyOutputRequired
			//		&& Q1 > energyOutputRequired)
			//	{
			//		Del_m /= 5.;
			//		m_dot_ms1 += Del_m;
			//	}
			//	else if (Q_to_water > energyOutputRequired
			//		&& Q1 > energyOutputRequired)
			//	{
			//		m_dot_ms1 -= Del_m*(2./3.);
			//	}
			//}

			//Determining C_min and C_max
			m_dot_w * c_w < m_dot_ms1*c_ms ? (C_min = m_dot_w*c_w, C_max = m_dot_ms1*c_ms) :
				(C_min = m_dot_ms1*c_ms, C_max = m_dot_w*c_w);
			/*C_max = m_dot_ms1*c_ms;
			C_min = m_dot_w*c_w;*/
			C_r = C_min / C_max;

			//Determining coefficient outside the tubes
			V_dot_ms = m_dot_ms1 / MS_DENSITY;
			S_T = _tubesSpacing;
			S_D = S_T*sqrt(5) / 2.;
			shellCrossArea = _tubesLength*((_totalRows + 1) * S_T)/(_nbOfBaffles + 1.);

			V = V_dot_ms / shellCrossArea;

			//Incropera P.439
			if (S_T*(sqrt(5) - 1.) < _tubesDout)
			{
				V_max = V*S_T / (2.*(S_D - _tubesDout));
			}
			else{
				V_max = V*S_T / (S_T / _tubesDout);
			}

			visc_ms = MoltenSalt::fComputeViscosity(0.5*(T_in_ms + T_out_ms));

			Re_ms = MS_DENSITY*V_max*_tubesDout / visc_ms;
			Pr_ms = HEAT_CAPACITY * visc_ms / MS_CONDUCTIVITY;
			C1 = fComputeC1(S_T, _tubesDout);
			m = fComputeM(S_T, _tubesDout);

			Nus_ms = 1.13 * C1 * pow(Re_ms, m) * pow(Pr_ms, 1. / 3.);

			if (_totalRows < 10)
			{
				C2 = fComputeC2();
				Nus_ms *= C2;
			}
			else{ C2 = 1.; }
			h_ms = Nus_ms*MS_CONDUCTIVITY / _tubesDout;

			//According to Incropera P.688, for a shell-tubes exchanger, NUT1
			//is assumed identical for every shell with NTU = n(NTU1)
			//Note that the tubes thickness and thermal resistance is neglected.
			shellTransferArea = _nbOfTubes *_nbOfPassesPerShell*
				2. * M_PI*(_tubesDout / 2.)*
				_tubesLength;

			UA_shell = shellTransferArea*pow((1. / h_ms) + (1. / h_w) +
				((_tubesDout / (2.*SS_COND))
				* log(_tubesDout / _tubesDin)), -1.);

			NTU1 = UA_shell / C_min;
			eps1 = 2. * pow(
				1. + C_r + sqrt(1. + C_r*C_r)
				*(1 + exp(-NTU1*sqrt(1. + C_r*C_r)))
				/ (1 - exp(-NTU1*sqrt(1 + C_r*C_r)))
				, -1.);

			_nbOfShells == 1 ? eps = eps1 :
				eps = (pow((1. - eps1*C_r) / (1. - eps1), _nbOfShells) - 1)
				* pow(pow((1. - eps1*C_r) / (1. - eps1), _nbOfShells) - C_r, -1.);

			T_out_ms = T_in_ms - eps*C_min * (T_in_ms - T_in_w) / (m_dot_ms1*c_ms);

			Q_to_water = c_ms * m_dot_ms1 * (T_in_ms - T_out_ms);

			m_dot_ms2 = m_dot_ms1;
			if (Q_to_water < energyOutputRequired
				&& Q1 < energyOutputRequired)
			{
				m_dot_ms1 += Del_m;
			}
			else if (Q_to_water > energyOutputRequired
				&& Q1 < energyOutputRequired)
			{
				Del_m /= 5.;
				m_dot_ms1 -= Del_m;
			}
			else if (Q_to_water < energyOutputRequired
				&& Q1 > energyOutputRequired)
			{
				Del_m /= 5.;
				m_dot_ms1 += Del_m;
			}
			else if (Q_to_water > energyOutputRequired
				&& Q1 > energyOutputRequired)
			{
				m_dot_ms1 -= Del_m*(2. / 3.);
			}
			Q1 = Q_to_water;

			++count;

			if (count > 20 && m_dot_ms1 > 4.*maximumFlow){ break; }
		}

		if (count >= 500){
			std::range_error noConvergence("Error : could not converge on steam generator outlet conditions.");
			throw noConvergence;
		}
	}
	catch (...)
	{
		T_out_ms = T_in_ms - Q_to_water / (h_w_o - h_w_i);
		m_dot_ms1 = energyOutputRequired / (c_ms*(T_in_ms - T_out_ms));
	}

	if (m_dot_ms1 > maximumFlow || m_dot_ms1 < 0.) {
		_heatTransfered.push_back(0.);
		_output->set_massFlow(0);
		_input->set_massFlow(0);
		return 0.;
	}
	_output->set_temperature(T_out_ms);
	_output->set_massFlow(m_dot_ms1);
	_input->set_massFlow(m_dot_ms1);
	_heatTransfered.push_back(Q_to_water);

	return m_dot_ms1;
}

double HeatExchanger::fComputeC1(double s_t, double d)
{
	double R = s_t / d;
	double C1 = -0.2476*pow(R, 3.) + 1.5442*pow(R, 2.) - 3.0702*R + 2.4266;

	return 	C1;
}

double HeatExchanger::fComputeM(double s_t, double d)
{
	double R = s_t / d;
	double M = 0.0389*pow(R,3.) - 0.2326*pow(R,2.) + 0.4426*R + 0.2903;

	return 	M;
}

double HeatExchanger::fComputeC2()
{
	double C2 = 0.0;
	switch ((int)_totalRows)
	{
	case 1: C2 = 0.68; break;
	case 2: C2 = 0.75; break;
	case 3: C2 = 0.83; break;
	case 4: C2 = 0.89; break;
	case 5: C2 = 0.92; break;
	case 6: C2 = 0.95; break;
	case 7: C2 = 0.97; break;
	case 8: C2 = 0.98; break;
	case 9: C2 = 0.99; break;
	}
	return C2;
}


//valid for P = 6.8 MPA, 600K to 1000K
double HeatExchanger::fComputeWaterEnthalpy(double T_o)
{
	double h_w = 0.146*pow(T_o, 2.) + 2182.4*T_o + 2.*pow(10.,6);

	return h_w;
}

double HeatExchanger::computeYieldPressure()
{
	double YT = SS316_YIELD_PRESSURE;
	double yieldPressure = (_tubesDout - _tubesDin)* YT / (0.5*(_tubesDin + _tubesDout));

	return yieldPressure;
}

double HeatExchanger::computePressureInTubes(double steamRate)
{
	if (steamRate == 0.){ return 0.; }
	double A_tubes = M_PI*pow(_tubesDin / 2., 2.); //m^2
	double waterM_dot = steamRate;
	double V = waterM_dot / (WATER_DENSITY * A_tubes * _nbOfTubes); //m/s
	double mu = WATER_300K_1ATM_VISCOSITY; //change value to get appropriate conditions
	double Re = V * _tubesDin / mu;
	double Lambda, S, Pm, Dh, P, L;
	std::out_of_range description("Reynolds number out of range");


	if (Re < 2300){
		Lambda = 64. / Re;
	}
	else if (Re < 4000){
		Lambda = 0.5*(64. / Re + 0.3164*pow(Re, -0.25));
	}
	else if (Re < 100000){
		Lambda = 0.3164*pow(Re, -0.25);
	}
	else{
		throw description;
	}

	S = A_tubes;
	Pm = M_PI*_tubesDin;
	Dh = 4 * S / Pm;
	L = _tubesLength * _nbOfPassesPerShell * _nbOfShells;
	P = Lambda * L * WATER_DENSITY * V*V / (2 * Dh);

	return P;
}

//----------------------------------------------------------------------------
//Function returns the pressure drop across the heat exchanger shell for the molten salt
//Pressure drop model taken from Edward S. Gaddis
//Notation is according to the paper
double HeatExchanger::computePressureInShells()
{
	double A_E, A_F, A_FR, A_FG, A_B;
	double beta;
	double S_baf, H_baf;
	double d_a, d_g, D_i, D_bun, d_S;
	double dP, dP_Q, dP_QE, dP_F, dP_Fl, dP_Ft, dP_S;
	double dP_Qo;
	double eps, eps_l, eps_t;
	double f_z, f_B, f_zl, f_zt, f_alv, f_atv;
	double L_e;
	double m_dot_w, m_dot_ms;
	double Re;
	double Tin_w, To_w, Tin_ms, To_ms;
	double eta_ms, eta_w;
	double rho;
	double R_B;
	double w_e, w_p, w_z;
	double n_wF;
	int n_w, n_wE, N_tubes, N_pass, N_tot, N_c, N_baf;

	m_dot_ms = _input->get_massFlow();
	if (m_dot_ms == 0.){ return 0.; }

	//initializing variables
	A_FG = _windowArea_FG;
	A_FR = _windowArea_FR;
	A_F = _windowArea_F;
	d_S = _tubesSpacing;
	d_a = _tubesDout;
	d_g = _window_EqDiameter;
	D_i = _shellWidth;
	D_bun = _bundle_EqDiameter;
	N_c = _totalRows;
	N_tubes = _nbOfTubes;
	N_pass = _nbOfPassesPerShell;
	N_tot = N_tubes * N_pass;
	rho = MS_DENSITY;
	Tin_ms = _input->get_temperature();
	To_ms = _output->get_temperature();
	eta_ms = MoltenSalt::fComputeViscosity(0.5*(Tin_ms + To_ms));
	H_baf = _baffleCut * D_i;
	S_baf = _baffleSpacing;
	N_baf = _nbOfBaffles;
	Tin_w = 300;
	To_w = _powerblock->get_temperature();
	eta_w = WATER_300K_1ATM_VISCOSITY; //find average viscosity;
	m_dot_w = _powerblock->get_steamRate();

	//computing dP_Qo -> crossflow sections
	n_w = _crossFlowRows;
	L_e = _L_E;
	A_E = S_baf * L_e;
	w_e = (m_dot_ms / rho) / A_E;
	Re = w_e * d_a * rho / eta_ms;
	f_zt = pow(eta_w / eta_ms, 0.14);

	f_zl = pow(eta_w / eta_ms, 0.57 / pow(((4 * _a*_b / M_PI) - 1)*Re, 0.25));

	f_atv = 2.5 + (1.2 / pow(_a - 0.85, 1.08)) + 0.4*pow(_b / _a - 1, 3.) - 0.01*pow(_a / _b - 1, 3.);
	eps_t = f_atv / pow(Re, 0.25);

	if (_c == 0)
	{
		f_alv = 280 * M_PI * (pow(pow(_b, 0.5) - 0.6, 2.) + 0.75) / ((4 * _a*_b - M_PI)*pow(_a, 1.6));
	}
	else{
		f_alv = 280 * M_PI*(pow(pow(_b, 0.5) - 0.6, 2.) + 0.75) / ((4 * _a*_b - M_PI)*pow(_c, 1.6));
	}

	eps_l = f_alv / Re;
	eps = eps_l*f_zl + eps_t*f_zt*(1 - exp(-(Re + 200) / 1000.));

	dP_Qo = eps * n_w * rho*pow(w_e, 2.) / 2.;
	
	//fL is supposed to be always 1 because we suppose no built in imperfections
	//fB is considered because for a low number of tubes the difference between the bundle
	//area and the shell's circular area can be significant. This also requires no additional
	//design variables.

	//------------------------------------------------------------------------
	//fB
	if (_e < D_i - D_bun)
	{
		A_B = S_baf * (D_i - D_bun - _e);
	}
	else { A_B = 0.; }
	
	if (Re < 100) { beta = 4.5; }
	else { beta = 3.7; }
	R_B = A_B / A_E;
	f_B = exp(-beta * R_B);
	//------------------------------------------------------------------------

	dP_Q = dP_Qo*f_B;

	//----------------------------------------------------
	/*dP_QE -> end sections are identical to other crossflow sections*/
	n_wE = ceil(N_c * (1 - H_baf / D_i));
	dP_QE = dP_Qo* f_B * n_wE / n_w;

	//----------------------------------------------------
	/*dP_F -> window sections*/
	
	w_p = (m_dot_ms / rho) / A_F;
	w_z = sqrt(w_e * w_p);
	n_wF = 0.8*_window_Nrows;

	dP_Fl = (56. * n_wF / (_e * rho * w_z / eta_ms)
		+ 52. * _baffleSpacing / (d_g*d_g*w_z*rho / eta_ms) + 2.)
		*rho * w_z * w_z / 2.;

	dP_Ft = (0.6 * n_wF + 2) *rho*w_z*w_z / 2.;

	if (Re < 100){ f_z = f_zl; }
	else { f_z = f_zt; }

	dP_F = f_z * sqrt(dP_Fl*dP_Fl + dP_Ft*dP_Ft);

	/*review usage of H_baf definition has changed from baffle size (height) to baffle cut (portion WITHOUT baffle)*/

	//dP_S for nozzle inlets and outlets
	double eps_s, w_s, d_s;
	eps_s = 2.;
	w_s = (m_dot_ms / rho) / _nozzlesArea;
	d_s = _nozzlesDiameter;
	dP_S = eps_s * rho * w_s*w_s / 2.;

	//final value for dP
	dP = _nbOfShells*((_nbOfBaffles - 1)*dP_Q + 2*dP_QE + _nbOfBaffles*dP_F + dP_S);

	return dP;
}
