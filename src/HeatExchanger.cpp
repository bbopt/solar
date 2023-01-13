/*-------------------------------------------------------------------------------*/
/*  SOLAR - The solar thermal power plant simulator                              */
/*  https://github.com/bbopt/solar                                               */
/*                                                                               */
/*  Miguel Diago, Sebastien Le Digabel, Mathieu Lemyre-Garneau, Bastien Talgorn  */
/*                                                                               */
/*  Polytechnique Montreal / GERAD                                               */
/*  sebastien.le-digabel@polymtl.ca                                              */
/*                                                                               */
/*  This program is free software: you can redistribute it and/or modify it      */
/*  under the terms of the GNU Lesser General Public License as published by     */
/*  the Free Software Foundation, either version 3 of the License, or (at your   */
/*  option) any later version.                                                   */
/*                                                                               */
/*  This program is distributed in the hope that it will be useful, but WITHOUT  */
/*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        */
/*  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License  */
/*  for more details.                                                            */
/*                                                                               */
/*  You should have received a copy of the GNU Lesser General Public License     */
/*  along with this program. If not, see <http://www.gnu.org/licenses/>.         */
/*                                                                               */
/*-------------------------------------------------------------------------------*/
#include "HeatExchanger.hpp"

/*-------------------------------------------------------------------------*/
/*                             constructor #1                              */
/*-------------------------------------------------------------------------*/
HeatExchanger::HeatExchanger ( MoltenSalt * input  ,
			       MoltenSalt * output ,
			       Powerblock * powerblock ) :
  _input      ( input      ) ,
  _output     ( output     ) ,
  _powerblock ( powerblock )   {

  //exchangerModel ought to be 1 or 2
  //1- simple energy balance model with specified outlet conditions
  //2- efficiency model with floating conditions for the molten salt outlet
  _exchangerModel         = 1;
  _inletWaterTemperature  = T_ATM;
  _inletWaterPressure     = P_ATM;
  _outletSteamTemperature = powerblock->get_temperature();
  _outletSteamPressure    = powerblock->get_pressure();
  _heatTransferred.reserve(86400); // 24x60x60
}

/*-------------------------------------------------------------------------*/
/*                             constructor #2                              */
/*-------------------------------------------------------------------------*/
HeatExchanger::HeatExchanger ( MoltenSalt * input          ,
			       MoltenSalt * output         ,
			       Powerblock * powerblock     ,
			       double       tubesLength    ,
			       double       tubesDin       ,
			       double       tubesDout      ,
			       double       tubesSpacing   ,
			       double       baffleCut      ,
			       int          nbOfBaffles    ,
			       int          nbOfTubes      ,
			       int          passesPerShell ,
			       int          nbOfShells       ) :
  _input              ( input          ) , 
  _output             ( output         ) ,
  _powerblock         ( powerblock     ) ,
  _tubesLength        ( tubesLength    ) ,
  _tubesDin           ( tubesDin       ) ,
  _tubesDout          ( tubesDout      ) ,
  _tubesSpacing       ( tubesSpacing   ) ,
  _baffleCut          ( baffleCut      ) ,
  _nbOfBaffles        ( nbOfBaffles    ) ,
  _nbOfTubes          ( nbOfTubes      ) ,
  _nbOfPassesPerShell ( passesPerShell ) ,
  _nbOfShells         ( nbOfShells     )   {
  
  _exchangerModel         = 2;
  _inletWaterTemperature  = T_ATM;
  _inletWaterPressure     = P_ATM;
  _outletSteamTemperature = powerblock->get_temperature();
  _outletSteamPressure    = powerblock->get_pressure();

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
  int    N_tot, N_H;

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
  r_h = (_tubesSpacing/2) * sqrt(4.0 / 3.0);
  A_h = r_h *r_h * 3 * sqrt(3.0 / 4.0);
  
  //find the number of hexagonal layers to get enough tubes
  bool flag = false;
  N_H = 1;
  while (!flag) {
    if (1 + 6 * N_H*(N_H - 1) / 2 >= N_tot)
      flag = true;
    else
      ++N_H;
  }
  R_H = _tubesSpacing * (2 * N_H - 1);
  r_H = sqrt(N_tot * A_h / PI);

  u = 1 / sqrt(N_tot);
  _shellWidth = 2 * (u*R_H + (1 - u)*r_H);

  if (_baffleSpacing / _shellWidth < 0.2 || _baffleSpacing / _shellWidth > 1.0)
    throw std::invalid_argument ("inappropriate value for shell diameter");

  _shellCrossSection = PI*pow(_shellWidth / 2.0, 2.0);
  _longitudinalPitch = _tubesSpacing * sin(PI / 3.0);
  _totalRows = myround(floor(_shellWidth / _longitudinalPitch));

  _bundleArea = N_tot*A_h;
  _bypassArea = PI*pow(_shellWidth / 2.0, 2.0) - _bundleArea;
  if (_bypassArea < 0.0)
    _bypassArea = 0.0;
  _bundle_EqDiameter = sqrt(_shellWidth*_shellWidth - 4.0 * _bypassArea / PI);

  //cross flow section properties
  _crossFlowRows = myround(ceil((1.0 - 2.0 * _baffleCut)*_totalRows));
  _L_E = (_tubesSpacing - _tubesDout) * (_crossFlowRows + 1.0);

  // Hidden constraint: if _tubesSpacing == _tubesDout, then _L_E==0.0 will be a problem
  // It is checked later and it throws an exception  
    
  //computing window section properties
  _windowAngle = 2 * acos(1 - 2 * _baffleCut);
  _windowArea_FG = (PI*pow(_shellWidth, 2.0) / 4.0)*(_windowAngle / (2 * PI))
    - (_shellWidth - 2 * _shellWidth*_baffleCut)*_shellWidth * sin(_windowAngle / 2.0) / 4;
  _window_Nrows = ceil(_shellWidth*_baffleCut / _longitudinalPitch);
  _window_Ntubes = ceil(_windowArea_FG *N_tot / _shellCrossSection);
  _windowArea_FR = _window_Ntubes * PI*pow(_tubesDout / 2.0, 2.0);
  _windowArea_F = _windowArea_FG - _windowArea_FR;
  _window_EqPerimeter = _shellWidth * (_windowAngle / 2.0) + PI*_tubesDout*_window_Ntubes;
  _window_EqDiameter = 4.0 * _windowArea_F / _window_EqPerimeter;

  //Other geometric properties
  _a = _tubesSpacing / _tubesDout;
  _b = _longitudinalPitch / _tubesDout;
  _c = 0;
  _e = (_a - 1) * _tubesDout;
  if ( _b < 0.5*sqrt(2 * _a + 1) ) {
    _c = sqrt(pow(a / 2.0, 2.0) + _b*_b);
    _e = (_c - 1)*_tubesDout;
  }
  
  //We assume nozzle diameters to be so that the outlet/inlet area is the same as the
  //smallest between a window section and the crossflow sections
  _baffleSpacing * _L_E < _windowArea_F ?
			  _nozzlesArea = _baffleSpacing*_L_E/5.0 : _nozzlesArea = _windowArea_F/5.0;
  
  _nozzlesDiameter = sqrt(4 * _nozzlesArea / PI);
}

/*-------------------------------------------------------------------------*/
double HeatExchanger::fComputeRequiredMoltenSaltMassFlow ( double energyOutputRequired, double inputMSTemp ) const {
/*-------------------------------------------------------------------------*/
  if ( inputMSTemp > MELTING_POINT && energyOutputRequired > 0.0 )
    return -energyOutputRequired / ( (_output->get_temperature()-inputMSTemp) * HEAT_CAPACITY );
  return 0.0;
}

/*-------------------------------------------------------------------------*/
double HeatExchanger::fEnergyToPowerBlock ( int timeInSeconds ) {
/*-------------------------------------------------------------------------*/
  double Q_transferred =
    HEAT_CAPACITY * timeInSeconds * _input->get_massFlow() * (_input->get_temperature() - _output->get_temperature());
  if ( _exchangerModel == 1 )
    _heatTransferred.push_back ( Q_transferred );
  return Q_transferred;
}

/*-------------------------------------------------------------------------*/
double HeatExchanger::fComputeRequiredMoltenSaltMassFlow ( double energyOutputRequired ,
							   double inputMSTemp          ,
							   double maximumFlow            ) {
/*-------------------------------------------------------------------------*/
  if ( energyOutputRequired == 0.0 ) { 
    _heatTransferred.push_back ( 0.0 );
    return 0.0;
  }
  
  //Subscripts ms and w go for "molten salt" and "water" respectively
  double m_dot_ms1, m_dot_ms2, m_dot_w, Del_m; //kg/s
  double V_dot_ms, V_dot_w; //m^3/s
  double vel_w, S_D, S_T, V, V_max;
  double shellCrossArea, C1, C2, m;
  double Re_w, Re_ms;
  double Pr_w, Pr_ms;
  double h_w_i, h_w_o; //specific enthalpy
  double h_ms, h_w; //convection coefficients
  double T_in_ms, T_in_w, T_out_ms, T_out_w;
  double C_min, C_max, C_r; //m_dot * c   C_ms  C_w
  double c_ms, c_w;
  double eps, eps1, NTU1, eps_requis; // NTU
  double shellTransferArea, UA_shell;  // transferArea
  double Nus_w, Nus_ms, f;
  double visc_ms; // visc_w;
  double Q_to_water, Q1;

  int count = 0;
  
  T_in_ms  = inputMSTemp;
  T_out_ms = _output->get_temperature();
  T_in_w   = _inletWaterTemperature;
  T_out_w  = _powerblock->get_temperature();
  
  //ADD enthalpy calculation. total energy must comprise phase change energy
  //Using enthalpy value, fine m_dot_w
  //values for enthalpy must be provided by turbine model.
  h_w_i = WATER_300K_1ATM_ENTHALPY;
  h_w_o = _powerblock->get_hotEnthalpy();
  c_ms  = HEAT_CAPACITY;

  // h_ms_i = HEAT_CAPACITY*T_in_ms;
  m_dot_w = _powerblock->get_steamRate();
  Q_to_water = m_dot_w * (h_w_o - h_w_i);

  //effective c_w
  c_w = (h_w_o - h_w_i) / (T_out_w - T_in_w);

  //Finding first draft for m_dot
  m_dot_ms1 = Q_to_water / (c_ms*(T_in_ms - T_out_ms));

  if ( m_dot_ms1 > maximumFlow ) { 
    _heatTransferred.push_back(0.0);
    return 0.0;
  }

  m_dot_ms2 = 0.0;

  //Determining convection coefficient in the tube
  V_dot_w = m_dot_w / WATER_DENSITY;
  vel_w = V_dot_w / (_nbOfTubes*PI*pow(_tubesDin / 2.0, 2.0));
  Re_w = WATER_DENSITY*vel_w*_tubesDin / WATER_300K_1ATM_VISCOSITY;
  Pr_w = WATER_HEAT_CAPACITY*WATER_300K_1ATM_VISCOSITY / WATER_300K_1ATM_CONDUCTIVITY;
  if (Re_w <= 3000.0) {
    //for low Re the flow will be laminar and we can't use this
    //equation. For laminar flow in circular tubes with constant
    //heat flux on the surface Nus =  4.36
    Nus_w = 4.36;
  }
  else {
    //for high Re the flow is assumed to be turbulent
    f = pow(0.790*log(Re_w) - 1.64, -2.0);
    Nus_w = ((f / 8.)*(Re_w - 1000)*Pr_w) / (1.0 + 12.7*sqrt(f / 8.0)*(pow(Pr_w, 2.0 / 3.0) - 1.0));
  }

  h_w = Nus_w * WATER_300K_1ATM_CONDUCTIVITY / _tubesDin;
  
  //Verify if it is possible to reach the desired heat transfer with maximum flow
  m_dot_w * c_w < maximumFlow*c_ms ? (C_min = m_dot_w*c_w, C_max = maximumFlow*c_ms) :
    (C_min = maximumFlow*c_ms, C_max = m_dot_w*c_w);

  //epsilon that is required
  eps_requis = Q_to_water / (C_min * (T_in_ms - T_in_w));

  if (eps_requis > 0.95) {
    _output->set_massFlow(0);
    _input->set_massFlow(0);
    _heatTransferred.push_back(0);
    return 0.0;
  }

  C_r = C_min / C_max;
  //Determining coefficient outside the tubes
  V_dot_ms = maximumFlow / MS_DENSITY;
  S_T = _tubesSpacing;
  S_D = S_T*sqrt(5) / 2.;
  shellCrossArea = _tubesLength*((_totalRows + 1.0) * S_T)/(_nbOfBaffles + 1.);
  
  V = V_dot_ms / shellCrossArea;

  //Incropera P.439
  if (S_T*(sqrt(5) - 1.) < _tubesDout)
    V_max = V*S_T / (2.*(S_D - _tubesDout));
  else
    V_max = V*S_T / (S_T / _tubesDout);
  
  visc_ms = MoltenSalt::fComputeViscosity(0.5*(T_in_ms + T_out_ms));

  Re_ms = MS_DENSITY*V_max*_tubesDout / visc_ms;
  Pr_ms = HEAT_CAPACITY * visc_ms / MS_CONDUCTIVITY;
  C1 = fComputeC1(S_T, _tubesDout);
  m = fComputeM(S_T, _tubesDout);
  
  Nus_ms = 1.13 * C1 * pow(Re_ms, m) * pow(Pr_ms, 1. / 3.);

  if ( _totalRows < 10 ) {
    C2 = fComputeC2();
    Nus_ms *= C2;
  }
  else
    C2 = 1.0;

  h_ms = Nus_ms*MS_CONDUCTIVITY / _tubesDout;

  //According to Incropera P.688, for a shell-tubes exchanger, NUT1
  //is assumed identical for every shell with NTU = n(NTU1)
  //Note that the tubes thickness and thermal resistance is neglected.
  shellTransferArea = _nbOfTubes *_nbOfPassesPerShell* 2.0 * PI*(_tubesDin / 2.0)*_tubesLength;
  
  UA_shell = shellTransferArea *
    pow( (_tubesDin / (_tubesDout * h_ms)) + (1.0 / h_w) +
	 ((_tubesDin / (2.0*SS_COND))* log(_tubesDout / _tubesDin)) , -1.0 );

  NTU1 = UA_shell / C_min;
  eps1 = 2.0 * pow( 1.0 + C_r + sqrt(1. + C_r*C_r)*
		    (1 + exp(-NTU1*sqrt(1. + C_r*C_r))) /
		    (1 - exp(-NTU1*sqrt(1 + C_r*C_r))) , -1.0);

  _nbOfShells == 1 ? eps = eps1 :
    eps = (pow((1.0 - eps1*C_r) / (1.0 - eps1), _nbOfShells) - 1) *
    pow(pow((1.0 - eps1*C_r) / (1.0 - eps1), _nbOfShells) - C_r, -1.0);
	
  Q_to_water = 0.0;
  Q1 = 0.0;
  Del_m = m_dot_ms1;


  
  try {

    while ( ( fabs(Q_to_water - energyOutputRequired) > 1e2 || fabs(m_dot_ms1 - m_dot_ms2)/m_dot_ms1 > 0.001 )
	    && count < 500) {
    
      //Determining C_min and C_max
      m_dot_w * c_w < m_dot_ms1*c_ms ? (C_min = m_dot_w*c_w, C_max = m_dot_ms1*c_ms) :
	(C_min = m_dot_ms1*c_ms, C_max = m_dot_w*c_w);

      C_r = C_min / C_max;

      //Determining coefficient outside the tubes
      V_dot_ms = m_dot_ms1 / MS_DENSITY;
      S_T = _tubesSpacing;
      S_D = S_T*sqrt(5) / 2.;
      shellCrossArea = _tubesLength*((_totalRows + 1.0) * S_T)/(_nbOfBaffles + 1.0);
      
      V = V_dot_ms / shellCrossArea;

      //Incropera P.439
      if (S_T*(sqrt(5) - 1.0) < _tubesDout)
	V_max = V*S_T / (2.0*(S_D - _tubesDout));
      else
	V_max = V*S_T / (S_T / _tubesDout);
      
      visc_ms = MoltenSalt::fComputeViscosity(0.5*(T_in_ms + T_out_ms));

      Re_ms = MS_DENSITY*V_max*_tubesDout / visc_ms;
      Pr_ms = HEAT_CAPACITY * visc_ms / MS_CONDUCTIVITY;
      C1 = fComputeC1(S_T, _tubesDout);
      m = fComputeM(S_T, _tubesDout);

      Nus_ms = 1.13 * C1 * pow(Re_ms, m) * pow(Pr_ms, 1.0 / 3.0);

      if ( _totalRows < 10 ) {
	C2 = fComputeC2();
	Nus_ms *= C2;
      }
      else
	C2 = 1.0;

      h_ms = Nus_ms*MS_CONDUCTIVITY / _tubesDout;

      //According to Incropera P.688, for a shell-tubes exchanger, NUT1
      //is assumed identical for every shell with NTU = n(NTU1)
      //Note that the tubes thickness and thermal resistance is neglected.
      shellTransferArea = _nbOfTubes *_nbOfPassesPerShell*
	2.0 * PI*(_tubesDout / 2.0)*
	_tubesLength;

      UA_shell = shellTransferArea*pow((1.0 / h_ms) + (1.0 / h_w) +
				       ((_tubesDout / (2.0*SS_COND))
					* log(_tubesDout / _tubesDin)), -1.0);
      
      NTU1 = UA_shell / C_min;
      eps1 = 2.0 * pow(1.0 + C_r + sqrt(1.0 + C_r*C_r) *
		       (1 + exp(-NTU1*sqrt(1. + C_r*C_r))) /
		       (1 - exp(-NTU1*sqrt(1 + C_r*C_r))) , -1.0);
      
      _nbOfShells == 1 ? eps = eps1 :
	eps = (pow((1. - eps1*C_r) / (1. - eps1), _nbOfShells) - 1) *
	pow(pow((1. - eps1*C_r) / (1. - eps1), _nbOfShells) - C_r, -1.0);

      T_out_ms = T_in_ms - eps*C_min * (T_in_ms - T_in_w) / (m_dot_ms1*c_ms);

      Q_to_water = c_ms * m_dot_ms1 * (T_in_ms - T_out_ms);

      m_dot_ms2 = m_dot_ms1;
      if ( Q_to_water < energyOutputRequired && Q1 < energyOutputRequired )
	m_dot_ms1 += Del_m;
      else if ( Q_to_water > energyOutputRequired && Q1 < energyOutputRequired ) {
	Del_m /= 5.0;
	m_dot_ms1 -= Del_m;
      }
      else if ( Q_to_water < energyOutputRequired && Q1 > energyOutputRequired ) {
	Del_m /= 5.0;
	m_dot_ms1 += Del_m;
      }
      else if ( Q_to_water > energyOutputRequired && Q1 > energyOutputRequired )
	m_dot_ms1 -= Del_m*(2.0 / 3.0);

      Q1 = Q_to_water;

      ++count;
      
      if ( count > 20 && m_dot_ms1 > 4.0*maximumFlow )
	break;
    }

    if ( count >= 500 )
      throw std::range_error ( "could not converge on steam generator outlet conditions" );
  }
  catch (...) {
    T_out_ms  = T_in_ms - Q_to_water / (h_w_o - h_w_i);
    m_dot_ms1 = energyOutputRequired / (c_ms*(T_in_ms - T_out_ms));
  }

  if ( m_dot_ms1 > maximumFlow || m_dot_ms1 < 0.0 ) {
    _heatTransferred.push_back(0.);
    _output->set_massFlow(0);
    _input->set_massFlow(0);
    return 0.0;
  }
  _output->set_temperature(T_out_ms);
  _output->set_massFlow(m_dot_ms1);
  _input->set_massFlow(m_dot_ms1);
  _heatTransferred.push_back(Q_to_water);
  
  return m_dot_ms1;
}

double HeatExchanger::fComputeC1 ( double s_t, double d ) const {
  double R = s_t / d;
  return -0.2476*pow(R, 3.0) + 1.5442*pow(R, 2.0) - 3.0702*R + 2.4266;
}

double HeatExchanger::fComputeM ( double s_t, double d ) const {
  double R = s_t / d;
  return 0.0389*pow(R,3.0) - 0.2326*pow(R,2.0) + 0.4426*R + 0.2903;
}

double HeatExchanger::fComputeC2 ( void ) const {
  switch ( _totalRows ) {
  case 1: return 0.68;
  case 2: return 0.75;
  case 3: return 0.83;
  case 4: return 0.89;
  case 5: return 0.92;
  case 6: return 0.95;
  case 7: return 0.97;
  case 8: return 0.98;
  case 9: return 0.99;
  }
  return 0.0;
}

/*--------------------------------------------------------------------------*/
double HeatExchanger::computePressureInTubes ( double steamRate ) const {
/*--------------------------------------------------------------------------*/

  if ( steamRate == 0.0 )
    return 0.0;

  double A_tubes = PI*pow(_tubesDin / 2.0, 2.0); //m^2
  double V  = steamRate / (WATER_DENSITY * A_tubes * _nbOfTubes); //m/s
  double Re = V * _tubesDin / (WATER_300K_1ATM_VISCOSITY);
  double Lambda;
  
  if ( Re < 2300 )
    Lambda = 64.0 / Re;
  else if ( Re < 4000 )
    Lambda = 0.5*(64.0 / Re + 0.3164*pow(Re, -0.25));
  else if ( Re < 100000 )
    Lambda = 0.3164*pow(Re, -0.25);
  else
    throw std::out_of_range ( "Reynolds number out of range" );
  
  return Lambda * _tubesLength * _nbOfPassesPerShell * _nbOfShells * WATER_DENSITY * V*V / (8.0 * A_tubes / (PI*_tubesDin));
}

/*------------------------------------------------------------------------------------------*/
/*  Function returns the pressure drop across the heat exchanger shell for the molten salt  */
/*  Pressure drop model taken from Edward S. Gaddis; Notation is according to the paper     */
/*------------------------------------------------------------------------------------------*/
double HeatExchanger::computePressureInShells ( void ) const {

  // It is here that we check the hidden constraint _tubesSpacing = _tubesDout
  // which gives _L_E == 0.0; If not checked, w_e will be a NaN
  if ( _L_E == 0.0 )
    throw std::invalid_argument ("tubes spacing is equal to tubes outer diameter");
  
  double A_E, A_F, A_B;
  double beta;
  double S_baf, H_baf;
  double d_a, d_g, D_i, D_bun;
  double dP, dP_Q, dP_QE, dP_F, dP_Fl, dP_Ft, dP_S;
  double dP_Qo;
  double eps, eps_l, eps_t;
  double f_z, f_B, f_zl, f_zt, f_alv, f_atv;
  double m_dot_ms;
  double Re;
  double Tin_ms, To_ms;
  double eta_ms, eta_w;
  double rho;
  double R_B;
  double w_e, w_p, w_z;
  double n_wF;
  int    n_w, n_wE, N_c;

  m_dot_ms = _input->get_massFlow();
  if ( m_dot_ms == 0.0)
    return 0.0;

  //initializing variables
  // A_FG = _windowArea_FG;
  // A_FR = _windowArea_FR;
  A_F    = _windowArea_F;
  d_a    = _tubesDout;
  d_g    = _window_EqDiameter;
  D_i    = _shellWidth;
  D_bun  = _bundle_EqDiameter;
  N_c    = _totalRows;
  rho    = MS_DENSITY;
  Tin_ms = _input->get_temperature();
  To_ms  = _output->get_temperature();
  eta_ms = MoltenSalt::fComputeViscosity(0.5*(Tin_ms + To_ms));
  H_baf  = _baffleCut * D_i;
  S_baf  = _baffleSpacing;
  eta_w  = WATER_300K_1ATM_VISCOSITY; //find average viscosity;
  
  //computing dP_Qo -> crossflow sections
  n_w   = _crossFlowRows;
  A_E   = S_baf * _L_E;
  
  w_e   = (m_dot_ms / rho) / A_E;
  Re    = w_e * d_a * rho / eta_ms;
  f_zt  = pow(eta_w / eta_ms, 0.14);
  f_zl  = pow(eta_w / eta_ms, 0.57 / pow(((4 * _a*_b / PI) - 1)*Re, 0.25));
  f_atv = 2.5 + (1.2 / pow(_a - 0.85, 1.08)) + 0.4*pow(_b / _a - 1, 3.) - 0.01*pow(_a / _b - 1, 3.0);
  eps_t = f_atv / pow(Re, 0.25);

  if (_c == 0)
    f_alv = 280 * PI * (pow(pow(_b, 0.5) - 0.6, 2.0) + 0.75) / ((4 * _a*_b - PI)*pow(_a, 1.6));
  else
    f_alv = 280 * PI*(pow(pow(_b, 0.5) - 0.6, 2.0) + 0.75) / ((4 * _a*_b - PI)*pow(_c, 1.6));

  eps_l = f_alv / Re;
  eps = eps_l*f_zl + eps_t*f_zt*(1 - exp(-(Re + 200) / 1000.0));

  dP_Qo = eps * n_w * rho*pow(w_e, 2.0) / 2.0;

  //fL is supposed to be always 1 because we suppose no built in imperfections
  //fB is considered because for a low number of tubes the difference between the bundle
  //area and the shell's circular area can be significant. This also requires no additional
  //design variables.
  
  //fB
  if (_e < D_i - D_bun)
    A_B = S_baf * (D_i - D_bun - _e);
  else
    A_B = 0.0;
	
  if ( Re < 100 )
    beta = 4.5;
  else
    beta = 3.7;
  R_B = A_B / A_E;
  f_B = exp(-beta * R_B);

  dP_Q = dP_Qo*f_B;
    
  // dP_QE -> end sections are identical to other crossflow sections
  n_wE = myround(ceil(N_c * (1.0 - H_baf / D_i)));
  dP_QE = dP_Qo* f_B * n_wE / n_w;
  
  // dP_F -> window sections
	
  w_p  = (m_dot_ms / rho) / A_F;
  w_z  = sqrt(w_e * w_p);
  n_wF = 0.8*_window_Nrows;

  dP_Fl = (56.0 * n_wF / (_e * rho * w_z / eta_ms)
	   + 52.0 * _baffleSpacing / (d_g*d_g*w_z*rho / eta_ms) + 2.0) * rho * w_z * w_z / 2.;

  dP_Ft = (0.6 * n_wF + 2) *rho*w_z*w_z / 2.0;

  f_z = ( Re < 100 ) ? f_zl: f_zt;

  dP_F = f_z * sqrt(dP_Fl*dP_Fl + dP_Ft*dP_Ft);

  // review usage of H_baf definition has changed from baffle size (height) to baffle cut
  // (portion WITHOUT baffle)
  
  //dP_S for nozzle inlets and outlets
  double eps_s, w_s;
  eps_s = 2.0;
  w_s = (m_dot_ms / rho) / _nozzlesArea;
  dP_S = eps_s * rho * w_s*w_s / 2.0;

  //final value for dP
  dP = _nbOfShells*((_nbOfBaffles - 1)*dP_Q + 2*dP_QE + _nbOfBaffles*dP_F + dP_S);
 
  return dP;
}
