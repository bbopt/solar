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
#include "CentralReceiver.hpp"

CentralReceiver::CentralReceiver ( MoltenSalt* input               ,
				   MoltenSalt* output              ,
				   double      apertureHeight      ,
				   double      apertureWidth       ,
				   double      insulationThickness ,
				   double      tubeDi              ,
				   double      tubeThickness       ,
				   int         nbTubes               ) :
  _input                ( input                  ) ,
  _output               ( output                 ) ,
  _apertureHeight       ( apertureHeight         ) ,
  _apertureWidth        ( apertureWidth          ) ,
  _insulationThickness  ( insulationThickness    ) ,
  _tubesInsideDiameter  ( tubeDi                 ) ,
  _tubesOutsideDiameter ( tubeDi + tubeThickness ) ,
  _numberOfTubes        ( nbTubes                )   {

  if (_numberOfTubes * _tubesOutsideDiameter > PI*_apertureWidth/2.0) {
    _numberOfTubes = (int)floor(PI*_apertureWidth / _tubesOutsideDiameter);
    _numberOfPasses = 1;
  }
  else {
    _numberOfPasses = (int)floor((PI*_apertureWidth / 2.) / 
				 (_numberOfTubes*_tubesOutsideDiameter));
  }
  _receiverSurfaceArea = PI*_apertureWidth*_apertureHeight/2;
  _receiverEfficiency = 1.;
  _losses.reserve(24);
  _efficiency.reserve(24);
  _surfaceTemperature.reserve(24);
}

/*----------------------------------------------------------------------------*/
/*  Returns the amount of molten salt that can be heated to the design point  */
/*  conditions with the provided amount of energy                             */
/*----------------------------------------------------------------------------*/
double CentralReceiver::computeEnergyToFluid ( double Q_in ) {

  //The procedure will go through if energy is inputed from the field
  if (Q_in <= 0.0){ 
    _losses.push_back(0.0);
    _efficiency.push_back(0.0);
    _surfaceTemperature.push_back(0.0);
    _msRate.push_back(0.0);
    return 0.0;
  }

  double To = _output->get_temperature();
  double Ti  = _input->get_temperature();
  double eff = 0.7; //default efficiency assumed at 70%
  double m_dot_1 = 0.;
  double T_ms = (To + Ti) / 2.;
  double Del_T;
  double mu = MoltenSalt::fComputeViscosity(T_ms);
  double m_dot_2, h_ms, T_re_sur,Nus, Re, Pr, V, Q_loss, f;
  double Q_cond, Q_em, Q_ref, Q_tot1, Q_tot2, Q_abs;
  int count  = 0;
  int count2 = 0;
  
  m_dot_2 = eff * Q_in / (HEAT_CAPACITY * (To - Ti));
  
  try {

    while ( fabs(m_dot_2 - m_dot_1) > fabs(0.001*m_dot_1) && count < 150 ) {
          
      Q_tot1  = 0;
      m_dot_1 = m_dot_2;
      /*  remove */
      V  = m_dot_1 / (1.0* _numberOfTubes * MS_DENSITY * PI * pow(_tubesInsideDiameter, 2.0) / 4.0);
      Re = MS_DENSITY * V * _tubesInsideDiameter / mu;
      Pr = HEAT_CAPACITY * mu / MS_CONDUCTIVITY;
      
      if (Re <= 3000.0) {
	//for low Re the flow will be laminar and we can't use this
	//equation. For laminar flow in circular tubes with constant
	//heat flux on the surface Nus =  4.36
	Nus = 4.36;
      }
      else {
	//for high Re the flow is assumed to be turbulent
	f = pow(0.790*log(Re) - 1.64, -2.0);
	Nus = ((f / 8.0)*(Re - 1000)*Pr) / (1.0 + 12.7*sqrt(f / 8.0)*(pow(Pr, 2.0 / 3.0) - 1.0));
      }
      
      h_ms = Nus * MS_CONDUCTIVITY / _tubesInsideDiameter;
      T_re_sur = (eff * Q_in / (_numberOfTubes*_numberOfPasses * _apertureHeight) )
	* ( (_tubesOutsideDiameter / (h_ms * _tubesInsideDiameter))													+ ((_tubesOutsideDiameter / (2.*SS_COND)) * log(_tubesOutsideDiameter / _tubesInsideDiameter)) )
	+ T_ms;

      Del_T  = T_re_sur/2.0;
      count2 = 0;

      while ( fabs(Q_tot1 - Q_in) > 100 ) {
	
	if ( count2 > 0 ) {
	  if ( Q_tot1 < Q_in && Q_tot2 < Q_in )
	    T_re_sur += Del_T;
	  else if ( Q_tot1 > Q_in && Q_tot2 < Q_in ) {
	    Del_T    /= 5.0;
	    T_re_sur -= Del_T;
	  }
	  else if ( Q_tot1 < Q_in && Q_tot2 > Q_in ) {
	    Del_T    /= 5.0;
	    T_re_sur += Del_T;
	  }
	  else if ( Q_tot1 > Q_in && Q_tot2 > Q_in )
	    T_re_sur -= (Del_T*(2.0/3.0));
	}
	
	Q_tot2 = Q_tot1;
	Q_abs = (T_re_sur - T_ms)*(_numberOfTubes*_numberOfPasses * _apertureHeight)
	  / ( (_tubesOutsideDiameter / (h_ms * _tubesInsideDiameter))
	      + ((_tubesOutsideDiameter / (2.0*SS_COND))
		 * log(_tubesOutsideDiameter / _tubesInsideDiameter)) );
	  
	// Compute all losses type:
	Q_ref  = computeReflectionLosses(Q_in);
	Q_em   = computeEmissionLosses(T_re_sur);
	Q_cond = computeConductionLosses(T_re_sur);
	Q_loss = Q_ref + Q_em + Q_cond;
	
	Q_tot1 = Q_loss + Q_abs;
	
	++count2;

	if ( count2 > 1E6 ) // Added 2022-11-11
	  throw std::runtime_error ( "could not find converging value for central receiver absorbed energy (case 1)" );
      }
      
      eff     = 1 - (Q_loss / Q_in);
      m_dot_2 = eff * Q_in / (HEAT_CAPACITY * (To - Ti));
      
      ++count;
    }
    
    if ( count >= 150 )
      throw std::runtime_error ( "could not find converging value for central receiver absorbed energy (case 2)" );
    if ( eff <= 0.0 ) {
      eff     = 0.0;
      m_dot_2 = 0.0;
    }
  }
  catch ( const std::runtime_error & ) {
    m_dot_2  = 0.0;
    eff      = 0.0;
    Q_loss   = 0.0;
    T_re_sur = 0.0;
  }
  
  _input->set_massFlow  ( m_dot_2 );
  _output->set_massFlow ( m_dot_2 );
  _receiverEfficiency = eff;
  
  // Gathering data:
  _losses.push_back             ( Q_loss   );
  _efficiency.push_back         ( eff      );
  _surfaceTemperature.push_back ( T_re_sur );
  _msRate.push_back             ( m_dot_2  );
  
  return m_dot_2;
}

/*-------------------------------------------------------------------------*/
double CentralReceiver::computeEmissionLosses ( double T_sur ) const {
/*-------------------------------------------------------------------------*/

  double losses, Fr, eps_avg;
  Fr = (_apertureHeight * _apertureWidth) / (_receiverSurfaceArea + PI*pow(_apertureWidth/2.0,2.0));
  eps_avg = EPSILON_RECEIVER_SURF / (EPSILON_RECEIVER_SURF + (1 - EPSILON_RECEIVER_SURF)*Fr);

  losses = eps_avg * BOLTZMANN * (pow(T_sur, 4.0) - pow(T_ATM, 4.0)) * Fr * _receiverSurfaceArea;

  return losses;
}

/*-------------------------------------------------------------------------*/
double CentralReceiver::computeConductionLosses ( double T ) const {
/*-------------------------------------------------------------------------*/
  double k_0 = 0.043; //W/mK
  double k_1 = 1.3*pow(10.0, -4); //W/mK^2
  double k_insul = k_0 + k_1*(((T + T_ATM) / 2.0) - 273); //W/mK
  double A_out = _apertureHeight*PI*_apertureWidth/2.;
  double t_insul = _insulationThickness;
  double Re = WIND_VELOCITY * (_apertureWidth + 2.0*t_insul) / AIR_VISCOSITY;
  double Pr = 0.707; //1 atm 30 deg C
  double Nu;

  int count = 0;

  //Hilpert's correlation for Nu, finding h
  double C, m;
  if (Re >= 4.0 && Re < 40.0) {
    C = 0.911;
    m = 0.385;
  }
  if (Re >= 40.0 && Re < 4000.0) {
    C = 0.683;
    m = 0.466;
  }
  if (Re >= 4000.0 && Re < 40000.0) {
    C = 0.193;
    m = 0.618;
  }
  if (Re >= 40000.0 /*  should be upperbound  */) {
    C = 0.027;
    m = 0.805;
  }
  Nu = C*pow(Re, m)*pow(Pr, 1.0 / 3.0);
  double h_out = AIR_CONDUCTIVITY*Nu / (_apertureWidth  + 2.0*t_insul);
  
  //The overall heat transfer coefficient including outside surface convection
  double UA = (PI*0.5*_apertureWidth*_apertureHeight)
    / ( + (0.5*_apertureWidth / k_insul)* log((0.5*_apertureWidth + t_insul) / (0.5*_apertureWidth))
	+ ((0.5*_apertureWidth) / ((0.5*_apertureWidth + t_insul)*h_out)) );
  
  //The overall heat transfer coefficient for conduction through the tank wall
  double UA_cond = (PI*0.5*_apertureWidth*_apertureHeight) / 
    ((0.5*_apertureWidth / k_insul)*log((0.5*_apertureWidth + t_insul) / (0.5*_apertureWidth)) );
  
  //The heat transfer resistance for convection 
  double R_conv = 1.0 / (h_out*A_out);
  double k_conv = 1.0 / R_conv;
  double k_rad_wet = BOLTZMANN*EPSILON_OUT*A_out;
	
  //First approximation of heat transfer to outer air excluding radiation losses
  double q_out = UA*(T - T_ATM);
  double T_surf_out = T - q_out / UA_cond;
  
  double q_1, q_2;
  q_1 = 0.;
  q_2 = q_out;

  try {
    
    while ( fabs(q_2 - q_1) / q_2 > 0.001  && count < 150 ) {
    
      q_1 = q_2;
	
      T_surf_out = fSolveForT(k_rad_wet, k_conv, T, q_1, 0.01);
      
      k_insul = k_0 + k_1*(((T_surf_out + T) / 2.0) - 273.0); //W/mK
      
      UA_cond = (PI*0.5*_apertureWidth*_apertureHeight)
	/ ( +(0.5*_apertureWidth / k_insul)*log((0.5*_apertureWidth + t_insul) / (0.5*_apertureWidth)) );
      q_2 = (T - T_surf_out)*UA_cond;
      ++count;
    }
    
    if ( count >= 150 )
      throw std::runtime_error ( "could not find converging value for receiver conduction losses rate" );

  }
  catch ( const std::runtime_error & ) {
    q_2 = q_out;
  }

  q_out = q_2;
  return q_out;
}

/*----------------------------------------------------------*/
/*  This function solves the typical k1*T^4 + k2*T - q = 0  */
/*  equation using Newton's method                          */
/*----------------------------------------------------------*/
double CentralReceiver::fSolveForT ( double coef_T4, double coef_T, double T_max, double q, double eps ) const {
  double T_1, T_2;
  double g_k, Dg_k;
  int    count = 0;

  T_1 = 0.;
  T_2 = T_max;
  
  try {
    while (fabs(T_2 - T_1) > eps && count < 150) {
      T_1 = T_2;
      g_k = coef_T4*pow(T_1, 4.) + coef_T*T_1 - (q + coef_T4*pow(T_ATM, 4.0) + coef_T*T_ATM);
      Dg_k = 4.*coef_T4*pow(T_1, 3.0) + coef_T;

      T_2 = T_1 - (g_k / Dg_k);
      ++count;
			
      if (T_2 < T_ATM && T_1 < T_ATM) 
	throw std::range_error
		("Newton method gives impossible result for external surface temperature (Receiver) Setting to T_max");

    }

    if ( count >= 150 )
      throw std::runtime_error ("Newton method could not converge to an external surface temperature (Receiver)"); 
  }
  catch ( const std::runtime_error & ) {
    T_2 = T_max;
  }
  return T_2;
}

/*--------------------------------------------------------------*/
double CentralReceiver::computePressureInTubes ( void ) const {
/*--------------------------------------------------------------*/
 
  // 2022-10-21:
  // -----------
  // With SOLAR3 and x6.txt:
  // 22.63675024781991 32.246473558685054 142.86672891680968 26.92517304730389 21.093252906684114 22.0 29.857785323972358 0.36449306092140193 2.7899958004033376 870.2245260346899 49.556166148312144 20.253259690932857 0.03720503051690485 3.489708072839301 565.4948792302774 36.0 1.5450597157331276 0.023546568862376955 0.047943056313100675 8.0
  // when fid=1, we get Re=325864 and the exeption is thrown -- Hidden constraint.
  // but for fid < 1, the evaluation can complete (but is imprecise).
  
  
  double A_tubes = PI * pow ( _tubesInsideDiameter / 2.0, 2.0 ); //m^2
  double V       = _input->get_massFlow() / (MS_DENSITY * A_tubes * _numberOfTubes); //m/s
  double mu      = MoltenSalt::fComputeViscosity(_input->get_temperature());
  double Re      = V * _tubesInsideDiameter / mu;
  double Lambda;
	
  if ( V > 0 ) {
    if ( Re < 2300 )
      Lambda = 64.0 / Re;
    else if ( Re < 4000 )
      Lambda = 0.5*(64.0 / Re + 0.3164*pow(Re, -0.25));
    else if ( Re < 100000 )
      Lambda = 0.3164*pow(Re, -0.25);
    else {      
      throw std::out_of_range ( "Reynolds number out of range" );
    }
    return Lambda * _apertureHeight * _numberOfPasses * MS_DENSITY * V*V / (8.0 * A_tubes / (PI*_tubesInsideDiameter));
  }
  return 0.0;
}
