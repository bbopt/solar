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
#include "MoltenSalt.hpp"

//Coefficients for viscosity determination
//Viscosity equation :
//eta = a + b*T + c*T^2 + d*T^3
//a, b, c and d vary as a function of the molar concentration of NaNO3 and KNO3
//here x is the concentration of NaNO3.
//a = 0.0431x2 - 7.9905x + 382.68
#define a_A 0.0431
#define a_B (-7.9905)
#define a_C 382.68
//b x 10 ^2 = -0.019x2 + 3.5305x - 164.19
#define b_A (-0.019)
#define b_B 3.5305
#define b_C (-164.19)
//c x 10^5 = 0.028x2 - 5.2045x + 236.84
#define c_A 0.028
#define c_B (-5.2045)
#define c_C 236.84
//d x 10^8 = -0.0138x2 + 2.5561x - 114.23
#define d_A (-0.0138)
#define d_B 2.5561
#define d_C (-114.23)

/*-------------------------------------------------------------------------------*/
/*                                  constructor #1                               */
/*-------------------------------------------------------------------------------*/
MoltenSalt::MoltenSalt ( double temp, double pres, double masf ) :
  _temperature(temp), _pressure(pres), _massFlow(masf) {
  _enthalpy = HEAT_CAPACITY * _temperature;
  fComputeViscosity();
}

/*-------------------------------------------------------------------------------*/
/*                                  constructor #2                               */
/*-------------------------------------------------------------------------------*/
MoltenSalt::MoltenSalt ( double temp, double pres ) :
  _temperature(temp), _pressure(pres), _massFlow(0.0) {
  _enthalpy = HEAT_CAPACITY * _temperature;
  fComputeViscosity();
}

/*-------------------------------------------------------------------------------*/
/*                                 copy constructor                              */
/*-------------------------------------------------------------------------------*/
MoltenSalt::MoltenSalt ( MoltenSalt& moltenSalt ) {
  _temperature = moltenSalt._temperature;
  _enthalpy    = moltenSalt._enthalpy;
  _pressure    = moltenSalt._pressure;
  _massFlow    = moltenSalt._massFlow;
  _viscosity   = moltenSalt._viscosity;
}

/*-------------------------------------------------------------*/
void MoltenSalt::set_temperature ( double temp ) {
/*-------------------------------------------------------------*/
  _temperature = temp;
  _enthalpy = _temperature * HEAT_CAPACITY;
  fComputeViscosity();
}

/*-------------------------------------------------------------*/
void MoltenSalt::set_enthalpy ( double enth ) {
/*-------------------------------------------------------------*/
  _enthalpy = enth;
  _temperature = _enthalpy / HEAT_CAPACITY;
  fComputeViscosity();
}

/*-------------------------------------------------------------*/
void MoltenSalt::fComputeViscosity ( void ) {
/*-------------------------------------------------------------*/
  double T = _temperature;
  double C = PERCENT_MASS_NANO3;
  
  // calculate mol percent from percent mass
  double N_NaNO3 = C / MOL_MASS_NANO3;
  double N_KNO3  = (1 - C) / MOL_MASS_KNO3;
  double N_total = N_NaNO3 + N_KNO3;
  double C_Mol   = (N_NaNO3 / N_total)*100.0;
  
  double a = a_A * pow(C_Mol, 2.0) + a_B * C_Mol + a_C;
  double b = b_A * pow(C_Mol, 2.0) + b_B * C_Mol + b_C;
  double c = c_A * pow(C_Mol, 2.0) + c_B * C_Mol + c_C;
  double d = d_A * pow(C_Mol, 2.0) + d_B * C_Mol + d_C;

  _viscosity = (a + b*T*pow(10.0, -2.)
		+ c*pow(T, 2.)*pow(10.0, -5.0)
		+ d*pow(T, 3.)*pow(10.0, -8.0))/1000.0;

  // The expression originally gives the viscosity in cp
  // 1 cp = 0.001 kg/ms so we divide by 1000 to get viscosity
  // in standard SI units
}

/*-------------------------------------------------------------*/
double MoltenSalt::fComputeViscosity ( double T ) {
/*-------------------------------------------------------------*/

  double C = PERCENT_MASS_NANO3;

  // calculate mol percent from percent mass
  double N_NaNO3 = C / MOL_MASS_NANO3;
  double N_KNO3  = (1 - C) / MOL_MASS_KNO3;
  double N_total = N_NaNO3 + N_KNO3;
  double C_Mol   = (N_NaNO3 / N_total)*100.0;
  
  double a = a_A * pow(C_Mol, 2.) + a_B * C_Mol + a_C;
  double b = b_A * pow(C_Mol, 2.) + b_B * C_Mol + b_C;
  double c = c_A * pow(C_Mol, 2.) + c_B * C_Mol + c_C;
  double d = d_A * pow(C_Mol, 2.) + d_B * C_Mol + d_C;
  double mu = (a + b*T*pow(10.0, -2.0)
	       + c*pow(T, 2.0)*pow(10.0, -5.0)
	       + d*pow(T, 3.0)*pow(10.0, -8.0)) / 1000.0;
  
  return mu;
}
