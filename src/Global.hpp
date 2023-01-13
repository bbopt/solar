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
#ifndef __GLOBALS_H__
#define __GLOBALS_H__

#include <vector>
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <stdlib.h>
#include <string.h>
#include "Constants.hpp"

double kernelSmoothing ( std::vector<double>& , std::vector<double>&, double );
double kernelSmoothing ( std::vector<double>&, double );
double kernelSmoothing ( std::vector<int>&, std::vector<double>&, int );

// custom round function:
int myround ( const double x );

// is_int check function:
bool is_int  ( const double x );

// put a sting in upper cases:
std::string toupper ( std::string );

// convert a string to an integer or a double:
bool string_to_int    ( const std::string & s , int    & i );
bool string_to_double ( const std::string & s , double & x );

extern double Q_heat_hot;
extern double Q_heat_cold;

//coefficients for heliostat field prices
extern double M_Kh, M_Kf, D_Kh, D_Kf, P_Kh, P_Kf, C_Kh, C_Kf, W_Kh, W_Kf, Ma_Kh, Ma_Kf, I_Kh, I_Kf;
//Base costs for heliostat price/m[2
extern double C_M, C_D, C_P, C_C, C_W, C_Ma, C_I;
extern double A_h_ref, A_f_ref;

//Demand profile coefficients (Winter data)(Ontario)
extern double demandProfile_W[24];
extern double demandProfile_S[24];

//Turbine efficiency correction for partial loads
extern double a, b, c, d, A1, B1, C1, D1, A2, B2, C2, D2, A3, B3, C3, D3, A4, B4, C4, D4;

//Turbine efficiency for condensing multistage
extern double a_e, b_e, c_e, d_e, A1_e, B1_e, C1_e, D1_e, A2_e, B2_e, C2_e, D2_e, A3_e, B3_e, C3_e, D3_e, A4_e, B4_e, C4_e, D4_e;

//Turbine efficiency for non-condensing multistage
extern double A1_en, B1_en, C1_en, D1_en, A2_en, B2_en, C2_en, D2_en, A3_en, B3_en, C3_en, D3_en, A4_en, B4_en, C4_en, D4_en;

extern std::string separatorString();

/*-------------------------------------------------------------------------*/
/*  custom exception class to catch simulation (controlled) interruptions  */
/*-------------------------------------------------------------------------*/
class Simulation_Interruption : public std::exception {
private:
  mutable std::string _what;
public:
  Simulation_Interruption ( const std::string s ) : _what(s) {}

#ifdef _MSC_VER
  virtual ~Simulation_Interruption ( void ) noexcept {}
#else
  virtual ~Simulation_Interruption ( void ) throw() {}
#endif

  const char * what ( void ) const throw() { return _what.c_str(); }
};

#endif
