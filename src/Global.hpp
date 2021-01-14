#ifndef _GLOBALS_H
#define _GLOBALS_H

#include <vector>
#include <iostream>
#include <cmath>
#include "Constants.hpp"

// TOTO: plein de trucs a virer ici

double kernelSmoothing(std::vector<double>& , std::vector<double>&, double);
double kernelSmoothing(std::vector<double>&, double);
double kernelSmoothing(std::vector<int>&, std::vector<double>&, int);

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

// TOTO A VIRER
// #if defined(WIN32) || defined(_WIN32) 
// 	#define PATH_SEPARATOR "\\" 
// #else 
// 	#define PATH_SEPARATOR "/" 
// #endif

#endif
