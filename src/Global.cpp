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
#include "Global.hpp"

std::string separatorString() {
#ifdef _WIN32
	return "\\";
#else
	return "/";
#endif
}

/*-------------------------------------------------------------------------------*/
double kernelSmoothing ( std::vector<int>& xData, std::vector<double>& yData, int xValue ) {
/*-------------------------------------------------------------------------------*/

  double kernel;
  double yOutput   = 0.;
  double sumKernel = 0.;

  for (unsigned int i = 0; i < yData.size(); ++i) {
    kernel = exp(-(xValue*1.0 - xData[i]*1.0)*(xValue*1.0 - xData[i]*1.0) / 2.0);
    yOutput   += yData[i] * kernel;
    sumKernel += kernel;
  }
  yOutput /= sumKernel;
  return yOutput;
}

/*-------------------------------------------------------------------------------*/
double kernelSmoothing ( std::vector<double>& xData, std::vector<double>& yData, double xValue ) {
/*-------------------------------------------------------------------------------*/

  double kernel;
  double yOutput   = 0.0;
  double sumKernel = 0.0;

  for ( unsigned int i = 0; i < yData.size(); ++i ) {
    kernel     = exp(-(xValue - xData[i])*(xValue - xData[i] ) / 2.0);
    yOutput   += yData[i] * kernel;
    sumKernel += kernel;
  }
  yOutput /= sumKernel;
  return yOutput;
}

/*------------------------------------------*/
/*           custom round function          */
/*------------------------------------------*/
int  myround ( const double x ) {
  return static_cast<int> ((x < 0.0 ? -std::floor(.5-x) : std::floor(.5+x)));
}

/*------------------------------------------*/
/*                  is_int                  */
/*------------------------------------------*/
bool is_int ( const double x ) {
  return (myround(x)==x);
}

/*------------------------------------------*/
/*                  toupper                 */
/*------------------------------------------*/
std::string toupper ( std::string s ) {
  for ( size_t i = 0 ; i < s.size() ; ++i )
    s[i] = std::toupper(s[i]);
  return s;
}

/*-----------------------------------------------------------*/
/*              convert a string to an integer               */
/*-----------------------------------------------------------*/
bool string_to_int ( const std::string & s , int & i ) {
  std::stringstream ss(s);
  ss >> i;
  std::ostringstream oss;
  oss << i;
  std::string s2 = oss.str();
  if ( s != s2 )
    return false;
  return true;
}

/*-----------------------------------------------------------*/
/*              convert a string to a double                 */
/*-----------------------------------------------------------*/
/* For testing the function:                                 */
/*                                                           */
// std::cout << setprecision(12);
// std::string s;
// std::cout << "s: ";
// std::cin  >> s; 
// double x;
// if ( string_to_double(s,x) )
//   std::cout << "x=" << x << std::endl;
// else
//   std::cout << "Error" << std::endl;
/*                                                           */
/*-----------------------------------------------------------*/
bool string_to_double ( const std::string & s , double & x ) {
  x = 0.0;
  if ( s.size() == 0 )
    return false;
  char * err;
  x = strtod(s.c_str(),&err);
  return strcmp(err,"")==0;
  
  // C++ 11 version:
  // ---------------
  // x = 0.0;
  // if ( s.size() == 0 )
  //   return false;
  // size_t idx = 0;
  // try {
  //   x = std::stod(s,&idx);
  //   idx = s.size();  
  // }
  // catch ( std::invalid_argument ) {
  //   x = 0.0;
  //   return false;
  // }
  // if ( idx != s.size() ) {
  //   x = 0.0;
  //   return false;
  // }
  // return true;
}

double Q_heat_hot = 0.;
double Q_heat_cold = 0.;
double A_h_ref = REFERENCE_HELIOSTAT_AREA;
double A_f_ref = REFERENCE_FIELD_AREA;
double C_M = BASE_COST_MIRROR_MODULES;
double C_D = BASE_COST_DRIVES;
double C_P = BASE_COST_PEDESTALS;
double C_C = BASE_COST_CONTROLS;
double C_W = BASE_COST_WIRING;
double C_Ma = BASE_COST_MANUFACTURING;
double C_I = BASE_COST_INSTALLATION;
double M_Kh = K_H_MIRROR_MODULES;
double M_Kf = K_F_MIRROR_MODULES;
double D_Kh = K_H_DRIVES;
double D_Kf = K_F_DRIVES;
double P_Kh = K_H_PEDESTALS;
double P_Kf = K_F_PEDESTALS;
double C_Kh = K_H_CONTROLS;
double C_Kf = K_F_CONTROLS;
double W_Kh = K_H_WIRING;
double W_Kf = K_F_WIRING;
double Ma_Kh = K_H_MANUFACTURING;
double Ma_Kf = K_F_MANUFACTURING;
double I_Kh = K_H_INSTALLATION;
double I_Kf = K_F_INSTALLATION;


//From Ontario data (hourly avg for January)
double demandProfile_W[24] = {
	0.788644304,
	0.772938009,
	0.765316641,
	0.762298778,
	0.769051442,
	0.791984229,
	0.839497754,
	0.895343566,
	0.913218848,
	0.920895446,
	0.9262841,
	0.92505675,
	0.91838587,
	0.913501846,
	0.907996712,
	0.912633496,
	0.943412468,
	0.993129676,
	1.0,
	0.981337917,
	0.959682636,
	0.924242451,
	0.875145881,
	0.82175819
};

//From Ontario Data (hourly avg for July)
double demandProfile_S[24] = {
	0.748110459,
	0.7299731,
	0.717280819,
	0.713135609,
	0.722752429,
	0.745615021,
	0.797838866,
	0.853923489,
	0.898128934,
	0.937161092,
	0.964329556,
	0.978766229,
	0.988973205,
	0.991086972,
	0.993815489,
	0.999634953,
	1.0,
	0.982788467,
	0.9592214,
	0.94178601,
	0.940785717,
	0.899158148,
	0.835514451,
	0.780717057,
};

//turbine basic efficiency coefficients //  Condensing
double A1_e = 4.37877350555386;
double B1_e = -5.886246936292e-6;
double C1_e = 4.37824615985392e-10;
double D1_e = 4.9052617628984e-14;
double A2_e = -1.387851243025345e2;
double B2_e = -6.589550191385539e-2;
double C2_e = 9.729602056575966e-6;
double D2_e = -5.27789047337105e-10;
double A3_e = 3.98069674568167e4;
double B3_e = 2.980006167118549e1;
double C3_e = -5.19648879302579e-3;
double D3_e = 3.2495031828337e-7;
double A4_e = -4.522126353628195e6;
double B4_e = -5.198744177935268e3;
double C4_e = 1.06667544357578;
double D4_e = -7.3315284072605e-5;
double a_e = 0.0, b_e = 0.0, c_e = 0.0, d_e = 0.0;

//turbine basic efficiency coefficients //  non-Condensing
double A1_en = 4.3579522488849;
double B1_en = 1.442866040264652e-5;
double C1_en = -4.49835638024881e-9;
double D1_en = 2.79985661346127e-13;
double A2_en = -3.45898953966307e1;
double B2_en = -1.58632573100788e-1;
double C2_en = 3.29243029645988e-5;
double D2_en = -2.12357823149299e-9;
double A3_en = -3.7896243395658e4;
double B3_en = 1.11304136766883e2;
double C3_en = -2.86203302643892e-2;
double D3_en = 2.033476179976455e-6;
double A4_en = 1.6118197770998e7;
double B4_en = -3.03480606125437e4;
double C4_en = 8.7404857491443;
double D4_en = -6.61654062984245e-4;

//Powerblock efficiency correction for partial load
double A1 = -9.67554468989956e-2;
double B1 = -3.864271558704358e-1;
double C1 = 1.044184496346528e-1;
double D1 = -8.181357318296115e-3;
double A2 = -1.150178844252596e-2;
double B2 = 1.528358655177229e-2;
double C2 = -3.893637295630038e-3;
double D2 = 3.169814158242376e-4;
double A3 = 2.44354207528106e-4;
double B3 = -1.89257100557497e-4;
double C3 = 4.747677297164477e-5;
double D3 = -4.00379604790941e-6;
double A4 = -1.19638058919324e-6;
double B4 = 7.50620062251362e-7;
double C4 = -1.899443544787457e-7;
double D4 = 1.653733711580176e-8;
double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
