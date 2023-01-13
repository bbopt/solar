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
#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

const double PI = 3.14159265358979323846264338327950;

#define DEG_TO_RAD PI/180.0
#define RAD_TO_DEG 180.0/PI
#define DECLINATION 23.45 //degrees
#define EARTH_OMEGA 0.25 //degrees per minute

#define EXTRATERRESTRIAL_INSOLATION 1368.0 //W/m^2
#define ATMOSPHERE_ATTENUATION 30.0        //percent losses of sun radiation energy before it reaches sea level on a clear day

//This file contains the values for global constants used by all components
//atmospheric properties
//air properties at 30 degC and 1 atm
#define AIR_CONDUCTIVITY 0.0263 //W/mK
#define AIR_VISCOSITY 0.00001589 // m^2/s
#define T_ATM 303.0 //air temperature K
#define WIND_VELOCITY 7.0 // m/s
#define P_ATM 103250.0 //Pascal

//Properties of liquid water
#define WATER_DENSITY 9997.0 //kg/m^3
#define WATER_HEAT_CAPACITY 4180.0 //J/kg*K
#define STEAM_785K_6_8MPA_ENTHALPY 4090437.0 //J/kg
#define WATER_300K_1ATM_ENTHALPY 1254000.0 //J/kg
#define WATER_300K_1ATM_VISCOSITY 8.32e-4
#define WATER_300K_1ATM_CONDUCTIVITY 610e-3

//Constants pertaining to radiative transfer
#define BOLTZMANN 5.67e-8 // W/m^2 K^4
#define EPSILON_MS 0.95 //emissivity of molten salt (quasi-black body assumption)
#define EPSILON_SS 0.35 //emissivity of stainless steel coating
#define EPSILON_OUT 0.3 //emissivity of the insulation outter coating
#define EPSILON_RECEIVER_SURF 0.95 // quasi black body assumption
#define RECEIVER_SURF_REFLECTIVITY 0.04

//Stainless steel properties
#define SS_COND 23.9 //stainless steel conductivity W/mK
#define SS_THICKNESS 0.04 //4 cm stainless steel wall and floor

//Molten Salt properties
//See CSP_REVIEW_MEETING 042413 presentation
//University of Alabama
#define MELTING_POINT 495.0 // 222. + 273. //Kelvin
#define MS_DENSITY 1840.0 //kg/m^3
#define HEAT_CAPACITY 1530.0 // J/(kg*K)
#define ENERGY_DENSITY 756.0 // MJ/m^3
#define MS_CONDUCTIVITY 1.16 //W/(m K)
#define PERCENT_MASS_NANO3 0.6 // kg_NaNO3/kg_moltenSalt
#define MOL_MASS_NANO3 0.08499 // kg/mol
#define MOL_MASS_KNO3  0.10110 // kg/mol

//Economics
#define COST_NANO3_KNO3 1.2 //US$/kg
#define REFERENCE_HELIOSTAT_AREA 30.0 //m^2
#define REFERENCE_FIELD_AREA 235000.0 //m^2
#define BASE_COST_MIRROR_MODULES 39.0 // US$/m^2
#define BASE_COST_DRIVES 71.0 //US$/m^2
#define BASE_COST_PEDESTALS 17.0
#define BASE_COST_CONTROLS 27.0
#define BASE_COST_WIRING 18.0
#define BASE_COST_MANUFACTURING 54.0
#define BASE_COST_INSTALLATION 11.0
#define K_H_MIRROR_MODULES 0.000827445
#define K_F_MIRROR_MODULES -7.80008e-7
#define K_H_DRIVES -0.002639289
#define K_F_DRIVES 0.0
#define K_H_PEDESTALS 0.006816719
#define K_F_PEDESTALS 4.62604e-7
#define K_H_CONTROLS -0.010308435
#define K_F_CONTROLS -9.307e-7
#define K_H_WIRING -0.00687229
#define K_F_WIRING 2.33234e-7
#define K_H_MANUFACTURING -0.001545098
#define K_F_MANUFACTURING -1.08627e-6
#define K_H_INSTALLATION -0.008572889
#define K_F_INSTALLATION 1.37257e-6
#define COST_OF_TOWER_CONSTANT 3.0e7
#define COST_OF_TOWER_LIN -285868.0
#define COST_OF_TOWER_QUAD 1835.7
#define CERAMIC_FIBER_INSULATION_COST 4960.60 //$/m^3 (using conversion to USD in 2009 and USD inflation since 2009 (around 11.5%)
#define STORAGE_TANK_FOUNDATION_COST_COEF 0.050869 //$/kg of ms
#define STORAGE_TANK_FOUNDATION_COST_CONST 241696.0
#define RECEIVER_REFERENCE_SURFACE 1571.0 //m^2
#define RECEIVER_REFERENCE_COST 51000000.0
#define RECEIVER_EXPONENT_COEFFICIENT 0.7
#define RECEIVER_MTL_COST_ESC_RATE 1.12
#define POWERBLOCK_REFERENCE_COST 86960000.0 //USD
#define POWERBLOCK_REFERENCE_POWER 115000000.0 //W
#define POWERBLOCK_EXPONENT_COEF 0.8
#define POWERBLOCK_MTL_COST_ESC_RATE 1.0
#define STEAM_GENERATOR_EXPONENT -1.112
#define STEAM_GENERATOR_CONSTANT 3.0e6
#define STEAM_GEN_REF_POWER 115.0 //MWe turbine gross)
#define STEAM_GEN_REF_EQUIP 7.3e6 //USD
#define STEAM_GEN_REF_PUMP 3.01e6
#define STEAM_GEN_REF_PIPE 1.2e6
#define STEAM_GEN_REF_CONTROL 3.9e6
#define STEAM_GEN_REF_SUPP 9.6e6
#define STEAM_GEN_SCALE    0.8

//Properties of stainless steel 316
#define SS316_YIELD_PRESSURE 290.0e6

//Siemens Turbines info
//CHECK ALL OUTLET CONDITIONS THEY ARE ERRONEOUS
//STEAM OUTLET OUGHT NOT TO BE ATMOSPHERIC CONDITIONS
//OUGHT TO BE LOW PRESSURE (0.07bar)
#define SST110_PRESSURE 130.0e5
#define SST110_TEMPERATURE 803.0  // (530. + 273.)
#define SST110_MIN 0.0
#define SST110_MAX 7.0e6
#define SST110_ENTHALPY 3415.992e3 //j/kg
#define SST110_STAGES 1 //info not found
#define SST110_OUTLET_PRESSURE 6894.757
#define SST110_OUTLET_TEMP 381.0  // (108 + 273)
#define SST110_OUTLET_ENTHALPY 2726.342e3 //j/kg
#define SST110_CONDENSE false
#define SST110_REHEAT false
#define SST120_PRESSURE 130.0e5
#define SST120_TEMPERATURE 803.0  // (530. + 273.)
#define SST120_MIN 0.0
#define SST120_MAX 10.0e6
#define SST120_ENTHALPY 3415.992e3 //j/kg
#define SST120_STAGES 1 //info not found
#define SST120_OUTLET_PRESSURE 103.25e3
#define SST120_OUTLET_TEMP 381.0 // (108 + 273)
#define SST120_OUTLET_ENTHALPY 2726.342e3 //j/kg
#define SST120_CONDENSE false
#define SST120_REHEAT false
#define SST300_PRESSURE 120.0e5
#define SST300_TEMPERATURE 793.0 // (520. + 273.)
#define SST300_MIN 10.0e6
#define SST300_MAX 50.0e6
#define SST300_ENTHALPY 3401.04e3
#define SST300_STAGES 2 //info not found
#define SST300_OUTLET_PRESSURE 103.25e3
#define SST300_OUTLET_TEMP 392.88 // (119.88 + 273)
#define SST300_OUTLET_ENTHALPY 2661.995e3 //j/kg
#define SST300_CONDENSE true
#define SST300_REHEAT false
#define SST400_PRESSURE 140.0e5
#define SST400_TEMPERATURE 813.0 // (540. + 273.)
#define SST400_MIN 30.0e6
#define SST400_MAX 65.0e6
#define SST400_ENTHALPY 3431.768e3
#define SST400_STAGES 3 //info not found
#define SST400_OUTLET_PRESSURE 103.25e3
#define SST400_OUTLET_TEMP 392.88 // (119.88 + 273)
#define SST400_OUTLET_ENTHALPY 2659.863e3 //j/kg
#define SST400_CONDENSE true
#define SST400_REHEAT false
#define SST600_PRESSURE 140.0e5
#define SST600_TEMPERATURE 813.0 // (540. + 273.)
#define SST600_MIN 5.0e6
#define SST600_MAX 100.0e6
#define SST600_ENTHALPY 3431.768e3
#define SST600_STAGES 4 //info not found
#define SST600_OUTLET_PRESSURE 103.25e3
#define SST600_OUTLET_TEMP 394.78 // (121.78 + 273)
#define SST600_OUTLET_ENTHALPY 2707.978e3 //j/kg
#define SST600_CONDENSE false
#define SST600_REHEAT false
#define SST700_PRESSURE 165.0e5
#define SST700_TEMPERATURE 858.0 // (585. + 273.)
#define SST700_MIN 20.0e6
#define SST700_MAX 175.0e6
#define SST700_ENTHALPY 3527.97e3
#define SST700_STAGES 5 //info not found
#define SST700_OUTLET_PRESSURE 103.253e3
#define SST700_OUTLET_TEMP 412.56  // (139.56 + 273)
#define SST700_OUTLET_ENTHALPY 2745.867e3
#define SST700_CONDENSE false
#define SST700_REHEAT true
#define SST800_PRESSURE 140.0e5
#define SST800_TEMPERATURE 813.0 // (540. + 273.)
#define SST800_MIN 50.0e6
#define SST800_MAX 150.0e6
#define SST800_ENTHALPY 3431.768e3
#define SST800_STAGES 5 //info not found
#define SST800_OUTLET_PRESSURE 103.25e3
#define SST800_OUTLET_TEMP 394.78 // (121.78 + 273)
#define SST800_OUTLET_ENTHALPY 2707.978e3 //j/kg
#define SST800_CONDENSE false
#define SST800_REHEAT true
#define SST900_PRESSURE 165e5
#define SST900_TEMPERATURE 858.0  // (585. + 273.)
#define SST900_MIN 49.0e6
#define SST900_MAX 250.0e6
#define SST900_ENTHALPY 3527.97e3
#define SST900_STAGES 6 //info not found
#define SST900_OUTLET_PRESSURE 103.253e3
#define SST900_OUTLET_TEMP 412.56 // (139.56 + 273)
#define SST900_OUTLET_ENTHALPY 2745.867e3
#define SST900_CONDENSE false
#define SST900_REHEAT true

#endif
