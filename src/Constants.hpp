#ifndef _CONSTANTS_H_
#define _CONSTANTS_H_

#define PI 3.14159265358979323846 //rads
#define DEG_TO_RAD PI/180
#define RAD_TO_DEG 180./PI
#define DECLINATION 23.45 //degrees
#define EARTH_OMEGA 0.25 //degrees per minute

#define EXTRATERRESTRIAL_INSOLATION 1368. //W/m^2
#define ATMOSPHERE_ATTENUATION 30 //percent losses of sun radiation energy before it reaches sea level on a clear day

//This file contains the values for global constants used by all components
//atmospheric properties
//air properties at 30 degC and 1 atm
#define AIR_CONDUCTIVITY (26.3/1000.) //W/mK
#define AIR_VISCOSITY (15.89/1000000.) // m^2/s
#define T_ATM 303. //air temperature K
#define WIND_VELOCITY 7. // m/s
#define P_ATM 103.25*1000 //Pascal

//Properties of liquid water
#define WATER_DENSITY (9997.) //kg/m^3
#define WATER_HEAT_CAPACITY (4180.) //J/kg*K
#define STEAM_785K_6_8MPA_ENTHALPY (4090437.) //J/kg
#define WATER_300K_1ATM_ENTHALPY (1254000.) //J/kg
#define WATER_300K_1ATM_VISCOSITY (8.32*pow(10.,-4.))
#define WATER_300K_1ATM_CONDUCTIVITY (610*pow(10.,-3.))

//Constants pertaining to radiative transfer
#define BOLTZMANN (5.67* pow(10.,-8.)) // W/m^2 K^4
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
#define MELTING_POINT (222. + 273.) //Kelvin
#define MS_DENSITY (1.84 *1000.) //kg/m^3
#define HEAT_CAPACITY (1.53 * 1000.) // J/(kg*K)
#define ENERGY_DENSITY (756.) // MJ/m^3
#define MS_CONDUCTIVITY 1.16 //W/(m K)
#define PERCENT_MASS_NANO3 0.6 // kg_NaNO3/kg_moltenSalt
#define MOL_MASS_NANO3 0.08499 // kg/mol
#define MOL_MASS_KNO3  0.10110 // kg/mol

//Economics
#define COST_NANO3_KNO3 (1.2) //US$/kg
#define REFERENCE_HELIOSTAT_AREA (30.) //m^2
#define REFERENCE_FIELD_AREA (235000.) //m^2
#define BASE_COST_MIRROR_MODULES (39.) // US$/m^2
#define BASE_COST_DRIVES (71.) //US$/m^2
#define BASE_COST_PEDESTALS (17.)
#define BASE_COST_CONTROLS (27.)
#define BASE_COST_WIRING (18.)
#define BASE_COST_MANUFACTURING (54.)
#define BASE_COST_INSTALLATION (11.)
#define K_H_MIRROR_MODULES (0.000827445)
#define K_F_MIRROR_MODULES (-7.80008*pow(10.,-7.))
#define K_H_DRIVES (-0.002639289)
#define K_F_DRIVES (0.)
#define K_H_PEDESTALS (0.006816719)
#define K_F_PEDESTALS (4.62604*pow(10.,-7.))
#define K_H_CONTROLS (-0.010308435)
#define K_F_CONTROLS (-9.307*pow(10.,-7.))
#define K_H_WIRING (-0.00687229)
#define K_F_WIRING (2.33234*pow(10.,-7.))
#define K_H_MANUFACTURING (-0.001545098)
#define K_F_MANUFACTURING (-1.08627*pow(10.,-6.))
#define K_H_INSTALLATION (-0.008572889)
#define K_F_INSTALLATION (1.37257*pow(10.,-6.))
#define COST_OF_TOWER_CONSTANT (3.*pow(10.,7.))
#define COST_OF_TOWER_LIN (-285868.)
#define COST_OF_TOWER_QUAD (1835.7)
#define CERAMIC_FIBER_INSULATION_COST (4960.60) //$/m^3 (using conversion to USD in 2009 and USD inflation since 2009 (around 11.5%)
#define STORAGE_TANK_FOUNDATION_COST_COEF (50.869/1000.) //$/kg of ms
#define STORAGE_TANK_FOUNDATION_COST_CONST (241696)
#define RECEIVER_REFERENCE_SURFACE (1571.) //m^2
#define RECEIVER_REFERENCE_COST (51000000.)
#define RECEIVER_EXPONENT_COEFFICIENT (0.7)
#define RECEIVER_MTL_COST_ESC_RATE (1.12)
#define POWERBLOCK_REFERENCE_COST (86960000.) //USD
#define POWERBLOCK_REFERENCE_POWER (115000000.) //W
#define POWERBLOCK_EXPONENT_COEF (0.8)
#define POWERBLOCK_MTL_COST_ESC_RATE (1.0)
#define STEAM_GENERATOR_EXPONENT (-1.112)
#define STEAM_GENERATOR_CONSTANT (3.*pow(10.,6.))
#define STEAM_GEN_REF_POWER (115) //MWe turbine gross)
#define STEAM_GEN_REF_EQUIP (7.3e6) //USD
#define STEAM_GEN_REF_PUMP (3.01e6)
#define STEAM_GEN_REF_PIPE (1.2e6)
#define STEAM_GEN_REF_CONTROL (3.9e6)
#define STEAM_GEN_REF_SUPP (9.6e6)
#define STEAM_GEN_SCALE (0.8)

//Properties of stainless steel 316
#define SS316_YIELD_PRESSURE (290e6)

//Siemens Turbines info
//CHECK ALL OUTLET CONDITIONS THEY ARE ERRONEOUS
//STEAM OUTLET OUGHT NOT TO BE ATMOSPHERIC CONDITIONS
//OUGHT TO BE LOW PRESSURE (0.07bar)
#define SST110_PRESSURE (130e5)
#define SST110_TEMPERATURE (530. + 273.)
#define SST110_MIN (0.)
#define SST110_MAX (7e6)
#define SST110_ENTHALPY (3415.992e3) //j/kg
#define SST110_STAGES 1 //info not found
#define SST110_OUTLET_PRESSURE (6894.757)
#define SST110_OUTLET_TEMP (108 + 273)
#define SST110_OUTLET_ENTHALPY (2726.342e3) //j/kg
#define SST110_CONDENSE false
#define SST110_REHEAT false
#define SST120_PRESSURE (130e5)
#define SST120_TEMPERATURE (530. + 273.)
#define SST120_MIN (0.)
#define SST120_MAX (10e6)
#define SST120_ENTHALPY (3415.992e3) //j/kg
#define SST120_STAGES 1 //info not found
#define SST120_OUTLET_PRESSURE (103.25e3)
#define SST120_OUTLET_TEMP (108 + 273)
#define SST120_OUTLET_ENTHALPY (2726.342e3) //j/kg
#define SST120_CONDENSE false
#define SST120_REHEAT false
#define SST300_PRESSURE (120e5)
#define SST300_TEMPERATURE (520. + 273.)
#define SST300_MIN (10e6)
#define SST300_MAX (50e6)
#define SST300_ENTHALPY (3401.04e3)
#define SST300_STAGES 2 //info not found
#define SST300_OUTLET_PRESSURE (103.25e3)
#define SST300_OUTLET_TEMP (119.88 + 273)
#define SST300_OUTLET_ENTHALPY ( 2661.995e3) //j/kg
#define SST300_CONDENSE true
#define SST300_REHEAT false
#define SST400_PRESSURE (140e5)
#define SST400_TEMPERATURE (540. + 273.)
#define SST400_MIN (30e6)
#define SST400_MAX (65e6)
#define SST400_ENTHALPY (3431.768e3)
#define SST400_STAGES 3 //info not found
#define SST400_OUTLET_PRESSURE (103.25e3)
#define SST400_OUTLET_TEMP (119.88 + 273)
#define SST400_OUTLET_ENTHALPY (2659.863e3) //j/kg
#define SST400_CONDENSE true
#define SST400_REHEAT false
#define SST600_PRESSURE (140e5)
#define SST600_TEMPERATURE (540. + 273.)
#define SST600_MIN (5e6)
#define SST600_MAX (100e6)
#define SST600_ENTHALPY (3431.768e3)
#define SST600_STAGES 4 //info not found
#define SST600_OUTLET_PRESSURE (103.25e3)
#define SST600_OUTLET_TEMP (121.78 + 273)
#define SST600_OUTLET_ENTHALPY (2707.978e3) //j/kg
#define SST600_CONDENSE false
#define SST600_REHEAT false
#define SST700_PRESSURE (165e5)
#define SST700_TEMPERATURE (585. + 273.)
#define SST700_MIN (20e6)
#define SST700_MAX (175e6)
#define SST700_ENTHALPY (3527.97e3)
#define SST700_STAGES 5 //info not found
#define SST700_OUTLET_PRESSURE (103.253e3)
#define SST700_OUTLET_TEMP (139.56 + 273)
#define SST700_OUTLET_ENTHALPY (2745.867e3)
#define SST700_CONDENSE false
#define SST700_REHEAT true
#define SST800_PRESSURE (140e5)
#define SST800_TEMPERATURE (540. + 273.)
#define SST800_MIN (50e6)
#define SST800_MAX (150e6)
#define SST800_ENTHALPY (3431.768e3)
#define SST800_STAGES 5 //info not found
#define SST800_OUTLET_PRESSURE (103.25e3)
#define SST800_OUTLET_TEMP (121.78 + 273)
#define SST800_OUTLET_ENTHALPY (2707.978e3) //j/kg
#define SST800_CONDENSE false
#define SST800_REHEAT true
#define SST900_PRESSURE (165e5)
#define SST900_TEMPERATURE (585. + 273.)
#define SST900_MIN (49e6)
#define SST900_MAX (250e6)
#define SST900_ENTHALPY (3527.97e3)
#define SST900_STAGES 6 //info not found
#define SST900_OUTLET_PRESSURE (103.253e3)
#define SST900_OUTLET_TEMP (139.56 + 273)
#define SST900_OUTLET_ENTHALPY (2745.867e3)
#define SST900_CONDENSE false
#define SST900_REHEAT true


#endif
