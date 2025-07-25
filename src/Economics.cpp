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
#include "Economics.hpp"

/*--------------------------------------------------------*/
/*                       constructor                      */
/*--------------------------------------------------------*/
Economics::Economics ( int    numberOfHeliostats          ,
		       double hotStorageInsul             ,
		       double coldStorageInsul            ,
		       double hotStorageHeight            ,
		       double coldStorageHeight           , // new in V2 (P.B., 2025-07)
		       double hotStorageDiameter          ,
		       double coldStorageDiameter         , // new in V2 (P.B., SLD, 2025-07-24)
		       double heightOfTower               ,
		       double lengthOfHeliostats          ,
		       double widthOfHeliostats           ,
		       double turbinePower                ,
		       double heightOfReceiverAperture    ,
		       double widthOfAperture             ,
		       int    receiverNbOfTubes           ,
		       double receiverTubesOutsideDiam    ,
		       int    exchangerModel              , // new in V2
		       double exchangerTubesDiameter      ,
		       double exchangerTubesLength        , 
		       int    exchangerNumberOfTubes      ,
		       int    exchangerTubePassesPerShell ,
		       int    exchangerNumberOfShell        ) :
  _nbOfHeliostats                 ( numberOfHeliostats          ) , 
  _hotStorageInsulationThickness  ( hotStorageInsul             ) ,
  _coldStorageInsulationThickness ( coldStorageInsul            ) ,
  _hotStorageHeight               ( hotStorageHeight            ) ,
  _coldStorageHeight              ( coldStorageHeight           ) , // (P.B., 2025-07)
  _hotStorageDiameter             ( hotStorageDiameter          ) ,
  _coldStorageDiameter            ( coldStorageDiameter         ) , // (P.B., SLD, 2025-07-24)
  _heightOfTower                  ( heightOfTower               ) ,
  _lengthOfHeliostats             ( lengthOfHeliostats          ) ,
  _widthOfHeliostats              ( widthOfHeliostats           ) ,
  _turbineNominalPowerOutput      ( turbinePower                ) ,
  _heightOfReceiverAperture       ( heightOfReceiverAperture    ) ,
  _widthOfReceiverAperture        ( widthOfAperture             ) ,
  _receiverNumberOfTubes          ( receiverNbOfTubes           ) ,
  _receiverTubesDout              ( receiverTubesOutsideDiam    ) , 
  _exchangerModel                 ( exchangerModel              ) ,
  _exchangerTubesOutterDiameter   ( exchangerTubesDiameter      ) ,
  _exchangerTubesLength           ( exchangerTubesLength        ) ,
  _exchangerNumberOfTubes         ( exchangerNumberOfTubes      ) ,
  _exchangerTubePassesPerShell    ( exchangerTubePassesPerShell ) ,
  _exchangerNumberOfShell         ( exchangerNumberOfShell      )   {

  _costOfField           = 0.0;
  _costPerHeliostat      = 0.0;
  _costOfTower           = 0.0;
  _costOfStorage         = 0.0;
  _costOfPowerblock      = 0.0;
  _costOfSteamGenerator  = 0.0;
  _costOfReceiver        = 0.0;
  _totalCost             = 0.0;


  // below are some initializations from a previous version
  // (we keep them here just in case):
  
  // module costs:
  // _costOfField          = 0.0;
  // _costPerHeliostat     = 0.0;
  // _costOfTower          = 0.0;
  // _costOfStorage        = 0.0;
  // _costOfPowerblock     = 0.0;
  // _costOfSteamGenerator = 0.0;
  // _costOfReceiver       = 0.0;
  // _totalCost            = 0.0;

  // // parameters:
  // _hotStorageInsulationThickness  = 0.01;
  // _coldStorageInsulationThickness = 0.01;
  // _hotStorageHeight               = 1.0;
  // _coldStorageHeight              = 1.0;  // new attribute to compute the cost of storage without bugs (P.B., 2025-07)
  // _hotStorageDiameter             = 1.0;
  // _coldStorageDiameter            = 1.0;  // new in V2 (P.B., SLD, 2025-07-24) 
  // _receiverInsulationThickness    = 0.01;
  // _heightOfTower                  = 20.0;
  // _heightOfReceiverAperture       = 1.0;
  // _widthOfReceiverAperture        = 1.0;
  // _receiverNumberOfTubes          = 0;
  // _receiverTubesDout              = 0.0;
  // _lengthOfHeliostats             = 1.0;
  // _widthOfHeliostats              = 1.0;
  // _nbOfHeliostats                 = 1;
  // _reflectiveArea                 = 0.0;
  // _totalMoltenSaltMass            = 0.0;
  // _turbineNominalPowerOutput      = 0.0;
  // _exchangerModel                 = 1;
  // _exchangerTubesOutterDiameter   = 0.0005;
  // _exchangerTubesLength           = 1.0;
  // _exchangerNumberOfTubes         = 1;
  // _exchangerTubePassesPerShell    = 1;
  // _exchangerNumberOfShell         = 1;
}

/*--------------------------------------------------------*/
double Economics::evaluateCostOfHeliostat ( void ) {
/*--------------------------------------------------------*/

  double A_h = _widthOfHeliostats*_lengthOfHeliostats; //reflective area of 1 heliostat
  double A_f = _reflectiveArea;
  double X   = A_h - A_h_ref;
  double Y   = A_f - A_f_ref;
  double mirrorCostPerM2        = C_M*exp(M_Kh*X + M_Kf*Y);
  double drivesCostPerM2        = C_D*exp(C_Kh*X + C_Kf*Y);
  double pedestalsCostPerM2     = C_P*exp(P_Kh*X + P_Kf*Y);
  double controlsCostPerM2      = C_C*exp(C_Kh*X + C_Kf*Y);
  double wiringCostPerM2        = C_W*exp(W_Kh*X + W_Kf*Y);
  double manufacturingCostPerM2 = C_Ma*exp(Ma_Kh*X + Ma_Kf*Y);
  double installationCostPerM2  = C_I*exp(I_Kh*X + I_Kf*Y);

  _costPerHeliostat = A_h * ( mirrorCostPerM2        +
			      drivesCostPerM2        +
			      pedestalsCostPerM2     +
			      controlsCostPerM2      +
			      wiringCostPerM2        +
			      manufacturingCostPerM2 +
			      installationCostPerM2    );

  return _costPerHeliostat;
}

/*--------------------------------------------------------*/
double Economics::evaluateCostOfField ( void ) {
/*--------------------------------------------------------*/
  evaluateCostOfHeliostat();
  _costOfField = _nbOfHeliostats * _costPerHeliostat;
  return _costOfField;
}

/*--------------------------------------------------------*/
double Economics::evaluateCostOfTower ( void) {
/*--------------------------------------------------------*/
  _costOfTower = COST_OF_TOWER_CONSTANT + COST_OF_TOWER_LIN*_heightOfTower + COST_OF_TOWER_QUAD*pow(_heightOfTower, 2.0);
  return _costOfTower;
}

/*--------------------------------------------------------*/
double Economics::evaluateCostOfReceiver ( void ) {
/*--------------------------------------------------------*/
  double A_r; //absorbing surface area
  if (_receiverNumberOfTubes != 0 && _receiverTubesDout != 0) {
    int numberOfPasses = static_cast<int>(floor((PI*_widthOfReceiverAperture / 2.0) /
					    (_receiverNumberOfTubes*_receiverTubesDout)));
    A_r = PI*_receiverTubesDout*_heightOfReceiverAperture*_receiverNumberOfTubes*numberOfPasses;
  }
  else
    A_r = PI*_widthOfReceiverAperture*_heightOfReceiverAperture / 2.0;

  _costOfReceiver = 1.0*RECEIVER_REFERENCE_COST*RECEIVER_MTL_COST_ESC_RATE *
    pow(A_r / (1.0*RECEIVER_REFERENCE_SURFACE), 1.0*RECEIVER_EXPONENT_COEFFICIENT);

  return _costOfReceiver;
}

/*--------------------------------------------------------*/
double Economics::evaluateCostOfStorage ( void ) {
/*--------------------------------------------------------*/
  
  double moltenSaltPerKg = COST_NANO3_KNO3; // $/kg
  
  // Total molten salt inventory is assumed to be that of the volume of the full cold tank

  // double moltenSaltVolume  = _hotStorageHeight*1.1*PI*pow(_hotStorageDiameter / 2.0, 2.0); // OLD VERSION (v1)
  double moltenSaltVolume     = _coldStorageHeight   *PI*pow(_hotStorageDiameter / 2.0, 2.0); // NEW VERSION (v2, P.B., SLD, 2025-07-23)
  double totalMoltenSaltMass = moltenSaltVolume*MS_DENSITY;
  double moltenSaltCost       = moltenSaltPerKg * totalMoltenSaltMass;

  // OLD VERSION (v1):
  // double insulationVolume = PI*_hotStorageHeight*(pow(_hotStorageDiameter / 2.0 +
  // 						      _hotStorageInsulationThickness + 0.04, 2.0) -
  // 						  pow(_hotStorageDiameter / 2.0 + 0.04, 2.0)) +
  //   PI*_hotStorageHeight*1.1*(pow(_hotStorageDiameter / 2.0 + _coldStorageInsulationThickness + 0.04, 2.0) -
  // 			      pow(_hotStorageDiameter / 2.0 + 0.04, 2.0)) +
  //   PI*pow(_hotStorageDiameter / 2.0 + 0.04, 2.0)*(_coldStorageInsulationThickness + _hotStorageInsulationThickness);

  // UPDATED VERSION (P.B., 2025-07): Replaced _hotStorageHeight*1.1 by the new attribute _coldStorageHeight:
  // double insulationVolume = PI*_hotStorageHeight *(pow(_hotStorageDiameter / 2.0 + _hotStorageInsulationThickness  + 0.04, 2.0) - pow(_hotStorageDiameter / 2.0 + 0.04, 2.0)) +
  //                           PI*_coldStorageHeight*(pow(_hotStorageDiameter / 2.0 + _coldStorageInsulationThickness + 0.04, 2.0) - pow(_hotStorageDiameter / 2.0 + 0.04, 2.0)) +
  //                           PI*pow(_hotStorageDiameter / 2.0 + 0.04, 2.0)*(_coldStorageInsulationThickness + _hotStorageInsulationThickness);
  
  // VERSION for V2 (SLD and P.B., 2025-07-25): We use the other new attribute _coldStorageDiameter:
  double insulationVolume
    = PI*_hotStorageHeight *(pow(_hotStorageDiameter  / 2.0 + _hotStorageInsulationThickness  + 0.04, 2.0) - pow(_hotStorageDiameter  / 2.0 + 0.04, 2.0))
    + PI*_coldStorageHeight*(pow(_coldStorageDiameter / 2.0 + _coldStorageInsulationThickness + 0.04, 2.0) - pow(_coldStorageDiameter / 2.0 + 0.04, 2.0))
    + PI*pow(_hotStorageDiameter  / 2.0 + 0.04, 2.0) * _hotStorageInsulationThickness
    + PI*pow(_coldStorageDiameter / 2.0 + 0.04, 2.0) * _coldStorageInsulationThickness;  
  
  // The 4cm thick stainless steel tank is considered when evaluating the total volume of insulation needed.
  // The design tank diameter is the inner diameter

  double ceramicFiberCost = CERAMIC_FIBER_INSULATION_COST * insulationVolume;
	
  double foundationCost   = totalMoltenSaltMass * STORAGE_TANK_FOUNDATION_COST_COEF + 1.*STORAGE_TANK_FOUNDATION_COST_CONST;
  
  _costOfStorage = ceramicFiberCost + moltenSaltCost + foundationCost;

  return _costOfStorage;
}

/*--------------------------------------------------------*/
double Economics::evaluateCostOfSteamGenerator ( void ) {
/*--------------------------------------------------------*/

  double areaOfAShell, priceOfASquareMeter;
  double equip, pump, pipe, control, supp;

  //Cost is provided in table 5 of the Roadmap report for cost reduction
  //We use the Roadmap baseline cost $/kWe for consistency

  if ( _exchangerModel == 2 ) {
    areaOfAShell = _exchangerTubesLength*PI*2.0*(_exchangerTubesOutterDiameter / 2.0) *
      _exchangerNumberOfTubes*_exchangerTubePassesPerShell;
    priceOfASquareMeter   = 1.0*STEAM_GENERATOR_CONSTANT*pow(areaOfAShell, 1.0*STEAM_GENERATOR_EXPONENT);
    _costOfSteamGenerator = areaOfAShell*priceOfASquareMeter*_exchangerNumberOfShell;
  }
  else {
    equip   = 1.0*STEAM_GEN_REF_EQUIP   * pow((_turbineNominalPowerOutput*(1e-6)) / (1.0*STEAM_GEN_REF_POWER), 1.0*STEAM_GEN_SCALE);
    pump    = 1.0*STEAM_GEN_REF_PUMP    * pow((_turbineNominalPowerOutput*(1e-6)) / (1.0*STEAM_GEN_REF_POWER), 1.0*STEAM_GEN_SCALE);
    pipe    = 1.0*STEAM_GEN_REF_PIPE    * pow((_turbineNominalPowerOutput*(1e-6)) / (1.0*STEAM_GEN_REF_POWER), 1.0*STEAM_GEN_SCALE);
    control = 1.0*STEAM_GEN_REF_CONTROL * pow((_turbineNominalPowerOutput*(1e-6)) / (1.0*STEAM_GEN_REF_POWER), 1.0*STEAM_GEN_SCALE);
    supp    = 1.0*STEAM_GEN_REF_SUPP    * pow((_turbineNominalPowerOutput*(1e-6)) / (1.0*STEAM_GEN_REF_POWER), 1.0*STEAM_GEN_SCALE);
    
    _costOfSteamGenerator = equip + pump + pipe + control + supp;
  }
  return _costOfSteamGenerator;
}

/*--------------------------------------------------------*/
double Economics::evaluateCostOfPowerblock ( void ) {
/*--------------------------------------------------------*/

  //cost is provided in table
  _costOfPowerblock = (1.0*POWERBLOCK_REFERENCE_COST*POWERBLOCK_MTL_COST_ESC_RATE) *
    pow(_turbineNominalPowerOutput / (1.0*POWERBLOCK_REFERENCE_POWER), 1.0*POWERBLOCK_EXPONENT_COEF);
  return _costOfPowerblock;
}

/*--------------------------------------------------------*/
double Economics::evaluateTotalInvestmentCost ( void ) {
/*--------------------------------------------------------*/
  _totalCost = _costOfField + _costOfTower + _costOfStorage +
    _costOfPowerblock + _costOfSteamGenerator + _costOfReceiver;
  return _totalCost;
}
