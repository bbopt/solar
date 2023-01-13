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
#include "HtfCycle.hpp"

/*----------------------------------------------------------------*/
/*                          constructor #1                        */
/*----------------------------------------------------------------*/
HtfCycle::HtfCycle ( double       receiverTemp           ,
		     double       storageHeight          ,
		     double       storageDiameter        ,
		     double       insulationThickness    ,
		     double       exchangerExitTemp      ,
		     Powerblock * powerblock             ,
		     double       apertureHeight         ,
		     double       apertureWidth          ,
		     double       receiverTubesDin       ,
		     double       receiverTubesThickness ,
		     int          receiverNbTubes        ,
		     int          timeInterval             )  :

  // Mechanical components:
  _centralReceiver ( &_centralReceiverInlet ,
		     &_centralReceiverOutlet,
		     apertureHeight         ,
		     apertureWidth          ,
		     insulationThickness    ,
		     receiverTubesDin       ,
		     receiverTubesThickness ,
		     receiverNbTubes          ) ,

  _hotStorage  ( &_centralReceiverOutlet ,
		 &_steamGeneratorInlet   ,
		 storageHeight           ,
		 storageDiameter         ,
		 insulationThickness       ) ,

  _coldStorage ( &_steamGeneratorOutlet ,
		 &_centralReceiverInlet ,
		 1.2*storageHeight      ,
		 storageDiameter        ,
		 insulationThickness      ) ,

  _steamGenerator ( &_steamGeneratorInlet, &_steamGeneratorOutlet, powerblock ) ,

  // MoltenSalt conditions:
  _centralReceiverInlet  ( exchangerExitTemp, P_ATM, 0.0 ) ,
  _centralReceiverOutlet ( receiverTemp     , P_ATM, 0.0 ) ,
  _steamGeneratorInlet   ( receiverTemp     , P_ATM, 0.0 ) ,
  _steamGeneratorOutlet  ( exchangerExitTemp, P_ATM, 0.0 ) , 

  _timeInterval       ( timeInterval )   {
  
  _steamGenOutletMsRate.reserve ( 1440 * _timeInterval ); // 1440=24x60
  _steamGenOutletTemp.reserve   ( 1440 * _timeInterval );
  _storageHeat.reserve          ( 1440 * _timeInterval );
  
  _minColdStorageTemp    = exchangerExitTemp;
  _minHotStorageTemp     = receiverTemp;
  _minSteamGenOutletTemp = exchangerExitTemp;
  }

/*----------------------------------------------------------------*/
/*                          constructor #2                        */
/*----------------------------------------------------------------*/
HtfCycle::HtfCycle ( double       receiverTemp                ,
		     double       storageHeight               ,
		     double       storageDiameter             ,
		     double       insulationThickness         ,
		     double       exchangerExitTemp           ,
		     Powerblock * powerblock                  ,
		     double       apertureHeight              ,
		     double       apertureWidth               ,
		     double       receiverTubesDin            ,
		     double       receiverTubesThickness      ,
		     int          receiverNbTubes             ,
		     int          timeInterval                ,
		     double       exchangerTubesDin           ,
		     double       exchangerTubesDout          ,
		     double       exchangerTubesLength        ,
		     double       exchangerTubesSpacing       ,
		     double       baffleCut                   ,
		     int          nbOfBaffles                 ,
		     int          exchangerNbOfTubes          ,
		     int          exchangerNbOfPassesPerShell ,
		     int          exchangerNbOfShells           ) :

  // Mechanical components:
  _centralReceiver ( &_centralReceiverInlet ,
		     &_centralReceiverOutlet,
		     apertureHeight         ,
		     apertureWidth          ,
		     insulationThickness    ,
		     receiverTubesDin       ,
		     receiverTubesThickness ,
		     receiverNbTubes          ) ,

  _hotStorage ( &_centralReceiverOutlet,
		&_steamGeneratorInlet  ,
		storageHeight          ,
		storageDiameter        ,
		insulationThickness      ) ,

  _coldStorage ( &_steamGeneratorOutlet,
		 &_centralReceiverInlet,
		 1.2*storageHeight     ,
		 storageDiameter       ,
		 insulationThickness     ) ,

  _steamGenerator ( &_steamGeneratorInlet      ,
		    &_steamGeneratorOutlet     ,
		    powerblock                 ,
		    exchangerTubesLength       ,
		    exchangerTubesDin          ,
		    exchangerTubesDout         ,
		    exchangerTubesSpacing      ,
		    baffleCut                  ,
		    nbOfBaffles                ,
		    exchangerNbOfTubes         ,
		    exchangerNbOfPassesPerShell,
		    exchangerNbOfShells          ) ,

  // MoltenSalt conditions:
  _centralReceiverInlet  ( exchangerExitTemp, P_ATM, 0.0 ) ,
  _centralReceiverOutlet ( receiverTemp     , P_ATM, 0.0 ) ,
  _steamGeneratorInlet   ( receiverTemp     , P_ATM, 0.0 ) ,
  _steamGeneratorOutlet  ( exchangerExitTemp, P_ATM, 0.0 ) ,

  _timeInterval       ( timeInterval )  {
 
  _steamGenOutletMsRate.reserve ( 1440 * _timeInterval ); // 1440=24x60
  _steamGenOutletTemp.reserve   ( 1440 * _timeInterval );
  _storageHeat.reserve          ( 1440 * _timeInterval );
  _minColdStorageTemp    = exchangerExitTemp;
  _minHotStorageTemp     = receiverTemp;
  _minSteamGenOutletTemp = exchangerExitTemp;
}

/*----------------------------------------------------------------------*/
/*  This function generates one full iteration of the HTF cycle for     */
/*  a whole time interval for a given energy input. The cycle goes as:  */
/*  1- Determine receiver's outlet properties (mass flow/temperature)   */
/*  2- Determine inflow requirements for the steam generator            */
/*  3- Determine hot storage level and temperature at interval end      */
/*  4- Determine cold storage level and temperature at interval end     */
/*----------------------------------------------------------------------*/
void HtfCycle::fOperateCycle ( int    timeInSeconds       ,
			       double energyFromField     ,
			       double requiredPowerOutput ,
			       bool   low_fid               ) {
  
  double hotMass_i, hotMass_f, hotMass_avg;
  double hotTemp_i, hotTemp_f, hotTemp_avg, hotTemp_avg2;
  double hotLevel_i, hotLevel_f, hotLevel_avg;
  double U_i, U_f; // U_avg;
  double Q_dot_i, Q_dot_avg; // Q_dot_f
  double availableMass;
  double hotMassFlow_in, hotMassFlow_out, hotMassFlow_out2;
  double hotTemp_in, hotTemp_out;
  double heatedMass;
  double temperaturePrecision = 0.01 / timeInSeconds;
  size_t cas, count;
  int steamGeneratorModel = _steamGenerator.get_exchangerModel();
  double Cp = HEAT_CAPACITY;
  double minimumHotMass = 0.05*MS_DENSITY*_hotStorage.get_heightOfStorage()
    *PI*pow(_hotStorage.get_diameterOfStorage() / 2.0, 2.0);

  hotMassFlow_in = _centralReceiver.computeEnergyToFluid ( energyFromField );
  hotMass_i      = _hotStorage.get_storedMass();
  
  Q_heat_hot     = 0.0;
  Q_heat_cold    = 0.0;
  
  if      ( hotMassFlow_in == 0.0 && hotMass_i == 0.0 ) { cas = 1; }
  else if ( hotMassFlow_in == 0.0 && hotMass_i  > 0.0 ) { cas = 2; }
  else if ( hotMassFlow_in  > 0.0 && hotMass_i == 0.0 ) { cas = 3; }
  else if ( hotMassFlow_in  > 0.0 && hotMass_i  > 0.0 ) { cas = 4; }
  else if ( hotMassFlow_in == 0.0 && hotMass_i  < 0.0 ) { hotMass_i = 0.0; cas = 1; }
  else if ( hotMassFlow_in  > 0.0 && hotMass_i  < 0.0 ) { hotMass_i = 0.0; cas = 3; }
  else {
    std::ostringstream oss;
    oss << "error in the full iteration of the HTF cycle: hotMassFlow_in = "
	<< hotMassFlow_in << " hotMass_i = " << hotMass_i;
    throw Simulation_Interruption ( oss.str() );
  }
  
  switch ( cas ) {
      
  case 1:
    hotTemp_i       = 0.0;
    hotLevel_i      = 0.0;
    hotMassFlow_out = 0.0;
    availableMass   = 0.0;
    heatedMass      = 0.0;
    break;
      
  case 2:
    heatedMass    = 0.0;
    availableMass = hotMass_i;
    hotTemp_i     = _hotStorage.get_storedTemperature();
    if (steamGeneratorModel == 1)
      hotMassFlow_out = _steamGenerator.fComputeRequiredMoltenSaltMassFlow ( requiredPowerOutput, hotTemp_i);
    else
      hotMassFlow_out
	= _steamGenerator.fComputeRequiredMoltenSaltMassFlow (requiredPowerOutput, hotTemp_i, availableMass / timeInSeconds);
    hotLevel_i = _hotStorage.fComputeStorageLevel(hotMass_i);
    U_i = hotMass_i* Cp * hotTemp_i;
    Q_dot_i = _hotStorage.fComputeEnergyLosses(hotTemp_i, hotLevel_i);
    hotTemp_out = hotTemp_i;
    hotMassFlow_out2 = 0.0;
      
    if ( availableMass < timeInSeconds*hotMassFlow_out ||
	 availableMass - timeInSeconds*hotMassFlow_out < minimumHotMass ||
	 hotMassFlow_out == 0.0 ) {
      hotMass_f = hotMass_i;
      hotLevel_f = hotLevel_i;
      hotMassFlow_out  = 0.0;
      hotMassFlow_out2 = 0.0;
      U_f = U_i - timeInSeconds*Q_dot_i;
      hotTemp_f = U_f / (Cp*hotMass_f);
      if (hotTemp_f < MELTING_POINT) {
	Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	hotTemp_f = MELTING_POINT;
      }
      hotTemp_avg = (hotTemp_i + hotTemp_f) / 2.0;
      hotTemp_avg2 = 0.0;
      count = 0;

      if ( !low_fid ) {
	while ( fabs(hotTemp_avg2 - hotTemp_avg) > 0.001 && count < 15 ) {
	  hotTemp_avg2 = hotTemp_avg;
	  Q_dot_avg = _hotStorage.fComputeEnergyLosses(hotTemp_avg, hotLevel_i);
	  U_f = U_i - timeInSeconds*Q_dot_avg;
	  hotTemp_f = U_f / (Cp * hotMass_f);
	  if (hotTemp_f < MELTING_POINT) {
	    Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	    hotTemp_f = MELTING_POINT;
	  }
	  hotTemp_avg = (hotTemp_i + hotTemp_f) / 2.0;
	  ++count;
	}
      }
    }
    else {
      hotMass_f    = hotMass_i - timeInSeconds * hotMassFlow_out;
      hotMass_avg  = (hotMass_i + hotMass_f) / 2.0;
      hotLevel_avg = _hotStorage.fComputeStorageLevel(hotMass_avg);
      U_f = U_i - timeInSeconds*(hotMassFlow_out*Cp*hotTemp_out - Q_dot_i);
      hotTemp_f = U_f / (Cp * hotMass_f);
      if (hotTemp_f < MELTING_POINT) {
	Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	hotTemp_f = MELTING_POINT;
      }
      hotTemp_avg = (hotTemp_i + hotTemp_f) / 2.0;
      hotTemp_avg2 = 0.0;
      if ( steamGeneratorModel == 1 )
	hotMassFlow_out = _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_avg);
      else
	hotMassFlow_out
	  = _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_avg, availableMass / timeInSeconds);

      count = 0;
      
      if ( !low_fid ) {
	while ((fabs(hotTemp_avg2 - hotTemp_avg) > temperaturePrecision ||
		fabs(hotMassFlow_out2 - hotMassFlow_out) / hotMassFlow_out > 0.001) &&
	       count < 15) {
	  hotTemp_avg2 = hotTemp_avg;
	  hotMassFlow_out2 = hotMassFlow_out;
	  hotMass_f = hotMass_i - timeInSeconds*hotMassFlow_out;
	  hotMass_avg = (hotMass_f + hotMass_i) / 2.0;
	  hotLevel_avg = _hotStorage.fComputeStorageLevel(hotMass_avg);
	  Q_dot_avg = _hotStorage.fComputeEnergyLosses(hotTemp_avg, hotLevel_avg);
	  U_f = U_i - timeInSeconds*(Cp*hotMassFlow_out*hotTemp_avg + Q_dot_avg);
	  hotTemp_f = U_f / (Cp*hotMass_f);
	  if (hotTemp_f < MELTING_POINT) {
	    Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	    hotTemp_f = MELTING_POINT;
	  }
	  hotTemp_avg = (hotTemp_i + hotTemp_f) / 2.0;
	  if ( steamGeneratorModel == 1 )
	    hotMassFlow_out = _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_avg);
	  else
	    hotMassFlow_out =
	      _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_avg, availableMass / timeInSeconds);
	  ++count;
	}
      }

      // Removed 2022-12-15:
      // hotMassFlow_out*timeInSeconds > availableMass ? hotMassFlow_out = 0.0 : 0;

      // and replaced by:
      if ( hotMassFlow_out*timeInSeconds > availableMass )
	hotMassFlow_out = 0.0;

      hotMass_f = hotMass_i - timeInSeconds*hotMassFlow_out;
  }
      
  // Set final conditions for the flows and storage
  _hotStorage.set_storage(hotMass_f, hotTemp_f);
  _steamGeneratorInlet.set_massFlow(hotMassFlow_out);
  _steamGeneratorInlet.set_temperature(hotTemp_avg);
  _steamGeneratorOutlet.set_massFlow(hotMassFlow_out);
      
  break;
      
  case 3:
    hotTemp_in = _centralReceiverOutlet.get_temperature();
    hotTemp_i = hotTemp_in;
    hotLevel_i = 0.0;
    hotMass_i = 0.0;
    availableMass = hotMassFlow_in*timeInSeconds - minimumHotMass;
    if (steamGeneratorModel == 1)
      hotMassFlow_out = _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_i);
    else
      hotMassFlow_out =
	_steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_i, availableMass / timeInSeconds);     
    hotMass_f = hotMass_i + timeInSeconds*(hotMassFlow_in - hotMassFlow_out);
    hotLevel_f = _hotStorage.fComputeStorageLevel(hotMass_f);
      
    if ( timeInSeconds*hotMassFlow_in < minimumHotMass ) {
      hotMass_f = 0.0;
      Q_dot_avg = 0.0;
      hotMassFlow_out = 0.0;
      hotMass_f = timeInSeconds*hotMassFlow_in;
	 
      //thermal losses are neglected for these intervals
      //Setting final conditions
      _hotStorage.set_storage(hotMass_f, hotTemp_i);
      _centralReceiverInlet.set_massFlow(hotMassFlow_in);
      _centralReceiverOutlet.set_massFlow(hotMassFlow_in);
      _steamGeneratorInlet.set_massFlow(hotMassFlow_out);
      _steamGeneratorOutlet.set_massFlow(hotMassFlow_out);
      heatedMass = timeInSeconds*hotMassFlow_in;
      break;
    }
    else if ( hotLevel_f <= _hotStorage.get_heightOfStorage() && timeInSeconds*hotMassFlow_in > minimumHotMass ) {
      availableMass = timeInSeconds*hotMassFlow_in - minimumHotMass;
      if (timeInSeconds*hotMassFlow_out > availableMass)
	hotMassFlow_out = 0.0;
      hotMass_f = timeInSeconds*(hotMassFlow_in - hotMassFlow_out);
      hotMass_avg = hotMass_f / 2.0;
      hotLevel_avg = _hotStorage.fComputeStorageLevel(hotMass_avg);
	
      Q_dot_i = _hotStorage.fComputeEnergyLosses(hotTemp_i, hotLevel_avg);
      U_f = timeInSeconds*(hotMassFlow_in*Cp*hotTemp_in - hotMassFlow_out*Cp*hotTemp_in - Q_dot_i);
      hotTemp_f = U_f / (Cp*hotMass_f);
      if ( hotTemp_f < MELTING_POINT ) {
	Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	hotTemp_f = MELTING_POINT;
      }
      hotTemp_avg = (hotTemp_i + hotTemp_f) / 2.0;
      hotTemp_avg2 = 0.0;
      count = 0;

      if ( !low_fid ) {
	while ( fabs(hotTemp_avg2 - hotTemp_avg) > temperaturePrecision
		&& count < 15 ) {
	  hotTemp_avg2 = hotTemp_avg;
	    
	  Q_dot_avg = _hotStorage.fComputeEnergyLosses(hotTemp_avg, hotLevel_avg);
	  
	  U_f = timeInSeconds*(hotMassFlow_in*Cp*hotTemp_in - hotMassFlow_out*Cp*hotTemp_avg - Q_dot_avg);
	  hotTemp_f = U_f / (Cp*hotMass_f);
	  if (hotTemp_f < MELTING_POINT) {
	    Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	    hotTemp_f = MELTING_POINT;
	  }
	  hotTemp_avg = (hotTemp_i - hotTemp_f) / 2.0;
	  ++count;
	}
      }
	
      //setting final conditions
      _hotStorage.set_storage(hotMass_f, hotTemp_f);
      _steamGeneratorInlet.set_massFlow(hotMassFlow_out);
      _steamGeneratorInlet.set_temperature(hotTemp_avg);
      _steamGeneratorOutlet.set_massFlow(hotMassFlow_out);
	heatedMass = timeInSeconds*hotMassFlow_in;
	
    }
    else if (hotLevel_f > _hotStorage.get_heightOfStorage() && hotMassFlow_out == 0.0) {
      hotLevel_f = _hotStorage.get_heightOfStorage();
      hotMass_f = MS_DENSITY*hotLevel_f * PI*pow(_hotStorage.get_diameterOfStorage(), 2.0) / 4.0;
      U_i = hotMass_f * hotTemp_i * Cp;
      Q_dot_avg = _hotStorage.fComputeEnergyLosses(hotTemp_i, hotLevel_f);
      U_f = U_i - timeInSeconds*Q_dot_avg;
      hotTemp_f = U_f / (Cp*hotMass_f);
      if (hotTemp_f < MELTING_POINT) {
	Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	hotTemp_f = MELTING_POINT;
      }
      
      hotMassFlow_in = 0.0;
      
      //Setting final conditions
      heatedMass = hotMass_f;
      _hotStorage.set_storage(hotMass_f, hotTemp_f);
      _steamGeneratorInlet.set_massFlow(0.0);
      _steamGeneratorOutlet.set_massFlow(0.0);
      break;
    }
    else if (hotLevel_f > _hotStorage.get_heightOfStorage()) {
      hotLevel_f = _hotStorage.get_heightOfStorage();
      hotMass_f = MS_DENSITY*hotLevel_f * PI*pow(_hotStorage.get_diameterOfStorage(), 2.0) / 4.0;
      hotTemp_f = hotTemp_in;
      hotMassFlow_in = hotMassFlow_out;
      
      //Setting final conditions
      heatedMass = hotMass_f + timeInSeconds*hotMassFlow_out;
      _hotStorage.set_storage(hotMass_f, hotTemp_f);
      _steamGeneratorInlet.set_massFlow(hotMassFlow_out);
      _steamGeneratorOutlet.set_massFlow(hotMassFlow_out);
      _steamGeneratorInlet.set_temperature(hotTemp_in);
      break;
    }
    else if (hotMassFlow_out == 0.0) {
      hotMassFlow_out = 0.0;
      hotMass_f = hotMassFlow_in*timeInSeconds;
      hotMass_avg = hotMass_f / 2.0;
      hotLevel_f = _hotStorage.fComputeStorageLevel(hotMass_f);
      hotLevel_avg = _hotStorage.fComputeStorageLevel(hotMass_avg);
      
      U_i = 0.0;
      Q_dot_i = _hotStorage.fComputeEnergyLosses(hotTemp_i, hotLevel_avg);
      U_f = timeInSeconds*(hotMassFlow_in*Cp*hotTemp_i - Q_dot_i);
      hotTemp_f = U_f / (Cp*hotMass_f);
      if (hotTemp_f < MELTING_POINT) {
	Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	hotTemp_f = MELTING_POINT;
      }
      
      hotTemp_avg = (hotTemp_i + hotTemp_f) / 2.0;
      hotTemp_avg2 = 0.0;
      count = 0;
      if ( !low_fid ) {
	while (fabs(hotTemp_avg2 - hotTemp_avg) > temperaturePrecision && count < 15) {
	  hotTemp_avg2 = hotTemp_avg;
	  Q_dot_avg = _hotStorage.fComputeEnergyLosses(hotTemp_avg, hotLevel_avg);
	  U_f = timeInSeconds*(hotMassFlow_in*Cp*hotTemp_in - Q_dot_avg);
	  
	  hotTemp_f = U_f / (Cp*hotMass_f);
	  if (hotTemp_f < MELTING_POINT) {
	    Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	    hotTemp_f = MELTING_POINT;
	  }
	  hotTemp_avg = (hotTemp_i + hotTemp_f) / 2.0;
	  ++count;
	}
      }
      heatedMass = hotMassFlow_in*timeInSeconds;
      _hotStorage.set_storage(hotMass_f, hotTemp_f);
      _steamGeneratorInlet.set_massFlow(0.0);
      _steamGeneratorOutlet.set_massFlow(0.0);
      _centralReceiverOutlet.set_massFlow(hotMassFlow_in);
      _centralReceiverInlet.set_massFlow(hotMassFlow_in);
    }
    else {
      hotMass_f = timeInSeconds*(hotMassFlow_in - hotMassFlow_out);
      hotTemp_f = hotTemp_i;
      //losses are neglected for this time frame in this case just like in the case
      //where the tank is full. The temperature of the fluid being outputed is considered
      //equal to the entrance temperature.
      
      //setting final conditions
      heatedMass = timeInSeconds * hotMassFlow_in;
      _hotStorage.set_storage(hotMass_f, hotTemp_f);
      _steamGeneratorInlet.set_massFlow(hotMassFlow_out);
      _steamGeneratorOutlet.set_massFlow(hotMassFlow_out);
      _steamGeneratorInlet.set_temperature(hotTemp_f);
      _centralReceiverOutlet.set_massFlow(hotMassFlow_in);
      _centralReceiverInlet.set_massFlow(hotMassFlow_in);
      break;
    }
      
    break;
      
    case 4:
      hotMass_i = _hotStorage.get_storedMass();
      hotTemp_i = _hotStorage.get_storedTemperature();
      hotTemp_in = _centralReceiverOutlet.get_temperature();
      hotLevel_i = _hotStorage.fComputeStorageLevel(hotMass_i);
      
      availableMass = hotMass_i + timeInSeconds*hotMassFlow_in;
      if (steamGeneratorModel == 1)
	hotMassFlow_out = _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_i);
      else
	hotMassFlow_out =
	  _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_i, availableMass / timeInSeconds);

      if ( availableMass - timeInSeconds*hotMassFlow_out < minimumHotMass ) {
	/*hotMassFlow_out = availableMass / timeInSeconds;*/
	if (availableMass > minimumHotMass) {
	  hotMassFlow_out = 0.0;
	  hotMass_f = availableMass - timeInSeconds*hotMassFlow_out;
	}
	else {
	  hotMassFlow_out = 0.0;
	  hotMass_f = availableMass;
	}

	U_i = hotMass_i * Cp * hotTemp_i;
	hotMass_avg = (hotMass_i + hotMass_f) / 2.0;
	hotLevel_avg = _hotStorage.fComputeStorageLevel(hotMass_avg);
	Q_dot_i = _hotStorage.fComputeEnergyLosses(hotTemp_i, hotLevel_avg);
	U_f = U_i + timeInSeconds*(hotMassFlow_in * Cp * hotTemp_in - hotMassFlow_out * Cp * hotTemp_i - Q_dot_i);
	hotTemp_f = U_f / (Cp*hotMass_f);
	if (hotTemp_f < MELTING_POINT) {
	  Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	  hotTemp_f = MELTING_POINT;
	}

	hotTemp_avg = (hotTemp_i + hotTemp_f) / 2.0;
	hotTemp_avg2 = 0.0;
	count = 0;
	if ( !low_fid ) {
	  while (fabs(hotTemp_avg2 - hotTemp_avg) > temperaturePrecision && count < 15) {
	    hotTemp_avg2 = hotTemp_avg;
	    Q_dot_avg = _hotStorage.fComputeEnergyLosses(hotTemp_avg, hotLevel_avg);
	    
	    U_f = U_i + timeInSeconds*(hotMassFlow_in * Cp * hotTemp_in - hotMassFlow_out * Cp * hotTemp_avg - Q_dot_avg);
	    hotTemp_f = U_f / (Cp*hotMass_f);
	    hotTemp_avg = (hotTemp_i + hotTemp_f) / 2.0;
	    ++count;
	  }
	}
	//in this case, thermal losses will be neglected

	heatedMass = timeInSeconds*hotMassFlow_in;
	_hotStorage.set_storage(hotMass_f, hotTemp_f);
	_centralReceiverOutlet.set_massFlow(hotMassFlow_in);
	_centralReceiverInlet.set_massFlow(hotMassFlow_in);
	_steamGeneratorInlet.set_massFlow(hotMassFlow_out);
	_steamGeneratorInlet.set_temperature((hotMass_i*hotTemp_i + hotMassFlow_in*hotTemp_in*timeInSeconds) / availableMass);
	_steamGeneratorOutlet.set_massFlow(hotMassFlow_out);
	break;
      }

      hotMass_f = hotMass_i + timeInSeconds*(hotMassFlow_in - hotMassFlow_out);
      hotMass_avg = (hotMass_i + hotMass_f) / 2.0;
      hotLevel_f = _hotStorage.fComputeStorageLevel(hotMass_f);
      hotLevel_avg = _hotStorage.fComputeStorageLevel(hotMass_avg);

      if (hotLevel_f > _hotStorage.get_heightOfStorage()) {
	hotLevel_f = _hotStorage.get_heightOfStorage();
	hotMass_f = MS_DENSITY*hotLevel_f*PI*pow(_hotStorage.get_diameterOfStorage(), 2.0) / 4.0;

	hotLevel_avg = hotLevel_f;
	hotMass_avg = hotMass_f;
	
	hotMassFlow_in = (hotMass_f - hotMass_i + hotMassFlow_out*timeInSeconds) / timeInSeconds;
	hotTemp_in = _centralReceiverOutlet.get_temperature();
	hotTemp_f = (hotTemp_i*hotMass_i + hotMassFlow_in*hotTemp_in - hotMassFlow_out*hotTemp_i) / hotMass_f;
	if (hotTemp_f < MELTING_POINT) {
	  Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	  hotTemp_f = MELTING_POINT;
	}
	
	hotTemp_avg = (hotTemp_i + hotTemp_f) / 2.0;
	
	Q_dot_avg = _hotStorage.fComputeEnergyLosses(hotTemp_avg, hotLevel_avg);
	if (steamGeneratorModel == 1)
	  hotMassFlow_out = _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_avg);
	else
	  hotMassFlow_out =
	    _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_avg, availableMass / timeInSeconds);
	hotTemp_avg2 = 0.0;
	hotMassFlow_out2 = 0.0;
	count = 0;
	if ( !low_fid ) {
	  while ((fabs(hotTemp_avg2 - hotTemp_avg) > temperaturePrecision
		  || fabs(hotMassFlow_out2 - hotMassFlow_out) / hotMassFlow_out > 0.001)
		 && count < 15) {
	    hotTemp_avg2 = hotTemp_avg;
	    hotMassFlow_out2 = hotMassFlow_out;
	    
	    hotMassFlow_in = (hotMass_f - hotMass_i + hotMassFlow_out*timeInSeconds) / timeInSeconds;
	    
	    U_i = hotMass_i*Cp*hotTemp_i;
	    U_f = U_i + timeInSeconds*(Cp*hotMassFlow_in*hotTemp_in - Cp*hotMassFlow_out*hotTemp_avg - Q_dot_avg);
	    hotTemp_f = U_f / (Cp*hotMass_f);
	    if (hotTemp_f < MELTING_POINT) {
	      Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	      hotTemp_f = MELTING_POINT;
	    }
	    
	    hotTemp_avg = (hotTemp_i + hotTemp_f) / 2.0;
	    
	    if (steamGeneratorModel == 1)
	      hotMassFlow_out = _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_avg);
	    else
	      hotMassFlow_out =
		_steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_avg, availableMass / timeInSeconds);
	    Q_dot_avg = _hotStorage.fComputeEnergyLosses(hotTemp_avg, hotLevel_avg);
	    ++count;
	  }
	}

	hotMassFlow_out*timeInSeconds > availableMass ? hotMassFlow_out = 0. : 0;
	hotMassFlow_in = (hotMass_f - hotMass_i + hotMassFlow_out*timeInSeconds) / timeInSeconds;

	//setting final conditions
	_hotStorage.set_storage(hotMass_f, hotTemp_f);
	_centralReceiverInlet.set_massFlow(hotMassFlow_in);
	_centralReceiverOutlet.set_massFlow(hotMassFlow_in);
	_steamGeneratorInlet.set_massFlow(hotMassFlow_out);
	_steamGeneratorInlet.set_temperature(hotTemp_avg);
	_steamGeneratorOutlet.set_massFlow(hotMassFlow_out);
	heatedMass = hotMassFlow_in*timeInSeconds;
	break;
      }

      hotTemp_in = _centralReceiverOutlet.get_temperature();
      Q_dot_i = _hotStorage.fComputeEnergyLosses(hotTemp_i, hotLevel_avg);

      U_i = hotMass_i * Cp*hotTemp_i;
      U_f = U_i + timeInSeconds*(Cp*hotMassFlow_in*hotTemp_in - Cp*hotMassFlow_out*hotTemp_i - Q_dot_i);
      hotTemp_f = U_f / (Cp*hotMass_f);
      if (hotTemp_f < MELTING_POINT) {
	Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	hotTemp_f = MELTING_POINT;
      }

      hotTemp_avg = (hotTemp_i + hotTemp_f) / 2.0;
      hotTemp_avg2 = 0.0;

      if (steamGeneratorModel == 1) 
	hotMassFlow_out = _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_avg);
      else
	hotMassFlow_out =
	  _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_avg, availableMass / timeInSeconds);
      hotMassFlow_out2 = 0.0;
      count = 0;
      if ( !low_fid ) {
	while ((fabs(hotTemp_avg2 - hotTemp_avg) > temperaturePrecision
		|| fabs(hotMassFlow_out2 - hotMassFlow_out) > 0.001)
	       && count < 15) {
	  hotTemp_avg2 = hotTemp_avg;
	  hotMassFlow_out2 = hotMassFlow_out;
	  
	  hotMass_f = hotMass_i + timeInSeconds*(hotMassFlow_in - hotMassFlow_out);
	  hotMass_avg = (hotMass_i + hotMass_f) / 2.0;
	  hotLevel_avg = _hotStorage.fComputeStorageLevel(hotMass_avg);
	
	  Q_dot_avg = _hotStorage.fComputeEnergyLosses(hotTemp_avg, hotLevel_avg);
	  
	  U_f = U_i + timeInSeconds*(Cp*hotMassFlow_in*hotTemp_in - Cp*hotMassFlow_out*hotTemp_avg - Q_dot_avg);
	  hotTemp_f = U_f / (Cp*hotMass_f);
	  if (hotTemp_f < MELTING_POINT) {
	    Q_heat_hot = (MELTING_POINT - hotTemp_f) * Cp * hotMass_f;
	    hotTemp_f = MELTING_POINT;
	  }
	  hotTemp_avg = (hotTemp_i + hotTemp_f) / 2.0;
	  
	  if (steamGeneratorModel == 1)
	    hotMassFlow_out = _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_avg);
	  else
	    hotMassFlow_out =
	      _steamGenerator.fComputeRequiredMoltenSaltMassFlow(requiredPowerOutput, hotTemp_avg, availableMass / timeInSeconds);
	  
	  ++count;
	}
      }
      hotMass_f = hotMass_i + timeInSeconds*(hotMassFlow_in - hotMassFlow_out);
      
      //Final conditions
      heatedMass = hotMassFlow_in * timeInSeconds;
      _hotStorage.set_storage(hotMass_f, hotTemp_f);
      _centralReceiverInlet.set_massFlow(hotMassFlow_in);
      _centralReceiverOutlet.set_massFlow(hotMassFlow_in);
      _steamGeneratorInlet.set_massFlow(hotMassFlow_out);
      _steamGeneratorInlet.set_temperature(hotTemp_avg);
      _steamGeneratorOutlet.set_massFlow(hotMassFlow_out);
      break;
  }

  //Now evaluate the cold storage
  double coldTemp_i, coldTemp_f, coldTemp_avg, coldTemp_avg2;
  double coldMass_i, coldMass_f, coldMass_avg;
  double coldLevl_i, coldLevl_avg;
  double coldMassFlow_in, coldMassFlow_out;
  double coldTemp_in;
  double maxColdMass;

  coldMass_i = _coldStorage.get_storedMass();
  coldTemp_i = _coldStorage.get_storedTemperature();
  coldLevl_i = _coldStorage.fComputeStorageLevel(coldMass_i);
  coldMassFlow_in = _steamGeneratorOutlet.get_massFlow();
  coldMassFlow_out = heatedMass / timeInSeconds;
  coldTemp_in = _steamGeneratorOutlet.get_temperature();

  //preliminary approximation
  coldMass_f = coldMass_i + timeInSeconds*(coldMassFlow_in - coldMassFlow_out);

  // coldLevl_f = _coldStorage.fComputeStorageLevel(coldMass_f);
  _coldStorage.fComputeStorageLevel(coldMass_f);
	
  coldMass_avg = (coldMass_i + coldMass_f) / 2.0;
  coldLevl_avg = _coldStorage.fComputeStorageLevel(coldMass_avg);

  Q_dot_i = _coldStorage.fComputeEnergyLosses(coldTemp_i, coldLevl_i);
  if (Q_dot_i < 0.0)
    Q_dot_i = 0.0;
	
  U_i = coldMass_i * Cp * coldTemp_i;
  U_f = U_i + timeInSeconds*(coldMassFlow_in*Cp*coldTemp_in - coldMassFlow_out*Cp*coldTemp_i - Q_dot_i);

  coldTemp_f = U_f / (Cp*coldMass_f);
  coldTemp_avg = (coldTemp_i + coldTemp_f) / 2.0;
  coldTemp_avg2 = 0.0;

  count = 0;

  if ( !low_fid ) {
    while ( fabs(coldTemp_avg2 - coldTemp_avg) > temperaturePrecision && count < 15) {
    
      coldTemp_avg2 = coldTemp_avg;

      Q_dot_avg = _coldStorage.fComputeEnergyLosses(coldTemp_avg, coldLevl_avg);

      U_f = U_i + timeInSeconds*(coldMassFlow_in*Cp*coldTemp_in - coldMassFlow_out*Cp*coldTemp_avg - Q_dot_avg);
    
      coldTemp_f = U_f / (Cp*coldMass_f);
      if (coldTemp_f < MELTING_POINT) {
	Q_heat_cold = (MELTING_POINT - coldTemp_f) * Cp * coldMass_f;
	coldTemp_f = MELTING_POINT;
      }
      coldTemp_avg = (coldTemp_i + coldTemp_f) / 2.0;
      ++count;
    }
  }
  
  maxColdMass = _coldStorage.get_heightOfStorage() *
    PI * pow(_coldStorage.get_diameterOfStorage()/2.0,2.0) * MS_DENSITY;
  
  coldMass_f = maxColdMass - _hotStorage.get_storedMass();
  _coldStorage.set_storage(coldMass_f, coldTemp_f);
  
  // Data gathering:
  
  _steamGenOutletMsRate.push_back( _steamGeneratorOutlet.get_massFlow()    );
  _steamGenOutletTemp.push_back  ( _steamGeneratorOutlet.get_temperature() );
  
  if (_hotStorage.get_storedTemperature() != 0.0
      && _hotStorage.get_storedTemperature() < _minHotStorageTemp)
    _minHotStorageTemp = _hotStorage.get_storedTemperature();

  if (_coldStorage.get_storedTemperature() != 0.0
      && _coldStorage.get_storedTemperature() < _minColdStorageTemp)
    _minColdStorageTemp = _coldStorage.get_storedTemperature();

  if (_steamGeneratorOutlet.get_temperature() < _minSteamGenOutletTemp)
    _minSteamGenOutletTemp = _steamGeneratorOutlet.get_temperature();

  _storageHeat.push_back(Q_heat_cold + Q_heat_hot);
}

/*--------------------------------------------*/
void HtfCycle::initiateColdStorage ( void ) {
/*--------------------------------------------*/
  double volume            = PI*pow(_coldStorage.get_diameterOfStorage() / 2.0, 2.0)*_coldStorage.get_heightOfStorage();
  double storedMass        = volume * MS_DENSITY;
  double storedTemperature = _minColdStorageTemp;
  _coldStorage.set_storage(storedMass, storedTemperature);
}

/*---------------------------------------------------------------------------------*/
void HtfCycle::setStorage ( double percentHot, double tempHot, double tempCold ) {
/*---------------------------------------------------------------------------------*/
  double H      = _hotStorage.get_heightOfStorage();
  double cold_H = _coldStorage.get_heightOfStorage();
  _hotStorage.set_storage2(H*(percentHot / 100.0), tempHot);
  _coldStorage.set_storage2(cold_H - (H*percentHot / 100.0), tempCold);
}
