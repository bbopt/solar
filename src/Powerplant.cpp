#include "Powerplant.hpp"

Powerplant::Powerplant ( Clock          & time       ,
			 Sun            & sun        ,
			 int              model_type ,
			 HeliostatField * field      ,
			 HtfCycle       * htfCycle   ,
			 Powerblock     * powerblock ,
			 Economics      * economics    ) :
  _time            ( time       ) ,
  _sun             ( sun        ) ,
  _model_type      ( model_type ) ,
  _moltenSaltLoop  ( htfCycle   ) ,
  _heliostatsField ( field      ) ,
  _powerblock      ( powerblock ) ,
  _investmentCost  ( economics  )   {}

class foncteurSum
{
private:
	double* _sum;

public:
	foncteurSum(double* sum) { _sum = sum; }

	void operator()(double d){ (*_sum) += d; }
};

//This function simply runs the simulation
void Powerplant::fSimulatePowerplant()
{

	//Options:
	//1- Simulate the heliostats field or replace it with an input file
	//2- Simulate the heat exchanger or use the energy balance basic model
	//3- Use a dT dependant model for poweblock efficiency or use a simple
	//   constant efficiency model to convert thermal energy in expected electrical power
	//4- Simulate only the heliostat field
	//5- Use a demand model based on an input file, use a constant demand with 
	//   a fixated value, or use a daily demand model (summer or winter options) with data from ontario.

	if ( _model_type == 1 ) // heliostats field
	{
		fSimulateHeliostatField();
	}
	else if ( _model_type == 2 ) // whole plant
	{
		fSimulateHeliostatField();

		//Pre-setting the reserves
		_hotStorageLevel.reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
		_coldStorageLevel.reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
		_hotStorageTemp.reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
		_coldStorageTemp.reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
		_energyToPowerBlockWatts.reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
		_receiverPumpHead.reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
		_powerplantPowerOutput.reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
		_demandCompliance.reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
		_powerblock->get_powerOutput().reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
		_steamRate.reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
		_msRateSteamGen.reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
		_receiverOutletFlow.reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
		_heliostatsFieldPar.reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
		//Determining yield pressure for the receiver

		if (_moltenSaltLoop->get_steamGenerator().get_exchangerModel() == 2)
		{
			_pressureTubesSide.reserve(_time.get_numberOfIncrements()* _time.get_sizeOfIncrements() * 60);
			_maximumPressureInExchanger = 0.;
		}
		_maximumPressureInReceiver = 0.;
		_yieldPressureReceiver = _moltenSaltLoop->get_centralReceiver().computeYieldPressure();
		_yieldPressureExchanger = _moltenSaltLoop->get_steamGenerator().computeYieldPressure();
		// double pressureInReceiver, pressureInExchangerTubes;
		double thermalPowerNeeded, turbineThermalPower, steamRate;

		for (int j = 0; j < _time.get_numberOfIncrements(); ++j)
		{
			for (int i = 0; i < _time.get_sizeOfIncrements() * 60; ++i)
			{
				turbineThermalPower = _powerblock->fComputeRequiredThermalEnergy(_demand[j]);
				steamRate = fComputeSteamRate(turbineThermalPower);
				thermalPowerNeeded = fComputeThermalEnergy(steamRate);
				_powerblock->get_requiredThermalPower().push_back(thermalPowerNeeded);
				_moltenSaltLoop->fOperateCycle(1, _heliostatFieldPowerOutput[j], thermalPowerNeeded);

				/*Data gathering*/
				_hotStorageLevel.push_back(_moltenSaltLoop->get_hotStorage().get_heightOfVolumeStored());
				_coldStorageLevel.push_back(_moltenSaltLoop->get_coldStorage().get_heightOfVolumeStored());
				_hotStorageTemp.push_back(_moltenSaltLoop->get_hotStorage().get_storedTemperature());
				_coldStorageTemp.push_back(_moltenSaltLoop->get_coldStorage().get_storedTemperature());
				_energyToPowerBlockWatts.push_back(_moltenSaltLoop->get_steamGenerator().fEnergyToPowerBlock(1/**/));


				if (_heliostatFieldPowerOutput[j] > 0 && _heliostatsFieldModel == 1)
				{ _heliostatsFieldPar.push_back(55. * _heliostatsField->get_listOfHeliostats().size()); }
				else { _heliostatsFieldPar.push_back(0.); }

				if (_moltenSaltLoop->get_steamGeneratorOutlet().get_massFlow() > 0.)
				{
					//Filling powerplant power output vector
					_powerplantPowerOutput.push_back(_powerblock->get_Pout());

					//Filling demand compliance vector
					if (_powerplantPowerOutput.back() >= _demand[j]
						&& _demand[j] > 0.){
						_demandCompliance.push_back(1);
					}
					else if (_demand[j] > 0.){ _demandCompliance.push_back(_powerplantPowerOutput.back() / _demand[j]); }
				}
				else 
				{
					_powerplantPowerOutput.push_back(0.);
					//Filling demand compliance vector
					if (_demand[j] > 0.){ _demandCompliance.push_back(0.); }
				}

				//recording pressure in receiver tubes
				_receiverPumpHead.push_back(_moltenSaltLoop->get_centralReceiver().computePressureInTubes());
				_receiverPumpHead.back() > _maximumPressureInReceiver ? _maximumPressureInReceiver = _receiverPumpHead.back() : 1;

				if (_moltenSaltLoop->get_steamGenerator().get_exchangerModel() == 2)
				{
					//recording pressure in steam generator tubes
					_pressureTubesSide.push_back(_moltenSaltLoop->get_steamGenerator().computePressureInTubes(steamRate));

					//recording pressure in steam generator shell
					_pressureShellSide.push_back(_moltenSaltLoop->get_steamGenerator().computePressureInShells());

					_pressureTubesSide.back() + _powerblock->get_pressure() > _maximumPressureInExchanger ? _maximumPressureInExchanger = _pressureTubesSide.back() : 1;
				}

				//recording molten salt flows to receiver and exchanger
				_receiverOutletFlow.push_back(_moltenSaltLoop->get_centralReceiverOutlet().get_massFlow());
				_msRateSteamGen.push_back(_moltenSaltLoop->get_steamGeneratorOutlet().get_massFlow());

				//recording steam rate
				_steamRate.push_back(steamRate);
			}
		}

		double sumCompliance = 0;
		foncteurSum somme(&sumCompliance);
		std::for_each(_demandCompliance.begin(), _demandCompliance.end(), somme);
		unsigned int compliance = 0;
		if (_demandCompliance.size() > 0){ compliance = _demandCompliance.size(); }
		else{ compliance = 1; }
		_overallComplianceToDemand = sumCompliance*100. / compliance;
	}
}



void Powerplant::fSimulateHeliostatField()
{
	if (_heliostatsFieldModel == 1)
	{
		//Generating heliostats field
		_heliostatsField->fGenerateField();
		_heliostatFieldPowerOutput.reserve(_time.get_numberOfIncrements());

		//Economics
		int numberOfHeliostats = static_cast<int>(_heliostatsField->get_listOfHeliostats().size());
		_investmentCost->set_nbOfHeliostats(numberOfHeliostats);

		_totalEnergyConcentrated = 0.;

		//Simulating heliostats field
		_heliostatsField->fGenerateSunrays();
		for (int i = 0; i < _time.get_numberOfIncrements(); ++i)
		{
			_heliostatFieldPowerOutput.push_back(_heliostatsField->fCalculateTotalEnergyOutput());
			_heliostatsField->get_sun()->get_time().fTimeIncrement();
			_totalEnergyConcentrated += _heliostatFieldPowerOutput[i] * _time.get_sizeOfIncrements() * 60;
		}

		_reflectiveSurface = _heliostatsField->get_listOfHeliostats().size() *
			_heliostatsField->get_heliostatLength() * _heliostatsField->get_heliostatWidth();
		_investmentCost->set_reflectiveArea(_reflectiveSurface);

		_fieldSurface = _heliostatsField->get_maxAngularDeviation() * 2 * (M_PI / 180.)*
			(pow(_heliostatsField->get_towerHeight()*_heliostatsField->get_maxDistanceFromTower(), 2.)
			- pow(_heliostatsField->get_towerHeight()*_heliostatsField->get_minDistanceFromTower(), 2.));

		_heliostatFieldEfficiency = _heliostatsField->get_fieldEfficiency();

		//converting to kWh
		_totalEnergyConcentrated /= (1000. * 3600); // TOTO calculer direct
	}
	else if (_heliostatsFieldModel == 2)
	{

	  // _heliostatFieldPowerOutput has been filled with pre-fixed values
	  // (see methods set_heliostatFieldPowerOutput_XXX at the end of this file)
	  
	  
	  // TOTO TUTU  Virer tout ca ci-dessous

	  
	  // //std::string pathToSunData = pathToBlackbox;
	  // //int pos = pathToSunData.find_last_of(separatorString());
	  // //pathToSunData = pathToBlackbox.substr(0, pos);
	  // //	pathToSunData.append(_pathToFieldData);

	  // //Reading input file
	  // ifstream sunData("TOTO TATA", ios::in);
	  // if (sunData)
	  // 	{
	  // 		for (int i = 0; i < _time.get_numberOfIncrements(); ++i)
	  // 		{
	  // 			if (!sunData.eof())
	  // 			{
	  // 				sunData >> pwrOutput;
	  // 				_heliostatFieldPowerOutput.push_back(pwrOutput);
	  // 			}
	  // 			else{
	  // 				std::invalid_argument data("couldn't read solar energy output from file \n");
	  // 				throw data;
	  // 			}
	  // 		}
	  // 		sunData.close();
	  // 	}
	  // 	else
	  // 	{
	  // 		std::invalid_argument data("couldn't read solar energy output from file \n");
	  // 		throw data;
	  // 	}

	}
}

double Powerplant::fComputeSteamRate(double& turbineEnergy)
{
	double steamRate = turbineEnergy /
		(_powerblock->get_hotEnthalpy() - _powerblock->get_turbineOutletEnthalpy());

	_powerblock->set_steamRate(steamRate);

	return steamRate;
}

double Powerplant::fComputeThermalEnergy(double& steamRate)
{
	double thermalEnergy = steamRate * (_powerblock->get_hotEnthalpy() - WATER_300K_1ATM_ENTHALPY);

	return thermalEnergy;
}

//Computes energy lost in main components in order to either maintain temperature or
//drive the main pumps.
class foncteurSum9
{
private:
	double* _sum;

public:
	foncteurSum9(double* s){ _sum = s; }

	void operator()(double d){ (*_sum) = (*_sum) + d; }
};

double Powerplant::fComputeParasiticLosses()
{
	double W_rec, W_shell, W_steam, Q_helios, Q_heat;
	double total = 0.;
	W_rec = 0;
	W_shell = 0;
	W_steam = 0;
	Q_helios = 0;
	Q_heat = 0;

	W_rec = std::inner_product(_receiverPumpHead.begin(), _receiverPumpHead.end(),
		_receiverOutletFlow.begin(), 0.0);

	W_rec /= MS_DENSITY;

	W_shell = std::inner_product(_pressureShellSide.begin(), _pressureShellSide.end(),
		_msRateSteamGen.begin(), 0.0);
	W_shell /= MS_DENSITY;

	W_steam = std::inner_product(_pressureTubesSide.begin(), _pressureTubesSide.end(),
		_steamRate.begin(), 0.0);
	W_steam /= WATER_DENSITY;

	foncteurSum9 sumQ_heat(&Q_heat);
	std::for_each(_moltenSaltLoop->get_storageHeatV().begin(),
		_moltenSaltLoop->get_storageHeatV().end(),
		sumQ_heat);

	foncteurSum9 sumQ_helios(&Q_helios);
	std::for_each(_heliostatsFieldPar.begin(), _heliostatsFieldPar.end(),
		sumQ_helios);

	total = W_rec + W_shell + W_steam + Q_heat + Q_helios;

	return total;
}

double Powerplant::fComputeParasiticsForPb3()
{
	double W_rec, W_shell, W_steam, Q_helios, Q_heat;
	double total = 0.;
	W_rec = 0;
	W_shell = 0;
	W_steam = 0;
	Q_helios = 0;
	Q_heat = 0;

	W_rec = std::inner_product(_receiverPumpHead.begin(), _receiverPumpHead.end(),
		_receiverOutletFlow.begin(), 0.0);

	W_rec /= MS_DENSITY;

	foncteurSum9 sumQ_heat(&Q_heat);
	std::for_each(_moltenSaltLoop->get_storageHeatV().begin(),
		_moltenSaltLoop->get_storageHeatV().end(),
		sumQ_heat);

	foncteurSum9 sumQ_helios(&Q_helios);
	std::for_each(_heliostatsFieldPar.begin(), _heliostatsFieldPar.end(),
		sumQ_helios);

	total = 2.*W_rec + Q_heat + Q_helios;

	return total;
}

double Powerplant::fComputeParasiticsForPb7()
{
	double W_rec;

	W_rec = 0.;
	W_rec = std::inner_product(_receiverPumpHead.begin(), _receiverPumpHead.end(),
		_receiverOutletFlow.begin(), 0.0);

	W_rec /= MS_DENSITY;

	return W_rec;
}


double Powerplant::fComputeParasiticsForPb9()
{
	double W_rec, Q_helios;
	double total = 0.;

	W_rec = 0.;
	Q_helios = 0.;
	
	W_rec = std::inner_product(_receiverPumpHead.begin(), _receiverPumpHead.end(),
		_receiverOutletFlow.begin(), 0.0);
	W_rec /= MS_DENSITY;

	foncteurSum9 somme(&Q_helios);
	std::for_each(_heliostatsFieldPar.begin(), _heliostatsFieldPar.end(),
		somme);
	

	total = W_rec + Q_helios;
	return total;
}


// TOTO TUTU

// For solar6:
bool Powerplant::set_heliostatFieldPowerOutput_MINCOST_TS ( void ) {
  
  if ( _time.get_numberOfIncrements() != 24 )
    return false;
  
  _heliostatFieldPowerOutput.resize(24);
  for ( int i = 0 ; i < 24 ; ++i )
    _heliostatFieldPowerOutput[i] = 0.0;

  _heliostatFieldPowerOutput[ 6] = 682099024;
  _heliostatFieldPowerOutput[ 7] = 1008520474;
  _heliostatFieldPowerOutput[ 8] = 1191661173;
  _heliostatFieldPowerOutput[ 9] = 1288329801;
  _heliostatFieldPowerOutput[10] = 1336037074;
  _heliostatFieldPowerOutput[11] = 1336737373;
  _heliostatFieldPowerOutput[12] = 1331348545;
  _heliostatFieldPowerOutput[13] = 1335783774;
  _heliostatFieldPowerOutput[14] = 1331944544;
  _heliostatFieldPowerOutput[15] = 1276484312;
  _heliostatFieldPowerOutput[16] = 1190983224;
  _heliostatFieldPowerOutput[17] = 1006943559;
  _heliostatFieldPowerOutput[18] = 683127123;

  return true;
}

// For solar5:
bool Powerplant::set_heliostatFieldPowerOutput_MAXCOMP_HTF1 ( void ) {

  if ( _time.get_numberOfIncrements() != 720 )
    return false;

  _heliostatFieldPowerOutput.resize(720);
  for ( int i = 0 ; i < 720 ; ++i )
    _heliostatFieldPowerOutput[i] = 0.0;

  _heliostatFieldPowerOutput[7]=6380504;
  _heliostatFieldPowerOutput[8]=31279320;
  _heliostatFieldPowerOutput[9]=37648690;
  _heliostatFieldPowerOutput[10]=40425920;
  _heliostatFieldPowerOutput[11]=41574850;
  _heliostatFieldPowerOutput[12]=42259420;
  _heliostatFieldPowerOutput[13]=42785530;
  _heliostatFieldPowerOutput[14]=41391120;
  _heliostatFieldPowerOutput[15]=40368350;
  _heliostatFieldPowerOutput[16]=32116270;
  _heliostatFieldPowerOutput[17]=6488462;
  _heliostatFieldPowerOutput[31]=5742453.6;
  _heliostatFieldPowerOutput[32]=28151388;
  _heliostatFieldPowerOutput[33]=33883821;
  _heliostatFieldPowerOutput[34]=36383328;
  _heliostatFieldPowerOutput[35]=37417365;
  _heliostatFieldPowerOutput[36]=38033478;
  _heliostatFieldPowerOutput[37]=38506977;
  _heliostatFieldPowerOutput[38]=37252008;
  _heliostatFieldPowerOutput[39]=36331515;
  _heliostatFieldPowerOutput[40]=28904643;
  _heliostatFieldPowerOutput[41]=5839615.8;
  _heliostatFieldPowerOutput[55]=6061478.8;
  _heliostatFieldPowerOutput[56]=29715354;
  _heliostatFieldPowerOutput[57]=35766255.5;
  _heliostatFieldPowerOutput[58]=38404624;
  _heliostatFieldPowerOutput[59]=39496107.5;
  _heliostatFieldPowerOutput[60]=40146449;
  _heliostatFieldPowerOutput[61]=40646253.5;
  _heliostatFieldPowerOutput[62]=39321564;
  _heliostatFieldPowerOutput[63]=38349932.5;
  _heliostatFieldPowerOutput[64]=30510456.5;
  _heliostatFieldPowerOutput[65]=6164038.9;
  _heliostatFieldPowerOutput[79]=5104403.2;
  _heliostatFieldPowerOutput[80]=25023456;
  _heliostatFieldPowerOutput[81]=30118952;
  _heliostatFieldPowerOutput[82]=30319440;
  _heliostatFieldPowerOutput[83]=31181137.5;
  _heliostatFieldPowerOutput[84]=29581594;
  _heliostatFieldPowerOutput[85]=29949871;
  _heliostatFieldPowerOutput[86]=28973784;
  _heliostatFieldPowerOutput[87]=28257845;
  _heliostatFieldPowerOutput[88]=22481389;
  _heliostatFieldPowerOutput[89]=4541923.4;
  _heliostatFieldPowerOutput[103]=4466352.8;
  _heliostatFieldPowerOutput[104]=21895524;
  _heliostatFieldPowerOutput[105]=26354083;
  _heliostatFieldPowerOutput[106]=32340736;
  _heliostatFieldPowerOutput[107]=33259880;
  _heliostatFieldPowerOutput[108]=33807536;
  _heliostatFieldPowerOutput[109]=34228424;
  _heliostatFieldPowerOutput[110]=33112896;
  _heliostatFieldPowerOutput[111]=32294680;
  _heliostatFieldPowerOutput[112]=25693016;
  _heliostatFieldPowerOutput[113]=5190769.6;
  _heliostatFieldPowerOutput[127]=6061478.8;
  _heliostatFieldPowerOutput[128]=29715354;
  _heliostatFieldPowerOutput[129]=35766255.5;
  _heliostatFieldPowerOutput[130]=38404624;
  _heliostatFieldPowerOutput[131]=39496107.5;
  _heliostatFieldPowerOutput[132]=40146449;
  _heliostatFieldPowerOutput[133]=40646253.5;
  _heliostatFieldPowerOutput[134]=39321564;
  _heliostatFieldPowerOutput[135]=38349932.5;
  _heliostatFieldPowerOutput[136]=30510456.5;
  _heliostatFieldPowerOutput[137]=6164038.9;
  _heliostatFieldPowerOutput[151]=6061478.8;
  _heliostatFieldPowerOutput[152]=29715354;
  _heliostatFieldPowerOutput[153]=35766255.5;
  _heliostatFieldPowerOutput[154]=38404624;
  _heliostatFieldPowerOutput[155]=39496107.5;
  _heliostatFieldPowerOutput[156]=40146449;
  _heliostatFieldPowerOutput[157]=36367700.5;
  _heliostatFieldPowerOutput[158]=31043340;
  _heliostatFieldPowerOutput[159]=26239427.5;
  _heliostatFieldPowerOutput[160]=19269762;
  _heliostatFieldPowerOutput[161]=3893077.2;
  _heliostatFieldPowerOutput[175]=5104403.2;
  _heliostatFieldPowerOutput[176]=25023456;
  _heliostatFieldPowerOutput[177]=30118952;
  _heliostatFieldPowerOutput[178]=32340736;
  _heliostatFieldPowerOutput[179]=33259880;
  _heliostatFieldPowerOutput[180]=33807536;
  _heliostatFieldPowerOutput[181]=38506977;
  _heliostatFieldPowerOutput[182]=37252008;
  _heliostatFieldPowerOutput[183]=40368350;
  _heliostatFieldPowerOutput[184]=32116270;
  _heliostatFieldPowerOutput[185]=6488462;
  _heliostatFieldPowerOutput[199]=6380504;
  _heliostatFieldPowerOutput[200]=31279320;
  _heliostatFieldPowerOutput[201]=37648690;
  _heliostatFieldPowerOutput[202]=40425920;
  _heliostatFieldPowerOutput[203]=41574850;
  _heliostatFieldPowerOutput[204]=42259420;
  _heliostatFieldPowerOutput[205]=42785530;
  _heliostatFieldPowerOutput[206]=41391120;
  _heliostatFieldPowerOutput[207]=40368350;
  _heliostatFieldPowerOutput[208]=32116270;
  _heliostatFieldPowerOutput[209]=6488462;
  _heliostatFieldPowerOutput[223]=6380504;
  _heliostatFieldPowerOutput[224]=31279320;
  _heliostatFieldPowerOutput[225]=37648690;
  _heliostatFieldPowerOutput[226]=40425920;
  _heliostatFieldPowerOutput[227]=41574850;
  _heliostatFieldPowerOutput[228]=38033478;
  _heliostatFieldPowerOutput[229]=38506977;
  _heliostatFieldPowerOutput[230]=37252008;
  _heliostatFieldPowerOutput[231]=36331515;
  _heliostatFieldPowerOutput[232]=28904643;
  _heliostatFieldPowerOutput[233]=5839615.8;
  _heliostatFieldPowerOutput[247]=5742453.6;
  _heliostatFieldPowerOutput[248]=25023456;
  _heliostatFieldPowerOutput[249]=26354083;
  _heliostatFieldPowerOutput[250]=24255552;
  _heliostatFieldPowerOutput[251]=20787425;
  _heliostatFieldPowerOutput[252]=16903768;
  _heliostatFieldPowerOutput[253]=17114212;
  _heliostatFieldPowerOutput[254]=12417336;
  _heliostatFieldPowerOutput[255]=28257845;
  _heliostatFieldPowerOutput[256]=28904643;
  _heliostatFieldPowerOutput[257]=5839615.8;
  _heliostatFieldPowerOutput[271]=5742453.6;
  _heliostatFieldPowerOutput[272]=28151388;
  _heliostatFieldPowerOutput[273]=41413559;
  _heliostatFieldPowerOutput[274]=44468512;
  _heliostatFieldPowerOutput[275]=47811077.5;
  _heliostatFieldPowerOutput[276]=49020927.2;
  _heliostatFieldPowerOutput[277]=50059070.1;
  _heliostatFieldPowerOutput[278]=48841521.6;
  _heliostatFieldPowerOutput[279]=48038336.5;
  _heliostatFieldPowerOutput[280]=38539524;
  _heliostatFieldPowerOutput[281]=6488462;
  _heliostatFieldPowerOutput[295]=6380504;
  _heliostatFieldPowerOutput[296]=31279320;
  _heliostatFieldPowerOutput[297]=26354083;
  _heliostatFieldPowerOutput[298]=24255552;
  _heliostatFieldPowerOutput[299]=41574850;
  _heliostatFieldPowerOutput[300]=42259420;
  _heliostatFieldPowerOutput[301]=42785530;
  _heliostatFieldPowerOutput[302]=33112896;
  _heliostatFieldPowerOutput[303]=24221010;
  _heliostatFieldPowerOutput[304]=19269762;
  _heliostatFieldPowerOutput[305]=3893077.2;
  _heliostatFieldPowerOutput[319]=3828302.4;
  _heliostatFieldPowerOutput[320]=18767592;
  _heliostatFieldPowerOutput[321]=22589214;
  _heliostatFieldPowerOutput[322]=26276848;
  _heliostatFieldPowerOutput[323]=27855149.5;
  _heliostatFieldPowerOutput[324]=29581594;
  _heliostatFieldPowerOutput[325]=34228424;
  _heliostatFieldPowerOutput[326]=28973784;
  _heliostatFieldPowerOutput[327]=24221010;
  _heliostatFieldPowerOutput[328]=16058135;
  _heliostatFieldPowerOutput[329]=3244231;
  _heliostatFieldPowerOutput[343]=6380504;
  _heliostatFieldPowerOutput[344]=31279320;
  _heliostatFieldPowerOutput[345]=41413559;
  _heliostatFieldPowerOutput[346]=40425920;
  _heliostatFieldPowerOutput[347]=33259880;
  _heliostatFieldPowerOutput[348]=29581594;
  _heliostatFieldPowerOutput[349]=25671318;
  _heliostatFieldPowerOutput[350]=28973784;
  _heliostatFieldPowerOutput[351]=40368350;
  _heliostatFieldPowerOutput[352]=32116270;
  _heliostatFieldPowerOutput[353]=6488462;
  _heliostatFieldPowerOutput[367]=5742453.6;
  _heliostatFieldPowerOutput[368]=28151388;
  _heliostatFieldPowerOutput[369]=33883821;
  _heliostatFieldPowerOutput[370]=36383328;
  _heliostatFieldPowerOutput[371]=37417365;
  _heliostatFieldPowerOutput[372]=38033478;
  _heliostatFieldPowerOutput[373]=38506977;
  _heliostatFieldPowerOutput[374]=37252008;
  _heliostatFieldPowerOutput[375]=36331515;
  _heliostatFieldPowerOutput[376]=28904643;
  _heliostatFieldPowerOutput[377]=5839615.8;
  _heliostatFieldPowerOutput[391]=6061478.8;
  _heliostatFieldPowerOutput[392]=29715354;
  _heliostatFieldPowerOutput[393]=35766255.5;
  _heliostatFieldPowerOutput[394]=38404624;
  _heliostatFieldPowerOutput[395]=39496107.5;
  _heliostatFieldPowerOutput[396]=40146449;
  _heliostatFieldPowerOutput[397]=40646253.5;
  _heliostatFieldPowerOutput[398]=39321564;
  _heliostatFieldPowerOutput[399]=38349932.5;
  _heliostatFieldPowerOutput[400]=30510456.5;
  _heliostatFieldPowerOutput[401]=6164038.9;
  _heliostatFieldPowerOutput[415]=4466352.8;
  _heliostatFieldPowerOutput[416]=21895524;
  _heliostatFieldPowerOutput[417]=26354083;
  _heliostatFieldPowerOutput[418]=28298144;
  _heliostatFieldPowerOutput[419]=29102395;
  _heliostatFieldPowerOutput[420]=29581594;
  _heliostatFieldPowerOutput[421]=29949871;
  _heliostatFieldPowerOutput[422]=28973784;
  _heliostatFieldPowerOutput[423]=28257845;
  _heliostatFieldPowerOutput[424]=22481389;
  _heliostatFieldPowerOutput[425]=4541923.4;
  _heliostatFieldPowerOutput[439]=5104403.2;
  _heliostatFieldPowerOutput[440]=25023456;
  _heliostatFieldPowerOutput[441]=30118952;
  _heliostatFieldPowerOutput[442]=32340736;
  _heliostatFieldPowerOutput[443]=33259880;
  _heliostatFieldPowerOutput[444]=33807536;
  _heliostatFieldPowerOutput[445]=34228424;
  _heliostatFieldPowerOutput[446]=33112896;
  _heliostatFieldPowerOutput[447]=32294680;
  _heliostatFieldPowerOutput[448]=25693016;
  _heliostatFieldPowerOutput[449]=5190769.6;
  _heliostatFieldPowerOutput[463]=6061478.8;
  _heliostatFieldPowerOutput[464]=29715354;
  _heliostatFieldPowerOutput[465]=35766255.5;
  _heliostatFieldPowerOutput[466]=38404624;
  _heliostatFieldPowerOutput[467]=39496107.5;
  _heliostatFieldPowerOutput[468]=40146449;
  _heliostatFieldPowerOutput[469]=40646253.5;
  _heliostatFieldPowerOutput[470]=39321564;
  _heliostatFieldPowerOutput[471]=38349932.5;
  _heliostatFieldPowerOutput[472]=30510456.5;
  _heliostatFieldPowerOutput[473]=6164038.9;
  _heliostatFieldPowerOutput[487]=5423428.4;
  _heliostatFieldPowerOutput[488]=23459490;
  _heliostatFieldPowerOutput[489]=24471648.5;
  _heliostatFieldPowerOutput[490]=24255552;
  _heliostatFieldPowerOutput[491]=24944910;
  _heliostatFieldPowerOutput[492]=25355652;
  _heliostatFieldPowerOutput[493]=25671318;
  _heliostatFieldPowerOutput[494]=24834672;
  _heliostatFieldPowerOutput[495]=24221010;
  _heliostatFieldPowerOutput[496]=19269762;
  _heliostatFieldPowerOutput[497]=3893077.2;
  _heliostatFieldPowerOutput[511]=5742453.6;
  _heliostatFieldPowerOutput[512]=28151388;
  _heliostatFieldPowerOutput[513]=37648690;
  _heliostatFieldPowerOutput[514]=40425920;
  _heliostatFieldPowerOutput[515]=41574850;
  _heliostatFieldPowerOutput[516]=42259420;
  _heliostatFieldPowerOutput[517]=42785530;
  _heliostatFieldPowerOutput[518]=41391120;
  _heliostatFieldPowerOutput[519]=40368350;
  _heliostatFieldPowerOutput[520]=32116270;
  _heliostatFieldPowerOutput[521]=6488462;
  _heliostatFieldPowerOutput[535]=5104403.2;
  _heliostatFieldPowerOutput[536]=21895524;
  _heliostatFieldPowerOutput[537]=18824345;
  _heliostatFieldPowerOutput[538]=12127776;
  _heliostatFieldPowerOutput[539]=8314970;
  _heliostatFieldPowerOutput[540]=8451884;
  _heliostatFieldPowerOutput[541]=12835659;
  _heliostatFieldPowerOutput[542]=12417336;
  _heliostatFieldPowerOutput[543]=20184175;
  _heliostatFieldPowerOutput[544]=16058135;
  _heliostatFieldPowerOutput[545]=2595384.8;
  _heliostatFieldPowerOutput[559]=2552201.6;
  _heliostatFieldPowerOutput[560]=18767592;
  _heliostatFieldPowerOutput[561]=33883821;
  _heliostatFieldPowerOutput[562]=36383328;
  _heliostatFieldPowerOutput[563]=37417365;
  _heliostatFieldPowerOutput[564]=29581594;
  _heliostatFieldPowerOutput[565]=38506977;
  _heliostatFieldPowerOutput[566]=28973784;
  _heliostatFieldPowerOutput[567]=36331515;
  _heliostatFieldPowerOutput[568]=25693016;
  _heliostatFieldPowerOutput[569]=5839615.8;
  _heliostatFieldPowerOutput[583]=5742453.6;
  _heliostatFieldPowerOutput[584]=28151388;
  _heliostatFieldPowerOutput[585]=33883821;
  _heliostatFieldPowerOutput[586]=36383328;
  _heliostatFieldPowerOutput[587]=37417365;
  _heliostatFieldPowerOutput[588]=38033478;
  _heliostatFieldPowerOutput[589]=38506977;
  _heliostatFieldPowerOutput[590]=37252008;
  _heliostatFieldPowerOutput[591]=36331515;
  _heliostatFieldPowerOutput[592]=28904643;
  _heliostatFieldPowerOutput[593]=5839615.8;
  _heliostatFieldPowerOutput[607]=7465189.68;
  _heliostatFieldPowerOutput[608]=36909597.6;
  _heliostatFieldPowerOutput[609]=44801941.1;
  _heliostatFieldPowerOutput[610]=48511104;
  _heliostatFieldPowerOutput[611]=41574850;
  _heliostatFieldPowerOutput[612]=42259420;
  _heliostatFieldPowerOutput[613]=42785530;
  _heliostatFieldPowerOutput[614]=41391120;
  _heliostatFieldPowerOutput[615]=40368350;
  _heliostatFieldPowerOutput[616]=32116270;
  _heliostatFieldPowerOutput[617]=6488462;
  _heliostatFieldPowerOutput[631]=6380504;
  _heliostatFieldPowerOutput[632]=25023456;
  _heliostatFieldPowerOutput[633]=22589214;
  _heliostatFieldPowerOutput[634]=24255552;
  _heliostatFieldPowerOutput[635]=24944910;
  _heliostatFieldPowerOutput[636]=25355652;
  _heliostatFieldPowerOutput[637]=25671318;
  _heliostatFieldPowerOutput[638]=24834672;
  _heliostatFieldPowerOutput[639]=24221010;
  _heliostatFieldPowerOutput[640]=19269762;
  _heliostatFieldPowerOutput[641]=3893077.2;
  _heliostatFieldPowerOutput[655]=5104403.2;
  _heliostatFieldPowerOutput[656]=21895524;
  _heliostatFieldPowerOutput[657]=22589214;
  _heliostatFieldPowerOutput[658]=20212960;
  _heliostatFieldPowerOutput[659]=20787425;
  _heliostatFieldPowerOutput[660]=25355652;
  _heliostatFieldPowerOutput[661]=29949871;
  _heliostatFieldPowerOutput[662]=33112896;
  _heliostatFieldPowerOutput[663]=36331515;
  _heliostatFieldPowerOutput[664]=32116270;
  _heliostatFieldPowerOutput[665]=6488462;
  _heliostatFieldPowerOutput[679]=6380504;
  _heliostatFieldPowerOutput[680]=31279320;
  _heliostatFieldPowerOutput[681]=37648690;
  _heliostatFieldPowerOutput[682]=40425920;
  _heliostatFieldPowerOutput[683]=41574850;
  _heliostatFieldPowerOutput[684]=42259420;
  _heliostatFieldPowerOutput[685]=42785530;
  _heliostatFieldPowerOutput[686]=41391120;
  _heliostatFieldPowerOutput[687]=40368350;
  _heliostatFieldPowerOutput[688]=32116270;
  _heliostatFieldPowerOutput[689]=6488462;
  _heliostatFieldPowerOutput[703]=6380504;
  _heliostatFieldPowerOutput[704]=31279320;
  _heliostatFieldPowerOutput[705]=37648690;
  _heliostatFieldPowerOutput[706]=40425920;
  _heliostatFieldPowerOutput[707]=41574850;
  _heliostatFieldPowerOutput[708]=42259420;
  _heliostatFieldPowerOutput[709]=42785530;
  _heliostatFieldPowerOutput[710]=41391120;
  _heliostatFieldPowerOutput[711]=28257845;
  _heliostatFieldPowerOutput[712]=16058135;
  _heliostatFieldPowerOutput[713]=2595384.8;
    
  return true;
}
