#include "Economics.hpp"

Economics::Economics()
{
	//module costs
	_costOfField = 0.;
	_costPerHeliostat = 0.;
	_costOfTower = 0.;
	_costOfStorage = 0.;
	_costOfPowerblock = 0.;
	_costOfSteamGenerator = 0.;
	_costOfReceiver = 0.;
	_totalCost = 0.;

	//parameters
	_hotStorageInsulationThickness = 0.01;
	_coldStorageInsulationThickness = 0.01;
	_hotStorageHeight = 1.;
	_storageDiameter = 1.;
	_receiverInsulationThickness = 0.01;
	_heightOfTower = 20.;
	_heightOfReceiverAperture = 1.;
	_widthOfReceiverAperture = 1.;
	_receiverNumberOfTubes = 0;
	_receiverTubesDout = 0.;
	_lengthOfHeliostats = 1.;
	_widthOfHeliostats = 1.;
	_nbOfHeliostats = 1;
	_reflectiveArea = 0.;
	_totalMoltenSaltMass = 0.;
	_turbineNominalPowerOutput = 0.;
	_exchangerModel = 1;
	_exchangerTubesOutterDiameter = 0.0005;
	_exchangerTubesLength = 1.;
	_exchangerNumberOfTubes = 1;
	_exchangerTubePassesPerShell = 1;
	_exchangerNumberOfShell = 1;
}

Economics::Economics(int numberOfHeliostats, double hotStorageInsul, double coldStorageInsul,
	double hotStorageHeight, double storageDiameter, double receiverInsul, double heightOfTower, double lengthOfHeliostats, double widthOfHeliostats,
	double turbinePower, double heightOfReceiverAperture, double widthOfAperture,
	double exchangerTubesDiameter, double exchangerTubesLength, 
	int exchangerNumberOfTubes, int exchangerTubePassesPerShell, int exchangerNumberOfShell)
	:
  _hotStorageInsulationThickness(hotStorageInsul),
  _coldStorageInsulationThickness(coldStorageInsul),
  _hotStorageHeight(hotStorageHeight),
  _storageDiameter(storageDiameter),
  _receiverInsulationThickness(receiverInsul),
  _heightOfTower(heightOfTower),
  _heightOfReceiverAperture(heightOfReceiverAperture),
  _widthOfReceiverAperture(widthOfAperture),
  _lengthOfHeliostats(lengthOfHeliostats),
  _widthOfHeliostats(widthOfHeliostats),
  _nbOfHeliostats(numberOfHeliostats),
  _turbineNominalPowerOutput(turbinePower),
  _exchangerTubesOutterDiameter(exchangerTubesDiameter),
  _exchangerTubesLength(exchangerTubesLength),
  _exchangerNumberOfTubes(exchangerNumberOfTubes),
  _exchangerTubePassesPerShell(exchangerTubePassesPerShell),
  _exchangerNumberOfShell(exchangerNumberOfShell)
{
	//_reflectiveArea = widthOfHeliostats*lengthOfHeliostats*numberOfHeliostats;
	//evaluateCostOfHeliostat();
	//evaluateCostOfField();
	//evaluateCostOfTower();
	//evaluateCostOfStorage();
	//evaluateCostOfPowerblock();
	//evaluateCostOfSteamGenerator();
	//evaluateTotalInvestmentCost();
	_costOfField = 0.;
	_costPerHeliostat = 0.;
	_costOfTower = 0.;
	_costOfStorage = 0.;
	_costOfPowerblock = 0.;
	_costOfSteamGenerator = 0.;
	_costOfReceiver = 0.;
	_totalCost = 0.;
	_receiverTubesDout = 0.;
	_receiverNumberOfTubes = 0;
}

double Economics::evaluateCostOfHeliostat()
{
	double A_h = _widthOfHeliostats*_lengthOfHeliostats; //reflective area of 1 heliostat
	double A_f = _reflectiveArea;
	double X = A_h - A_h_ref;
	double Y = A_f - A_f_ref;
	double mirrorCostPerM2 = C_M*exp(M_Kh*X + M_Kf*Y);
	double drivesCostPerM2 = C_D*exp(C_Kh*X + C_Kf*Y);
	double pedestalsCostPerM2 = C_P*exp(P_Kh*X + P_Kf*Y);
	double controlsCostPerM2 = C_C*exp(C_Kh*X + C_Kf*Y);
	double wiringCostPerM2 = C_W*exp(W_Kh*X + W_Kf*Y);
	double manufacturingCostPerM2 = C_Ma*exp(Ma_Kh*X + Ma_Kf*Y);
	double installationCostPerM2 = C_I*exp(I_Kh*X + I_Kf*Y);

	_costPerHeliostat = A_h * (mirrorCostPerM2
		+ drivesCostPerM2
		+ pedestalsCostPerM2
		+ controlsCostPerM2
		+ wiringCostPerM2
		+ manufacturingCostPerM2
		+ installationCostPerM2);

	return _costPerHeliostat;
}

double Economics::evaluateCostOfField()
{
	evaluateCostOfHeliostat();
	_costOfField = _nbOfHeliostats * _costPerHeliostat;
	return _costOfField;
}

double Economics::evaluateCostOfTower()
{
	double c0 = 1.*COST_OF_TOWER_CONSTANT;
	double c1 = 1.*COST_OF_TOWER_LIN;
	double c2 = 1.*COST_OF_TOWER_QUAD;

	_costOfTower = c0 + c1*_heightOfTower + c2*pow(_heightOfTower, 2.);

	return _costOfTower;
}

double Economics::evaluateCostOfReceiver()
{

	int numberOfPasses; 
	double A_r;//absorbing surface area
	if (_receiverNumberOfTubes != 0 && _receiverTubesDout != 0)
	{
		numberOfPasses = (int)floor((M_PI*_widthOfReceiverAperture / 2.) /
			(_receiverNumberOfTubes*_receiverTubesDout));
		A_r = M_PI*_receiverTubesDout*_heightOfReceiverAperture
			*_receiverNumberOfTubes*numberOfPasses;
	}
	else {
		A_r = M_PI*_widthOfReceiverAperture*_heightOfReceiverAperture / 2.;
	}

	_costOfReceiver = 1.*RECEIVER_REFERENCE_COST*RECEIVER_MTL_COST_ESC_RATE*pow(A_r / (1.*RECEIVER_REFERENCE_SURFACE), 1.*RECEIVER_EXPONENT_COEFFICIENT);
	return _costOfReceiver;
}

double Economics::evaluateCostOfStorage()
{
	double moltenSaltPerKg = COST_NANO3_KNO3; //$/kg
	//total molten salt inventory is assumed to be that of the volume of the full cold
	//tank.
	double moltenSaltVolume = _hotStorageHeight*1.1*M_PI*pow(_storageDiameter / 2., 2.);
	double _totalMoltenSaltMass = moltenSaltVolume*MS_DENSITY;
	double moltenSaltCost = moltenSaltPerKg * _totalMoltenSaltMass;

	double insulationVolume = M_PI*_hotStorageHeight*(pow(_storageDiameter / 2. + _hotStorageInsulationThickness + 0.04, 2.) - pow(_storageDiameter / 2. + 0.04, 2.))
		+ M_PI*_hotStorageHeight*1.1*(pow(_storageDiameter / 2. + _coldStorageInsulationThickness + 0.04, 2.) - pow(_storageDiameter / 2. + 0.04, 2.))
		+ M_PI*pow(_storageDiameter / 2. + 0.04, 2.)*(_coldStorageInsulationThickness + _hotStorageInsulationThickness);

	//the 4cm thick stainless steel tank is considered when evaluating the total volume of insulation needed.
	//The design tank diameter is the inner diameter.

	double ceramicFiberCost = CERAMIC_FIBER_INSULATION_COST*insulationVolume;
	
	double foundationCost = _totalMoltenSaltMass * STORAGE_TANK_FOUNDATION_COST_COEF + 1.*STORAGE_TANK_FOUNDATION_COST_CONST;

	_costOfStorage = ceramicFiberCost + moltenSaltCost + foundationCost;

	return _costOfStorage;
}

double Economics::evaluateCostOfSteamGenerator()
{
	double areaOfAShell, priceOfASquareMeter;
	double equip, pump, pipe, control, supp;
	//Cost is provided in table 5 of the Roadmap report for cost reduction
	//We use the Roadmap baseline cost $/kWe for consistency
	if (_exchangerModel == 2)
	{
		areaOfAShell = _exchangerTubesLength*M_PI*2.*(_exchangerTubesOutterDiameter / 2.)
			*_exchangerNumberOfTubes*_exchangerTubePassesPerShell;

		priceOfASquareMeter = 1.*STEAM_GENERATOR_CONSTANT*pow(areaOfAShell, 1.*STEAM_GENERATOR_EXPONENT);

		_costOfSteamGenerator = areaOfAShell*priceOfASquareMeter*_exchangerNumberOfShell;
	}
	else {
		equip = 1.*STEAM_GEN_REF_EQUIP * pow((_turbineNominalPowerOutput*(1e-6)) / (1.*STEAM_GEN_REF_POWER), 1.*STEAM_GEN_SCALE);
		pump = 1.*STEAM_GEN_REF_PUMP * pow((_turbineNominalPowerOutput*(1e-6)) / (1.*STEAM_GEN_REF_POWER), 1.*STEAM_GEN_SCALE);
		pipe = 1.*STEAM_GEN_REF_PIPE * pow((_turbineNominalPowerOutput*(1e-6)) / (1.*STEAM_GEN_REF_POWER), 1.*STEAM_GEN_SCALE);
		control = 1.*STEAM_GEN_REF_CONTROL * pow((_turbineNominalPowerOutput*(1e-6)) / (1.*STEAM_GEN_REF_POWER), 1.*STEAM_GEN_SCALE);
		supp = 1.*STEAM_GEN_REF_SUPP * pow((_turbineNominalPowerOutput*(1e-6)) / (1.*STEAM_GEN_REF_POWER), 1.*STEAM_GEN_SCALE);

		_costOfSteamGenerator = equip + pump + pipe + control + supp;
	}
	return _costOfSteamGenerator;
}

double Economics::evaluateCostOfPowerblock()
{
	//cost is provided in table
	_costOfPowerblock = (1.*POWERBLOCK_REFERENCE_COST*POWERBLOCK_MTL_COST_ESC_RATE)*pow(_turbineNominalPowerOutput / (1.*POWERBLOCK_REFERENCE_POWER), 1.*POWERBLOCK_EXPONENT_COEF);
	return _costOfPowerblock;
}

double Economics::evaluateTotalInvestmentCost()
{
	_totalCost = _costOfField
		+ _costOfTower
		+ _costOfStorage
		+ _costOfPowerblock
		+ _costOfSteamGenerator
		+ _costOfReceiver;

	return _totalCost;
}
