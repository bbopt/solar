#include "ThermalStorage.hpp"

ThermalStorage::ThermalStorage(MoltenSalt* input, MoltenSalt* output, double height, double diameter, double thickness)
:_inputHTF(input),
_outputHTF(output),
_heightOfStorage(height),
_diameterOfStorage(diameter),
_thicknessOfInsulation(thickness)
{
	_storedTemperature = 0.;
	_storedMass = 0.;
}

//
double ThermalStorage::fComputeStorageMass(int timeInterval)
{
	double storedMass;
	double timeInSeconds = timeInterval * 60.;
	storedMass = _storedMass + timeInSeconds * (_inputHTF->get_massFlow() - _outputHTF->get_massFlow());
	return storedMass;
}

double ThermalStorage::fComputeStorageLevel()
{
	_heightOfVolumeStored = (_storedMass / MS_DENSITY) / (pow(_diameterOfStorage, 2.)*M_PI / 4.);
	return _heightOfVolumeStored;
}

double ThermalStorage::fComputeStorageLevel(double storedMass)
{
	double storageLevel = (storedMass / MS_DENSITY) / (M_PI*pow(0.5*_diameterOfStorage, 2.));
	return storageLevel;
}

double ThermalStorage::fInitialStorageMass(int timeInterval)
{
	double timeInSeconds = timeInterval * 60;
	double availableMass = _storedMass + timeInSeconds* _inputHTF->get_massFlow();
	return availableMass;
}

void ThermalStorage::set_storage(double mass, double temperature)
{
	double volume = mass / MS_DENSITY;
	_storedTemperature = temperature;
	if (volume <= M_PI*pow(_diameterOfStorage / 2., 2.)*_heightOfStorage)
	{
		_storedMass = mass;
		_heightOfVolumeStored = volume / (M_PI*pow(_diameterOfStorage / 2., 2.));
	}
	else{
		_storedMass = M_PI*pow(_diameterOfStorage / 2., 2.)*_heightOfStorage*MS_DENSITY;
		_heightOfVolumeStored = _heightOfStorage;
	}

	_outputHTF->set_temperature(_storedTemperature);
}

void ThermalStorage::set_storage2(double level, double temperature)
{
	double area = pow(_diameterOfStorage / 2., 2.)*M_PI;
	double volume = area * level;
	double mass = volume * MS_DENSITY;
	_storedMass = mass;
	_storedTemperature = temperature;
	_heightOfVolumeStored = level;
	_outputHTF->set_temperature(_storedTemperature);
}


double ThermalStorage::fComputeEnergyLosses(double T, double H)
{
	double Q_loss, Q_out_wet, Q_out_floor, Q_out_rad;
	double T_bottom = 90. + 273; //K
	double k_0 = 0.043; //W/mK
	double k_1 = 1.3*pow(10., -4); //W/mK^2
	double k_insul = k_0 + k_1*(((T_bottom + T) / 2.) - 273); //W/mK
	double A_floor = M_PI*pow(_diameterOfStorage / 2., 2.);
	double k_ss = SS_COND; //W/mK (conductivity of stainless steel)
	double t_ss = SS_THICKNESS; //stainless steel thickness in m
	double t_insul = _thicknessOfInsulation;
	int count = 0, count2 = 0;

	//Calculating floor losses --------------------------------------------------------------
	//using default outside temperature of 90 degC
	double R_floor = (t_ss / (k_ss*A_floor)) + (t_insul / (k_insul*A_floor));

	Q_out_floor = (T - T_bottom) / R_floor;

	//calculating wetted wall losses --------------------------------------------------------
	//1- determine outflow by neglecting the outside radiations
	//2- determine the surface temperature needed to get the corresponding convection rate
	//3- compute the amount of radiation flux with this temperature
	//4- determine the surface temperature if the remaining heat flux goes to convection
	//5- cycle until convergence.

	double Re_wall = WIND_VELOCITY * (_diameterOfStorage + t_ss + t_insul) / AIR_VISCOSITY;
	double Pr = 0.707; //1 atm 30 deg C
	double Nu_wall;

	//Hilpert's correlation for Nu, finding h
	double C, m;
	if (Re_wall >= 4. && Re_wall < 40.)
        { C = 0.911; m = 0.385; }
	else if (Re_wall >= 40. && Re_wall < 4000.)
        { C = 0.683; m = 0.466;}
	else if (Re_wall >= 4000. && Re_wall < 40000.)
        { C = 0.193; m = 0.618; }
	else if (Re_wall >= 40000. /*  should be upperbound  */)
        { C = 0.027; m = 0.805; }
    else
    {
        throw 10;
    }
        
	Nu_wall = C*pow(Re_wall, m)*pow(Pr, 1. / 3.);
	double h_wall = AIR_CONDUCTIVITY*Nu_wall / (_diameterOfStorage + 2.*t_ss + 2.*t_insul);
	double A_wet = 2.*M_PI*(0.5*_diameterOfStorage + t_ss + t_insul)*H;

	k_insul = k_0 + k_1*(((T_ATM + T) / 2.) - 273); //W/mK

	//The overall heat transfer coefficient including outside surface convection
	double UA = (M_PI*_diameterOfStorage*H) / (
		(0.5*_diameterOfStorage / k_ss)* log((0.5*_diameterOfStorage + t_ss) / (0.5*_diameterOfStorage))
		+(0.5*_diameterOfStorage / k_insul)* log((0.5*_diameterOfStorage + t_ss + t_insul) / (0.5*_diameterOfStorage + t_ss))
		+((0.5*_diameterOfStorage) / ((0.5*_diameterOfStorage + t_ss + t_insul)*h_wall))
		);

	//The overall heat transfer coefficient for conduction through the tank wall
	double UA_cond =  (M_PI*_diameterOfStorage*H) / (
		(0.5*_diameterOfStorage / k_ss)*log((0.5*_diameterOfStorage + t_ss) / (0.5*_diameterOfStorage))
		+(0.5*_diameterOfStorage / k_insul)*log((0.5*_diameterOfStorage + t_ss + t_insul) / (0.5*_diameterOfStorage + t_ss))
		);

	//The heat transfer resistance for convection 
	double R_conv_wall = 1/(h_wall*A_wet);
	
	//First approximation of heat transfer to outer air excluding radiation losses
        Q_out_wet = UA*(T - T_ATM);

	//Defining the outside surface temperature needed for Q_out_wet to go out through convection
	double T_surf_wet = T -  Q_out_wet*UA_cond;

	double q_wet_1=0.;
	double q_wet_2 = Q_out_wet;
	double k_rad_wet = BOLTZMANN*EPSILON_OUT*A_wet;
	double k_conv_wet = 1. / R_conv_wall;

	try{
		while (fabs(q_wet_2 - q_wet_1) / q_wet_2 > 0.001 && count < 150)
		{
			q_wet_1 = q_wet_2;

			T_surf_wet = fSolveForT(k_rad_wet, k_conv_wet, T, T_ATM, q_wet_1, 0.01);

			k_insul = k_0 + k_1*(((T_surf_wet + T) / 2.) - 273); //W/mK

			UA_cond = (M_PI*_diameterOfStorage*H) / (
				(0.5*_diameterOfStorage / k_ss)*log((0.5*_diameterOfStorage + t_ss) / (0.5*_diameterOfStorage))
				+ (0.5*_diameterOfStorage / k_insul)*log((0.5*_diameterOfStorage + t_ss + t_insul) / (0.5*_diameterOfStorage + t_ss))
				);

			q_wet_2 = (T - T_surf_wet)*UA_cond;
			++count;
		}

		if (count >= 150){
			logic_error noConvergence("couldn't converge to storage external surface temperature for wetted section.");
			throw noConvergence;
		}
	}
	catch (logic_error& e){
		//Setting to preliminary value
		q_wet_2 = Q_out_wet;
	}

	Q_out_wet = q_wet_2;

	//Calculating dry wall losses through radiation from the molten salt surface to the inner tank
	//An iterative procedure is used
	//T_sur_t, inner top surface temperature
	//T_sur_w, inner dry wall surface temperature
	//T_out_t, outter top surface temperature
	//T_out_w, outter dry wall surface temperature
	//--------------------------------------------------
	//1- Assume T_sur_t = T_sur_w = T_ms when starting
	//2- find q_dot = (T_sur - T_a) / R_tot for both surface
	//   R_tot includes R_cond and R_conv
	//3- Knowing q_dot, find the exterior surface temperature necessary if all the heat flux
	//   was occuring through convection only.
	//   T_out = q_dot*R_conv + T_a
	//4- With this outter surface temperature, find radations emitted by the surface
	//   q_dot_rad = epsilon*BOLTZMAN*A*(T_out^4 - T_a^4)
	//5- With this surface temperature, compute the amount of energy transfered from the inner
	//   surface to the outer surface by conduction. q_dot_cond = (T_sur - T_out)/R_cond
	//6- Now the only energy left for convection would be q_dot_conv = q_dot_cond - q_dot_rad
	//7- go back to 3- using q_dot = q_dot_conv Do this until convergence of T_out
	//---------------------------------------------------
	//8- After convergence of the outer temperature and finding out heat flux to the outside,
	//   find the inner surface temperature necessary for this heat flux to happen via radiation
	//   from the inner tank. For this, use the resistance analogy for radiative heat transfer
	//   inside enclosed surface presented in example 13.3 of [incropera].
	//
	//   For this, solve the linear system for the J`s of each surface for each surface
	//  ----------------------------------------------------------------------------------------------------------------------------------
	//   Jt                      Jw                 Jm 
	//   (AmFmt + AtFtw)         (-AtFtw)           (-AmFmt)                                 = q_dot_top
	//   (-AtFtw)                (AmFmw + AtFtw)    (-AmFmw)                                 = q_dot_wall
	//   (-AmFmt)                (-AmFmw)           ((Eps_m A_m)/(1- Eps_m) + AmFmt + AmFmw) = (BOLTZMANN*T_ms^4)/((1- Eps_m)/Eps_m A_m)
	//  ----------------------------------------------------------------------------------------------------------------------------------
	//   After solving this system, use the values obtained for the J`s to compute inner surface temperatures
	//   q_dot = (BOLTZMANN* T^4 - J)/((1-Eps)/Eps*A) Solving for T
	//   
	//9- Using the temperatures found for the inner surface, Return to step 2- if temperature
	//   has changed by more than a given tolerance.

	double H_dry = _heightOfStorage - H;
	double Q_top_cond;
	double A_dry_wall = M_PI*_diameterOfStorage*H_dry;
	double A_top = M_PI*pow(0.5*_diameterOfStorage, 2.);
	double A_w_out = 2.*M_PI*(_heightOfStorage - H)*(0.5*_diameterOfStorage + t_ss + t_insul);

	double S = 2 +pow(H_dry/(0.5*_diameterOfStorage),2.);
	double F_ms_top = 0.5*(S - sqrt(pow(S,2.)-4));
	double F_ms_wall = 1-F_ms_top;
	double F_top_wall = F_ms_wall;
	// double F_wall_top = 0.5*(A_top/A_dry_wall);
	// double F_wall_ms = F_wall_top;
	// double F_wall_wall= 1 - (F_wall_top + F_wall_ms);
	double q_t_cond_1, q_t_cond_2, q_w_cond_1, q_w_cond_2;
	// double q_w_conv, q_t_conv, q_m_rad, q_t_rad;
	// double q_t_rad_1;  q_w_rad;
	// double q_t_rad_2;
	// double q_w_rad_1;
	// double q_w_rad_2;
	double J_m = 0.0, J_t, J_w;
	double h_top, Nu_top = 0.0, Re_top;
	// double x_turb;
	double T_o_t, T_o_w;
	double T_i_t_1, T_i_t_2, T_i_w_1, T_i_w_2;
	double R_w_tot, R_w_conv, R_w_cond, R_t_tot, R_t_conv, R_t_cond;
	double a11, a12, a13, b1;
	double a21, a22, a23, b2;
	double a31, a32, a33, b3;
	double a11_T, a12_T, a13_T;
	double a21_T, a22_T, a23_T;
	double a31_T, a32_T, a33_T;
	double adj11, adj12, adj13;
	double adj21, adj22, adj23;
	double adj31, adj32, adj33;
	double inv11, inv12, inv13;
	double inv21, inv22, inv23;
	double inv31, inv32, inv33;
	double det_A;
	double k_rad_t;
	double k_rad_w;
	double k_conv_t;
	double k_conv_w;
	bool enabler = false;

	if (H_dry > 0.01)
	{

		T_i_t_1 = T - 5.;
		T_i_w_1 = T - 5.;
		T_i_t_2 = 0.;
		T_i_w_2 = 0.;

		//calculating convection coefficient for top surface
		Re_top = WIND_VELOCITY * (_diameterOfStorage + t_ss + t_insul) / AIR_VISCOSITY;
		if (Re_top <= 5. * pow(10., 5.)) { Nu_top = 0.664*pow(Re_top, 1. / 2.)*pow(Pr, 1. / 3.); }
		if (Re_top > 5.*pow(10., 5.) && Re_top <= pow(10., 7.)){ Nu_top = 0.037*pow(Re_top, 4. / 5.)*pow(Pr, 1. / 3.); }
		h_top = Nu_top*AIR_CONDUCTIVITY / _diameterOfStorage;

		//Compute thermal resistances values
		R_w_tot = pow(
			(M_PI*_diameterOfStorage*(_heightOfStorage - H)) / (
			(0.5*_diameterOfStorage / k_ss)*log((0.5*_diameterOfStorage + t_ss) / (0.5*_diameterOfStorage))
			+ (0.5*_diameterOfStorage / k_insul)*log((0.5*_diameterOfStorage + t_ss + t_insul) / (0.5*_diameterOfStorage + t_ss))
			+ ((0.5*_diameterOfStorage) / ((0.5*_diameterOfStorage + t_ss + t_insul)*h_wall))
			)
			, -1.);

		R_w_cond = pow(
			(M_PI*_diameterOfStorage*(_heightOfStorage - H)) / (
			(0.5*_diameterOfStorage / k_ss)*log((0.5*_diameterOfStorage + t_ss) / (0.5*_diameterOfStorage)) +
			+(0.5*_diameterOfStorage / k_insul)*log((0.5*_diameterOfStorage + t_ss + t_insul) / (0.5*_diameterOfStorage + t_ss))
			)
			, -1.);

		R_w_conv = 1 / (h_wall * A_w_out);

		//***insert calculation for h_top

		R_t_tot = R_floor + 1. / (h_top*A_floor);
		R_t_cond = (t_ss / (k_ss*A_floor)) + (t_insul / (k_insul*A_floor));
		R_t_conv = 1. / (h_top*A_floor);

		k_rad_t = BOLTZMANN*EPSILON_SS*A_top;
		k_rad_w = BOLTZMANN*EPSILON_SS*A_w_out;
		k_conv_t = 1. / R_t_conv;
		k_conv_w = 1. / R_w_conv;

		//Find initial values for heat transfer
		/*q_t_cond_2 = (T_i_t_1 - T_ATM) / R_t_tot;
		q_w_cond_2 = (T_i_w_1 - T_ATM) / R_w_tot;*/
		q_t_cond_2 = (T_i_t_1 - T_ATM) / R_t_cond;
		q_w_cond_2 = (T_i_w_1 - T_ATM) / R_w_cond;
		q_t_cond_1 = 0.;
		q_w_cond_1 = 0.;

		//Suppose that all heat is dissipated through convection, what T_o_t is needed
		/*T_o_t = q_t_cond_2*R_t_conv + T_ATM;
		T_o_w = q_w_cond_2*R_w_conv + T_ATM;*/
		T_o_t = T_ATM + 10.;
		T_o_w = T_ATM + 10.;

		//definition of the linear system for radiosity J using the thermal resistance system
		//for three surface.
		//solve AJ = b
		//Define all elements of A according to equations listed in the document
		a11 = A_floor*F_ms_top + A_top*F_top_wall /* ((EPSILON*A_floor)/(1. - EPSILON))*/;
		a12 = -A_top*F_top_wall;
		a13 = -A_floor*F_ms_top;
		a21 = -A_top*F_top_wall;
		a22 = A_floor*F_ms_wall + A_top*F_top_wall /* ((EPSILON*A_dry_wall)/(1. - EPSILON))*/;
		a23 = -A_floor*F_ms_wall;
		a31 = -A_floor*F_ms_top;
		a32 = -A_floor*F_ms_wall;
		a33 = (EPSILON_MS * A_floor / (1. - EPSILON_MS)) + A_floor*F_ms_top + A_floor*F_ms_wall;
		b3 = BOLTZMANN * pow(T, 4.) / ((1. - EPSILON_MS) / (EPSILON_MS*A_floor));

		//find the determinant of A
		det_A = a11*(a22*a33 - a32*a23)
			- a12*(a21*a33 - a31*a23)
			+ a13*(a21*a32 - a31*a22);

		//Use the following procedure to find the inverted matrix A^-1
		// inv(A) = (1/det_A)*adj(A)
		a11_T = a11; a12_T = a21; a13_T = a31;
		a21_T = a12; a22_T = a22; a23_T = a32;
		a31_T = a13; a32_T = a23; a33_T = a33;

		adj11 = (a22_T*a33_T - a32_T*a23_T);
		adj12 = -(a21_T*a33_T - a31_T*a23_T);
		adj13 = a21_T*a32_T - a22_T*a31_T;
		adj21 = -(a12_T*a33_T - a31_T*a13_T);
		adj22 = a11_T*a33_T - a31_T*a13_T;
		adj23 = -(a11_T*a32_T - a31_T*a12_T);
		adj31 = a12_T*a23_T - a22_T*a13_T;
		adj32 = -(a11_T*a23_T - a21_T*a13_T);
		adj33 = a11_T*a22_T - a21_T*a12_T;

		inv11 = adj11 / det_A; inv12 = adj12 / det_A; inv13 = adj13 / det_A;
		inv21 = adj21 / det_A; inv22 = adj22 / det_A; inv23 = adj23 / det_A;
		inv31 = adj31 / det_A; inv32 = adj32 / det_A; inv33 = adj33 / det_A;

		try{
			count = 0;
			while (fabs(T_i_t_2 - T_i_t_1) / T_i_t_1 >= 0.001
				|| fabs(T_i_w_2 - T_i_w_1) / T_i_w_1 >= 0.001
				|| enabler != true)
			{
				enabler = true;
				T_i_t_2 = T_i_t_1;
				T_i_w_2 = T_i_w_1;

				////loop for top values
				count2 = 0;
				try{
					while (fabs(q_t_cond_2 - q_t_cond_1) / q_t_cond_2 > 0.001)
					{
						q_t_cond_1 = q_t_cond_2;

						T_o_t = fSolveForT(k_rad_t, k_conv_t, T_i_t_1, T_ATM, q_t_cond_1, 0.01);

						k_insul = k_0 + k_1*(((T_i_t_1 + T_o_t) / 2.) - 273); //W/mK

						R_t_cond = (t_ss / (k_ss*A_floor)) + (t_insul / (k_insul*A_floor));

						q_t_cond_2 = (T_i_t_1 - T_o_t) / R_t_cond;
						++count2;

						if (count2 >= 150){
							logic_error noConvergence("Can't converge on heat loss through top surface of storage");
							throw noConvergence;
						}
					}
				}
				catch(logic_error& e){
					//Setting temperature to worst case
					T_o_t = (T_ATM + T)/2.;
					q_t_cond_2 = (T_i_t_1 - T_ATM) / R_t_cond;
				}

				b1 = /*(BOLTZMANN*pow(T_i_t_1,4.)*(EPSILON*A_floor)/(1. - EPSILON))*/ -q_t_cond_2;

				//loop for wall values

				try{
					count2 = 0;
					while (fabs(q_w_cond_2 - q_w_cond_1) / q_w_cond_2 > 0.001)
					{
						q_w_cond_1 = q_w_cond_2;

						T_o_w = fSolveForT(k_rad_w, k_conv_w, T_i_w_1, T_ATM, q_w_cond_1, 0.01);

						k_insul = k_0 + k_1*(((T_i_w_1 + T_o_w) / 2.) - 273); //W/mK

						R_w_cond = pow(
							(M_PI*_diameterOfStorage*(_heightOfStorage - H)) / (
							(0.5*_diameterOfStorage / k_ss)*log((0.5*_diameterOfStorage + t_ss) / (0.5*_diameterOfStorage)) +
							+(0.5*_diameterOfStorage / k_insul)*log((0.5*_diameterOfStorage + t_ss + t_insul) / (0.5*_diameterOfStorage + t_ss))
							)
							, -1.);

						q_w_cond_2 = (T_i_w_1 - T_o_w) / R_w_cond;
						++count2;

						if (count2 >= 150) {
							logic_error noConvergence("Can't converge on heat loss through dry wall surface of storage");
							throw noConvergence;
						}
					}
				}
				catch (logic_error& e){

					//Setting temperature to worst case
					T_o_w = (T_ATM + T) / 2.;;
					q_w_cond_2 = (T_i_w_1 - T_ATM) / R_w_cond;
				}


				b2 = /*(BOLTZMANN*pow(T_i_w_1, 4.)*(EPSILON*A_dry_wall) / (1. - EPSILON))*/ -q_w_cond_2;

				J_t = inv11*b1 + inv12*b2 + inv13*b3;
				J_w = inv21*b1 + inv22*b2 + inv23*b3;
				J_m = inv31*b1 + inv32*b2 + inv33*b3;

				//T_i_t_1 = pow(((/*q_t_cond_2 +*/ (J_t-J_w)*A_top*F_top_wall + (J_t - J_m)*A_top*F_ms_top)*(1 - EPSILON) / (EPSILON*A_top) + J_t)/BOLTZMANN, 0.25);
				//T_i_w_1 = pow(((/*q_w_cond_2 +*/ (J_w - J_t)*A_dry_wall*F_wall_top + (J_w - J_m)*A_dry_wall*F_wall_ms)*(1 - EPSILON) / (EPSILON*A_dry_wall) + J_w)/BOLTZMANN, 0.25);

				if ((J_t*EPSILON_SS*A_top) / (1. - EPSILON_SS) > 0.){
					T_i_t_1 = fSolveForT_i((BOLTZMANN*EPSILON_SS*A_top) / (1. - EPSILON_SS),
						(1. / R_t_cond),
						T,
						T_o_t,
						(J_t*EPSILON_SS*A_top) / (1. - EPSILON_SS),
						0.001);
				}
				else { T_i_t_1 = (T + T_o_t) / 2.; }

				if ((J_w*EPSILON_SS*A_dry_wall) / (1. - EPSILON_SS) > 0.){
					T_i_w_1 = fSolveForT_i((BOLTZMANN*EPSILON_SS*A_dry_wall) / (1. - EPSILON_SS),
						(1. / R_w_cond),
						T,
						T_o_w,
						(J_w*EPSILON_SS*A_dry_wall) / (1. - EPSILON_SS),
						0.001);
				}
				else{ T_i_w_1 = (T + T_o_w) / 2.; }

				q_t_cond_2 = (T_i_t_1 - T_o_t) / R_t_cond;
				q_w_cond_2 = (T_i_w_1 - T_o_w) / R_w_cond;

				++count;

				if (count >= 150){ logic_error noConvergence("Could not converge to values for inside dry surfaces of storage");
				throw noConvergence;
				}
			}
			//Total radiosity from the molten salt surface is obtained with the 
			//end results of the procedure


			//The total radiative losses are obtained using the final value of Jm
			Q_out_rad = (BOLTZMANN*pow(T, 4.) - J_m) / ((1 - EPSILON_MS) / (EPSILON_MS*A_floor));
		}
		catch (logic_error& e)
		{
			//Assuming worst case losses
			J_m = BOLTZMANN*pow((T + T_ATM)/2.,4.);
			Q_out_rad = (BOLTZMANN*pow(T, 4.) - J_m) / (((1 - EPSILON_MS) / (EPSILON_MS*A_floor)) + (1/A_floor) + ((1 - EPSILON_SS)/(EPSILON_SS*A_floor)));
		}

		//Returning the total losses.
		Q_loss = Q_out_rad + Q_out_floor + Q_out_wet;
	}

	//In the case that H_dry = 0, that is, the storage is full, then
	//no radiation losses will be considered and losses through the top
	//will be computed the same way they were computed through the wetted
	//wall.
	else{

		Re_top = WIND_VELOCITY * (_diameterOfStorage + t_ss + t_insul) / AIR_VISCOSITY;
		if (Re_top <= 5. * pow(10., 5.)) { Nu_top = 0.664*pow(Re_top, 1. / 2.)*pow(Pr, 1. / 3.); }
		if (Re_top > 5.*pow(10., 5.) && Re_top <= pow(10., 7.)){ Nu_top = 0.037*pow(Re_top, 4. / 5.)*pow(Pr, 1. / 3.); }
		h_top = Nu_top*AIR_CONDUCTIVITY / _diameterOfStorage;

		A_top = M_PI*pow(_diameterOfStorage / 2., 2.);
		k_insul = k_0 + k_1*(((T_ATM + T) / 2.) - 273);

		R_t_tot = R_floor + 1. / (h_top*A_top);
		R_t_cond = (t_ss / (k_ss*A_top)) + (t_insul / (k_insul*A_top));
		R_t_conv = 1. / (h_top*A_top);

		k_rad_t = BOLTZMANN*EPSILON_SS*A_top;
		k_conv_t = 1. / R_t_conv;

		//First approximation of heat transfer to outer air excluding radiation losses
		Q_top_cond = (T - T_ATM)/R_t_tot;

		//Defining the outside surface temperature needed for Q_out_wet to go out through convection
		double T_o_t = T - Q_out_wet*UA_cond;

		double q_t_cond_1 = 0.;
		double q_t_cond_2 = Q_top_cond;

		try{
			count2 = 0;
			while (fabs(q_t_cond_2 - q_t_cond_1) / q_t_cond_2 > 0.001)
			{
				q_t_cond_1 = q_t_cond_2;

				T_o_t = fSolveForT(k_rad_t, k_conv_t, T, T_ATM, q_t_cond_1, 0.01);

				k_insul = k_0 + k_1*(((T_o_t + T) / 2.) - 273); //W/mK

				R_t_cond = (t_ss / (k_ss*A_top)) + (t_insul / (k_insul*A_top));

				q_t_cond_2 = (T - T_o_t) / R_t_cond;

				++count2;

				if (count2 >= 150){
					logic_error noConvergence("Could not find convergence for conduction losses through top wetted surface.");
					throw noConvergence;
				}
			}
		}
		catch (logic_error& e){

			//Setting to preliminary value
			q_t_cond_2 = Q_top_cond;
		}

		Q_top_cond = q_t_cond_2;

		//Returning the total losses.
		Q_loss = Q_top_cond + Q_out_floor + Q_out_wet;
	}

	return Q_loss;
}

double ThermalStorage::fComputeStorageTemperature(int timeInterval)
{
	double timeInSeconds = timeInterval * 60.;
	double rateOfLosses = 0.;
	double totalEnergy, storedTemperature;
	fComputeStorageLevel();

	//The new temperature is calculated sequentially : consider the mass at the beginning of the time interval
	//  and add all of the mass sent to the storage at once. Then take the output temperature to be that of the thermal 
	//  storage after this modification has been done.

	storedTemperature = (_storedMass * _storedTemperature + timeInSeconds* _inputHTF->get_massFlow() * _inputHTF->get_temperature())
		/ (_storedMass + timeInSeconds * _inputHTF->get_massFlow());

	rateOfLosses = fComputeEnergyLosses(storedTemperature, _heightOfVolumeStored);
	totalEnergy = HEAT_CAPACITY * _storedMass * _storedTemperature - timeInSeconds * rateOfLosses;

	_storedTemperature = totalEnergy / (HEAT_CAPACITY * _storedMass);

	//using current mass and temperature, add the new mass sent from the splitter and using design temperature
	//Compute new temperature
	//Then compute losses/time
	//Using losses/time and given time interval, compute the temperature at the end of the time interval after losses
	//Use end value for the exit temperature (to mixer)

	_outputHTF->set_temperature(_storedTemperature);

	return _storedTemperature;
}

double ThermalStorage::fInitialStorageTemperature(int timeInterval)
{
	int timeInSeconds = timeInterval * 60;
	double storedTemperature = (_storedMass * _storedTemperature + timeInSeconds* _inputHTF->get_massFlow() * _inputHTF->get_temperature())
		/ (_storedMass + timeInSeconds * _inputHTF->get_massFlow());
	return storedTemperature;
}

//This function solves the typical k1*T^4 + k2*T - q = 0 equation using Newton<s method.
double ThermalStorage::fSolveForT(double coef_T4, double coef_T, double T_max, double T_min, double q, double eps)
{
	double T_1, T_2;
	double g_k, Dg_k;
	int count = 0;

	T_1 = 0.;
	T_2 = T_max;
	try{
		while (fabs(T_2 - T_1) > eps && count < 150)
		{
			T_1 = T_2;
			g_k = coef_T4*pow(T_1, 4.) + coef_T*T_1 - (q + coef_T4*pow(T_min, 4.) + coef_T*T_min);
			Dg_k = 4.*coef_T4*pow(T_1, 3.) + coef_T;
			//**** insert throw if Dg_k = 0.
			T_2 = T_1 - (g_k / Dg_k);

			++count;

			if (T_2 < T_min && T_1 < T_min){
				logic_error invalidTemperature("Newton's method could not converge to a valid temperature for storage");
				throw invalidTemperature;
			}
		}

		if (count >= 150){
			logic_error noConvergence("Newton's method could not converge to a valid external temperature for storage.");
			throw noConvergence;
		}
	}
	catch (logic_error &e){
		T_2 = (T_max + T_min) / 2.;
	}
	return T_2;
}

double ThermalStorage::fSolveForT_i(double coef_T4, double coef_T, double T_max, double T_min, double q, double eps)
{
	double T_1, T_2;
	double g_k, Dg_k;
	int counter = 0;
	bool newtonFailed = false;
	double delta_T;
	double g1; // g2;

	T_1 = 0.;
	T_2 = T_max;
	/**DEAD CODE*/
	if (q<0.)
	{
		g1 = 0.;
	}

	while (fabs(T_2 - T_1) > eps)
	{
		T_1 = T_2;
		g_k = q - coef_T4*pow(T_1,4.) - coef_T*(T_1 - T_min);
		Dg_k = - 4.*coef_T4*pow(T_1, 3.) - coef_T;
		//**** insert throw if Dg_k = 0.
		T_2 = T_1 - (g_k / Dg_k);
		++counter;
		if (counter > 150)
		{
			newtonFailed = true;
			break;
		}
	}

	try{
		if (newtonFailed == true || T_2 <= T_min)
		{
			delta_T = (T_max - T_min) * eps;
			T_1 = T_max;
			g1 = coef_T4 * pow(T_1, 4.) - coef_T*(T_1 - T_min);
			while (g1 > q)
			{
				T_2 = T_1;
				T_1 -= delta_T;
				g1 = coef_T4 * pow(T_1, 4.) + coef_T*(T_1 - T_min);
				if (T_1 < T_min) { logic_error invalidTemperature("Could not find valid dry surface temperature for storage.");
				throw invalidTemperature;
				}
			}
			T_2 = (T_1 + T_2) / 2.;
		}
	}
	catch (logic_error& e){
		T_2 = T_min;
	}
	return T_2;
}
