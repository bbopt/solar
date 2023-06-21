EXE                = ../bin/solar
COMPILATOR         = g++
COMPILATOR_OPTIONS = -O3 -Wall -Wextra -pedantic -Wcast-align
LIBS               = -lc -lm
INCLUDE            = -I.
COMPILE            = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
OBJS               = main.o CentralReceiver.o Time_Manager.o Economics.o Global.o Turbine.o \
                     HeatExchanger.o Heliostat.o HeliostatField.o HtfCycle.o helpFunctions.o sampling.o Random_Pickup.o RNG.o \
                     MoltenSalt.o Powerblock.o Powerplant.o Problem.o Scenario.o Sunray.o Sun.o ThermalStorage.o \
                     Clock.o Evaluator.o

$(EXE): $(OBJS)
	$(COMPILATOR) -o $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS)

CentralReceiver.o: CentralReceiver.cpp CentralReceiver.hpp MoltenSalt.o Constants.hpp
	$(COMPILE) CentralReceiver.cpp

Time_Manager.o: Time_Manager.cpp Time_Manager.hpp
	$(COMPILE) Time_Manager.cpp

Economics.o: Economics.cpp Economics.hpp Global.o
	$(COMPILE) Economics.cpp

Global.o: Global.cpp Global.hpp
	$(COMPILE) Global.cpp

HeatExchanger.o: HeatExchanger.cpp HeatExchanger.hpp MoltenSalt.o Powerblock.o Powerblock.o Turbine.o Global.o Constants.hpp
	$(COMPILE) HeatExchanger.cpp

Heliostat.o: Heliostat.cpp Heliostat.hpp Sun.hpp Constants.hpp
	$(COMPILE) Heliostat.cpp

HeliostatField.o: HeliostatField.cpp HeliostatField.hpp Sun.o Heliostat.o Global.o
	$(COMPILE) HeliostatField.cpp

helpFunctions.o: helpFunctions.hpp helpFunctions.cpp Scenario.o
	$(COMPILE) helpFunctions.cpp

HtfCycle.o: HtfCycle.cpp HtfCycle.hpp CentralReceiver.o ThermalStorage.o HeatExchanger.o MoltenSalt.o Global.o Constants.hpp
	$(COMPILE) HtfCycle.cpp

main.o: main.cpp sampling.o Evaluator.o
	$(COMPILE) main.cpp

sampling.o: sampling.cpp sampling.hpp Evaluator.o Random_Pickup.o
	$(COMPILE) sampling.cpp

Evaluator.o: Evaluator.cpp Evaluator.hpp Problem.o helpFunctions.o Global.o Clock.o
	$(COMPILE) Evaluator.cpp

Clock.o: Clock.cpp Clock.hpp
	$(COMPILE) Clock.cpp

MoltenSalt.o: MoltenSalt.cpp MoltenSalt.hpp Constants.hpp
	$(COMPILE) MoltenSalt.cpp

Powerblock.o: Powerblock.cpp Powerblock.hpp Turbine.o Global.o
	$(COMPILE) Powerblock.cpp

Powerplant.o: Powerplant.cpp Powerplant.hpp HtfCycle.o HeliostatField.o Economics.o Powerblock.o Constants.hpp Global.o Time_Manager.o
	$(COMPILE) Powerplant.cpp

Problem.o: Problem.cpp Problem.hpp Global.o
	$(COMPILE) Problem.cpp

Scenario.o: Scenario.cpp Scenario.hpp Global.o Powerplant.o Problem.o
	$(COMPILE) Scenario.cpp

Random_Pickup.o: Random_Pickup.cpp Random_Pickup.hpp RNG.o
	$(COMPILE) Random_Pickup.cpp

RNG.o: RNG.cpp RNG.hpp
	$(COMPILE) RNG.cpp

Sun.o: Sun.cpp Sun.hpp Time_Manager.o Sunray.o RNG.o
	$(COMPILE) Sun.cpp

Sunray.o: Sunray.cpp Sunray.hpp Heliostat.o
	$(COMPILE) Sunray.cpp

ThermalStorage.o: ThermalStorage.cpp ThermalStorage.hpp Constants.hpp MoltenSalt.o Global.o
	$(COMPILE) ThermalStorage.cpp

Turbine.o: Turbine.cpp Turbine.hpp Constants.hpp Global.o
	$(COMPILE) Turbine.cpp

clean:
	@echo "   cleaning obj files"
	@rm -f $(OBJS)

del:
	@echo "   cleaning trash files"
	@rm -f core *~
	@echo "   cleaning obj files"
	@rm -f $(OBJS)
	@echo "   cleaning exe file"
	@rm -f $(EXE)
