CCXX = g++

F = gfortran -O0

CXXFLAGS = -std=c++11 -lgfortran -O3 -lstdc++ -march=native -mtune=native

INCLUDE_HEP = -I/home/artfin/Downloads/hep-mc-0.5/include/
INCLUDE_EIGEN = -I/usr/local/include/eigen3

LDFLAGS = -lm

###################################################################
# hamiltonians (kinetic parts)
HAM_INCLUDE = -I../HAM/include/
HAM_SRC = ../HAM/src/
HAM_OBJ = ../HAM/obj/

$(HAM_OBJ)co2ar_hamiltonian.o: $(HAM_SRC)co2ar_hamiltonian.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(LDFLAGS) -c $^ -o $@

$(HAM_OBJ)co2h2_hamiltonian.o: $(HAM_SRC)co2h2_hamiltonian.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(LDFLAGS) -c $^ -o $@ $(HAM_INCLUDE)

###################################################################

###################################################################
# potentials
POT_INCLUDE = -I../POT/include/
POT_SRC = ../POT/src/
POT_OBJ = ../POT/obj/

$(POT_OBJ)co2ar_potential_julia.o: $(POT_SRC)co2ar_potential_julia.cpp
	$(CCXX) $(CXXFLAGS) $(LDFLAGS) -c $^ -o $@

$(POT_OBJ)co2h2_pes.o: $(POT_SRC)co2h2_pes.f
	$(F) -c $^ -o $@
###################################################################

###################################################################
# dipolies
DIP_INCLUDE = -I../DIP/include/
DIP_SRC = ../DIP/src/
DIP_OBJ = ../DIP/obj/

$(DIP_OBJ)co2h2_dipole_lr.o: $(DIP_SRC)co2h2_dipole_lr.cpp
	$(CCXX) $(CXXFLAGS) $(LDFLAGS) -c $^ -o $@ $(DIP_INCLUDE)
###################################################################

###################################################################
# CO2-H2
# zeroth_moment.cpp -- calculation of Zeroth Moment for CO2-H2; full phase space integration

CO2H2_SRC_PATH = src/co2h2/
CO2H2_BIN_PATH = bin/co2h2/

zeroth_moment: $(HAM_OBJ)co2h2_hamiltonian.o $(DIP_OBJ)co2h2_dipole_lr.o $(POT_OBJ)co2h2_pes.o $(CO2H2_SRC_PATH)zeroth_moment.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_HEP) $^ -o $(CO2H2_BIN_PATH)$@ $(HAM_INCLUDE) $(POT_INCLUDE) $(DIP_INCLUDE)

###################################################################



