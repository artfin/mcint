# MPI compiler wrapper
MPICCXX = mpic++

# usual compiler
CCXX = g++

CXXFLAGS = -std=c++11 -O3 -lstdc++ -march=native -mtune=native

INCLUDE_EIGEN = -I /usr/local/include/eigen3
INCLUDE_GSL = -lgsl -lgslcblas

LDFLAGS = -lm

MPI_PROGRAMS_PATH = ./mpi/mpi_diatomics.cpp

###################################################################
# hamiltonians (kinetic parts)
HAM_INCLUDE = -I../HAM/include/
HAM_SRC = ../HAM/src/
HAM_OBJ = ../HAM/obj/

$(HAM_OBJ)co2ar_hamiltonian.o: $(HAM_SRC)co2ar_hamiltonian.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(LDFLAGS) -c $^ -o $@
###################################################################

###################################################################
# potentials
POT_INCLUDE = -I../POT/include/
POT_SRC = ../POT/src/
POT_OBJ = ../POT/obj/

$(POT_OBJ)co2ar_potential_julia.o: $(POT_SRC)co2ar_potential_julia.cpp
	$(CCXX) $(CXXFLAGS) $(LDFLAGS) -c $^ -o $@
###################################################################


##################################################################
# CO2-AR
# alpha ~ 10.0
# co2ar_rconst -- R = const, full-dimensional phase space MCMC 
# 	output: theta, pR, pTheta, Jx, Jy, Jz distibutions
# 		( angular momentum in Cartesian coordinates )
#	hamiltonian consists only of kinetic energy ( NO POTENTIAL PART )
#
# co2ar_rconst_jangles -- R = const, full-dimensional phase space MCMC
# 	output: theta, pR, pTheta, (angular)phi, (angular)theta, J
#		( angular momentum in spherical coordinates )
#	hamiltonian consists only of kinetic energy ( NO POTENTIAL PART )
#
# co2ar_rconst_jangles_icpertraj -- the same as co2ar_rconst_jangles
# 	but output is formatted for trajectory calculation
#
#
# alpha ~ 1.0
# co2ar_bound -- full-dimensional phase space MCMC, distributions for bound dimers ( H < 0 ) 
# 	output: R, theta, pR, pTheta, Jx, Jy, Jz distributions
#  		( angular momentum in Cartesian coordinates )
#  	FULL hamiltonian: kinetic + potential parts
#
#
# alpha ~ 1e-5
# co2ar_rconst_lagrange_binning -- R = const, full-dimensional LAGRANGE phase space MCMC
#		output: theta, rdot, thetadot, omega_x, omega_y, omega_z
# 	lagrangian consists only of kinetic energy ( NO POTENTIAL PART )
#	using gsl-histogram for insitu binning of all lagrange variables
#
# co2ar_rconst_lconst -- R = const, L = const ( rotational moment of CO2 )
#   full-dimensional phase space MCMC
#   output: R, theta, pR, pTheta, Jx, Jy, Jz distributions
#   	( angular momentum in Cartesian coordinates )
#   hamiltonian consists only of kinetic energy ( NO POTENTIAL PART )
#
#

CO2AR_SRC_PATH = src/co2ar/
CO2AR_BIN_PATH = bin/co2ar/

co2ar_rconst_lconst: $(HAM_OBJ)co2ar_hamiltonian.o $(CO2AR_SRC_PATH)co2ar_rconst_lconst.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_GSL) $(LDFLAGS) $^ -o $(CO2AR_BIN_PATH)$@ $(HAM_INCLUDE)

co2ar_rconst_lconst_odyssey: $(HAM_OBJ)co2ar_hamiltonian.o $(CO2AR_SRC_PATH)co2ar_rconst_lconst_odyssey.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_GSL) $(LDFLAGS) $^ -o $(CO2AR_BIN_PATH)$@ $(HAM_INCLUDE)

co2ar_rconst: $(HAM_OBJ)co2ar_hamiltonian.o $(CO2AR_SRC_PATH)co2ar_rconst_rep.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(LDFLAGS) $^ -o $(CO2AR_BIN_PATH)$@ $(HAM_INCLUDE) 

co2ar_rconst_autobinning: $(HAM_OBJ)co2ar_hamiltonian.o $(CO2AR_SRC_PATH)co2ar_rconst_autobinning.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_GSL) $(LDFLAGS) $^ -o $(CO2AR_BIN_PATH)$@ $(HAM_INCLUDE) 

co2ar_rconst_jangles: $(HAM_OBJ)co2ar_hamiltonian.o $(CO2AR_SRC_PATH)co2ar_rconst_jangles.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(LDFLAGS) $^ -o $(CO2AR_BIN_PATH)$@ $(HAM_INCLUDE)

co2ar_rconst_jangles_icpertraj: $(HAM_OBJ)co2ar_hamiltonian.o $(CO2AR_SRC_PATH)co2ar_rconst_jangles_icpertraj.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(LDFLAGS) $^ -o $(CO2AR_BIN_PATH)$@ $(HAM_INCLUDE)

co2ar_bound: $(POT_OBJ)co2ar_potential_julia.o $(HAM_OBJ)co2ar_hamiltonian.o $(CO2AR_SRC_PATH)co2ar_bound.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(LDFLAGS) $^ -o $(CO2AR_BIN_PATH)$@ $(HAM_INCLUDE) $(POT_INCLUDE)

co2ar_bound_binning: $(POT_OBJ)co2ar_potential_julia.o $(HAM_OBJ)co2ar_hamiltonian.o $(CO2AR_SRC_PATH)co2ar_bound_binning.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_GSL) $(LDFLAGS) $^ -o $(CO2AR_BIN_PATH)$@ $(HAM_INCLUDE) $(POT_INCLUDE)

co2ar_bound_icpertraj: $(POT_OBJ)co2ar_potential_julia.o $(HAM_OBJ)co2ar_hamiltonian.o $(CO2AR_SRC_PATH)co2ar_bound_icpertraj.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(LDFLAGS) $^ -o $(CO2AR_BIN_PATH)$@ $(HAM_INCLUDE) $(POT_INCLUDE)

co2ar_rconst_lagrange_binning: $(CO2AR_SRC_PATH)co2ar_rconst_lagrange_binning.cpp 
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_GSL) $(LDFLAGS) $^ -o $(CO2AR_BIN_PATH)$@

##################################################################

##################################################################
# DIATOMICS
#

DIATOMICS_SRC_PATH = src/diatomics/
DIATOMICS_BIN_PATH = bin/diatomics/

# optimal value of alpha: ~12.0 ( t = 300k )
#
# diatomics_rconst -- R = const, full-dimensional MCMC for two atom
# 	version of CO2-Ar (two particles with masses of co2 and ar)
# 	Metropolis-Hastings algorithm outputs points WITH repetition
#

diatomics_rconst: $(DIATOMICS_SRC_PATH)diatomics_rconst_rep.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $^ -o $(DIATOMICS_BIN_PATH)$@

diatomics_rconst_exact: $(DIATOMICS_SRC_PATH)diatomics_rconst_exact.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $(DIATOMICS_BIN_PATH)$@

##################################################################

##################################################################
# TESTS

TEST_SRC_PATH = src/test/
TEST_BIN_PATH = bin/test/

# 1d dimensional Markov chain
1d-gaussian: $(TEST_SRC_PATH)1d-gaussian.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $(TEST_BIN_PATH)$@
##################################################################


2d: 2d.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@


2d-exp-sin: 2d-exp-sin.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@

diatomics_exact: diatomics_exact.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

sine_sampling: sine_sampling.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

mpi_diatomics: $(MPI_PROGRAMS_PATH)mpi_diatomics.cpp
	$(MPICCXX) $(CXXFLAGS) $(STFLAGS) $^ -o $(MPI_PROGRAMS_PATH) $@

