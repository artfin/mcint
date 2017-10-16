# mpi gcc wrapper
MPICCXX = mpic++

# usual compiler
CCXX = g++

F = gfortran

CXXFLAGS = -std=c++11 -lm -Wall -O3 -march=native -mtune=native 
LINK_GFORTRAN = -lgfortran

INCLUDE_HEP = -I/home/artfin/Downloads/hep-mc-0.5/include -I/n/home11/artfin/hep-mc-0.5/include
INCLUDE_EIGEN = -I/usr/local/include/eigen3

SRC_PATH = ./src/

##########################################################################

HAM_SRC_PATH = ../../HAM/src/
HAM_OBJ_PATH = ../../HAM/obj/
HAM_INCLUDE_PATH = -I../../HAM/include/

$(HAM_OBJ_PATH)n2n2_hamiltonian.o: $(HAM_SRC_PATH)n2n2_hamiltonian.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@ 

$(HAM_OBJ_PATH)n2n2_hamiltonian_0_25.o: $(HAM_SRC_PATH)n2n2_hamiltonian_0_25.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@ 

$(HAM_OBJ_PATH)n2n2_hamiltonian_0_5.o: $(HAM_SRC_PATH)n2n2_hamiltonian_0_5.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@ 

$(HAM_OBJ_PATH)n2n2_hamiltonian_2_0.o: $(HAM_SRC_PATH)n2n2_hamiltonian_2_0.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@ 

$(HAM_OBJ_PATH)n2n2_hamiltonian_4_0.o: $(HAM_SRC_PATH)n2n2_hamiltonian_4_0.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@ 

##########################################################################


##########################################################################

POT_SRC_PATH = ../../POT/src/
POT_OBJ_PATH = ../../POT/obj/

$(POT_OBJ_PATH)pa2n2i.o: $(POT_SRC_PATH)pa2n2i.f
	@$(F) -c $^ -o $@

##########################################################################

full_parallel: $(SRC_PATH)full_parallel.cpp $(HAM_OBJ_PATH)n2n2_hamiltonian.o $(POT_OBJ_PATH)pa2n2i.o
	$(MPICCXX) $(CXXFLAGS) $(LINK_GFORTRAN) $(HAM_INCLUDE_PATH) $(INCLUDE_HEP) $(INCLUDE_EIGEN) $^ -o $@

##########################################################################


##########################################################################
# different lengths of N2 length

full_0_25: $(SRC_PATH)diff_lengths/full_parallel_0_25.cpp $(HAM_OBJ_PATH)n2n2_hamiltonian_0_25.o $(POT_OBJ_PATH)pa2n2i.o
	$(MPICCXX) $(CXXFLAGS) $(LINK_GFORTRAN) $(HAM_INCLUDE_PATH) $(INCLUDE_HEP) $(INCLUDE_EIGEN) $^ -o $@

full_0_5: $(SRC_PATH)diff_lengths/full_parallel_0_5.cpp $(HAM_OBJ_PATH)n2n2_hamiltonian_0_5.o $(POT_OBJ_PATH)pa2n2i.o
	$(MPICCXX) $(CXXFLAGS) $(LINK_GFORTRAN) $(HAM_INCLUDE_PATH) $(INCLUDE_HEP) $(INCLUDE_EIGEN) $^ -o $@

full_2_0: $(SRC_PATH)diff_lengths/full_parallel_2_0.cpp $(HAM_OBJ_PATH)n2n2_hamiltonian_2_0.o $(POT_OBJ_PATH)pa2n2i.o
	$(MPICCXX) $(CXXFLAGS) $(LINK_GFORTRAN) $(HAM_INCLUDE_PATH) $(INCLUDE_HEP) $(INCLUDE_EIGEN) $^ -o $@

full_4_0: $(SRC_PATH)diff_lengths/full_parallel_4_0.cpp $(HAM_OBJ_PATH)n2n2_hamiltonian_4_0.o $(POT_OBJ_PATH)pa2n2i.o
	$(MPICCXX) $(CXXFLAGS) $(LINK_GFORTRAN) $(HAM_INCLUDE_PATH) $(INCLUDE_HEP) $(INCLUDE_EIGEN) $^ -o $@

##########################################################################



