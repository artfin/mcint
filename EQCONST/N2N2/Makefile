# mpi gcc wrapper
MPICCXX = mpic++

# usual compiler
CCXX = g++

F = gfortran

CXXFLAGS = -std=c++11 -lgfortran -lm -O3 

INCLUDE_EIGEN = -I/usr/local/include/eigen3
INCLUDE_HEP = -I/home/artfin/Downloads/hep-mc-0.5/include

POT_SRC_PATH = ../../POT/src/
POT_OBJ_PATH = ../../POT/obj/

objects = pa2n2i.o
OBJ = $(patsubst %,$(POT_OBJ_PATH)%,$(objects))

$(POT_OBJ_PATH)%.o: $(POT_SRC_PATH)%.f 
	@echo [F] $< -- new
	@$(F) -c $< -o $@

$(TEST): $(TEST).cc
	@echo [INFO] Building test executable
	@$(CCXX) $(CXXFLAGS) $(TEST).cc -o $@

# Building target binary
simple: $(OBJ) 
	@echo [INFO] Creating Binary Executable
	@$(CCXX) $(CXXFLAGS) $(TARGET).cc -o $@ $^

hamiltonian.o: hamiltonian.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^-o $@

hamiltonian_wcoriolis.o: hamiltonian_wcoriolis.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) -c $^ -o $@

# Building full equilibrium binary
full: $(OBJ) hamiltonian.o
	@echo [INFO] Building Full Equilibrium Binary
	$(CCXX) $(CXXFLAGS) $(FULL_TARGET).cc $^ -o $@

potential_test: $(OBJ) potential_test.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

full_p: $(OBJ) hamiltonian.o full_p.cpp
	$(MPICCXX) $(CXXFLAGS) $(INCLUDE_HEP) $^ -o $@

full_pw: $(OBJ) hamiltonian_wcoriolis.o full_pw.cpp
	$(MPICCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_HEP) $^ -o $@

clean:
	rm simple || rm full || rm simps 
