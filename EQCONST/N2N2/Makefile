# mpi gcc wrapper
MPICCXX = mpic++

# usual compiler
CCXX = g++

F = gfortran

CXXFLAGS = -std=c++11 -lgfortran -I/home/artfin/Downloads/hep-mc-0.5/include -lm -I /usr/local/include/eigen3/

TARGET = simple
FULL_TARGET = full

POT_SRC_PATH = ../../POT/src/
POT_OBJ_PATH = ../../POT/obj/

objects = pa2n2i.o
OBJ = $(patsubst %,$(POT_OBJ_PATH)%,$(objects))

all: FULL_TARGET TARGET full_parallel

$(POT_OBJ_PATH)%.o: $(POT_SRC_PATH)%.f 
	@echo [F] $< -- new
	@$(F) -c $< -o $@

$(TEST): $(TEST).cc
	@echo [INFO] Building test executable
	@$(CCXX) $(CXXFLAGS) $(TEST).cc -o $@

# Building target binary
$(TARGET): $(OBJ) 
	@echo [INFO] Creating Binary Executable
	@$(CCXX) $(CXXFLAGS) $(TARGET).cc -o $@ $^

hamiltonian.o: hamiltonian.cpp
	$(CCXX) $(CXXFLAGS) -c hamiltonian.cpp -o $@

# Building full equilibrium binary
$(FULL_TARGET): $(OBJ) hamiltonian.o
	@echo [INFO] Building Full Equilibrium Binary
	$(CCXX) $(CXXFLAGS) $(FULL_TARGET).cc $^ -o $@

potential_test: $(OBJ) potential_test.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

full_parallel: $(OBJ) hamiltonian.o full_parallel.cpp
	$(MPICCXX) $(CXXFLAGS) $^ -o $@

clean:
	rm simple || rm full || rm simps 
