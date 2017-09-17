# MPI compiler wrapper
MPICCXX = mpic++

# usual compiler
CCXX = g++

CXXFLAGS = -std=c++11 -O3 -lstdc++ -march=native -mtune=native

STFLAGS = -I /usr/local/include/eigen3
LDFLAGS = -lm

MPI_PROGRAMS_PATH = ./mpi/

all: simple 2d diatomics co2ar 2d-exp-sin diatomics_exact co2ar_angles sine_sampling mpi_diatomics 

simple: simple.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

2d: 2d.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@

diatomics: diatomics.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@

co2ar: co2ar.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@

co2ar_angles: co2ar_angles.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@

2d-exp-sin: 2d-exp-sin.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@

diatomics_exact: diatomics_exact.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

sine_sampling: sine_sampling.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

mpi_diatomics: $(MPI_PROGRAMS_PATH)mpi_diatomics.cpp
	$(MPICCXX) $(CXXFLAGS) $(STFLAGS) $^ -o $(MPI_PROGRAMS_PATH)$@

