CCXX = g++
CXXFLAGS = -std=c++11

INCLUDE_EIGEN = -I /usr/local/include/eigen3
INCLUDE_GSL = -lgsl -lgslcblas

mpi_diatomics: mpi_diatomics.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_GSL) $^ -o $@
