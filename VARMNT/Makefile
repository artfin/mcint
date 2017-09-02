CCXX = g++
CXXFLAGS = -std=c++11 -O3 -DEIGEN_N_DEBUG -fopenmp 

STFLAGS = -I /home/artfin/Desktop/repos/gnuplot-iostream -I /usr/local/include/eigen3   
LDFLAGS = -lboost_iostreams -lboost_system -lm

TARGET = harmonic delta_variational guided_walk 2d_delta_var 

harmonic: harmonic.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@	

delta_variational: delta_variational.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@

guided_walk: guided_walk.cpp
	$(CCXX) $(CXXFLAGS) -lm $^ -o $@

2d_delta_var: 2d_delta_var.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@
