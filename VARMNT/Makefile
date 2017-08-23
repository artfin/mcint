CCXX = g++
CXXFLAGS = -std=c++11 

STFLAGS = -I /home/artfin/Desktop/repos/gnuplot-iostream 
LDFLAGS = -lboost_iostreams -lboost_system -lm

TARGET = harmonic 

$(TARGET): harmonic.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@	
