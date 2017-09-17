MPICC = mpic++

CCXX = g++

CXXFLAGS = -std=c++11 -march=native -mtune=native
STFLAGS = -I/home/artfin/Downloads/hep-mc-0.5/include

all: intro

intro: intro.cpp
	$(MPICC) $(CXXFLAGS) $(STFLAGS) $^ -o $@
