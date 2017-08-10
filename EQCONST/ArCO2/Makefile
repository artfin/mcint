.PHONY: all clean

EXECS = main precision

CCXX = g++
CXXFLAGS = -std=c++11 -O3 -I/home/artfin/Downloads/hep-mc-0.5/include -I /usr/local/include/eigen3/ -lm

all: ${EXECS}

precision: hamiltonian.o precision.cc
	$(CCXX) $(CXXFLAGS) $^ -o $@

main: hamiltonian.o main.cc
	$(CCXX) $(CXXFLAGS) $^ -o $@

hamiltonian.o: hamiltonian.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@

clean: 
	rm -rf *.o
