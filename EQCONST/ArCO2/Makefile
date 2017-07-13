.PHONY: all clean

EXECS = main

CCXX = g++
CXXFLAGS = -std=c++11 -I/home/artfin/Downloads/hep-mc-0.5/include -I /usr/local/include/eigen3/ -lm

all: ${EXECS}

main: hamiltonian.o main.cc
	$(CCXX) $(CXXFLAGS) $^ -o $@

hamiltonian.o: hamiltonian.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@

clean: 
	rm -rf *.o
