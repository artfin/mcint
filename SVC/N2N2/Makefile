CCXX = g++
F = gfortran

CXXFLAGS = -lgfortran

EXECS = main

main: print_hi.o main.cc
	$(CCXX) $(CXXFLAGS) $^ -o $@

print_hi.o: print_hi.f90
	$(F) -c $^ -o $@
