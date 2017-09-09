CCXX = g++
CXXFLAGS = -std=c++11 -O3

STFLAGS = -I /usr/local/include/eigen3
LDFLAGS = -lm

all: simple 2d diatomics co2ar 2d-exp-sin

simple: simple.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

2d: 2d.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@

diatomics: diatomics.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@

co2ar: co2ar.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@

2d-exp-sin: 2d-exp-sin.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@

