CCXX = g++
CXXFLAGS = -std=c++11

STFLAGS = -I /usr/local/include/eigen3
LDFLAGS = -lm

all: simple 2d

simple: simple.cpp
	$(CCXX) $(CXXFLAGS) $^ -o $@

2d: 2d.cpp
	$(CCXX) $(CXXFLAGS) $(STFLAGS) $(LDFLAGS) $^ -o $@
