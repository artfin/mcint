CCXX = g++
CXXFLAGS = -std=c++11

INCLUDE_GNUPLOT_IO = -I /home/artfin/Desktop/repos/gnuplot-iostream  

LINK_BOOST = -lboost_iostreams -lboost_system -lm

SRC_PATH = ./src/
OBJ_PATH = ./obj/
INC_PATH = -I ./include/

$(OBJ_PATH)ar_he_pes.o: $(SRC_PATH)ar_he_pes.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH)

$(OBJ_PATH)ar_he_dip.o: $(SRC_PATH)ar_he_dip.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH)

$(OBJ_PATH)ar_he_dip_der.o: $(SRC_PATH)ar_he_dip_der.cpp
	$(CCXX) $(CXXFLAGS) -c $^ -o $@ $(INC_PATH)

OBJ = ar_he_pes.o ar_he_dip.o ar_he_dip_der.o
objects = $(patsubst %,$(OBJ_PATH)%,$(OBJ))

main: $(objects) $(SRC_PATH)main.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_GNUPLOT_IO) $(LINK_BOOST) $^ -o $@ $(INC_PATH)

