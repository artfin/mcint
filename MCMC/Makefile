CCXX = g++

CXXFLAGS = -std=c++11 -O3  -mtune=native

INCLUDE_EIGEN = -I /usr/local/include/eigen3
INCLUDE_GSL = -lgsl -lgslcblas
LDFLAGS = -lm

SRC_DIR = ./src/
OBJ_DIR = ./obj/
INC_DIR = -I ./include/

#SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
#OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))

$(OBJ_DIR)mcmc.o: $(SRC_DIR)mcmc.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_GSL) $(INC_DIR)  -c $< -o $@

diatom_exact: $(SRC_DIR)diatom_exact.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_GSL) $(INC_DIR) $< -o $@

he_ar: $(SRC_DIR)he_ar.cpp $(OBJ_DIR)mcmc.o
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_GSL) $(INC_DIR) $^ -o $@
