CCXX = g++
MPICCXX = mpic++

CXXFLAGS = -std=c++11 -O3  -mtune=native

INCLUDE_HEP = -I /home/artfin/Downloads/hep-mc-0.5/include/
INCLUDE_EIGEN = -I /usr/local/include/eigen3
INCLUDE_GSL = -lgsl -lgslcblas
LDFLAGS = -lm

SRC_DIR = ./src/
OBJ_DIR = ./obj/
INC_DIR = -I ./include/
INC_DIRWO = ./include/

#SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
#OBJ_FILES = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRC_FILES))

$(OBJ_DIR)mcmc.o: $(SRC_DIR)mcmc.cpp $(INC_DIRWO)mcmc.hpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_GSL) $(INC_DIR)  -c $< -o $@

he_ar_integrator: $(SRC_DIR)he_ar_integrator.cpp $(INC_DIRWO)integrator.hpp $(INC_DIRWO)integrand.hpp
	$(MPICCXX) $(CXXFLAGS) $(INCLUDE_HEP) $(INC_DIR) $< -o $@

diatom_exact: $(SRC_DIR)diatom_exact.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_GSL) $(INC_DIR) $< -o $@

he_ar_mcmc: $(SRC_DIR)he_ar_mcmc.cpp $(OBJ_DIR)mcmc.o
	$(CCXX) $(CXXFLAGS) $(INCLUDE_EIGEN) $(INCLUDE_GSL) $(INC_DIR) $^ -o $@

TEST_DIR = ./tests/
tests = $(TEST_DIR)test_limit $(TEST_DIR)test_limits $(TEST_DIR)functor

test: $(tests)
	@echo [Compiling tests]

$(TEST_DIR)test_limit: $(TEST_DIR)limit.cpp ./include/limit.hpp 
	$(CCXX) $(CXXFLAGS) $(INC_DIR) $< -o $@	

$(TEST_DIR)test_limits: $(TEST_DIR)limits.cpp ./include/limits.hpp ./include/limits.hpp 
	$(CCXX) $(CXXFLAGS) $(INC_DIR) $< -o $@

$(TEST_DIR)functor: $(TEST_DIR)functor.cpp
	$(CCXX) $(CXXFLAGS) $(INCLUDE_HEP) $< -o $@
