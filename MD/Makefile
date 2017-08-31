CCXX = g++
CXXFLAGS = -std=c++11

LDFLAGS = -lm

SRC_PATH = ./src/
OBJ_PATH = ./obj/
INC_PATH = -I ./include

BIN_PATH = ./bin/
TEST_PATH = ./test/

# files to compile
objects = vector.o atom.o

OBJ = $(patsubst %,$(OBJ_PATH),$(objects))

$(OBJ_PATH)%.o: $(SRC_PATH)%.cpp
	@echo [CCXX] $<
	@$(CCXX) $(CXXFLAGS) $(LDFLAGS) $(INC_PATH) -o $@ -c $<

test: test_vector test_atom

test_vector: ./obj/vector.o 
	$(CCXX) $(CXXFLAGS) $< $(INC_PATH) $(TEST_PATH)$@.cpp -o $(TEST_PATH)$@

test_atom: ./obj/vector.o ./obj/atom.o
	$(CCXX) $(CXXFLAGS) $< $(INC_PATH) $(TEST_PATH)$@.cpp -o $(TEST_PATH)$@

