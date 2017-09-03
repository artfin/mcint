CCXX = g++
F = gfortran

CXXFLAGS = -std=c++11 -lgfortran -I/home/artfin/Downloads/hep-mc-0.5/include -lm -I/Users/mac/Downloads/hep-mc/include -L/usr/local/Cellar/gcc/6.3.0_1/lib/gcc/6/

TARGET = main

POT_SRC_PATH = ../../POT/src/
POT_OBJ_PATH = ../../POT/obj/

objects = pa2n2i.o
OBJ = $(patsubst %,$(POT_OBJ_PATH)%,$(objects))

$(POT_OBJ_PATH)%.o: $(POT_SRC_PATH)%.f 
	@echo [F] $< -- new
	@$(F) -c $< -o $@

# Building target binary
$(TARGET): $(OBJ)
	@echo [INFO] Creating Binary Executable
	@$(CCXX) $(CXXFLAGS) $(TARGET).cc -o $@ $^

