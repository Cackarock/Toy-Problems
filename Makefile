# -------------------------------------------------------
# Makefile for CS530 programming assignments
# Author: Ravi Shankar, Kenslee Moy, Floriana Ciaglia, Dylan Leman
# -------------------------------------------------------

# Usage
# -------------------------------------------------------
# make all         - makes all targets in bin/ 
# make test        - makes and runs all the unit tests
# make clean       - removes all files generated by make.


# Please tweak the following variable definitions as needed by your
# project, except GTEST_HEADERS, which you can use in your own targets
# but shouldn't modify.

# Points to the root of Google Test, relative to where this file is.
GTEST_DIR = deps/googletest/googletest

# Where to find user code, relative to where this file is.
SRC = src
TEST = test

# Other directories, relative to where this file is.
BIN = bin
OBJ = obj
LIB = lib
INCL = include

# Flags passed to the preprocessor
##################################
# Set Google Test's header directory as a system directory, so that
# the compiler doesn't generate warnings in Google Test headers.
#
# -isystem dir: Mark the directory 'dir'as a system directory, so that it
#               gets the same special treatment that is applied to the 
#               standard system directories
CPPFLAGS += -isystem $(GTEST_DIR)/include

# Flags passed to the C and C++ compiler
########################################
# -std=c++11: Set standard to C++11
#
# -g:         Produce debugging information in the operating system’s 
#             native format
#
# -Wall:      Enables all warning flags 
#             (https://gcc.gnu.org/onlinedocs/gcc/Warning-Options.html)
#
# -Wextra:    This enables extra warning flags that are not enabled by -Wall
#
# -pthread:   Link with the POSIX threads library
#
# -Idir:      Add the directory 'dir' to the list of directories to be searched 
#             for header files during preprocessing
#
# -llib:      Search the library named 'lib' when linking
CXXFLAGS += -std=c++11 -g -Wall -Wextra -Ideps -pthread
CFLAGS += -g -Wall -Wextra -Ideps -std=gnu99


# All targets
all: pi-mc mm 

debug: CFLAGS += -DDEBUG

debug: clean all

# All tests produced by this Makefile.  Remember to add new tests you
# created to the list.
TESTS = pi-mc-test 

# Runs all tests
test: $(TESTS)
	-$(BIN)/pi-mc-test
	-./test/matrix-test.sh matrix-test-results


# All Google Test headers.  Usually you shouldn't change this definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

# Builds gtest.a and gtest_main.a.
# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
$(OBJ)/gtest-all.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest-all.cc -o $@

$(OBJ)/gtest_main.o : $(GTEST_SRCS_)
	$(CXX) $(CPPFLAGS) -I$(GTEST_DIR) $(CXXFLAGS) -c \
            $(GTEST_DIR)/src/gtest_main.cc -o $@

$(LIB)/gtest.a : $(OBJ)/gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

$(LIB)/gtest_main.a : $(OBJ)/gtest-all.o $(OBJ)/gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

# Builds a/all test/s.  
# A test should link with either gtest.a or gtest_main.a, 
# depending on whether it defines its own main() function.
pi-mc-test: $(OBJ)/pi-mc.o \
		      $(OBJ)/pi-mc-test.o \
		      $(LIB)/gtest_main.a
	mpicxx $(CPPFLAGS) $(CXXFLAGS) $^ -o $(BIN)/$@


# Builds all unit tests
$(OBJ)/%-test.o: $(TEST)/%-test.cpp
	mpicxx $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $^

# pi-mc target
pi-mc: $(OBJ)/pi-mc-parallel-main.o $(OBJ)/pi-mc.o
	mpicc $(CFLAGS) -o $(BIN)/$@ $^ -lm

$(OBJ)/pi-mc-parallel-main.o: $(SRC)/pi-mc/pi-mc-parallel-main.c 
	mpicc $(CFLAGS) -c -o $@ $^ 

$(OBJ)/pi-mc.o: $(SRC)/pi-mc/pi-mc.c
	mpicc $(PFLAG) -c -o $@ $^ $(CFLAGS) 

# mm target
mm: $(OBJ)/parallel-mm.o $(OBJ)/mm.o $(OBJ)/mmio.o 
	mpicc $(CFLAGS) -o $(BIN)/$@ $^ -lm

$(OBJ)/parallel-mm.o: $(SRC)/mm/parallel-mm.c
	mpicc $(CFLAGS) -c -o $@ $^ 

$(OBJ)/mm.o: $(SRC)/mm/mm.c 
	mpicc $(CFLAGS) -c -o $@ $^  

$(OBJ)/mmio.o: $(INCL)/mmio.c
	$(CC) $(CFLAGS) -c -o $@ $^

# Builds all object files
#$(OBJ)/%.o: $(SRC)/%.c
#	$(CXX) $(PFLAG) -c -o $@ $^ $(CFLAGS)


# Clean up when done. 
# Removes all object, library and executable files
clean:
	rm -f $(BIN)/*
	rm -f $(OBJ)/*
	rm -f $(LIB)/*
	rm -f $(SRC)/*.o
	rm -f matrix-test-results
	rm -f test/test-output/*