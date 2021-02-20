# This file is run by the command `make`
# with `make` the file starts with the first target

# Directory variables
DIR		= ./
IDIR	= $(DIR)include/
MDIR	= $(DIR)modules/
TDIR	= $(DIR)test/
KDIR	= $(DIR)kissfft/

# variables of flags
CC = gcc
OPTIMIZATION = -O2
# -Wall gcc warnings during compilation
CFLAGS = -Wall $(OPTIMIZATION)
# Linking libraries
LIBS = -Lkissfft -lkissfft-double

# variables with executable files
TEST_LIN = inttest_linearity
TEST_ANAL = inttest_analytical
TEST_UNI = inttest_unitarity
# TEST_GAUSS = test_gaussian_wf
MAIN = main

# variable for the modules
# $(wildcar *) is the safe version of * and means all files in the modules
MODULES = $(wildcard $(MDIR)*.c)

# variables for including
INCLUDE = -I $(IDIR) -I $(KDIR)

# this is the first target all with the depency OBJECT_FILES
# it will look for the depency before running the below command(s)
all:
		# $(CC) $(CFLAGS) $(INCLUDE) $(MODULES) $(MAIN).c $(LIBS) -o $(MAIN).exe
		$(CC) $(CFLAGS) $(INCLUDE) $(MODULES) $(TDIR)$(TEST_LIN).c $(LIBS) -o $(TEST_LIN).exe
		$(CC) $(CFLAGS) $(INCLUDE) $(MODULES) $(TDIR)$(TEST_ANAL).c $(LIBS) -o $(TEST_ANAL).exe
		$(CC) $(CFLAGS) $(INCLUDE) $(MODULES) $(TDIR)$(TEST_UNI).c $(LIBS) -o $(TEST_UNI).exe
		# $(CC) $(CFLAGS) $(INCLUDE) $(MODULES) $(TDIR)$(TEST_GAUSS).c $(LIBS) -o $(TEST_GAUSS).exe

# maybe its better to compile and link seperately
# all: $(OBJECT_FILES)
# 		echo: "Linking: $@ ($(CC))"
# 		$(CC) -o $(TEST_LIN)
#
# $(OBJECT_FILES):
