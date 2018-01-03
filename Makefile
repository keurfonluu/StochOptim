# ======================================================================
#  Copyright (c) 2017 by MINES ParisTech, France
#  All rights reserved.
#
#  This software is furnished under a license and may be used and copied
#  only in  accordance with  the  terms  of  such  license  and with the
#  inclusion of the above copyright notice. This software or  any  other
#  copies thereof may not be provided or otherwise made available to any
#  other person.  No title to and ownership of  the  software is  hereby
#  transferred.
# ======================================================================


# ======================================================================
# Declarations for compiler
# ======================================================================

# ---- Compiler
COMPILER = gfortran

# ---- GNU compiler
ifeq ($(COMPILER), gfortran)
FC = gfortran
OMP = -fopenmp
MOD = -J
MPI_INC =
MPI_LIB =
GNU = 1
MPI = 0
else ifeq ($(COMPILER), mpif90)
FC = mpif90
OMP = -fopenmp
MOD = -J
MPI_INC =
MPI_LIB =
GNU = 1
MPI = 1
# ---- Intel compiler
else ifeq ($(COMPILER), ifort)
FC = ifort
OMP = -qopenmp
MOD = -module ./
MPI_INC	= #${I_MPI_ROOT}/intel64/include/
MPI_LIB	= #${I_MPI_ROOT}/intel64/lib
GNU = 0
MPI = 0
else ifeq ($(COMPILER), mpiifort)
FC = mpiifort
OMP = -qopenmp
MOD = -module ./
MPI_INC	= #${I_MPI_ROOT}/intel64/include/
MPI_LIB	= #${I_MPI_ROOT}/intel64/lib
GNU = 0
MPI = 1
endif


# ---- Options for GNU compiler
ifeq ($(GNU), 1)
FFLAGS = -O3 -cpp -funroll-loops -fno-protect-parens -specs override-specs.h
# ---- Options for Intel compiler
else ifeq ($(GNU), 0)
FFLAGS	= -O3 -fpp -assume byterecl
# FFLAGS = -O3 -fpp -assume byterecl -fp-model fast=2 -ftz -no-prec-div -no-prec-sqrt -fp-speculation fast -fno-fnalias -fma -ansi-alias -ipo -qopt-report5 -qopt-report-phase=VEC,LOOP -xCORE-AVX512 -qopt-zmm-usage=high -lifcore -lifport -align
endif

ifeq ($(MPI), 1)
FFLAGS += -Ddo_mpi
endif


# ======================================================================
# Declarations of executables to be compiled and various dependances
# ======================================================================
# Name of executable
TARGET1 = optimization.exe
TARGET2 = benchmark.exe
TARGET3 = sampling.exe
TARGET4 = sensitivity.exe
TEST = test.exe

# Directories
SRCDIR = src
OBJDIR = obj
LIBDIR = $(SRCDIR)/lib
TSTDIR = $(SRCDIR)/tests

# Link objects to create executable (tab required on second line)
OBJS1 = $(OBJDIR)/forlab.o \
				$(OBJDIR)/benchmark_functions.o \
				$(OBJDIR)/stochoptim.o \
        $(OBJDIR)/optimization.o

OBJS2 = $(OBJDIR)/forlab.o \
				$(OBJDIR)/benchmark_functions.o \
				$(OBJDIR)/stochoptim.o \
        $(OBJDIR)/benchmark.o

OBJS3 = $(OBJDIR)/forlab.o \
				$(OBJDIR)/benchmark_functions.o \
				$(OBJDIR)/stochoptim.o \
        $(OBJDIR)/sampling.o

OBJS4 = $(OBJDIR)/forlab.o \
				$(OBJDIR)/benchmark_functions.o \
				$(OBJDIR)/stochoptim.o \
        $(OBJDIR)/sensitivity.o

OBJST = $(OBJDIR)/forlab.o \
				$(OBJDIR)/stochoptim.o \
				$(OBJDIR)/test.o

# These routines depend on include file - recompile if include modified
ALL = $(TARGET1) $(TARGET2) $(TARGET3) $(TARGET4)
all: $(ALL)
optimization: $(TARGET1)
benchmark: $(TARGET2)
sampling: $(TARGET3)
sensitivity: $(TARGET4)
test: $(TEST)


# ======================================================================
# General rules, these should not require modification
# General rules for building ".o" objects from fortran programs or subroutines
# ======================================================================
vpath %.o obj
vpath %.mod obj
vpath %.f90 $(LIBDIR) $(TSTDIR)
COMPILE = $(FC) $(FFLAGS) $(OMP) $(MPI_INC) $(MPI_LIB)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: %.f90 | $(OBJDIR)
	$(COMPILE) -c $^ -o $@ $(MOD)$(OBJDIR)

$(TARGET1): $(OBJS1)
	$(COMPILE) -o $@ $(OBJS1)

$(TARGET2): $(OBJS2)
	$(COMPILE) -o $@ $(OBJS2)

$(TARGET3): $(OBJS3)
	$(COMPILE) -o $@ $(OBJS3)

$(TARGET4): $(OBJS4)
	$(COMPILE) -o $@ $(OBJS4)

$(TEST): $(OBJST)
	$(COMPILE) -o $@ $(OBJST)

# Utilities
.PHONY: $(ALL) test clean veryclean

clean:
	rm -rf $(ALL) $(TEST) $(OBJDIR)

veryclean:
	rm -rf $(TARGET1)* $(TARGET2)* $(TARGET3)* $(TARGET4)* $(TEST)* $(OBJDIR) output *.optrpt
