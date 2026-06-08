# Makefile for Full Diffusion Code (FDC)
#
# Compilation order:
#   1. Module files must be compiled before any file that uses them:
#        Module_Constants.f90  -> required by almost all files
#        Module_Variables.f90  -> required by Roots_finder, WeightFunction
#        Module_Dispersion.f90 -> required by WeightFunction, int_N_omega
#   2. All remaining .f90 files compiled together.
#
# Note: *.f90 in step 2 will re-compile the module files, but since their
# .mod files already exist this is harmless. A proper dependency-tracked
# Makefile would avoid this; see comment at bottom for a future improvement.

CC     = mpiifort
CCf90  = mpif90
TARGET = ../Execute/FDC_mpi

# Modules that must be compiled first, in order
MODULES = Module_Constants.f90 Module_Variables.f90 Module_Dispersion.f90

all:
	$(CC) $(MODULES) -lopenblas -O3 -c
	$(CC) *.f90 -lopenblas -O3 -o $(TARGET)

mpif90:
	$(CCf90) $(MODULES) -lopenblas -O3 -ffree-line-length-none -c
	$(CCf90) *.f90 -lopenblas -O3 -ffree-line-length-none -o $(TARGET)

# Lib MKL 2019.0 / OS Linux / Compiler Fortran / Arch Intel(R) 64
# Single Dynamic Library
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
mkl:
	$(CC) $(MODULES) -O3 -c \
		-I${MKLROOT}/include \
		-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl
	$(CC) *.f90 -O3 -o $(TARGET) \
		-I${MKLROOT}/include \
		-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl

clean:
	rm -f *.mod *.o $(TARGET)

