# Very simple make file without dependencies
CC=mpiifort # compiler
CCf90=mpif90 
TARGET=../Execute/FDC_mpi #make of thetarget

all:
	$(CC) Module_Constants.f90 Module_Variables.f90 -lopenblas -O3 -c 	
	$(CC) *.f90 -lopenblas -O3 -o $(TARGET)
mpif90:
	$(CCf90) Module_Constants.f90 Module_Variables.f90 -lopenblas -O3 -ffree-line-length-none -c 	
	$(CCf90) *.f90 -lopenblas -O3 -ffree-line-length-none -o $(TARGET)
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor
# Lib MKL 2019.0
# OS Linux
# Compiler Fortran
# Arch Inter(R) 64
# Single Dynamic Library
mkl:
	$(CC) Module_Constants.f90 Module_Variables.f90 -O3 -c 	\
		-I${MKLROOT}/include \
		-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl
	$(CC) *.f90 -O3 -o $(TARGET) \
		-I${MKLROOT}/include \
		-L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl

clean:	
	rm *.mod
	rm *.o
	rm $(TARGET)
