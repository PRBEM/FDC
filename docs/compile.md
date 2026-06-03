## GLIC

- Make sure you have an account in `GLIC`
- SSH into `GLIC`
- `cd` into `Src/` in the cloned FDC repository, then:
    ```bash
    module load iimpi
    module load OpenBLAS FlexiBLAS
    make all
    ```
## MP machines

- SSH into `MP` machine
- `cd` into `Src/` in the cloned FDC repository, then:
    ```bash
    mpif90 Module_Constants.f90 Module_Variables.f90 -llapack -lblas -O3 -c
    mpif90 *.f90 -o main -I/usr/lib/x86_64-linux-gnu/openmpi/include -ffree-line-length-none -L/usr/lib/ -l:libblas.so.3 -l:liblapack.so.3
    ```