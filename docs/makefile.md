## makefile

This project uses a simple `Makefile` to compile the MPI-based Fortran solver and generate the executable:

```
../Execute/FDC_mpi
```

The build system supports:

- Intel MPI Fortran (`mpiifort`) with OpenBLAS (default)
- Generic MPI Fortran (`mpif90`)
- Intel MKL (Math Kernel Library)

---

### Requirements

- MPI-enabled Fortran compiler:
  - `mpiifort` (Intel MPI) or
  - `mpif90`
- BLAS library:
  - OpenBLAS (default build), or
  - Intel MKL (MKL build)
- Linux environment with `make`

For MKL builds, the environment variable `MKLROOT` must be defined (e.g., via Intel oneAPI environment setup).

---

### Build Targets

#### Default Build (Intel MPI + OpenBLAS)

Compiles all `.f90` files using `mpiifort` and links against OpenBLAS.

```bash
make
```

Flags used:
- `-O3` (high optimization)
- `-lopenblas`

---

#### Alternative Build (mpif90)

Uses the generic MPI Fortran wrapper:

```bash
make mpif90
```

Additional flag:
- `-ffree-line-length-none` (removes free-form line length limits)

---

#### Intel MKL Build

Compiles and links against Intel MKL:

```bash
make mkl
```

Linked libraries:
- `-lmkl_rt`
- `-lpthread`
- `-lm`
- `-ldl`

This option is recommended for optimized numerical performance on Intel architectures.

---

### Cleaning the Build

Removes all generated files:

```bash
make clean
```

This deletes:
- Object files (`*.o`)
- Module files (`*.mod`)
- The executable (`../Execute/FDC_mpi`)

---

### Notes

- The Makefile does not implement automatic dependency tracking.
- All source files are recompiled at each build.
- The system is intentionally minimal for portability and simplicity.