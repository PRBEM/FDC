
#BSUB -a mvapich
#BSUB -n 64
#BSUB -J FDC_2018_UB_Chorus_Recalculate
#BSUB -q qintel
#BSUB -e err%J.mpi
#BSUB -o out%J.mpi

module load intel/intelmpi_2019u5
mpirun /home/dedong/bin/FDC_2019/full-diffusion-code-f90/Execute/FDC_mpi /home/dedong/bin/FDC_2019/full-diffusion-code-f90/Inputs/2018UBChorus.txt 
