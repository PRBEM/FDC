#!/bin/tcsh
#
# LSF batch script to run an MPI application
#
#BSUB -P UCLA0017        # project code
#BSUB -W 11:59              # wall-clock time (hrs:mins)
#BSUB -n 64                  # number of tasks in job         
#BSUB -R "span[ptile=16]"    # run 16 MPI tasks per node
#BSUB -J FDC            # job name
#BSUB -o FDC.%J.out     # output file name in which %J is replaced by the job ID
#BSUB -e FDC.%J.err     # error file name in which %J is replaced by the job ID
#BSUB -q regular         # queue
#BSUB -N
#BSUB -B
#run the executable
mpirun.lsf ./FDC dummyinput.txt
