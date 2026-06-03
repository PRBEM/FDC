#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=2-9:50
#SBATCH --array=0-23


module load iimpi
MACHINEFILE="mpinodes.$SLURM_JOB_ID"
srun -l /bin/hostname | sort -n | cut -f2 -d':' > $MACHINEFILE
 
mpirun -np $SLURM_NTASKS -machinefile $MACHINEFILE ./FDC_mpi ../Inputs/2026Chorus/LB_20/params_MLT$SLURM_ARRAY_TASK_ID.txt
