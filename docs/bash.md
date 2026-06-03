## Bash Script (GLIC)

The bash script that runs the entire FDC looks like this (only needed for `GLIC`):

```bash
#!/bin/bash                     
#SBATCH --nodes=1              
#SBATCH --ntasks-per-node=20   
#SBATCH --time=2-9:50          
#SBATCH --array=0-23           

module load iimpi                                              
MACHINEFILE="mpinodes.$SLURM_JOB_ID"                            
srun -l /bin/hostname | sort -n | cut -f2 -d':' > $MACHINEFILE  

mpirun -np $SLURM_NTASKS -machinefile $MACHINEFILE ./FDC_mpi ../Inputs/LB_MLTs/params_MLT$SLURM_ARRAY_TASK_ID.txt

```

| Code | What it does |
|------|--------------|
| #SBATCH --nodes=1 | Requests 1 compute node |
| #SBATCH --ntasks-per-node=20 | Runs 20 MPI tasks on that node |  *(removed "This asks to")*
| #SBATCH --time=2-9:50 | Limits runtime to 2 days, 9 hours, 50 minutes |  *(removed "This limits the max run time to")*
| #SBATCH --array=0-23 | Runs 24 jobs with `SLURM_ARRAY_TASK_ID` from 0-23 |  *(capitalized, removed "separate")*
| srun -l /bin/hostname \| sort -n \| cut -f2 -d':' > $MACHINEFILE | Extracts and stores hostnames | 
| ./FDC_mpi | Fortran MPI executable |

### Running on GLIC
```bash
cd Execute/  # Where the .sh file is stored 
sbatch name_of_file.sh
```
This will start the running of the code. To check status in `GLIC`:
```bash
sacct -S YYYY-MM-DD -u username
```
Replace `YYYY-MM-DD` with the date

**Output**  
After execution, you will these files in the `Execute/`.
- `*.plt`: The calculated diffusion coefficients.
- `*.txt`: The log file of execution.
- `*.out` and `slurm*`: Files containing details of the node and execution in `GLIC`.

## MP Machine
To run FDC on an MP machine, use this bash script from `Src/`:
```bash
mpirun -np 14 ./main ../Inputs/<input_file_name>.txt
```

| Code | What it does |
|------|------------|
| `mpirun` | MPI launcher used to start parallel programs. |
| `-np 14` | Run the program with 14 MPI processes. |
| `./main` | Executable of FDC. |
| `../Inputs/<input_file_name>.txt` |  Path to the input file |

**Output**  
After execution, you will find these files in the `Src/`:
- `*.plt`: The calculated diffusion coefficients.
- `*.txt`: The log file of execution.