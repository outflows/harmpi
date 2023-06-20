#!/bin/bash
#SBATCH -N 27
#SBATCH -C haswell
#SBATCH -q regular
#SBATCH -J harmpi1
#SBATCH --mail-user=gustavo.9891@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 300:00:00

#run the application:
srun -n 864 -c 2 --cpu_bind=cores /global/homes/g/gsoares/harmpi/harm 12 9 8
