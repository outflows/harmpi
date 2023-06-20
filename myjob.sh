#!/bin/bash
#SBATCH --partition=SP3
#SBATCH -J harmpi-teste1
#SBATCH --nodes=20
#SBATCH --ntasks-per-node=20
#SBATCH --time=168:0:00	
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gustavo.rodrigues.soares@usp.br

#run the application:
cd /scratch/6435163/harmpi/
mpirun -n 400 ./harm 10 10 4
