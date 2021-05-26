#!/bin/bash
#SBATCH --job-name=hello-world # Job name for tracking
#SBATCH --partition=cs402 # Partition you wish to use
#SBATCH --nodes=1 # Number of compute nodes used by your job
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2 # Number of CPU cores used on each node
#SBATCH --time=00:10:00
module load cs402-mpi
mpirun ./karman -o out.bin -v 2
gprof karman > out.gprof
./diffbin karman.original.bin out.bin