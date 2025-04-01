#!/bin/bash
#PBS -l select=1:ncpus=8:mem=2gb
#PBS -l walltime=6:00:00
#PBS -q short_cpuQ
module load mpich-3.2

cd $PBS_O_WORKDIR

mpirun.actual -n 8 ./tsp_epo_mpi