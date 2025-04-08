#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=3:00:00
#PBS -q short_cpuQ
module load mpich-3.2

cd $PBS_O_WORKDIR

mpirun.actual -n 1 ./tsp_epo_mpi cities.txt 10000 300