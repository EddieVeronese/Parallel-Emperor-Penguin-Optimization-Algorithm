#!/bin/bash
#PBS -l select=1:ncpus=32:mem=2gb
#PBS -l walltime=6:00:00
#PBS -q short_cpuQ
module load mpich-3.2
module load gcc

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=4

./tsp_epo_openmp