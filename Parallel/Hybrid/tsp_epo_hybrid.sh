#!/bin/bash
#PBS -l select=1:ncpus=32:mem=2gb
#PBS -l walltime=0:30:00
#PBS -q short_cpuQ
module load mpich-3.2
module load gcc91

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=2

mpirun.actual -n 1 ./tsp_epo_hybrid cities.txt 10000 300