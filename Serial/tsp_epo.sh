#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=0:60:00
#PBS -q short_cpuQ
module load mpich-3.2

cd $PBS_O_WORKDIR

mpirun.actual -n 1 ./tsp_epo