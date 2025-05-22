#!/bin/bash
#PBS -l select=1:ncpus=64:mem=2gb
#PBS -l walltime=5:00:00
#PBS -q short_cpuQ
module load mpich-3.2
module load gcc91

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=4

mpirun.actual -n 1 ./tsp_epo_hybrid cities.txt 15000 300

mpirun.actual -n 2 ./tsp_epo_hybrid cities.txt 15000 300

mpirun.actual -n 4 ./tsp_epo_hybrid cities.txt 15000 300

mpirun.actual -n 8 ./tsp_epo_hybrid cities.txt 15000 300

mpirun.actual -n 16 ./tsp_epo_hybrid cities.txt 15000 300