#!/bin/bash
#PBS -l select=1:ncpus=64:mem=2gb
#PBS -l walltime=5:00:00
#PBS -q short_cpuQ
module load mpich-3.2

cd $PBS_O_WORKDIR

mpirun.actual -n 1 ./tsp_epo_mpi cities.txt 15000 300

mpirun.actual -n 2 ./tsp_epo_mpi cities.txt 15000 300

mpirun.actual -n 4 ./tsp_epo_mpi cities.txt 15000 300

mpirun.actual -n 8 ./tsp_epo_mpi cities.txt 15000 300

mpirun.actual -n 16 ./tsp_epo_mpi cities.txt 15000 300

mpirun.actual -n 32 ./tsp_epo_mpi cities.txt 15000 300

mpirun.actual -n 64 ./tsp_epo_mpi cities.txt 15000 300