#!/bin/bash
#PBS -l select=1:ncpus=32:mem=2gb
#PBS -l walltime=5:00:00
#PBS -q short_cpuQ
module load mpich-3.2
module load gcc91

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=32

./tsp_epo_open cities.txt 5000 300

./tsp_epo_open cities.txt 7500 300

./tsp_epo_open cities.txt 10000 300

./tsp_epo_open cities.txt 12500 300

./tsp_epo_open cities.txt 15000 300


