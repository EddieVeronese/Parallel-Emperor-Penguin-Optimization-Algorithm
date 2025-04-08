#!/bin/bash
#PBS -l select=1:ncpus=32:mem=2gb
#PBS -l walltime=3:00:00
#PBS -q short_cpuQ
module load gcc91

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=32

./tsp_epo_open cities.txt 10000 300