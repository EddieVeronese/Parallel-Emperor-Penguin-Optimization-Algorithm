#!/bin/bash
#PBS -l select=1:ncpus=64:mem=2gb
#PBS -l walltime=0:30:00
#PBS -q short_cpuQ
module load mpich-3.2

cd $PBS_O_WORKDIR

./tsp_epo cities.txt 5000 300