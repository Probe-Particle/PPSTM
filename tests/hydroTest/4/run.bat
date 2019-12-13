#!/bin/bash
#PBS -N annulene
#PBS -q luna
#PBS -l select=1:ncpus=1:mem=48gb:scratch_local=10gb
#PBS -l walltime=23:59:59
#PBS -o stdout
#PBS -j oe
#PBS -m ae
. /packages/run/modules-2.0/init/sh
#module add openmpi-1.6-intel
module load intelcdk-17
ulimit -s unlimited
#cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=1
/storage/praha1/home/mendieta/progs/fireball.x > out.dat &
