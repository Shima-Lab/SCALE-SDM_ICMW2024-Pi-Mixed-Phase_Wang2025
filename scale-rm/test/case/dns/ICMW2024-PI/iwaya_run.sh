#!/bin/bash
#PBS -q M
#PBS -l select=4:ncpus=32:mpiprocs=32
#PBS -l walltime=24:00:00
#PBS -N scale-sdm

source /etc/profile.d/modules.sh
cd ${PBS_O_WORKDIR}
module load intel/2022.3.1 mpt hdf5/1.14.3 netcdf-c/4.9.2 netcdf-fortran/4.6.1

# run
mpiexec_mpt dplace -s1 ./scale-rm_init init.conf || exit
mpiexec_mpt dplace -s1 ./scale-rm run.conf || exit
