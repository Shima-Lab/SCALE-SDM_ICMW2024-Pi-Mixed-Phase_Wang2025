#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=la
#PJM -L node=7
#PJM --mpi proc=256
#PJM -L elapse=24:00:00
#PJM -g r22775
#PJM -j
#------- Program execution -------#
module load intel impi hdf5 netcdf netcdf-fortran

# run
mpiexec.hydra -n ${PJM_MPI_PROC} ./scale-rm_init init.conf || exit
mpiexec.hydra -n ${PJM_MPI_PROC} ./scale-rm run.conf || exit
