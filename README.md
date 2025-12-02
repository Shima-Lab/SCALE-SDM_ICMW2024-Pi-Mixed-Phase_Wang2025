Readme file for SCALE-SDM 5.2.6-2.2.2 for ICMW 2024 Pi Chamber Mixed-Phase Cloud Simulation Case to reproduce the results in Wang et al. (2025).

Corresponding Author: Shin-ichiro Shima (s_shima[at]sim.u-hyogo.ac.jp)

[//]:#(#################################################################################)
# General description
SCALE is a library of weather and climate models of the Earth and planets (Nishizawa et al., 2015; Sato et al., 2015, https://scale.riken.jp/).
SDM is a Lagrangian particle-based cloud microphysics model (Shima et al., 2009, 2020).
In order to run the ICMW 2024 Pi Chamber Mixed-Phase Cloud Simulation Case, we implemented the sidewalls and particle injection feature to SCALE-SDM 5.2.6-2.2.0. 

[//]:#(#################################################################################)
# Required software and supported environment
Fortran and C compiler are required to compile SCALE-SDM. 
MPI, NetCDF4, and HDF5 libraries are are also required.

The numerical experiments were conducted by using Intel Fortran/C compiler, Intel MPI, HDF5, and NetCDF.
For data analysis, R was used. 

[//]:#(#################################################################################)
# Set environment variable
```
$ export SCALE_SYS=Linux64-intel-impi
```

[//]:#(#################################################################################)
# Set compiler options
In the top directly, you can find a directory named sysdep.
Edit sysdep/Makedef.${SCALE_SYS}, and modify the compiler options according to your system.

For example, 
```
$ cd sysdep
$ vi Makedef.Linux64-intel-impi
----------------
...
-FFLAGS_FAST  = -fpp -m64 -O3 -xHost                 \
+FFLAGS_FAST  = -fpp -m64 -O3 -xCORE-AVX512 -parallel -qopenmp \
...
-CFLAGS_FAST  = -O3 -xHost -ip -ftz -mcmodel=medium -shared-intel
+CFLAGS_FAST  = -O3 -xCORE-AVX512 -parallel -qopenmp -ip -ftz -mcmodel=medium -shared-intel
...
----------------
```

[//]:#(#################################################################################)
# Clean
```
$ module load intel impi hdf5 netcdf netcdf-fortran
$ cd scale-rm/test/case/dns/ICMW-PI/
$ make allclean SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T
```

[//]:#(#################################################################################)
# Compile
```
$ module load intel impi hdf5 netcdf netcdf-fortran
$ cd scale-rm/test/case/dns/ICMW-PI/
$ time make -j5 SCALE_ENABLE_SDM=T SCALE_DISABLE_LOCALBIN=T
$ ln -fsv  `grep ^TOPDIR Makefile | sed s/\)//g | awk '{print $NF}'`/bin/scale-rm* .
```

[//]:#(#################################################################################)
# Run a batch job
```
$ vi grand_run.sh
#!/bin/sh
#------ pjsub option --------#
#PJM -L rscgrp=la
#PJM -L node=7
#PJM --mpi proc=256
#PJM -L elapse=48:00:00
#PJM -g r22775
#PJM -j
#------- Program execution -------#
module load intel impi hdf5 netcdf netcdf-fortran

# run
mpiexec.hydra -n ${PJM_MPI_PROC} ./scale-rm_init init.conf || exit
mpiexec.hydra -n ${PJM_MPI_PROC} ./scale-rm run.conf || exit
```
```
$ pjsub grand_run.sh
$ pjstat
$ pjstat --rsc -b
$ pjshowrsc --rscunit suba
$ tail -f LOG.pe000000
```

[//]:#(#################################################################################)
# Skew-T log-P Diagram of the Initial Atmospheric Sounding
```
$ cd analyse-env_R/
$ module purge
$ Rscript analyse-env.R
$ convert -rotate 90 -alpha off -density 400 skewTlogP.eps skewTlogP.png
$ ls
alldata.txt  analyse-env.R  skewTlogP.eps  skewTlogP.png  variable_info.txt
$ less variable_info.txt
$ less alldata.txt
$ evince skewTlogP.eps
$ cd ..
```

[//]:#(#################################################################################)
# Plot the Spatial Structure Using R
```
$ cd QHYD_2Dplot
$ module purge
$ Rscript QHYD_2Dplot.R
$ for i in *.pdf ; do convert -density 200 $i ${i/pdf/png} ; done
$ for i in QHYD_overlay ; do convert -delay 20 -loop 0 ${i}.*.png ${i}.gif ; done
$ cd ..
```

[//]:#(#################################################################################)
# Plot Particle Distribution: Aerosol Size Distribution (Number Density), Droplet Size Distribution (Mass Density), Terminal Velocities of Droplets, Distribution of SPs in x-y space.
```
$ cd ptl_dist_1D_R/
$ rm -f *.png *.pdf *.gif numbers_???.txt
$ Rscript ptl_dist_1D.history.R
$ for i in *.pdf; do convert -flatten -density 200 $i ${i/pdf/png}; done
$ for i in *.00000.png ; do var=${i//.00000.png} ; echo ${var} ; convert  -delay 20 -loop 0 ${var}.*.png ${var}.gif ; done
$ cd ..
```

[//]:#(#################################################################################)
# Save 1D, 2D, and 3D Data Using R
Run this after ptl_dist_1D_R analysis is done.
```
$ ln -s . run_00
$ cd rscripts
$ rm -f LES_Nxxx_*
$ for i in *.R; do Rscript $i || break ; done
$ for i in *.pdf; do convert -flatten -density 200 $i ${i/pdf/png}; done
$ cd ..
```
