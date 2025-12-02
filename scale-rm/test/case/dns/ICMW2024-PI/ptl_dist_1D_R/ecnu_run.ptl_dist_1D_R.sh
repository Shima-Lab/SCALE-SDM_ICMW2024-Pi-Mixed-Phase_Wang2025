#! /bin/sh -x
#PBS -q small
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -N ptl_dist_1D

cd ${PBS_O_WORKDIR}
module load intel/2015.6.233 zlib/1.2.8 xz/5.2.4 pcre/8.40 bzip2/1.0.6 openssl/1.1.1a curl/7.63.0 R/3.4.3

# run
Rscript ptl_dist_1D.netcdf.R
for i in *.pdf
do
    convert -density 400 $i ${i/pdf/png}
done
for i in mass_dens_drop num_dens_amsul term_vel_drop xz_SD
do
    convert -delay 20 -loop 0 ${i}.*.png ${i}.gif
done

################################################################################
