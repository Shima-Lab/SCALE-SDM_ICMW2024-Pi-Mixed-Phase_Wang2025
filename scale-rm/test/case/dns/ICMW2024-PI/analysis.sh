( cd analyse-env_R/
Rscript analyse-env.R
convert -rotate 90 -alpha off -density 400 skewTlogP.eps skewTlogP.png
cd ..

cd QHYD_2Dplot
rm -f *.png *.pdf *.gif
Rscript QHYD_2Dplot.R
for i in *.pdf ; do convert -density 200 $i ${i/pdf/png} ; done
for i in QHYD_overlay ; do convert -delay 20 -loop 0 ${i}.*.png ${i}.gif ; done
cd .. ) &

( cd ptl_dist_1D_R/
rm -f *.png *.pdf *.gif numbers_???.txt LES_Nxxx_DSD LES2_Nxxx_?SD_?????
Rscript ptl_dist_1D.history.R
for i in *.pdf; do convert -flatten -density 200 $i ${i/pdf/png}; done
for i in *.00000.png ; do var=${i//.00000.png} ; echo ${var} ; convert  -delay 20 -loop 0 ${var}.*.png ${var}.gif ; done
cd ..

ln -s . run_00
cd rscripts
rm -f LES_Nxxx_*
for i in *.R; do Rscript $i || break ; done
for i in *.pdf; do convert -flatten -density 200 $i ${i/pdf/png}; done
cd .. ) &

wait
