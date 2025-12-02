############################################################################
###
### Program to plot 2D graph for ICMW-PI
###
############################################################################
###
### History: 210304, S.Takahashi,   initital version
###
############################################################################

############################################################################
########## Loading Libraries
library(ncdf4)

############################################################################
########## Parameters
TIME_ITVL <-  1    # time interval of output file(s)

############################################################################
########## Functions

############################################################################
########## Read Data from Files
files <- dir("./", pattern="^LES_N[0-9x]+_2D$")
if(length(files) <= 0){
  cat("ERROR: Output LES_Nxxx_1D files not found.\n")
  quit()
}

parameters <- c("T", "Qv", "RH", "LWC", "IWC")

ncin <- nc_open(files[1])
time <- ncvar_get(ncin,"time")
time_units <- ncatt_get(ncin,"time","units")
z <- ncvar_get(ncin,"z")
z_units <- ncatt_get(ncin,"z","units")
nc_close(ncin)

#PLOT_TIME <- c("0", "30", "120")  # time for plot
#PLOT_TIME <- as.character(time)
tmp_increment <- 5
PLOT_TIME <- as.character(time[seq(1,length(time),tmp_increment)])

for(pt in PLOT_TIME){ 
  for(prm in parameters){
    cat(sprintf("Plotting 2D graph of %s at t=%s[s]\n", prm, pt))

    for(f in files){
      cat(sprintf(" Input file: %s\n", f))
      ncin <- nc_open(f)
      data <- ncvar_get(ncin, prm)
      units <- ncatt_get(ncin, prm, "units")
      dimnames(data)[[2]] <- time
      if(length(data) != length(z) * length(time)){
        cat("ERROR: Size of variables is not consistent.\n")
        nc_close(ncin)
        quit()
      }
      nc_close(ncin)

      if(units$hasatt){
        unit <- sprintf("[%s]",units$value)
      } else {
        unit <- ""
      }

      pdf(sprintf("%s.%s.%05d.pdf", f, prm, as.integer(pt)))
      par(mgp=c(3,1,0))
      plot(
        data[,pt],
        z,
        type="l",
        main=sprintf("Vertical profile of %s at t=%s", prm, pt),
        xlab=sprintf("%s%s", prm, unit),
        ylab="z[m]"
      )
      dev.off()
    }
    cat("done.\n")
  }
}

