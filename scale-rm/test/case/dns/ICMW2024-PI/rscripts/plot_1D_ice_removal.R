############################################################################
###
### Program to plot 1D graph for ICMW-PI
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
files <- dir("./", pattern="^LES_N[0-9x]+_1D_iremoval$")
if(length(files) <= 0){
  cat("ERROR: Output LES_Nxxx_1D files not found.\n")
  quit()
}

parameters <- c("N_ice_removal", "M_ice_removal")

ncin <- nc_open(files[1])
time <- ncvar_get(ncin,"time")
time_units <- ncatt_get(ncin,"time","units")
nc_close(ncin)

for(prm in parameters){
  cat(sprintf("Plotting 1D graph of %s.\n", prm))

  for(f in files){
    cat(sprintf(" Input file: %s\n", f))

    ncin <- nc_open(f)
    data <- ncvar_get(ncin, prm)
    units <- ncatt_get(ncin, prm, "units")
    dimnames(data)[[1]] <- time
    if(length(data) != length(time)){
      cat("ERROR: Size of variables is not consistent.\n")
      nc_close(ncin)
      quit()
    }
    if(units$hasatt){
      unit <- sprintf("[%s]", units$value)
    } else {
      unit <- ""
    }
    nc_close(ncin)

    pdf(sprintf("%s.%s.pdf", f, prm))
    par(mgp=c(3,1,0))
    plot(
        time,
        data,
        type="l",
        xlab="time[s]",
        ylab=sprintf("%s%s", prm, unit),
        main=sprintf("Domain average of %s", prm)
    )
    dev.off()
  }
  cat("done.\n")
}

