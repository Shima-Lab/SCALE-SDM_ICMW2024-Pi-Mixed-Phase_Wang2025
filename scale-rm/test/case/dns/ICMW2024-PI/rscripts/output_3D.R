############################################################################
###
### Program to output 3D data for ICMW-PI
### (for multiple run)
###
############################################################################
###
### History: 200320, S.Takahashi,   initital version
###          211223, S.Takahashi,   calculate by MPI, by one step
###
############################################################################

############################################################################
########## Loading Libraries
library(ncdf4)

############################################################################
########## Parameters
TIME_ITVL <- 30    # time interval of output file(s)
BORDER <- 0.125    # border of wall inside and outside(m)
EPSILON <- 0.001   # mean eddy dissipation rate(m^2/s^3)
RESULT_FILE <- "./LES_Nxxx_3D" # result file name

############################################################################
########## Read Directories
allrundirs <- dir("../",pattern="^run_[0-9][0-9]")
firstdir <- paste("../",allrundirs[1],"/",sep="")

############################################################################
########## Read Data from Files
##### Make a list of files,  mpiranks
allfiles = dir(firstdir,pattern="^history.")
tmp = strsplit(allfiles,"\\history.pe|\\.nc")
allmpiranks = unique(matrix(unlist(tmp),nrow=2)[2,])
names(allfiles) = allmpiranks
MPINUM <- length(allmpiranks)

##### Open the first file
ncin <- nc_open(paste(firstdir,allfiles[1],sep=""))

##### Grids
IMAX <- length(ncvar_get(ncin,"CX"))-4
JMAX <- length(ncvar_get(ncin,"CY"))-4
KMAX <- length(ncvar_get(ncin,"CZ"))-4
CDX <- ncvar_get(ncin,"CDX")
CDY <- ncvar_get(ncin,"CDY")
CDZ <- ncvar_get(ncin,"CDZ")
GRID_VOL <- CDX[3]*CDY[3]*CDZ[3]

##### Area informaion
CXG <- ncvar_get(ncin,"CXG")
CYG <- ncvar_get(ncin,"CYG")
CZ  <- ncvar_get(ncin,"CZ")
XNUM <- length(CXG)-4
YNUM <- length(CYG)-4
ZNUM <- length(CZ)-4
XMAX <- CXG[3] + CXG[XNUM+2]
YMAX <- CYG[3] + CYG[YNUM+2]
ZMAX <- CZ[3] + CZ[ZNUM+2]
XDAT <- CXG[3:(XNUM+2)]
YDAT <- CYG[3:(YNUM+2)]
ZDAT <-  CZ[3:(ZNUM+2)]

##### Border information
IIN <- c(1:XNUM)[(XDAT > BORDER) & ( XDAT < (XMAX-BORDER))]
JIN <- c(1:YNUM)[(YDAT > BORDER) & ( YDAT < (YMAX-BORDER))]
KIN <- c(1:ZNUM)[(ZDAT > BORDER) & ( ZDAT < (ZMAX-BORDER))]

##### Number of processors
XPRC <- XNUM %/% IMAX
YPRC <- YNUM %/% JMAX

##### Close
nc_close(ncin)

##### AllTimes
alltimes <- list(NULL, NULL)
names(alltimes) <- c('time','rundir')
for(rundir in allrundirs){
    ncin <- nc_open(paste("../",rundir,"/",allfiles[1],sep=""))
    time <- ncvar_get(ncin,"time")
    time_units <- ncatt_get(ncin,"time","units")
    nc_close(ncin)
    alltimes$time <- append(alltimes$time, time)
    for(t in time){
      alltimes$rundir <- append(alltimes$rundir, rundir)
    }
}
alltimes$time[alltimes$time %% TIME_ITVL != 0] <- NA
outtimes <- alltimes$time[!is.na(alltimes$time)]

########## Definition of NetCDF file
xdim    <- ncdim_def("x",    "m", XDAT)
ydim    <- ncdim_def("y",    "m", YDAT)
zdim    <- ncdim_def("z",    "m", ZDAT)
timedim <- ncdim_def("time", "seconds", outtimes, unlim=TRUE)

fillvalue <- -9.9999e+30
Tvar    <- ncvar_def("T", "K", list(xdim, ydim, zdim, timedim), fillvalue,
            longname="Temperature")
QVvar   <- ncvar_def("Qv", "kg/kg", list(xdim, ydim, zdim, timedim), fillvalue,
            longname="Water vapor mixing ratio")
Svar    <- ncvar_def("S", "", list(xdim, ydim, zdim, timedim), fillvalue,
            longname="Supersaturation ratio")
LWCvar   <- ncvar_def("LWC", "kgm^-3", list(xdim, ydim, zdim, timedim), fillvalue,
            longname="Liuquid water content")
RHvar   <- ncvar_def("RH", "%", list(xdim, ydim, zdim, timedim), fillvalue,
            longname="Relative humidity")

ncdflist <- list(Tvar, QVvar, Svar, LWCvar, RHvar)
outfile <- nc_create(filename=RESULT_FILE, force_v4=TRUE, ncdflist)

############################################################################
########## Functions
##### Open History files
openHistoryFiles <- function(rundir) {

  histfiles <- list(length(allmpiranks))
  for(mpirank in allmpiranks){
    histfiles[[as.numeric(mpirank)+1]] <-
      nc_open(paste("../", rundir, "/", allfiles[mpirank], sep=""))
  }
  return(histfiles)
}

##### Close History files
closeHistoryFiles <- function(histfiles) {

  for(mpirank in allmpiranks){
    nc_close(histfiles[[as.numeric(mpirank)+1]])
  }
}

##### Read History File
readHistoryFile <- function(histfiles, mpirank, idxtime, pname) {

  read <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], pname,
            c(1,1,1,idxtime), c(IMAX,JMAX,KMAX,1))
  if(!setequal(dim(read[,,]),c(IMAX,JMAX,KMAX))){
    message("ERROR: Size of variables is not consistent.mpirank=",mpirank)
    q()
  }
  return(read)
}

############################################################################
########## Main Process
cat("Start output_3D.R\n")

rundir_prv <- ""
i <- 0

for(t in 1:length(alltimes$time)) {

  if(is.na(alltimes$time[t])) {
    rundir_prv <- alltimes$rundir[t]
    next 
  }
  i <- i + 1

  # Control History Files
  if(rundir_prv != alltimes$rundir[t]) {
    if(rundir_prv != "") {
      closeHistoryFiles(histfiles)
    }
    cat(sprintf("Processing on rundir = %s\n",rundir))
    histfiles <- openHistoryFiles(rundir)
  }
  rundir_prv = alltimes$rundir[t]

  cat(sprintf("Processing on time = %d\n",outtimes[i]))

  # Read history files
  for(mpirank in allmpiranks){
    T <- readHistoryFile(histfiles, mpirank, t, "T")
    QV <- readHistoryFile(histfiles, mpirank, t, "QV")
    S <- readHistoryFile(histfiles, mpirank, t, "RH") / 100.0 - 1.0
    QC <- readHistoryFile(histfiles, mpirank, t, "QC_sd")
    QR <- readHistoryFile(histfiles, mpirank, t, "QR_sd")
    LWC <- (QC + QR) * readHistoryFile(histfiles, mpirank, t, "DENS")
    RH <- readHistoryFile(histfiles, mpirank, t, "RH")

    # Set position info
    is <- (as.numeric(mpirank) %% XPRC) * IMAX + 1
    js <- (as.numeric(mpirank) %/% XPRC) * JMAX + 1
    tryCatch(
    {
      ncvar_put(outfile, Tvar, T, c(is,js,1,i), c(IMAX,JMAX,KMAX,1))
      ncvar_put(outfile, QVvar, QV, c(is,js,1,i), c(IMAX,JMAX,KMAX,1))
      ncvar_put(outfile, Svar, S, c(is,js,1,i), c(IMAX,JMAX,KMAX,1))
      ncvar_put(outfile, LWCvar, LWC, c(is,js,1,i), c(IMAX,JMAX,KMAX,1))
      ncvar_put(outfile, RHvar, RH, c(is,js,1,i), c(IMAX,JMAX,KMAX,1))
    },
    error = function(e) {
      message(e) 
      nc_close(outfile)
      closeHistoryFiles(histfiles) 
      q()
    },
      silent = TRUE
    )
  }
}
closeHistoryFiles(histfiles)

# Write other data
tryCatch(
  {
    ncatt_put(outfile, 0, "author", "Shima-Lab.")
  },
  error = function(e) {
    message(e) 
  },
  finally = {
    nc_close(outfile)
  },
  silent = TRUE
)
cat("End output_3D.R\n")
