############################################################################
###
### Program to output 2D data for ICMW-PI
### (for multiple run)
###
############################################################################
###
### History: 200320, S.Takahashi,   initital version
###          210303, S.Takahashi,   add outputs
###                                 include top  and bottom
###          211222, S.Takahashi,   calculate by MPI, by one step
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
RESULT_FILE <- "./LES_Nxxx_2D" # result file name

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
zdim    <- ncdim_def("z",    "m", ZDAT[1:ZNUM])
timedim <- ncdim_def("time", "seconds", outtimes, unlim=TRUE)

fillvalue <- -9.9999e+30

Tvar    <- ncvar_def("T", "K", list(zdim, timedim), fillvalue,
            longname="Temperature")
QVvar   <- ncvar_def("Qv", "kg/kg", list(zdim, timedim), fillvalue,
            longname="Water vapor mixing ratio")
RHvar   <- ncvar_def("RH", "%", list(zdim, timedim), fillvalue,
            longname="Relative humidity")
LWCvar  <- ncvar_def("LWC","kgm^-3", list(zdim, timedim), fillvalue,
            longname="Liuquid water content")
IWCvar  <- ncvar_def("IWC","kgm^-3", list(zdim, timedim), fillvalue,
            longname="Ice water content")
# S2Svar   <- ncvar_def("Sigma2_S", "1", list(zdim, timedim), fillvalue,
#             longname="Profile of spatial variance of supersaturation")
# S_DROPvar   <- ncvar_def("S_drop", "1", list(zdim, timedim), fillvalue,
#             longname="Supersaturation at droplet locations")
# S2S_DROPvar  <- ncvar_def("Sigma2_S_drop", "1", list(zdim, timedim), fillvalue,
#             longname="The variance of the supersaturation at droplet locations")
# N_DROPvar   <- ncvar_def("N_drop", "cm^-3", list(zdim, timedim), fillvalue,
#             longname="Droplet number concentration")
# R_MEAN1var   <- ncvar_def("r_mean1", "um", list(zdim, timedim), fillvalue,
#             longname="Droplet mean radius")

ncdflist <- list(Tvar, QVvar, RHvar, LWCvar, IWCvar)
     #S2Svar, S_DROPvar, S2S_DROPvar,N_DROPvar, R_MEAN1var)
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

##### Get Domain Average(single parameter)
getDomainAverage <- function(histfiles, idxtime, param) {

  dav <- rep(0.0, length=ZNUM)
  for(mpirank in allmpiranks){
    rd <- readHistoryFile(histfiles, mpirank, idxtime, param)
    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")

    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))

    for(k in c(1:ZNUM)){
      # skip where i and j are both integer(0)
      if(length(i) != 0 && length(j) != 0){
        dav[k] <- dav[k] + sum(rd[i, j, k])
      }
    }
  }
  dav <- dav / (length(IIN) * length(JIN))
  return(dav)
}

##### Get Domain Average of LWC
getDomainAverageOfLWC <- function(histfiles, idxtime) {

  dav <- rep(0.0, length=ZNUM)
  for(mpirank in allmpiranks){
    QC <- readHistoryFile(histfiles, mpirank, idxtime, "QC_sd")
    QR <- readHistoryFile(histfiles, mpirank, idxtime, "QR_sd")
    DENS <- readHistoryFile(histfiles, mpirank, idxtime, "DENS")
    
    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")

    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))

    for(k in c(1:ZNUM)){
    if(length(i) != 0 && length(j) != 0 && length(k) != 0){
      LWC <- (QC + QR) * DENS
      dav[k] <- dav[k] + sum(LWC[i, j, k])* 1000.0
    }
    }
  }
  dav <- dav / (length(IIN) * length(JIN))
  return(dav)
}

#### Get Domain Average of IWC
getDomainAverageOfIWC <- function(histfiles, idxtime) {

  dav <- rep(0.0, length=ZNUM)
  for(mpirank in allmpiranks){
    QC <- readHistoryFile(histfiles, mpirank, idxtime, "QI_sd")
    QR <- readHistoryFile(histfiles, mpirank, idxtime, "QS_sd")
    QG <- readHistoryFile(histfiles, mpirank, idxtime, "QG_sd")
    DENS <- readHistoryFile(histfiles, mpirank, idxtime, "DENS")
    
    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")

    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))

    for(k in c(1:ZNUM)){
    if(length(i) != 0 && length(j) != 0 && length(k) != 0){
      IWC <- (QC + QR + QG) * DENS
      dav[k] <- dav[k] + sum(IWC[i, j, k])* 1000.0
    }
    }
  }
  dav <- dav / (length(IIN) * length(JIN))
  return(dav)
}

##### Get Sum(single parameter)
getSum <- function(histfiles, idxtime, param) {

  sm <- rep(0.0, length=ZNUM)
  for(mpirank in allmpiranks){
    rd <- readHistoryFile(histfiles, mpirank, idxtime, param)
    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")

    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))

    for(k in c(1:ZNUM)){
      # skip where i and j are both integer(0)
      if(length(i) != 0 && length(j) != 0){
        sm[k] <- sm[k] + sum(rd[i, j, k])
      }
    }
  }
  return(sm)
}

##### Get Domain Average of S_drop
# getDomainAverageOfSdrop <- function(histfiles, idxtime, nc_all) {

#   dav <- rep(0.0, length=ZNUM)
#   for(mpirank in allmpiranks){
#     RH <- readHistoryFile(histfiles, mpirank, idxtime, "RH")
#     NC <- readHistoryFile(histfiles, mpirank, idxtime, "NC")
#     x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
#     y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")

#     i <- na.omit(match(XDAT[IIN], x))
#     j <- na.omit(match(YDAT[JIN], y))

#     for(k in c(1:ZNUM)){
#       # skip where i and j are both integer(0)
#       if(length(i) != 0 && length(j) != 0){
#         dav[k] <- dav[k] + sum((RH[i, j, k] / 100.0 - 1.0) * NC[i, j, k])
#       }
#     }
#   }
#   dav <- ifelse(nc_all != 0, dav / nc_all, 0.0)
#   return(dav)
# }

##### Get Domain Average of r_mean1
# getDomainAverageOfRmean1 <- function(histfiles, idxtime) {

#   dav <- rep(0.0, length=ZNUM)
#   for(mpirank in allmpiranks){
#     NC <- readHistoryFile(histfiles, mpirank, idxtime, "NC")
#     DMOM1 <- readHistoryFile(histfiles, mpirank, idxtime, "DMOM1")
#     x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
#     y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")

#     i <- na.omit(match(XDAT[IIN], x))
#     j <- na.omit(match(YDAT[JIN], y))

#     R_MEAN1 <- ifelse(NC != 0, DMOM1 / NC, 0.0)

#     for(k in c(1:ZNUM)){
#       # skip where i and j are both integer(0)
#       if(length(i) != 0 && length(j) != 0){
#         dav[k] <- dav[k] + sum(R_MEAN1[i,j,k])
#       }
#     }
#   }
#   dav <- dav / (length(IIN) * length(JIN))
#   return(dav)
# }

##### Get Variance of S
# getVarianceOfS <- function(histfiles, idxtime, dav) {

#   var <- rep(0.0, length=ZNUM)
#   for(mpirank in allmpiranks){
#     S <- readHistoryFile(histfiles, mpirank, idxtime, "RH") / 100.0 - 1.0
#     x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
#     y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")

#     i <- na.omit(match(XDAT[IIN], x))
#     j <- na.omit(match(YDAT[JIN], y))

#     for(k in c(1:ZNUM)){
#       # skip where i and j are both integer(0)
#       if(length(i) != 0 && length(j) != 0){
#         for(n in i){
#           for(m in j){
#             var[k] <-  var[k] + (S[n, m, k] - dav[k])^2
#           }
#         }
#       }
#     }
#   }
#   var  <- var / (length(IIN) * length(JIN))
#   return(var)
# }

############################################################################
########## Main Process
cat("Start output_2D_new.R\n")

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

  # Calculate the domain average of each parameters
  T_DAV <- getDomainAverage(histfiles, t, "T")
  QV_DAV <- getDomainAverage(histfiles, t, "QV")
  RH_DAV <- getDomainAverage(histfiles, t, "RH")
  LWC_DAV <- getDomainAverageOfLWC(histfiles, t) * 1e-3
  IWC_DAV <- getDomainAverageOfIWC(histfiles, t) * 1e-3
  # NC_ALL <- getSum(histfiles, t, "NC")
  # ## convert m^-3 to cm^-3
  # N_DROP_DAV <- NC_ALL / (length(IIN) * length(JIN)) * 1e-6
  # S_DROP_DAV <- getDomainAverageOfSdrop(histfiles, t, NC_ALL)
  # ## convert m to um
  # R_MEAN1_DAV <- getDomainAverageOfRmean1(histfiles, t) * 1e+6

  # # Calculate the variance of each parameters
  # S_VAR <- getVarianceOfS(histfiles, t, RH_DAV / 100.0 - 1.0)
  # S_DROP_VAR <- ifelse(NC_ALL != 0, S_VAR / NC_ALL, 0.0)

  tryCatch(
  {
    ncvar_put(outfile, Tvar, T_DAV, c(1,i), c(ZNUM,1))
    ncvar_put(outfile, QVvar, QV_DAV, c(1,i), c(ZNUM,1))
    ncvar_put(outfile, RHvar, RH_DAV, c(1,i), c(ZNUM,1))
    ncvar_put(outfile, LWCvar, LWC_DAV, c(1,i), c(ZNUM,1))
    ncvar_put(outfile, IWCvar, IWC_DAV, c(1,i), c(ZNUM,1))
    # ncvar_put(outfile, S2Svar, S_VAR, c(1,i), c(ZNUM,1))
    # ncvar_put(outfile, S_DROPvar, S_DROP_DAV, c(1,i), c(ZNUM,1))
    # ncvar_put(outfile, S2S_DROPvar, S_DROP_VAR, c(1,i), c(ZNUM,1))
    # ncvar_put(outfile, N_DROPvar, N_DROP_DAV, c(1,i), c(ZNUM,1))
    # ncvar_put(outfile, R_MEAN1var, R_MEAN1_DAV, c(1,i), c(ZNUM,1))
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
cat("End output_2D_new.R\n")
