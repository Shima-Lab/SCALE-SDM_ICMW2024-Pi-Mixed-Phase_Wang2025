############################################################################
###
### Program to output 1D data for ICMW-PI
### (for multiple run)
###
############################################################################
###
### History: 200320, S.Takahashi,   initital version
###          210303, S.Takahashi,   add outputs
###                                 include the near wall region
###          211221, S.Takahashi,   calculate by MPI, by one step
###
############################################################################

############################################################################
########## Loading Libraries
library(ncdf4)


############################################################################
########## Parameters
TIME_ITVL <-  1    # time interval of output file(s)
BORDER <- 0.125    # border of wall inside and outside(m)
EPSILON <- 0.001   # mean eddy dissipation rate(m^2/s^3)
RESULT_FILE <- "./LES_Nxxx_1D" # result file name

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
timedim <- ncdim_def("time", "seconds", outtimes, unlim=TRUE)

fillvalue <- -9.9999e+30
Tvar    <- ncvar_def("T", "K", list(timedim), fillvalue,
            longname="Temperature")
QVvar   <- ncvar_def("Qv", "kg/kg", list(timedim), fillvalue,
            longname="Water vapor mixing ratio")
RHvar   <- ncvar_def("RH", "%", list(timedim), fillvalue,
            longname="Relative humidity")
S_DROPvar  <- ncvar_def("S_drop","1", list(timedim), fillvalue,
             longname="Supersaturation at droplet locations")
S_ICEvar  <- ncvar_def("S_ice","1", list(timedim), fillvalue,
             longname="Supersaturation at ice locations")
LWCvar  <- ncvar_def("LWC","kgm^-3", list(timedim), fillvalue,
            longname="Liuquid water content")
## convert g to kg
IWCvar  <- ncvar_def("IWC","kgm^-3", list(timedim), fillvalue,
            longname="ice water content")         
## convert cm^-3 to m^-3         
N_DROPvar  <- ncvar_def("N_drop","m^-3", list(timedim), fillvalue,
             longname="droplet number concentration")
## convert cm^-3 to m^-3
N_ICEvar  <- ncvar_def("N_ice","m^-3", list(timedim), fillvalue,
             longname="ice number concentration")
## convert cm^-3 to m^-3
N_ARSOLvar  <- ncvar_def("N_aerosol","m^-3", list(timedim), fillvalue,
             longname="aerosol number concentration")
DISP_Rvar  <- ncvar_def("disp_r","1", list(timedim), fillvalue,
             longname="relative dispersion of droplet size distribution")
R_DMEANvar  <- ncvar_def("r_dmean","um", list(timedim), fillvalue,
             longname="droplet mean radius")
R_IMEANvar  <- ncvar_def("r_imean","um", list(timedim), fillvalue,
             longname="ice mean radius")
SVvar  <- ncvar_def("Sigma2_S","1", list(timedim), fillvalue,
             longname="Variance of supersaturation")
TVvar   <- ncvar_def("Sigma2_T","K^2", list(timedim), fillvalue,
            longname="Variance of temperature")
QVVvar  <- ncvar_def("Sigma2_Qv","(kg/kg)^2", list(timedim), fillvalue,
            longname="Variance of water vapor mixing ratio")
S_DROPVvar  <- ncvar_def("Sigma2_S_drop","1", list(timedim), fillvalue,
             longname="Variance of supersaturation at droplet locations")
S_ICEVvar  <- ncvar_def("Sigma2_S_ice","1", list(timedim), fillvalue,
             longname="Variance of supersaturation at ice locations")
epsvar <- ncvar_def("epsilon","m^2s^-3", list(timedim), fillvalue,
             longname="eddy dissipation rate")
TKEvar  <- ncvar_def("TKE","m^2s^-2", list(timedim), fillvalue,
             longname="Turbulent kinetic energy")
ncdflist <- list(Tvar, QVvar, RHvar, S_DROPvar, S_ICEvar, LWCvar, IWCvar, N_DROPvar, N_ICEvar, N_ARSOLvar,
      DISP_Rvar, R_DMEANvar, R_IMEANvar, SVvar, TVvar, QVVvar, S_DROPVvar, S_ICEVvar, epsvar, TKEvar)
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

  dav <- 0.0
  for(mpirank in allmpiranks){
    rd <- readHistoryFile(histfiles, mpirank, idxtime, param)
    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    # inside of wall
    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))
    k <- na.omit(match(ZDAT[KIN], z))

    if(length(i) != 0 && length(j) != 0 && length(k) != 0){
      dav <- dav + sum(rd[i, j, k])
    }
  }
  
  dav <- dav / (length(IIN) * length(JIN) * length(KIN))
  return(dav)
}

##### Get Domain Average of LWC
getDomainAverageOfLWC <- function(histfiles, idxtime) {

  dav <- 0.0
  for(mpirank in allmpiranks){
    QC <- readHistoryFile(histfiles, mpirank, idxtime, "QC_sd")
    QR <- readHistoryFile(histfiles, mpirank, idxtime, "QR_sd")
    DENS <- readHistoryFile(histfiles, mpirank, idxtime, "DENS")
    
    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))
    k <- na.omit(match(ZDAT[KIN], z))

    if(length(i) != 0 && length(j) != 0 && length(k) != 0){
      LWC <- (QC + QR) * DENS
      dav <- dav + sum(LWC[i, j, k])* 1000.0
    }
  }
  dav <- dav / (length(IIN) * length(JIN) * length(KIN))
  return(dav)
}

##### Get Domain Average of IWC
getDomainAverageOfIWC <- function(histfiles, idxtime) {

  dav <- 0.0
  for(mpirank in allmpiranks){
    QI <- readHistoryFile(histfiles, mpirank, idxtime, "QI_sd")
    QS <- readHistoryFile(histfiles, mpirank, idxtime, "QS_sd")
    QG <- readHistoryFile(histfiles, mpirank, idxtime, "QG_sd")
    DENS <- readHistoryFile(histfiles, mpirank, idxtime, "DENS")
    
    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))
    k <- na.omit(match(ZDAT[KIN], z))

    if(length(i) != 0 && length(j) != 0 && length(k) != 0){
      IWC <- (QI + QS + QG) * DENS
      dav <- dav + sum(IWC[i, j, k])* 1000.0
    }
  }
  dav <- dav / (length(IIN) * length(JIN) * length(KIN))
  return(dav)
}

##### Get Sum(single parameter)
getSum <- function(histfiles, idxtime, param) {

  sm <- 0.0
  for(mpirank in allmpiranks){
    rd <- readHistoryFile(histfiles, mpirank, idxtime, param)
    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    # inside of wall
    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))
    k <- na.omit(match(ZDAT[KIN], z))

    if(length(i) != 0 && length(j) != 0 && length(k) != 0){
      sm <- sm + sum(rd[i, j, k])
    }
  }
  return(sm)
}

##### Get Domain Average of S_drop
getDomainAverageOfSdrop <- function(histfiles, idxtime, nc_all) {

  dav <- 0.0
  for(mpirank in allmpiranks){
    RH <- readHistoryFile(histfiles, mpirank, idxtime, "RH")
    NC <- readHistoryFile(histfiles, mpirank, idxtime, "NC")

    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))
    k <- na.omit(match(ZDAT[KIN], z))

    if(length(i) != 0 && length(j) != 0 && length(k) != 0){
      Sdrop <- (RH / 100.0 - 1.0) * NC
      dav <- dav + sum(Sdrop[i, j, k])
    }
  }
  if(nc_all != 0){
    dav <- dav / nc_all
  }else{
    dav <- 0.0
  }
  return(dav)
}

##### Get Domain Average of S_ice
getDomainAverageOfSice <- function(histfiles, idxtime, ni_all) {

  dav <- 0.0
  for(mpirank in allmpiranks){
    RH <- readHistoryFile(histfiles, mpirank, idxtime, "RH")
    NI <- readHistoryFile(histfiles, mpirank, idxtime, "NI")

    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))
    k <- na.omit(match(ZDAT[KIN], z))

    if(length(i) != 0 && length(j) != 0 && length(k) != 0){
      Sice <- (RH / 100.0 - 1.0) * NI
      dav <- dav + sum(Sice [i, j, k])
    }
  }
  #inside <- inside / (length(IIN) * length(JIN) * length(KIN))
  if(ni_all != 0){
    dav <- dav / ni_all
  }else{
    dav <- 0.0
  }
  return(dav)
}

##### Get Domain Average of disp_r
getDomainAverageOfDispr <- function(histfiles, idxtime) {

  dav <- 0.0
  for(mpirank in allmpiranks){
    NC <- readHistoryFile(histfiles, mpirank, idxtime, "NC")
    DMOM1 <- readHistoryFile(histfiles, mpirank, idxtime, "DMOM1")
    DMOM2 <- readHistoryFile(histfiles, mpirank, idxtime, "DMOM2")

    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))
    k <- na.omit(match(ZDAT[KIN], z))

    if(length(i) != 0 && length(j) != 0 && length(k) != 0){
      R_MEAN1 <- ifelse(NC != 0, DMOM1 / NC, 0.0)
      R_MEAN2 <- ifelse(NC != 0, DMOM2 / NC, 0.0)
      DISP_R <- ifelse(R_MEAN1 != 0,
                sqrt(pmax(R_MEAN2 - R_MEAN1^2, 0.0)) / R_MEAN1, 0.0)
      dav <- dav + sum(DISP_R[i, j, k])
    }
  }
  dav <- dav / (length(IIN) * length(JIN) * length(KIN))
  return(dav)
}

##### Get Variance(single parameter)
getVariance <- function(histfiles, idxtime, param, ave) {

  var <- 0.0
  for(mpirank in allmpiranks){
    rd <- readHistoryFile(histfiles, mpirank, idxtime, param)

    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))
    k <- na.omit(match(ZDAT[KIN], z))

    if(length(i) != 0 && length(j) != 0 && length(k) != 0){
        var <-  var + sum((rd[i, j, k] - ave)^2)
      }
  }
  var <- var / (length(IIN) * length(JIN) * length(KIN))
  return(var)
}

##### Get Variance of S
getVarianceOfS <- function(histfiles, idxtime, dav) {

  var <- 0.0
  for(mpirank in allmpiranks){
    S <- readHistoryFile(histfiles, mpirank, idxtime, "RH") / 100.0 - 1.0

    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))
    k <- na.omit(match(ZDAT[KIN], z))

    if(length(i) != 0 && length(j) != 0 && length(k) != 0){
      var <-  var + sum((S[i, j, k] - dav)^2)
    }
  }
  var  <- var / (length(IIN) * length(JIN) * length(KIN))
  return(var)
}

##### Get Variance of S_drop
getVarianceOfSdrop <- function(histfiles, idxtime, dav, nc_all) {

  var <- 0.0
  for(mpirank in allmpiranks){
    S <- readHistoryFile(histfiles, mpirank, idxtime, "RH") / 100.0 - 1.0
    NC <- readHistoryFile(histfiles, mpirank, idxtime, "NC")

    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))
    k <- na.omit(match(ZDAT[KIN], z))

    if(length(i) != 0 && length(j) != 0 && length(k) != 0){
      var <-  var + sum(NC[i, j, k] * (S[i, j, k] - dav)^2)
    }
  }
  var <- ifelse(nc_all != 0, var / nc_all, 0.0)
  return(var)
}

##### Get Variance of S_ice
getVarianceOfSice <- function(histfiles, idxtime, dav, ni_all) {

  var <- 0.0
  for(mpirank in allmpiranks){
    S <- readHistoryFile(histfiles, mpirank, idxtime, "RH") / 100.0 - 1.0
    NI <- readHistoryFile(histfiles, mpirank, idxtime, "NI")

    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))
    k <- na.omit(match(ZDAT[KIN], z))

    if(length(i) != 0 && length(j) != 0 && length(k) != 0){
      var <-  var + sum(NI[i, j, k] * (S[i, j, k] - dav)^2)
    }
  }
  var <- ifelse(ni_all != 0, var / ni_all, 0.0)
  return(var)
}

############################################################################
########## Main Process
cat("Start output_1D.R\n")

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
  TKE_SMG_DAV <- getDomainAverage(histfiles, t, "TKE_SMG")
  NC_ALL <- getSum(histfiles, t, "NC")
  NI_ALL <- getSum(histfiles, t, "NI")
  NA_ALL <- getSum(histfiles, t, "NASL")
  DMOM1_ALL <- getSum(histfiles, t, "DMOM1")
  IMOM1_ALL <- getSum(histfiles, t, "IMOM1")
  N_DROP_DAV <- NC_ALL / (length(IIN) * length(JIN) * length(KIN))
  N_ICE_DAV <- NI_ALL / (length(IIN) * length(JIN) * length(KIN))
  N_ARSOL_DAV <- NA_ALL / (length(IIN) * length(JIN) * length(KIN))
  S_DROP_DAV <- getDomainAverageOfSdrop(histfiles, t, NC_ALL)
  S_ICE_DAV <- getDomainAverageOfSice(histfiles, t, NI_ALL)
  ## convert g to kg
  LWC_DAV <- getDomainAverageOfLWC(histfiles, t) * 1e-3
  ## convert g to kg
  IWC_DAV <- getDomainAverageOfIWC(histfiles, t) * 1e-3
  DISP_R_DAV <- getDomainAverageOfDispr(histfiles, t)
  ## convert m to um
  R_DMEAN_DAV <- ifelse(NC_ALL != 0, DMOM1_ALL / NC_ALL  * 1e+6 , 0.0)
  R_IMEAN_DAV <- ifelse(NI_ALL != 0, IMOM1_ALL / NI_ALL  * 1e+6 , 0.0)

  # Calculate the variance of each parameters
  T_VAR <- getVariance(histfiles, t, "T", T_DAV)
  QV_VAR <- getVariance(histfiles, t, "QV", QV_DAV)
  S_VAR <- getVarianceOfS(histfiles, t, RH_DAV / 100.0 - 1.0)
  S_DROP_VAR <- getVarianceOfSdrop(histfiles, t, S_DROP_DAV, NC_ALL)
  S_ICE_VAR <- getVarianceOfSice(histfiles, t, S_ICE_DAV, NI_ALL)

  tryCatch(
  {
    ncvar_put(outfile, Tvar, T_DAV, i, 1)
    ncvar_put(outfile, QVvar, QV_DAV, i, 1)
    ncvar_put(outfile, RHvar, RH_DAV, i, 1)
    ncvar_put(outfile, S_DROPvar, S_DROP_DAV, i, 1)
    ncvar_put(outfile, S_ICEvar, S_ICE_DAV, i, 1)
    ncvar_put(outfile, LWCvar, LWC_DAV, i, 1)
    ncvar_put(outfile, IWCvar, IWC_DAV, i, 1)
    ncvar_put(outfile, N_DROPvar, N_DROP_DAV, i, 1)
    ncvar_put(outfile, N_ICEvar, N_ICE_DAV, i, 1)
    ncvar_put(outfile, N_ARSOLvar, N_ARSOL_DAV, i, 1)
    ncvar_put(outfile, DISP_Rvar, DISP_R_DAV, i, 1)
    ncvar_put(outfile, R_DMEANvar, R_DMEAN_DAV, i, 1)
    ncvar_put(outfile, R_IMEANvar, R_IMEAN_DAV, i, 1)
    ncvar_put(outfile, SVvar, S_VAR, i, 1)
    ncvar_put(outfile, TVvar, T_VAR, i, 1)
    ncvar_put(outfile, QVVvar, QV_VAR, i, 1)
    ncvar_put(outfile, S_DROPVvar, S_DROP_VAR, i, 1)
    ncvar_put(outfile, S_ICEVvar, S_ICE_VAR, i, 1)
    ncvar_put(outfile, epsvar, EPSILON, i, 1)
    ncvar_put(outfile, TKEvar, TKE_SMG_DAV, i, 1)
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

# # read numbers_XXX.txt 
# cat("start read ../run_??/ptl_dist_1D_R/numbers_XXX.txt\n")
# N_ARSOL_DAV   <- array(0, dim=c(length(outtimes)))
# names(N_ARSOL_DAV) <- outtimes

# allnumbsfiles <- NULL
# for(rundir in allrundirs){
#   tmp = dir(paste("../",rundir,"/ptl_dist_1D_R",sep=""),
#           pattern="^numbers_")
#   allnumbsfiles <- append(allnumbsfiles,
#                      paste("../",rundir,"/ptl_dist_1D_R/",tmp,sep=""))
# } 

# for(nmbfile in allnumbsfiles){
#   if (file.size(nmbfile) == 0) next
#   nmbs <- read.table(nmbfile, header=F, sep="=",comment.char="")
#   row <- 1
#   for(v1 in nmbs$V1){
#     if(v1 == "###### SD/RD numbers at time "){
#       tmp_time <- strsplit(as.character(nmbs$V2[row]), " ")[[1]][2]
#     }else if(v1 == "total aerosol number density "){
#       v2 <- strsplit(as.character(nmbs$V2[row]), " ")
#       if(is.element(tmp_time,outtimes)){
#         N_ARSOL_DAV[tmp_time] <- as.numeric(v2[[1]][2])
#       }
#     }
#     row <- row + 1
#   }
# }
## convert m^-3 to cm^-3
#N_ARSOL_DAV <- N_ARSOL_DAV * 1e-6

# Write other data
tryCatch(
  {
    #ncvar_put(outfile, N_ARSOLvar, N_ARSOL_DAV)
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
cat("End output_1D.R\n")
