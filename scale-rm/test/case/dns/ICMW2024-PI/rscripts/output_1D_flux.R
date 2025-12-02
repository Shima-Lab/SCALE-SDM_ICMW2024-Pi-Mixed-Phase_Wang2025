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
RESULT_FILE <- "./LES_Nxxx_1D_flux" # result file name

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
SFLX_B_SHvar  <- ncvar_def("H_flux_b","W/m2", list(timedim), fillvalue,
             longname="surface sensible heat flux at the bottom")
SFLX_B_QVvar  <- ncvar_def("qv_flux_b","kg/kg m/s", list(timedim), fillvalue,
             longname="surface moisuture flux at the bottom")
SFLX_T_SHvar  <- ncvar_def("H_flux_t","W/m2", list(timedim), fillvalue,
             longname="surface sensible heat flux at the top")
SFLX_T_QVvar  <- ncvar_def("qv_flux_t","kg/kg m/s", list(timedim), fillvalue,
             longname="surface moisuture flux at the top")
# SFLX_W_SHvar  <- ncvar_def("SFLX_W_SH","W/m2", list(timedim), fillvalue,
#              longname="surface sensible heat flux at the west side")
# SFLX_W_QVvar  <- ncvar_def("SFLX_W_QV","kg/kg m/s", list(timedim), fillvalue,
#              longname="surface moisuture flux at the west side")
# SFLX_E_SHvar  <- ncvar_def("SFLX_E_SH","W/m2", list(timedim), fillvalue,
#              longname="surface sensible heat flux at the east side")
# SFLX_E_QVvar  <- ncvar_def("SFLX_E_QV","kg/kg m/s", list(timedim), fillvalue,
#              longname="surface moisuture flux at the east side")
# SFLX_S_SHvar  <- ncvar_def("SFLX_S_SH","W/m2", list(timedim), fillvalue,
#              longname="surface sensible heat flux at the south side")
# SFLX_S_QVvar  <- ncvar_def("SFLX_S_QV","kg/kg m/s", list(timedim), fillvalue,
#              longname="surface moisuture flux at the south side")
# SFLX_N_SHvar  <- ncvar_def("SFLX_N_SH","W/m2", list(timedim), fillvalue,
#              longname="surface sensible heat flux at the north side")
# SFLX_N_QVvar  <- ncvar_def("SFLX_N_QV","kg/kg m/s", list(timedim), fillvalue,
#              longname="surface moisuture flux at the north side")
# ncdflist <- list(SFLX_B_SHvar, SFLX_B_QVvar, SFLX_T_SHvar, SFLX_T_QVvar, SFLX_W_SHvar, SFLX_W_QVvar,
#             SFLX_E_SHvar, SFLX_E_QVvar, SFLX_S_SHvar, SFLX_S_QVvar, SFLX_N_SHvar, SFLX_N_QVvar)
ncdflist <- list(SFLX_B_SHvar, SFLX_B_QVvar, SFLX_T_SHvar, SFLX_T_QVvar)
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

# ##### Read History File
# readHistoryFile <- function(histfiles, mpirank, idxtime, pname) {

#   read <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], pname,
#             c(1,1,1,idxtime), c(IMAX,JMAX,KMAX,1))
#   if(!setequal(dim(read[,,]),c(IMAX,JMAX,KMAX))){
#     message("ERROR: Size of variables is not consistent.mpirank=",mpirank)
#     q()
#   }
#   return(read)
# }

##### Read History File for flux (Top & Bottom)
readHistoryFileTB <- function(histfiles, mpirank, idxtime, pname) {

  read <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], pname,
            c(1,1,idxtime), c(IMAX,JMAX,1))
  if(!setequal(dim(read[,]),c(IMAX,JMAX))){
    message("ERROR: Size of variables is not consistent.mpirank=",mpirank)
    q()
  }
  return(read)
}

##### Read History File for flux (East & West)
readHistoryFileEW <- function(histfiles, mpirank, idxtime, pname) {

  read <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], pname,
            c(1,1,idxtime), c(JMAX,KMAX,1))
  if(!setequal(dim(read[,]),c(JMAX,KMAX))){
    message("ERROR: Size of variables is not consistent.mpirank=",mpirank)
    q()
  }
  return(read)
}

##### Read History File for flux (North & South)
readHistoryFileNS <- function(histfiles, mpirank, idxtime, pname) {

  read <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], pname,
            c(1,1,idxtime), c(IMAX,KMAX,1))
  if(!setequal(dim(read[,]),c(IMAX,KMAX))){
    message("ERROR: Size of variables is not consistent.mpirank=",mpirank)
    q()
  }
  return(read)
}

##### Get Domain Average(single parameter)
# getDomainAverage <- function(histfiles, idxtime, param) {

#   dav <- 0.0
#   for(mpirank in allmpiranks){
#     rd <- readHistoryFile(histfiles, mpirank, idxtime, param)
#     x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
#     y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
#     z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

#     # inside of wall
#     i <- na.omit(match(XDAT[IIN], x))
#     j <- na.omit(match(YDAT[JIN], y))
#     k <- na.omit(match(ZDAT[KIN], z))

#     if(length(i) != 0 && length(j) != 0 && length(k) != 0){
#       dav <- dav + sum(rd[i, j, k])
#     }
#   }
  
#   dav <- dav / (length(IIN) * length(JIN) * length(KIN))
#   return(dav)
# }

##### Get Domain Average of flux (Top & Bottom)
getDomainAverageOfTBflux <- function(histfiles, idxtime, param) {

  dav <- 0.0
  for(mpirank in allmpiranks){
    rd <- readHistoryFileTB(histfiles, mpirank, idxtime, param)
    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    #z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    # inside of wall
    i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))
    #k <- na.omit(match(ZDAT[KIN], z))

    if(length(i) != 0 && length(j) != 0){
      dav <- dav + sum(rd[i, j])
    }
  }
  
  dav <- dav / (length(IIN) * length(JIN))
  return(dav)
}

#### Get Domain Average of flux (East & West)
getDomainAverageOfEWflux <- function(histfiles, idxtime, param) {

  dav <- 0.0
  for(mpirank in allmpiranks){
    rd <- readHistoryFileEW(histfiles, mpirank, idxtime, param)
    #x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    # inside of wall
    #i <- na.omit(match(XDAT[IIN], x))
    j <- na.omit(match(YDAT[JIN], y))
    k <- na.omit(match(ZDAT[KIN], z))

    if(length(j) != 0 && length(k) != 0){
      dav <- dav + sum(rd[j, k])
    }
  }
  
  dav <- dav / (length(JIN) * length(KIN))
  return(dav)
}

#### Get Domain Average of flux (North & South)
getDomainAverageOfNSflux <- function(histfiles, idxtime, param) {

  dav <- 0.0
  for(mpirank in allmpiranks){
    rd <- readHistoryFileNS(histfiles, mpirank, idxtime, param)
    x <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "x")
    #y <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "y")
    z <- ncvar_get(histfiles[[as.numeric(mpirank)+1]], "z")

    # inside of wall
    i <- na.omit(match(XDAT[IIN], x))
    #j <- na.omit(match(YDAT[JIN], y))
    k <- na.omit(match(ZDAT[KIN], z))

    if(length(i) != 0 && length(k) != 0){
      dav <- dav + sum(rd[i, k])
    }
  }
  
  dav <- dav / (length(IIN) * length(KIN))
  return(dav)
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
  SFLX_B_SH_DAV <- getDomainAverageOfTBflux(histfiles, t, "SFLX_B_SH")
  SFLX_B_QV_DAV <- getDomainAverageOfTBflux(histfiles, t, "SFLX_B_QV")
  SFLX_T_SH_DAV <- getDomainAverageOfTBflux(histfiles, t, "SFLX_T_SH")
  SFLX_T_QV_DAV <- getDomainAverageOfTBflux(histfiles, t, "SFLX_T_QV")
  #SFLX_W_SH_DAV <- getDomainAverageOfTBflux(histfiles, t, "SFLX_W_SH")
  # SFLX_W_QV_DAV <- getDomainAverageOfEWflux(histfiles, t, "SFLX_W_QV")
  # SFLX_E_SH_DAV <- getDomainAverageOfEWflux(histfiles, t, "SFLX_E_SH")
  # SFLX_E_QV_DAV <- getDomainAverageOfEWflux(histfiles, t, "SFLX_E_QV")
  # SFLX_S_SH_DAV <- getDomainAverageOfNSflux(histfiles, t, "SFLX_S_SH")
  # SFLX_S_QV_DAV <- getDomainAverageOfNSflux(histfiles, t, "SFLX_S_QV")
  # SFLX_N_SH_DAV <- getDomainAverageOfNSflux(histfiles, t, "SFLX_N_SH")
  # SFLX_N_QV_DAV <- getDomainAverageOfNSflux(histfiles, t, "SFLX_N_QV")

  tryCatch(
  {
    ncvar_put(outfile, SFLX_B_SHvar, SFLX_B_SH_DAV, i, 1)
    ncvar_put(outfile, SFLX_B_QVvar, SFLX_B_QV_DAV, i, 1)
    ncvar_put(outfile, SFLX_T_SHvar, SFLX_T_SH_DAV, i, 1)
    ncvar_put(outfile, SFLX_T_QVvar, SFLX_T_QV_DAV, i, 1)
    #ncvar_put(outfile, SFLX_W_SHvar, SFLX_W_SH_DAV, i, 1)
    # ncvar_put(outfile, SFLX_W_QVvar, SFLX_W_QV_DAV, i, 1)
    # ncvar_put(outfile, SFLX_E_SHvar, SFLX_E_SH_DAV, i, 1)
    # ncvar_put(outfile, SFLX_E_QVvar, SFLX_E_QV_DAV, i, 1)
    # ncvar_put(outfile, SFLX_S_SHvar, SFLX_S_SH_DAV, i, 1)
    # ncvar_put(outfile, SFLX_S_QVvar, SFLX_S_QV_DAV, i, 1)
    # ncvar_put(outfile, SFLX_N_SHvar, SFLX_N_SH_DAV, i, 1)
    # ncvar_put(outfile, SFLX_N_QVvar, SFLX_N_QV_DAV, i, 1)
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
