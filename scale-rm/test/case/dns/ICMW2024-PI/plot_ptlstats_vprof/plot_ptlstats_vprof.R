###########################################################################
###
### Program to plot the x-z distribution of droplet size dist statistics
###
############################################################################

############################################################################
########## Loading Libraries
library(ncdf4)
library(fields)

############################################################################
########## Parameters
#pal <- colorRampPalette(c("blue", "orange", "red"), space = "Lab")
pal <- rainbow(100, start=0.4, end=1.0)

crt_qhyd <- 1.0e-5        # [kg/kg] critical mixing ratio to determine cloudy grid
CLOUD_BASE_WIDTH <- 1.0e3 # [m]

# read parameters
tmp <- read.table(textConnection(gsub("="," ",gsub(",", " ", readLines("../run.conf")))),header=FALSE,fill=TRUE,flush=T)
if( is.element("FLUX_CX", tmp$V1) ){
    FLUX_CX <- as.numeric(as.character(tmp[tmp$V1=="FLUX_CX",2]))
}else{
    FLUX_CX <- 5000.0
}

if( is.element("FLUX_CY", tmp$V1) ){
    FLUX_CY <- as.numeric(as.character(tmp[tmp$V1=="FLUX_CY",2]))
}else{
    FLUX_CY <- 5000.0
}

########## functions
#### calculate saturation vapor pressure from temperature
#### Teten's formula
temp.to.satvp <- function (temp)
{ # temp : temperature [K]
  para_a1 <- 6.11 #[hPa]
  para_a2 <- 17.27
  para_a3 <- 237.3 #[[deg C]]
  T0 <- 273.15      # 0 deg C in K [K]

  # saturation vapor pressure [hPa]
  satvp <- para_a1 * exp((para_a2 * (temp-T0) )/( (temp-T0) + para_a3))

  return(satvp)
}

#### calculate saturation vapor density from temperature
temp.to.satrhov <- function (temp)
{ # temp : temperature [K]
  gas_v <- 461.5 # specific gas constant for vapor [J/kg/K]

  # saturation vapor density [kg/m^3]
  satrhov <- 100.0*temp.to.satvp(temp)/gas_v/temp

  return(satrhov)
}

############################################################################
########## Read Data from Files
##### Make a list of files,  mpiranks
allfiles = dir("../",pattern="^history.")
tmp = strsplit(allfiles,"\\history.pe|\\.nc")
allmpiranks = unique(matrix(unlist(tmp),nrow=2)[2,])
names(allfiles) = allmpiranks

##### Open the first file
ncin <- nc_open(paste("../",allfiles[1],sep=""))

##### Times
alltimes <- ncvar_get(ncin,"time")
time_units <- ncatt_get(ncin,"time","units")

##### Close
nc_close(ncin)

###########################################################################
########## Region
##### read the first MPI rank for initialization
for(mpirank in allmpiranks[1]){
	##### Open
	ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

	##### Grids
	KMAX <- length(ncvar_get(ncin,"CZ"))-4
	CZ <- ncvar_get(ncin,"CZ")
        CDZ <- ncvar_get(ncin,"CDZ")

	##### Close
	nc_close(ncin)
}

############################################################################
#################### Plot Adiabatic Fraction

########## Cloud Base Height
CLOUD_BASE_HEIGHT <- rep(0.0,length(alltimes))
names(CLOUD_BASE_HEIGHT) <- alltimes
CLOUD_BASE_K <- rep(1000000,length(alltimes))
names(CLOUD_BASE_K) <- alltimes

##### loop of MPI rank
for(mpirank in allmpiranks){
        cat(sprintf("processing the rank = %s \n",mpirank))

        ##### Open
        ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

        ##### Grids
        IMAX <- length(ncvar_get(ncin,"CX"))-4
        JMAX <- length(ncvar_get(ncin,"CY"))-4
        KMAX <- length(ncvar_get(ncin,"CZ"))-4
        CDX <- ncvar_get(ncin,"CDX")
        CDY <- ncvar_get(ncin,"CDY")
        CDZ <- ncvar_get(ncin,"CDZ")
        CZ <- ncvar_get(ncin,"CZ")

        ##### Read QHYD_sd
        QHYD <- ncvar_get(ncin,"QHYD_sd")
        QHYD_units <- ncatt_get(ncin,"QHYD_sd","units")
        dimnames(QHYD)[[4]] <- alltimes
        if(!setequal(dim(QHYD[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Close
        nc_close(ncin)

        ##### loop of time
        for(time in alltimes){
                basek_mpi <- min(
                         apply(QHYD[1:IMAX,1:JMAX,1:KMAX,as.character(time)], c(1,2),
                               function(y){min(which(y>crt_qhyd),KMAX)})
                         )
                CLOUD_BASE_K[as.character(time)] <- min(CLOUD_BASE_K[as.character(time)], basek_mpi)
        }

}

CLOUD_BASE_K[which(CLOUD_BASE_K==KMAX)] <- 1
CLOUD_BASE_HEIGHT[] <- CZ[CLOUD_BASE_K[]+2]
#CLOUD_BASE_HEIGHT

########## Cloud Base Density, Temperature and Saturation Vapor Density
CLOUD_BASE_DENS <- rep(0.0,length(alltimes))
names(CLOUD_BASE_DENS) <- alltimes
CLOUD_BASE_TEMP <- rep(0.0,length(alltimes))
names(CLOUD_BASE_TEMP) <- alltimes
CLOUD_BASE_RHOVSAT <- rep(0.0,length(alltimes))
names(CLOUD_BASE_RHOVSAT) <- alltimes

count_cbgrd <- rep(0,length(alltimes))
names(count_cbgrd) <- alltimes

##### loop of MPI rank
for(mpirank in allmpiranks){
        cat(sprintf("processing the rank = %s \n",mpirank))

        ##### Open
        ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

        ##### Grids
        IMAX <- length(ncvar_get(ncin,"CX"))-4
        JMAX <- length(ncvar_get(ncin,"CY"))-4
        KMAX <- length(ncvar_get(ncin,"CZ"))-4
        CDX <- ncvar_get(ncin,"CDX")
        CDY <- ncvar_get(ncin,"CDY")
        CDZ <- ncvar_get(ncin,"CDZ")
        CX <- ncvar_get(ncin,"CX")
        CY <- ncvar_get(ncin,"CY")
        CZ <- ncvar_get(ncin,"CZ")

        ##### Read DENS
        DENS <- ncvar_get(ncin,"DENS")
        DENS_units <- ncatt_get(ncin,"DENS","units")
        dimnames(DENS)[[4]] <- alltimes
        if(!setequal(dim(DENS[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read T (Temperature)
        TEMP <- ncvar_get(ncin,"T")
        TEMP_units <- ncatt_get(ncin,"T","units")
        dimnames(TEMP)[[4]] <- alltimes
        if(!setequal(dim(TEMP[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Close
        nc_close(ncin)

        ##### loop of time
        for(time in alltimes){
                for(i in seq(1:IMAX)){
                for(j in seq(1:JMAX)){
                      	 if( (CX[i+2]>=(FLUX_CX-CLOUD_BASE_WIDTH) ) && (CX[i+2]<=(FLUX_CX+CLOUD_BASE_WIDTH) ) ){
                      	 if( (CY[j+2]>=(FLUX_CY-CLOUD_BASE_WIDTH) ) && (CY[j+2]<=(FLUX_CY+CLOUD_BASE_WIDTH) ) ){
                                CLOUD_BASE_DENS[as.character(time)] <- CLOUD_BASE_DENS[as.character(time)] + DENS[i,j,CLOUD_BASE_K[as.character(time)],as.character(time)]
                                CLOUD_BASE_TEMP[as.character(time)] <- CLOUD_BASE_TEMP[as.character(time)] + TEMP[i,j,CLOUD_BASE_K[as.character(time)],as.character(time)]
                                count_cbgrd[as.character(time)] <- count_cbgrd[as.character(time)] + 1
                      	 }
		      	 }
                }
                }
        }
}

CLOUD_BASE_DENS[] <- CLOUD_BASE_DENS[]/count_cbgrd[]
CLOUD_BASE_TEMP[] <- CLOUD_BASE_TEMP[]/count_cbgrd[]
CLOUD_BASE_RHOVSAT <- temp.to.satrhov(CLOUD_BASE_TEMP)

########## Plot Adiabatic Fraction
# Setting of bins
nbins = 100
min_AF <-  0.0 # []
max_AF <-  1.5 # []

AF.bin = seq(min_AF, max_AF, length=nbins)
d_AF.bin = AF.bin[2] - AF.bin[1]

col_scale_min <- log10(1.0e5)
col_scale_max <- log10(1.0e8)
 
##### initialize histogram
AF_dist <- array(rep(1.0e-100,nbins*KMAX*length(alltimes)),c(nbins,KMAX,length(alltimes)))
dimnames(AF_dist)[[3]] <- alltimes

##### initialize a list to store AFs of all cloudy grid 
AF_list <- list()
for(t in 1:length(alltimes) ){
      AF_list[[t]] <- list()
      for (idxk in 1:KMAX){
      	  AF_list[[t]][[idxk]] <- numeric(0)
      }
}
names(AF_list) <- alltimes

##### loop of MPI rank
    for(mpirank in allmpiranks){
	cat(sprintf("processing the rank = %s \n",mpirank))

	##### Open
	ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

	##### Grids
	IMAX <- length(ncvar_get(ncin,"CX"))-4
	JMAX <- length(ncvar_get(ncin,"CY"))-4
	KMAX <- length(ncvar_get(ncin,"CZ"))-4
	CX <- ncvar_get(ncin,"CX")
	CY <- ncvar_get(ncin,"CY")
	CZ <- ncvar_get(ncin,"CZ")
      	CDX <- ncvar_get(ncin,"CDX")
        CDY <- ncvar_get(ncin,"CDY")
        CDZ <- ncvar_get(ncin,"CDZ")

        ##### Read QR
        QR <- ncvar_get(ncin,"QR_sd")
        QR_units <- ncatt_get(ncin,"QR_sd","units")
        dimnames(QR)[[4]] <- alltimes
        if(!setequal(dim(QR[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read QC
        QC <- ncvar_get(ncin,"QC_sd")
        QC_units <- ncatt_get(ncin,"QC_sd","units")
        dimnames(QC)[[4]] <- alltimes
        if(!setequal(dim(QC[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read DENS
        DENS <- ncvar_get(ncin,"DENS")
        DENS_units <- ncatt_get(ncin,"DENS","units")
        dimnames(DENS)[[4]] <- alltimes
        if(!setequal(dim(DENS[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read T (Temperature)
        TEMP <- ncvar_get(ncin,"T")
        TEMP_units <- ncatt_get(ncin,"T","units")
        dimnames(TEMP)[[4]] <- alltimes
        if(!setequal(dim(TEMP[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read QHYD_sd
        QHYD <- ncvar_get(ncin,"QHYD_sd")
        QHYD_units <- ncatt_get(ncin,"QHYD_sd","units")
        dimnames(QHYD)[[4]] <- alltimes
        if(!setequal(dim(QHYD[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Close
        nc_close(ncin)

	### Making the Histogram
#	for(time in alltimes[which(as.numeric(alltimes)==7200)]){
     	for(time in alltimes){
		 #### calculate
                 LWC <- DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
                           (QC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]+QR[1:IMAX,1:JMAX,1:KMAX,as.character(time)]) 

                 RHOVSAT <- apply(TEMP[1:IMAX,1:JMAX,1:KMAX,as.character(time)], c(1,2,3), temp.to.satrhov)
                 ALWC <- DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
                         CLOUD_BASE_RHOVSAT[as.character(time)]/CLOUD_BASE_DENS[as.character(time)] -
                         RHOVSAT[1:IMAX,1:JMAX,1:KMAX]

                 AF <- LWC[1:IMAX,1:JMAX,1:KMAX]/ALWC[1:IMAX,1:JMAX,1:KMAX]
                 AF[1:IMAX,1:JMAX,1:CLOUD_BASE_K[as.character(time)]] <- 0.0

   		 # update the histogram
		 iAF <- array(rep(0,IMAX*JMAX*KMAX),c(IMAX,JMAX,KMAX))
               	 iAF[1:IMAX,1:JMAX,1:KMAX] <- findInterval(AF[1:IMAX,1:JMAX,1:KMAX], AF.bin, all.inside=TRUE)
                 for (idxi in 1:IMAX){
                 for (idxj in 1:JMAX){
                 for (idxk in 1:KMAX){
                    AF_dist[iAF[idxi,idxj,idxk],idxk,as.character(time)] <- AF_dist[iAF[idxi,idxj,idxk],idxk,as.character(time)] +
                                                      	CDX[idxi+2]*CDY[idxj+2]*CDZ[idxk+2]/CDZ[idxk+2]/d_AF.bin

                 }
                 }
                 }
		 # append cloudy grid data to the list
                 for (idxk in 1:KMAX){
	            tmp_AF   <- AF[1:IMAX,1:JMAX,idxk]
		    tmp_QHYD <- QHYD[1:IMAX,1:JMAX,idxk,as.character(time)]
		    AF_list[[as.character(time)]][[idxk]] <- append(AF_list[[as.character(time)]][[idxk]], tmp_AF[tmp_QHYD>crt_qhyd])
                 }
	}
    }

#### Plot Data
#    for(time in alltimes[which(as.numeric(alltimes)==7200)]){
     for(time in alltimes){
	## calculate percentile
     	ptile05 <- numeric(KMAX)
     	ptile50 <- numeric(KMAX)
     	ptile95 <- numeric(KMAX)
     	for (idxk in 1:KMAX){
	    ptile05[idxk] <- quantile(AF_list[[as.character(time)]][[idxk]],0.05)
	    ptile50[idxk] <- quantile(AF_list[[as.character(time)]][[idxk]],0.50)
	    ptile95[idxk] <- quantile(AF_list[[as.character(time)]][[idxk]],0.95)
	}

	## plot data
        pdf(paste("AF.vprof.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
    	image.plot(AF.bin, CZ[3:(KMAX+2)], log10(AF_dist[,,as.character(time)]),
                   zlim=c(col_scale_min,col_scale_max),
                   col = pal,
                   main=paste("Adiabatic Fraction Distribution (T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)\n(Volume Density log10([m^3/(unit AF)/m]))"),
                   xlab="Adiabatic Fraction []",
                   ylab="Height z [m]",
		   #ann=F,
		   useRaster=TRUE)

#	lines(ptile05,CZ[3:(KMAX+2)])
	lines(ptile50,CZ[3:(KMAX+2)])
#	lines(ptile95,CZ[3:(KMAX+2)])

    	#### Close all the files
    	graphics.off()
    }

############################################################################
#################### Plot the liquid water content
# Setting of bins
nbins = 100
min_lwc <-  0.0    # [kg/m^3]
max_lwc <-  8.0e-3 # [kg/m^3]

lwc.bin = seq(min_lwc, max_lwc, length=nbins)
d_lwc.bin = lwc.bin[2] - lwc.bin[1]

col_scale_min <- log10(1.0e7)
col_scale_max <- log10(1.0e10)
 
##### initialize histogram
lwc_dist <- array(rep(1.0e-100,nbins*KMAX*length(alltimes)),c(nbins,KMAX,length(alltimes)))
dimnames(lwc_dist)[[3]] <- alltimes

##### initialize a list to store data of all cloudy grid
lwc_list <- list()
for(t in 1:length(alltimes) ){
      lwc_list[[t]] <- list()
      for (idxk in 1:KMAX){
          lwc_list[[t]][[idxk]] <- numeric(0)
      }
}
names(lwc_list) <- alltimes

##### loop of MPI rank
    for(mpirank in allmpiranks){
	cat(sprintf("processing the rank = %s \n",mpirank))

	##### Open
	ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

	##### Grids
	IMAX <- length(ncvar_get(ncin,"CX"))-4
	JMAX <- length(ncvar_get(ncin,"CY"))-4
	KMAX <- length(ncvar_get(ncin,"CZ"))-4
	CX <- ncvar_get(ncin,"CX")
	CY <- ncvar_get(ncin,"CY")
	CZ <- ncvar_get(ncin,"CZ")
      	CDX <- ncvar_get(ncin,"CDX")
        CDY <- ncvar_get(ncin,"CDY")
        CDZ <- ncvar_get(ncin,"CDZ")

        ##### Read QR
        QR <- ncvar_get(ncin,"QR_sd")
        QR_units <- ncatt_get(ncin,"QR_sd","units")
        dimnames(QR)[[4]] <- alltimes
        if(!setequal(dim(QR[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read QC
        QC <- ncvar_get(ncin,"QC_sd")
        QC_units <- ncatt_get(ncin,"QC_sd","units")
        dimnames(QC)[[4]] <- alltimes
        if(!setequal(dim(QC[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read DENS
        DENS <- ncvar_get(ncin,"DENS")
        DENS_units <- ncatt_get(ncin,"DENS","units")
        dimnames(DENS)[[4]] <- alltimes
        if(!setequal(dim(DENS[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read QHYD_sd
        QHYD <- ncvar_get(ncin,"QHYD_sd")
        QHYD_units <- ncatt_get(ncin,"QHYD_sd","units")
        dimnames(QHYD)[[4]] <- alltimes
        if(!setequal(dim(QHYD[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Close
        nc_close(ncin)

	### Making the Histogram
#	for(time in alltimes[which(as.numeric(alltimes)==7200)]){
     	for(time in alltimes){
		 #### calculate
                 LWC <- DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
                           (QC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]+QR[1:IMAX,1:JMAX,1:KMAX,as.character(time)]) 

   		 # update the histogram
		 ilwc <- array(rep(0,IMAX*JMAX*KMAX),c(IMAX,JMAX,KMAX))
               	 ilwc[1:IMAX,1:JMAX,1:KMAX] <- findInterval(LWC[1:IMAX,1:JMAX,1:KMAX], lwc.bin, all.inside=TRUE)
                 for (idxi in 1:IMAX){
                 for (idxj in 1:JMAX){
                 for (idxk in 1:KMAX){
                    lwc_dist[ilwc[idxi,idxj,idxk],idxk,as.character(time)] <- lwc_dist[ilwc[idxi,idxj,idxk],idxk,as.character(time)] +
                                                      	CDX[idxi+2]*CDY[idxj+2]*CDZ[idxk+2]/CDZ[idxk+2]/d_lwc.bin
                 }
                 }
                 }
                 # append cloudy grid data to the list
                 for (idxk in 1:KMAX){
                    tmp_lwc   <- LWC[1:IMAX,1:JMAX,idxk]
                    tmp_QHYD <- QHYD[1:IMAX,1:JMAX,idxk,as.character(time)]
                    lwc_list[[as.character(time)]][[idxk]] <- append(lwc_list[[as.character(time)]][[idxk]], tmp_lwc[tmp_QHYD>crt_qhyd])
                 }

	}
    }

#### Plot Data
#    for(time in alltimes[which(as.numeric(alltimes)==7200)]){
     for(time in alltimes){
        ## calculate percentile
        ptile05 <- numeric(KMAX)
        ptile50 <- numeric(KMAX)
        ptile95 <- numeric(KMAX)
        for (idxk in 1:KMAX){
            ptile05[idxk] <- quantile(lwc_list[[as.character(time)]][[idxk]],0.05)
            ptile50[idxk] <- quantile(lwc_list[[as.character(time)]][[idxk]],0.50)
            ptile95[idxk] <- quantile(lwc_list[[as.character(time)]][[idxk]],0.95)
        }

        pdf(paste("lwc.vprof.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
    	image.plot(lwc.bin, CZ[3:(KMAX+2)], log10(lwc_dist[,,as.character(time)]),
                   zlim=c(col_scale_min,col_scale_max),
                   col = pal,
                   main=paste("Liquid Water Content Distribution (T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)\n(Volume Density log10([m^3/(kg/m^3)/m]))"),
                   xlab="Liquid Water Content [kg/m^3]",
                   ylab="Height z [m]",
		   #ann=F,
		   useRaster=TRUE)

#       lines(ptile05,CZ[3:(KMAX+2)])
        lines(ptile50,CZ[3:(KMAX+2)])
#       lines(ptile95,CZ[3:(KMAX+2)])

    	#### Close all the files
    	graphics.off()
    }

############################################################################
#################### Plot the standard deviation of droplets
# Setting of bins
nbins = 100
min_stddev <-   0.0    # [m]
max_stddev <-   8.0e-6 # [m]

stddev.bin = seq(min_stddev, max_stddev, length=nbins)
d_stddev.bin = stddev.bin[2] - stddev.bin[1]

col_scale_min <- log10(1.0e10)
col_scale_max <- log10(1.0e13)

##### initialize histogram
stddev_dist <- array(rep(1.0e-100,nbins*KMAX*length(alltimes)),c(nbins,KMAX,length(alltimes)))
dimnames(stddev_dist)[[3]] <- alltimes

##### initialize a list to store data of all cloudy grid
stddev_list <- list()
for(t in 1:length(alltimes) ){
      stddev_list[[t]] <- list()
      for (idxk in 1:KMAX){
          stddev_list[[t]][[idxk]] <- numeric(0)
      }
}
names(stddev_list) <- alltimes

##### loop of MPI rank
    for(mpirank in allmpiranks){
	cat(sprintf("processing the rank = %s \n",mpirank))

	##### Open
	ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

	##### Grids
	IMAX <- length(ncvar_get(ncin,"CX"))-4
	JMAX <- length(ncvar_get(ncin,"CY"))-4
	KMAX <- length(ncvar_get(ncin,"CZ"))-4
	CX <- ncvar_get(ncin,"CX")
	CY <- ncvar_get(ncin,"CY")
	CZ <- ncvar_get(ncin,"CZ")
      	CDX <- ncvar_get(ncin,"CDX")
        CDY <- ncvar_get(ncin,"CDY")
        CDZ <- ncvar_get(ncin,"CDZ")

       ##### Read SNC (Tentatively, this is number concentration of real droplets)
        DROP_NUM_CONC <- ncvar_get(ncin,"SNC")
        DROP_NUM_CONC_units <- ncatt_get(ncin,"SNC","units")
        dimnames(DROP_NUM_CONC)[[4]] <- alltimes
        if(!setequal(dim(DROP_NUM_CONC[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read Re_QC (Tentatively, this is the 1st momentum)
        MOM_1st <- ncvar_get(ncin,"Re_QC")
        MOM_1st_units <- ncatt_get(ncin,"Re_QC","units")
        dimnames(MOM_1st)[[4]] <- alltimes
        if(!setequal(dim(MOM_1st[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read Re_QR (Tentatively, this is the 2nd momentum)
        MOM_2nd <- ncvar_get(ncin,"Re_QR")
        MOM_2nd_units <- ncatt_get(ncin,"Re_QR","units")
        dimnames(MOM_2nd)[[4]] <- alltimes
        if(!setequal(dim(MOM_2nd[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read QHYD_sd
        QHYD <- ncvar_get(ncin,"QHYD_sd")
        QHYD_units <- ncatt_get(ncin,"QHYD_sd","units")
        dimnames(QHYD)[[4]] <- alltimes
        if(!setequal(dim(QHYD[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Close
        nc_close(ncin)

	### Making the Histogram
#	for(time in alltimes[which(as.numeric(alltimes)==7200)]){
     	for(time in alltimes){
		 #### calculate
                 MOM_0 <- DROP_NUM_CONC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
                          (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)]) # volume of each grid box
                 STD_DEV <- sqrt(pmax( MOM_2nd[1:IMAX,1:JMAX,1:KMAX,as.character(time)]/MOM_0[1:IMAX,1:JMAX,1:KMAX] -
                                      (MOM_1st[1:IMAX,1:JMAX,1:KMAX,as.character(time)]/MOM_0[1:IMAX,1:JMAX,1:KMAX])**2,
                                      0.0))

   		 # update the histogram
		 istddev <- array(rep(0,IMAX*JMAX*KMAX),c(IMAX,JMAX,KMAX))
               	 istddev[1:IMAX,1:JMAX,1:KMAX] <- findInterval(STD_DEV[1:IMAX,1:JMAX,1:KMAX], stddev.bin, all.inside=TRUE)
                 for (idxi in 1:IMAX){
                 for (idxj in 1:JMAX){
                 for (idxk in 1:KMAX){
                    stddev_dist[istddev[idxi,idxj,idxk],idxk,as.character(time)] <- stddev_dist[istddev[idxi,idxj,idxk],idxk,as.character(time)] +
                                                      	CDX[idxi+2]*CDY[idxj+2]*CDZ[idxk+2]/CDZ[idxk+2]/d_stddev.bin
                 }
                 }
                 }
                 # append cloudy grid data to the list
                 for (idxk in 1:KMAX){
                    tmp_stddev   <- STD_DEV[1:IMAX,1:JMAX,idxk]
                    tmp_QHYD <- QHYD[1:IMAX,1:JMAX,idxk,as.character(time)]
                    stddev_list[[as.character(time)]][[idxk]] <- append(stddev_list[[as.character(time)]][[idxk]], tmp_stddev[tmp_QHYD>crt_qhyd])
                 }

	}
    }

#### Plot Data
#    for(time in alltimes[which(as.numeric(alltimes)==7200)]){
     for(time in alltimes){
        ## calculate percentile
        ptile05 <- numeric(KMAX)
        ptile50 <- numeric(KMAX)
        ptile95 <- numeric(KMAX)
        for (idxk in 1:KMAX){
            ptile05[idxk] <- quantile(stddev_list[[as.character(time)]][[idxk]],0.05,na.rm=TRUE)
            ptile50[idxk] <- quantile(stddev_list[[as.character(time)]][[idxk]],0.50,na.rm=TRUE)
            ptile95[idxk] <- quantile(stddev_list[[as.character(time)]][[idxk]],0.95,na.rm=TRUE)
        }

        pdf(paste("stddev.vprof.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
    	image.plot(stddev.bin, CZ[3:(KMAX+2)], log10(stddev_dist[,,as.character(time)]),
                   zlim=c(col_scale_min,col_scale_max),
                   col = pal,
                   main=paste("Standard Deviation Distribution (T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)\n(Volume Density log10([m^3/(unit k)/m]))"),
                   xlab="Standard Deviation [m]",
                   ylab="Height z [m]",
		   #ann=F,
		   useRaster=TRUE)

#       lines(ptile05,CZ[3:(KMAX+2)])
        lines(ptile50,CZ[3:(KMAX+2)])
#       lines(ptile95,CZ[3:(KMAX+2)])

    	#### Close all the files
    	graphics.off()
    }

############################################################################
#################### Plot number concentration of droplets
# Setting of bins
nbins = 100
min_nconc <-  1.0e05 # [/m^3]
max_nconc <-  1.0e10 # [/m^3]

nconc.bin = seq(log10(min_nconc), log10(max_nconc), length=nbins)
d_nconc.bin = nconc.bin[2] - nconc.bin[1]

col_scale_min <- log10(1.0e5)
col_scale_max <- log10(1.0e7)

##### initialize histogram
nconc_dist <- array(rep(1.0e-100,nbins*KMAX*length(alltimes)),c(nbins,KMAX,length(alltimes)))
dimnames(nconc_dist)[[3]] <- alltimes

##### initialize a list to store data of all cloudy grid
nconc_list <- list()
for(t in 1:length(alltimes) ){
      nconc_list[[t]] <- list()
      for (idxk in 1:KMAX){
          nconc_list[[t]][[idxk]] <- numeric(0)
      }
}
names(nconc_list) <- alltimes

##### loop of MPI rank
    for(mpirank in allmpiranks){
	cat(sprintf("processing the rank = %s \n",mpirank))

	##### Open
	ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

	##### Grids
	IMAX <- length(ncvar_get(ncin,"CX"))-4
	JMAX <- length(ncvar_get(ncin,"CY"))-4
	KMAX <- length(ncvar_get(ncin,"CZ"))-4
	CX <- ncvar_get(ncin,"CX")
	CY <- ncvar_get(ncin,"CY")
	CZ <- ncvar_get(ncin,"CZ")
      	CDX <- ncvar_get(ncin,"CDX")
        CDY <- ncvar_get(ncin,"CDY")
        CDZ <- ncvar_get(ncin,"CDZ")

        ##### Read SNC (Tentatively, this is number concentration of real droplets)
        DROP_NUM_CONC <- ncvar_get(ncin,"SNC")
        DROP_NUM_CONC_units <- ncatt_get(ncin,"SNC","units")
        dimnames(DROP_NUM_CONC)[[4]] <- alltimes
        if(!setequal(dim(DROP_NUM_CONC[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read QHYD_sd
        QHYD <- ncvar_get(ncin,"QHYD_sd")
        QHYD_units <- ncatt_get(ncin,"QHYD_sd","units")
        dimnames(QHYD)[[4]] <- alltimes
        if(!setequal(dim(QHYD[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Close
        nc_close(ncin)

	### Making the Histogram
#	for(time in alltimes[which(as.numeric(alltimes)==7200)]){
     	for(time in alltimes){
   		 # update the histogram
		 inconc <- array(rep(0,IMAX*JMAX*KMAX),c(IMAX,JMAX,KMAX))
               	 inconc[1:IMAX,1:JMAX,1:KMAX] <- findInterval(log10(DROP_NUM_CONC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]), nconc.bin, all.inside=TRUE)
                 for (idxi in 1:IMAX){
                 for (idxj in 1:JMAX){
                 for (idxk in 1:KMAX){
                    nconc_dist[inconc[idxi,idxj,idxk],idxk,as.character(time)] <- nconc_dist[inconc[idxi,idxj,idxk],idxk,as.character(time)] +
                                                      	CDX[idxi+2]*CDY[idxj+2]*CDZ[idxk+2]/CDZ[idxk+2]/d_nconc.bin
                 }
                 }
                 }
                 # append cloudy grid data to the list
                 for (idxk in 1:KMAX){
                    tmp_nconc   <- log10(DROP_NUM_CONC[1:IMAX,1:JMAX,idxk,as.character(time)])
                    tmp_QHYD <- QHYD[1:IMAX,1:JMAX,idxk,as.character(time)]
                    nconc_list[[as.character(time)]][[idxk]] <- append(nconc_list[[as.character(time)]][[idxk]], tmp_nconc[tmp_QHYD>crt_qhyd])
                 }

	}
    }

#### Plot Data
#    for(time in alltimes[which(as.numeric(alltimes)==7200)]){
     for(time in alltimes){
        ## calculate percentile
        ptile05 <- numeric(KMAX)
        ptile50 <- numeric(KMAX)
        ptile95 <- numeric(KMAX)
        for (idxk in 1:KMAX){
            ptile05[idxk] <- quantile(nconc_list[[as.character(time)]][[idxk]],0.05)
            ptile50[idxk] <- quantile(nconc_list[[as.character(time)]][[idxk]],0.50)
            ptile95[idxk] <- quantile(nconc_list[[as.character(time)]][[idxk]],0.95)
        }

        pdf(paste("nconc.vprof.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
    	image.plot(nconc.bin, CZ[3:(KMAX+2)], log10(nconc_dist[,,as.character(time)]),
                   zlim=c(col_scale_min,col_scale_max),
                   col = pal,
                   main=paste("Droplet Number Concentration Distribution (T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)\n(Volume Density log10([m^3/(/m^3)/m]))"),
                   xlab="Droplet Number Concentration [/m^3]",
                   ylab="Height z [m]",
		   #ann=F,
		   useRaster=TRUE)

#       lines(ptile05,CZ[3:(KMAX+2)])
        lines(ptile50,CZ[3:(KMAX+2)])
#       lines(ptile95,CZ[3:(KMAX+2)])

    	#### Close all the files
    	graphics.off()
    }

############################################################################
#################### Plot the effective radius of droplets
# Setting of bins
nbins = 100
min_re <-   0.5e-6 # [m]
max_re <-  20.0e-6 # [m]

re.bin = seq(min_re, max_re, length=nbins)
d_re.bin = re.bin[2] - re.bin[1]

col_scale_min <- log10(1.0e10)
col_scale_max <- log10(1.0e13)
 
##### initialize histogram
re_dist <- array(rep(1.0e-100,nbins*KMAX*length(alltimes)),c(nbins,KMAX,length(alltimes)))
dimnames(re_dist)[[3]] <- alltimes

##### initialize a list to store data of all cloudy grid
re_list <- list()
for(t in 1:length(alltimes) ){
      re_list[[t]] <- list()
      for (idxk in 1:KMAX){
          re_list[[t]][[idxk]] <- numeric(0)
      }
}
names(re_list) <- alltimes

##### loop of MPI rank
    for(mpirank in allmpiranks){
	cat(sprintf("processing the rank = %s \n",mpirank))

	##### Open
	ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

	##### Grids
	IMAX <- length(ncvar_get(ncin,"CX"))-4
	JMAX <- length(ncvar_get(ncin,"CY"))-4
	KMAX <- length(ncvar_get(ncin,"CZ"))-4
	CX <- ncvar_get(ncin,"CX")
	CY <- ncvar_get(ncin,"CY")
	CZ <- ncvar_get(ncin,"CZ")
      	CDX <- ncvar_get(ncin,"CDX")
        CDY <- ncvar_get(ncin,"CDY")
        CDZ <- ncvar_get(ncin,"CDZ")

        ##### Read Re_QR (Tentatively, this is the 2nd momentum)
        MOM_2nd <- ncvar_get(ncin,"Re_QR")
        MOM_2nd_units <- ncatt_get(ncin,"Re_QR","units")
        dimnames(MOM_2nd)[[4]] <- alltimes
        if(!setequal(dim(MOM_2nd[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read QR
        QR <- ncvar_get(ncin,"QR_sd")
        QR_units <- ncatt_get(ncin,"QR_sd","units")
        dimnames(QR)[[4]] <- alltimes
        if(!setequal(dim(QR[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read QC
        QC <- ncvar_get(ncin,"QC_sd")
        QC_units <- ncatt_get(ncin,"QC_sd","units")
        dimnames(QC)[[4]] <- alltimes
        if(!setequal(dim(QC[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read DENS
        DENS <- ncvar_get(ncin,"DENS")
        DENS_units <- ncatt_get(ncin,"DENS","units")
        dimnames(DENS)[[4]] <- alltimes
        if(!setequal(dim(DENS[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read QHYD_sd
        QHYD <- ncvar_get(ncin,"QHYD_sd")
        QHYD_units <- ncatt_get(ncin,"QHYD_sd","units")
        dimnames(QHYD)[[4]] <- alltimes
        if(!setequal(dim(QHYD[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Close
        nc_close(ncin)

	### Making the Histogram
#	for(time in alltimes[which(as.numeric(alltimes)==7200)]){
     	for(time in alltimes){
		 #### calculate
                 coeff <- 1000.0*4.0*pi/3.0
                 MOM_3rd <- DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
                           (QC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]+QR[1:IMAX,1:JMAX,1:KMAX,as.character(time)]) *
                           (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)]) / coeff
                 Eff_Radi<- MOM_3rd[1:IMAX,1:JMAX,1:KMAX] / MOM_2nd[1:IMAX,1:JMAX,1:KMAX,as.character(time)]

   		 # update the histogram
		 ire <- array(rep(0,IMAX*JMAX*KMAX),c(IMAX,JMAX,KMAX))
               	 ire[1:IMAX,1:JMAX,1:KMAX] <- findInterval(Eff_Radi[1:IMAX,1:JMAX,1:KMAX], re.bin, all.inside=TRUE)
                 for (idxi in 1:IMAX){
                 for (idxj in 1:JMAX){
                 for (idxk in 1:KMAX){
                    re_dist[ire[idxi,idxj,idxk],idxk,as.character(time)] <- re_dist[ire[idxi,idxj,idxk],idxk,as.character(time)] +
                                                      	CDX[idxi+2]*CDY[idxj+2]*CDZ[idxk+2]/CDZ[idxk+2]/d_re.bin
                 }
                 }
                 }
                 # append cloudy grid data to the list
                 for (idxk in 1:KMAX){
                    tmp_re   <- Eff_Radi[1:IMAX,1:JMAX,idxk]
                    tmp_QHYD <- QHYD[1:IMAX,1:JMAX,idxk,as.character(time)]
                    re_list[[as.character(time)]][[idxk]] <- append(re_list[[as.character(time)]][[idxk]], tmp_re[tmp_QHYD>crt_qhyd])
                 }

	}
    }

#### Plot Data
#    for(time in alltimes[which(as.numeric(alltimes)==7200)]){
     for(time in alltimes){
        ## calculate percentile
        ptile05 <- numeric(KMAX)
        ptile50 <- numeric(KMAX)
        ptile95 <- numeric(KMAX)
        for (idxk in 1:KMAX){
            ptile05[idxk] <- quantile(re_list[[as.character(time)]][[idxk]],0.05)
            ptile50[idxk] <- quantile(re_list[[as.character(time)]][[idxk]],0.50)
            ptile95[idxk] <- quantile(re_list[[as.character(time)]][[idxk]],0.95)
        }

        pdf(paste("re.vprof.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
    	image.plot(re.bin, CZ[3:(KMAX+2)], log10(re_dist[,,as.character(time)]),
                   zlim=c(col_scale_min,col_scale_max),
                   col = pal,
                   main=paste("Effective Radius Distribution (T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)\n(Volume Density log10([m^3/(unit k)/m]))"),
                   xlab="Effective Radius [m]",
                   ylab="Height z [m]",
		   #ann=F,
		   useRaster=TRUE)

#       lines(ptile05,CZ[3:(KMAX+2)])
        lines(ptile50,CZ[3:(KMAX+2)])
#       lines(ptile95,CZ[3:(KMAX+2)])

    	#### Close all the files
    	graphics.off()
    }

############################################################################
#################### Plot k:=<r^3>/(effective radius)^3 = <r^2>^3/<r^3>^2 = M2^3/M0/M3^2
# Setting of bins
nbins = 100
min_k <- 0.0      # []
max_k <- 1.0

k_value.bin = seq(min_k, max_k, length=nbins)
d_k_value.bin = k_value.bin[2] - k_value.bin[1]

col_scale_min <- log10(1.0e5)
col_scale_max <- log10(1.0e8)
 
##### initialize histogram
k_value_dist <- array(rep(1.0e-100,nbins*KMAX*length(alltimes)),c(nbins,KMAX,length(alltimes)))
dimnames(k_value_dist)[[3]] <- alltimes

##### initialize a list to store data of all cloudy grid
k_value_list <- list()
for(t in 1:length(alltimes) ){
      k_value_list[[t]] <- list()
      for (idxk in 1:KMAX){
          k_value_list[[t]][[idxk]] <- numeric(0)
      }
}
names(k_value_list) <- alltimes

##### loop of MPI rank
    for(mpirank in allmpiranks){
	cat(sprintf("processing the rank = %s \n",mpirank))

	##### Open
	ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

	##### Grids
	IMAX <- length(ncvar_get(ncin,"CX"))-4
	JMAX <- length(ncvar_get(ncin,"CY"))-4
	KMAX <- length(ncvar_get(ncin,"CZ"))-4
	CX <- ncvar_get(ncin,"CX")
	CY <- ncvar_get(ncin,"CY")
	CZ <- ncvar_get(ncin,"CZ")
      	CDX <- ncvar_get(ncin,"CDX")
        CDY <- ncvar_get(ncin,"CDY")
        CDZ <- ncvar_get(ncin,"CDZ")

       	##### Read SNC (Tentatively, this is number concentration of real droplets)
        DROP_NUM_CONC <- ncvar_get(ncin,"SNC")
        DROP_NUM_CONC_units <- ncatt_get(ncin,"SNC","units")
        dimnames(DROP_NUM_CONC)[[4]] <- alltimes
        if(!setequal(dim(DROP_NUM_CONC[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

       	##### Read Re_QR (Tentatively, this is the 2nd momentum)
        MOM_2nd <- ncvar_get(ncin,"Re_QR")
        MOM_2nd_units <- ncatt_get(ncin,"Re_QR","units")
        dimnames(MOM_2nd)[[4]] <- alltimes
        if(!setequal(dim(MOM_2nd[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read QR
        QR <- ncvar_get(ncin,"QR_sd")
        QR_units <- ncatt_get(ncin,"QR_sd","units")
        dimnames(QR)[[4]] <- alltimes
        if(!setequal(dim(QR[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read QC
        QC <- ncvar_get(ncin,"QC_sd")
        QC_units <- ncatt_get(ncin,"QC_sd","units")
        dimnames(QC)[[4]] <- alltimes
        if(!setequal(dim(QC[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read DENS
        DENS <- ncvar_get(ncin,"DENS")
        DENS_units <- ncatt_get(ncin,"DENS","units")
        dimnames(DENS)[[4]] <- alltimes
        if(!setequal(dim(DENS[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

        ##### Read QHYD_sd
        QHYD <- ncvar_get(ncin,"QHYD_sd")
        QHYD_units <- ncatt_get(ncin,"QHYD_sd","units")
        dimnames(QHYD)[[4]] <- alltimes
        if(!setequal(dim(QHYD[,,,1]),c(IMAX,JMAX,KMAX))){
                cat(sprintf("ERROR: Size of variables is not consistent. Check HALO treatment,\n"))
                return()
        }

	##### Close
	nc_close(ncin)

	### Making the Histogram
#	for(time in alltimes[which(as.numeric(alltimes)==7200)]){
     	for(time in alltimes){
		 #### calculate
		 MOM_0 <- DROP_NUM_CONC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
		       	  (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)]) # volume of each grid box
		 coeff <- 1000.0*4.0*pi/3.0
		 MOM_3rd <- DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]* 
		 	   (QC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]+QR[1:IMAX,1:JMAX,1:KMAX,as.character(time)]) *
		       	   (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)]) / coeff
		 k_value <- MOM_2nd[1:IMAX,1:JMAX,1:KMAX,as.character(time)]**3 / 
		 	   (MOM_0[1:IMAX,1:JMAX,1:KMAX]*MOM_3rd[1:IMAX,1:JMAX,1:KMAX]**2)

   		 # update the histogram
		 ik_value <- array(rep(0,IMAX*JMAX*KMAX),c(IMAX,JMAX,KMAX))
               	 ik_value[1:IMAX,1:JMAX,1:KMAX] <- findInterval(k_value[1:IMAX,1:JMAX,1:KMAX], k_value.bin, all.inside=TRUE)
                 for (idxi in 1:IMAX){
                 for (idxj in 1:JMAX){
                 for (idxk in 1:KMAX){
                    k_value_dist[ik_value[idxi,idxj,idxk],idxk,as.character(time)] <- k_value_dist[ik_value[idxi,idxj,idxk],idxk,as.character(time)] +
                                                      	CDX[idxi+2]*CDY[idxj+2]*CDZ[idxk+2]/CDZ[idxk+2]/d_k_value.bin
                 }
                 }
                 }
                 # append cloudy grid data to the list
                 for (idxk in 1:KMAX){
                    tmp_k_value   <- k_value[1:IMAX,1:JMAX,idxk]
                    tmp_QHYD <- QHYD[1:IMAX,1:JMAX,idxk,as.character(time)]
                    k_value_list[[as.character(time)]][[idxk]] <- append(k_value_list[[as.character(time)]][[idxk]], tmp_k_value[tmp_QHYD>crt_qhyd])
                 }
	}
    }

#### Plot Data
#    for(time in alltimes[which(as.numeric(alltimes)==7200)]){
     for(time in alltimes){
        ## calculate percentile
        ptile05 <- numeric(KMAX)
        ptile50 <- numeric(KMAX)
        ptile95 <- numeric(KMAX)
        for (idxk in 1:KMAX){
            ptile05[idxk] <- quantile(k_value_list[[as.character(time)]][[idxk]],0.05,na.rm=TRUE)
            ptile50[idxk] <- quantile(k_value_list[[as.character(time)]][[idxk]],0.50,na.rm=TRUE)
            ptile95[idxk] <- quantile(k_value_list[[as.character(time)]][[idxk]],0.95,na.rm=TRUE)
        }

        pdf(paste("k_value.vprof.",sprintf("%05d",as.numeric(time)),".pdf",sep=""))
    	image.plot(k_value.bin, CZ[3:(KMAX+2)], log10(k_value_dist[,,as.character(time)]),
                   zlim=c(col_scale_min,col_scale_max),
                   col = pal,
                   main=paste("k value Distribution (T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)\n(Volume Density log10([m^3/(unit k)/m]))"),
                   xlab="k []",
                   ylab="Height z [m]",
		   #ann=F,
		   useRaster=TRUE)

#       lines(ptile05,CZ[3:(KMAX+2)])
        lines(ptile50,CZ[3:(KMAX+2)])
#       lines(ptile95,CZ[3:(KMAX+2)])

    	#### Close all the files
    	graphics.off()
    }

