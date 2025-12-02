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
# cross section point
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

CROSS_SEC_Y <- FLUX_CY # [m]

#pal <- colorRampPalette(c("blue", "orange", "red"), space = "Lab")
pal <- rainbow(1000, start=0.4, end=1.0)

CLOUD_BASE_WIDTH <- 1.0e3 # [m]

# read parameters
tmp <- read.table(textConnection(gsub("="," ",gsub(",", " ", readLines("../run.conf")))),header=FALSE,fill=TRUE,flush=T)
if( is.element("FLUX_CX", tmp$V1) ){
    FLUX_CX <- as.numeric(as.character(tmp[tmp$V1=="FLUX_CX",2]))
}else{
    FLUX_CX <- 5000.0
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
	FX_IMAX <- length(ncvar_get(ncin,"FX"))-4
	FZ_KMAX <- length(ncvar_get(ncin,"FZ"))-4
	FX <- ncvar_get(ncin,"FX")
	FZ <- ncvar_get(ncin,"FZ")

	xlim_min <- FX[3]         # horizontal range [m]
	xlim_max <- FX[2+FX_IMAX] # [m]
	ylim_min <- FZ[3]         # vertical range [m]
	ylim_max <- FZ[2+FZ_KMAX] # [m]
}

##### read the rest of MPI ranks
tmp_allmpiranks <- character(0)
for(mpirank in allmpiranks){
        ##### Open
        ncin <- nc_open(paste("../",allfiles[mpirank],sep=""))

        ##### Grids
        FX_IMAX <- length(ncvar_get(ncin,"FX"))-4
        FZ_KMAX <- length(ncvar_get(ncin,"FZ"))-4
        FY_JMAX <- length(ncvar_get(ncin,"FY"))-4
        FX <- ncvar_get(ncin,"FX")
        FZ <- ncvar_get(ncin,"FZ")
        FY <- ncvar_get(ncin,"FY")

        xlim_min <- min(xlim_min, FX[3])
        xlim_max <- max(xlim_max, FX[2+FX_IMAX])
        ylim_min <- min(ylim_min, FZ[3])
        ylim_max <- max(ylim_max, FZ[2+FZ_KMAX])

        if( (FY[3]<CROSS_SEC_Y) && (CROSS_SEC_Y<=FY[2+FY_JMAX]) ){
            tmp_allmpiranks <- c(tmp_allmpiranks, mpirank)
        }
}

if(length(tmp_allmpiranks)==0){
        cat(sprintf("ERROR: cross section cannot be defined\n"))
        stop()
}else{
        csect_mpiranks <- tmp_allmpiranks
}

w_plot_area <- (xlim_max-xlim_min)/1.0e4
h_plot_area <- (ylim_max-ylim_min)/1.0e4

############################################################################
########## Prepare a Vector for Devices
dev_num <- rep(0,length(alltimes))
names(dev_num) <- alltimes

############################################################################
#################### Plot adiabatic fraction

########## Cloud Base Height
crt_qhyd <- 1.0e-5 # [kg/kg] critical mixing ratio to determine cloudy grid

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

#CLOUD_BASE_DENS
#CLOUD_BASE_TEMP
#CLOUD_BASE_RHOVSAT

########## plot AF
zlim_min <-  1.0e-2
zlim_max <-  1.5

########## Plot 50 files at once
#skip <- as.integer(length(alltimes)/5000)+1
skip <- as.integer(length(alltimes)/50)+1
for(i0 in 1:skip){
    #### Prepare Devices (=Files)
#    for(time in alltimes[which(as.numeric(alltimes)==7200)]){
    for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
	cat(sprintf("processing time = %s [s]	\n",time))

	pdf(paste("AF.",sprintf("%05d",as.numeric(time)),".pdf",sep=""),width=w_plot_area*2+1.6,height=h_plot_area*2+0.8)
    	dev_num[as.character(time)] <- dev.cur()
	par(pin=c(w_plot_area*2,h_plot_area*2))
        par(cex=0.4)
	image.plot(c(0,1), c(0,1), diag(2)*0,
                   xlim=c(xlim_min,xlim_max),
                   ylim=c(ylim_min,ylim_max),
		   xaxs = "i", yaxs = "i",
		   asp = 1,
                   zlim=c(zlim_min,zlim_max),
		   col=pal,
                   main=paste("Adiabatic Fraction [] (Y=",
                              sprintf("%d",as.integer(CROSS_SEC_Y)),
                              "m, T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)"),
                   xlab="x [m]",
                   ylab="z [m]")
    }

    ##### loop of MPI rank
    for(mpirank in csect_mpiranks){
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

	##### Close
	nc_close(ncin)

	#### Plot Data
#	for(time in alltimes[which(as.numeric(alltimes)==7200)]){
	for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
		 #### calculate
		 LWC <- DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]* 
		 	   (QC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]+QR[1:IMAX,1:JMAX,1:KMAX,as.character(time)])

		 RHOVSAT <- apply(TEMP[1:IMAX,1:JMAX,1:KMAX,as.character(time)], c(1,2,3), temp.to.satrhov)
		 ALWC <- DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
		      	 CLOUD_BASE_RHOVSAT[as.character(time)]/CLOUD_BASE_DENS[as.character(time)] -
                         RHOVSAT[1:IMAX,1:JMAX,1:KMAX]

		 AF <- LWC[1:IMAX,1:JMAX,1:KMAX]/ALWC[1:IMAX,1:JMAX,1:KMAX]
		 AF[1:IMAX,1:JMAX,1:CLOUD_BASE_K[as.character(time)]] <- 0.0

		 #### set the device number
		 dev.set(dev_num[as.character(time)])

	 	 par(pin=c(w_plot_area*2,h_plot_area*2))
        	 par(cex=0.4)
		 par(new=T)
		 image.plot(CX[3:(IMAX+2)], CZ[3:(KMAX+2)],
		 	AF[1:IMAX,1:1,1:KMAX],
                 	xlim=c(xlim_min,xlim_max),
                 	ylim=c(ylim_min,ylim_max),
		   	xaxs = "i", yaxs = "i",
		 	asp = 1,
			zlim=c(zlim_min,zlim_max),
                   	col=pal,
                   	main=paste("Adiabatic Fraction [] (Y=",
                              sprintf("%d",as.integer(CROSS_SEC_Y)),
                              "m, T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)"),
                   	xlab="x [m]",
                   	ylab="z [m]",
			ann=F,
			useRaster=TRUE)
	}
    }

    #### Close all the files
    graphics.off()

}

############################################################################
#################### Plot liquid water content
zlim_min <-  1.0e-5
zlim_max <-  8.0e-3

########## Plot 50 files at once
#skip <- as.integer(length(alltimes)/5000)+1
skip <- as.integer(length(alltimes)/50)+1
for(i0 in 1:skip){
    #### Prepare Devices (=Files)
#    for(time in alltimes[which(as.numeric(alltimes)==7200)]){
    for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
	cat(sprintf("processing time = %s [s]	\n",time))

	pdf(paste("LWC.",sprintf("%05d",as.numeric(time)),".pdf",sep=""),width=w_plot_area*2+1.6,height=h_plot_area*2+0.8)
    	dev_num[as.character(time)] <- dev.cur()
	par(pin=c(w_plot_area*2,h_plot_area*2))
        par(cex=0.4)
	image.plot(c(0,1), c(0,1), diag(2)*0,
                   xlim=c(xlim_min,xlim_max),
                   ylim=c(ylim_min,ylim_max),
		   xaxs = "i", yaxs = "i",
		   asp = 1,
                   zlim=c(zlim_min,zlim_max),
		   col=pal,
                   main=paste("Liquid Water Content [kg/m^3] (Y=",
                              sprintf("%d",as.integer(CROSS_SEC_Y)),
                              "m, T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)"),
                   xlab="x [m]",
                   ylab="z [m]")
    }

    ##### loop of MPI rank
    for(mpirank in csect_mpiranks){
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

	##### Close
	nc_close(ncin)

	#### Plot Data
#	for(time in alltimes[which(as.numeric(alltimes)==7200)]){
	for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
		 #### calculate
		 LWC <- DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]* 
		 	   (QC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]+QR[1:IMAX,1:JMAX,1:KMAX,as.character(time)])

		 #### set the device number
		 dev.set(dev_num[as.character(time)])

	 	 par(pin=c(w_plot_area*2,h_plot_area*2))
        	 par(cex=0.4)
		 par(new=T)
		 image.plot(CX[3:(IMAX+2)], CZ[3:(KMAX+2)],
		 	pmin(LWC[1:IMAX,1:1,1:KMAX],zlim_max),
                 	xlim=c(xlim_min,xlim_max),
                 	ylim=c(ylim_min,ylim_max),
		   	xaxs = "i", yaxs = "i",
		 	asp = 1,
			zlim=c(zlim_min,zlim_max),
                   	col=pal,
                   	main=paste("Liquid Water Content [kg/m^3] (Y=",
                              sprintf("%d",as.integer(CROSS_SEC_Y)),
                              "m, T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)"),
                   	xlab="x [m]",
                   	ylab="z [m]",
			ann=F,
			useRaster=TRUE)
	}
    }

    #### Close all the files
    graphics.off()

}

############################################################################
#################### Plot k:=<r^3>/(effective radius)^3 = <r^2>^3/<r^3>^2 = M2^3/M0/M3^2
zlim_min <-  0.6
zlim_max <-  1.0

########## Plot 50 files at once
#skip <- as.integer(length(alltimes)/5000)+1
skip <- as.integer(length(alltimes)/50)+1
for(i0 in 1:skip){
    #### Prepare Devices (=Files)
#    for(time in alltimes[which(as.numeric(alltimes)==7200)]){
    for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
	cat(sprintf("processing time = %s [s]	\n",time))

	pdf(paste("k_value.",sprintf("%05d",as.numeric(time)),".pdf",sep=""),width=w_plot_area*2+1.6,height=h_plot_area*2+0.8)
    	dev_num[as.character(time)] <- dev.cur()
	par(pin=c(w_plot_area*2,h_plot_area*2))
        par(cex=0.4)
	image.plot(c(0,1), c(0,1), diag(2)*0,
                   xlim=c(xlim_min,xlim_max),
                   ylim=c(ylim_min,ylim_max),
		   xaxs = "i", yaxs = "i",
		   asp = 1,
                   zlim=c(zlim_min,zlim_max),
		   col=pal,
                   main=paste("k [] (Y=",
                              sprintf("%d",as.integer(CROSS_SEC_Y)),
                              "m, T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)"),
                   xlab="x [m]",
                   ylab="z [m]")
    }

    ##### loop of MPI rank
    for(mpirank in csect_mpiranks){
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

	##### Close
	nc_close(ncin)

	#### Plot Data
#	for(time in alltimes[which(as.numeric(alltimes)==7200)]){
	for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
		 #### calculate
		 MOM_0 <- DROP_NUM_CONC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
		       	  (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)]) # volume of each grid box
		 coeff <- 1000.0*4.0*pi/3.0
		 MOM_3rd <- DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]* 
		 	   (QC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]+QR[1:IMAX,1:JMAX,1:KMAX,as.character(time)]) *
		       	   (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)]) / coeff
		 k_value <- MOM_2nd[1:IMAX,1:1,1:KMAX,as.character(time)]**3 / 
		 	   (MOM_0[1:IMAX,1:1,1:KMAX]*MOM_3rd[1:IMAX,1:1,1:KMAX]**2)

		 #### set the device number
		 dev.set(dev_num[as.character(time)])

	 	 par(pin=c(w_plot_area*2,h_plot_area*2))
        	 par(cex=0.4)
		 par(new=T)
		 image.plot(CX[3:(IMAX+2)], CZ[3:(KMAX+2)],
		 	pmax(k_value,zlim_min),
                 	xlim=c(xlim_min,xlim_max),
                 	ylim=c(ylim_min,ylim_max),
		   	xaxs = "i", yaxs = "i",
		 	asp = 1,
			zlim=c(zlim_min,zlim_max),
                   	col=pal,
                   	main=paste("k [] (Y=",
                              sprintf("%d",as.integer(CROSS_SEC_Y)),
                              "m, T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)"),
                   	xlab="x [m]",
                   	ylab="z [m]",
			ann=F,
			useRaster=TRUE)
	}
    }

    #### Close all the files
    graphics.off()

}

############################################################################
#################### Plot the effective radius of droplets
zlim_min <-  0.5e-6
zlim_max <-  20.0e-6

########## Plot 50 files at once
#skip <- as.integer(length(alltimes)/5000)+1
skip <- as.integer(length(alltimes)/50)+1
for(i0 in 1:skip){
    #### Prepare Devices (=Files)
#    for(time in alltimes[which(as.numeric(alltimes)==7200)]){
    for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
	cat(sprintf("processing time = %s [s]	\n",time))

	pdf(paste("eff_radi.",sprintf("%05d",as.numeric(time)),".pdf",sep=""),width=w_plot_area*2+1.6,height=h_plot_area*2+0.8)
    	dev_num[as.character(time)] <- dev.cur()
	par(pin=c(w_plot_area*2,h_plot_area*2))
        par(cex=0.4)
	image.plot(c(0,1), c(0,1), diag(2)*0,
                   xlim=c(xlim_min,xlim_max),
                   ylim=c(ylim_min,ylim_max),
		   xaxs = "i", yaxs = "i",
		   asp = 1,
                   zlim=c(zlim_min,zlim_max),
		   col=pal,
                   main=paste("Effective Radius [m] (Y=",
                              sprintf("%d",as.integer(CROSS_SEC_Y)),
                              "m, T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)"),
                   xlab="x [m]",
                   ylab="z [m]")
    }

    ##### loop of MPI rank
    for(mpirank in csect_mpiranks){
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

	##### Close
	nc_close(ncin)

       	##### Calculate the 3rd moment

	#### Plot Data
#	for(time in alltimes[which(as.numeric(alltimes)==7200)]){
	for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
		 #### calculate
		 coeff <- 1000.0*4.0*pi/3.0
		 MOM_3rd <- DENS[1:IMAX,1:JMAX,1:KMAX,as.character(time)]* 
		 	   (QC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]+QR[1:IMAX,1:JMAX,1:KMAX,as.character(time)]) *
		       	   (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)]) / coeff
		 Eff_Radi<- MOM_3rd[1:IMAX,1:1,1:KMAX] / MOM_2nd[1:IMAX,1:1,1:KMAX,as.character(time)]

		 #### set the device number
		 dev.set(dev_num[as.character(time)])

	 	 par(pin=c(w_plot_area*2,h_plot_area*2))
        	 par(cex=0.4)
		 par(new=T)
		 image.plot(CX[3:(IMAX+2)], CZ[3:(KMAX+2)],
		 	Eff_Radi,
                 	xlim=c(xlim_min,xlim_max),
                 	ylim=c(ylim_min,ylim_max),
		   	xaxs = "i", yaxs = "i",
		 	asp = 1,
			zlim=c(zlim_min,zlim_max),
                   	col=pal,
                   	main=paste("Effective Radius [m] (Y=",
                              sprintf("%d",as.integer(CROSS_SEC_Y)),
                              "m, T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)"),
                   	xlab="x [m]",
                   	ylab="z [m]",
			ann=F,
			useRaster=TRUE)
	}
    }

    #### Close all the files
    graphics.off()

}

############################################################################
#################### Plot the standard deviation of droplets
zlim_min <-  0.1e-6 # [m]
zlim_max <-  8.0e-6 # [m]

########## Plot 50 files at once
#skip <- as.integer(length(alltimes)/5000)+1
skip <- as.integer(length(alltimes)/50)+1
for(i0 in 1:skip){
    #### Prepare Devices (=Files)
#    for(time in alltimes[which(as.numeric(alltimes)==7200)]){
    for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
	cat(sprintf("processing time = %s [s]	\n",time))

	pdf(paste("std_dev.",sprintf("%05d",as.numeric(time)),".pdf",sep=""),width=w_plot_area*2+1.6,height=h_plot_area*2+0.8)
    	dev_num[as.character(time)] <- dev.cur()
	par(pin=c(w_plot_area*2,h_plot_area*2))
        par(cex=0.4)
	image.plot(c(0,1), c(0,1), diag(2)*0,
                   xlim=c(xlim_min,xlim_max),
                   ylim=c(ylim_min,ylim_max),
		   xaxs = "i", yaxs = "i",
		   asp = 1,
                   zlim=c(zlim_min,zlim_max),
		   col=pal,
                   main=paste("Starndard Deviation of Droplet Radius [m] (Y=",
                              sprintf("%d",as.integer(CROSS_SEC_Y)),
                              "m, T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)"),
                   xlab="x [m]",
                   ylab="z [m]")
    }

    ##### loop of MPI rank
    for(mpirank in csect_mpiranks){
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

	##### Close
	nc_close(ncin)

	#### Plot Data
#	for(time in alltimes[which(as.numeric(alltimes)==7200)]){
	for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
		 #### calculate
		 MOM_0 <- DROP_NUM_CONC[1:IMAX,1:JMAX,1:KMAX,as.character(time)]*
		       	  (CDX[3:(IMAX+2)]%o%CDY[3:(JMAX+2)]%o%CDZ[3:(KMAX+2)]) # volume of each grid box
		 STD_DEV <- sqrt(pmax( MOM_2nd[1:IMAX,1:1,1:KMAX,as.character(time)]/MOM_0[1:IMAX,1:1,1:KMAX] -
		 	    	      (MOM_1st[1:IMAX,1:1,1:KMAX,as.character(time)]/MOM_0[1:IMAX,1:1,1:KMAX])**2,
				      0.0))
		 #### set the device number
		 dev.set(dev_num[as.character(time)])

	 	 par(pin=c(w_plot_area*2,h_plot_area*2))
        	 par(cex=0.4)
		 par(new=T)
		 image.plot(CX[3:(IMAX+2)], CZ[3:(KMAX+2)],
		 	pmin(STD_DEV,zlim_max),
                 	xlim=c(xlim_min,xlim_max),
                 	ylim=c(ylim_min,ylim_max),
		   	xaxs = "i", yaxs = "i",
		 	asp = 1,
			zlim=c(zlim_min,zlim_max),
                   	col=pal,
                   	main=paste("Starndard Deviation of Droplet Radius [m] (Y=",
                              sprintf("%d",as.integer(CROSS_SEC_Y)),
                              "m, T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)"),
                   	xlab="x [m]",
                   	ylab="z [m]",
			ann=F,
			useRaster=TRUE)
	}
    }

    #### Close all the files
    graphics.off()

}

############################################################################
#################### Plot number concentration of droplets
zlim_min <- 1.0e5
zlim_max <- 1.0e10

########## Plot 50 files at once
#skip <- as.integer(length(alltimes)/5000)+1
skip <- as.integer(length(alltimes)/50)+1
for(i0 in 1:skip){
    #### Prepare Devices (=Files)
#    for(time in alltimes[which(as.numeric(alltimes)==7200)]){
    for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
	cat(sprintf("processing time = %s [s]	\n",time))

	pdf(paste("num_conc.",sprintf("%05d",as.numeric(time)),".pdf",sep=""),width=w_plot_area*2+1.6,height=h_plot_area*2+0.8)
    	dev_num[as.character(time)] <- dev.cur()
	par(pin=c(w_plot_area*2,h_plot_area*2))
        par(cex=0.4)
	image.plot(c(0,1), c(0,1), diag(2)*0,
                   xlim=c(xlim_min,xlim_max),
                   ylim=c(ylim_min,ylim_max),
		   xaxs = "i", yaxs = "i",
		   asp = 1,
                   zlim=log10(c(zlim_min,zlim_max)),
		   col=pal,
                   main=paste("Number Concentration of Cloud and Rain Droplets log10([#/m^3]) (Y=",
                              sprintf("%d",as.integer(CROSS_SEC_Y)),
                              "m, T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)"),
                   xlab="x [m]",
                   ylab="z [m]")
    }

    ##### loop of MPI rank
    for(mpirank in csect_mpiranks){
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
	
       	##### Read SNC (Tentatively, this is number concentration of real droplets)
        DROP_NUM_CONC <- ncvar_get(ncin,"SNC")
        DROP_NUM_CONC_units <- ncatt_get(ncin,"SNC","units")
        dimnames(DROP_NUM_CONC)[[4]] <- alltimes

	##### Close
	nc_close(ncin)

	#### Plot Data
#	for(time in alltimes[which(as.numeric(alltimes)==7200)]){
	for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
		 #### set the device number
		 dev.set(dev_num[as.character(time)])

	 	 par(pin=c(w_plot_area*2,h_plot_area*2))
        	 par(cex=0.4)
		 par(new=T)
		 image.plot(CX[3:(IMAX+2)], CZ[3:(KMAX+2)],
		 	log10(pmin(DROP_NUM_CONC[1:IMAX,1:1,1:KMAX,as.character(time)],zlim_max)),
                 	xlim=c(xlim_min,xlim_max),
                 	ylim=c(ylim_min,ylim_max),
		   	xaxs = "i", yaxs = "i",
		 	asp = 1,
			zlim=log10(c(zlim_min,zlim_max)),
                   	col=pal,
                   	main=paste("Number Concentration of Cloud and Rain Droplets log10([#/m^3]) (Y=",
                              sprintf("%d",as.integer(CROSS_SEC_Y)),
                              "m, T=",
                              sprintf("%05d",as.numeric(time)),
                              "s)"),
                   	xlab="x [m]",
                   	ylab="z [m]",
			ann=F,
			useRaster=TRUE)
	}
    }

    #### Close all the files
    graphics.off()

}

