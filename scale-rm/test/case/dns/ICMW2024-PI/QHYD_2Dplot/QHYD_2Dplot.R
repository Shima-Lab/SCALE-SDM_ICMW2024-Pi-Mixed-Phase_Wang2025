###########################################################################
###
### Program to plot the x-z distribution of hydrometeor mixing ratio
###
############################################################################

############################################################################
########## Loading Libraries
library(ncdf4)

############################################################################
########## Parameters
# cross section point
CROSS_SEC_Y <- 1.0 # [m]

# read parameters
tmp <- read.table(textConnection(gsub("="," ",gsub(",", " ", readLines("../run.conf")))),header=FALSE,fill=TRUE,flush=T)

if( is.element("sdm_cold", tmp$V1) ){
    sdm_cold <- tmp[tmp$V1=="sdm_cold",2]
}else{
    sdm_cold <- ".false."
}

pal <- c(colorRampPalette(c(rgb(     0,   0,      1, 0.0), rgb(     0,   0,      1,0.6)), alpha = TRUE),
         colorRampPalette(c(rgb(     0,   1,      0, 0.0), rgb(     0,   1,      0,0.6)), alpha = TRUE),
         colorRampPalette(c(rgb(     1,   0,      0, 0.0), rgb(     1,   0,      0,0.6)), alpha = TRUE),
         colorRampPalette(c(rgb(     1,   1,      1, 0.0), rgb(     1,   1,      1,0.6)), alpha = TRUE),
         colorRampPalette(c(rgb(     1,   1,      0, 0.0), rgb(     1,   1,      0,0.6)), alpha = TRUE)
	 )
names(pal)=c("ice","snow","graupel","cloud water","rain")

zlim_min <- 1e-7
zlim_max <- 1e-2

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

############################################################################
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

w_plot_area <- (xlim_max-xlim_min)
h_plot_area <- (ylim_max-ylim_min)

############################################################################
########## Prepare a Vector for Devices
dev_num <- rep(0,length(alltimes))
names(dev_num) <- alltimes

############################################################################
########## Plot 50 files at once
skip <- as.integer(length(alltimes)/50)+1
for(i0 in 1:skip){
    #### Prepare Devices (=Files)
    for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
	cat(sprintf("processing time = %s [s]	\n",time))

	pdf(paste("QHYD_overlay.",sprintf("%05d",as.numeric(time)),".pdf",sep=""),width=w_plot_area*2+0.8,height=h_plot_area*2+0.8)
    	dev_num[as.character(time)] <- dev.cur()
	par(pin=c(w_plot_area*2,h_plot_area*2))
	par(bg="black", col="white", col.axis="white", 
	       col.lab="white", fg="white", col.main="white", col.sub="white", cex=0.4)
	image(c(0,1), c(0,1), diag(2)*0,
                   xlim=c(xlim_min,xlim_max),
                   ylim=c(ylim_min,ylim_max),
		   asp = 1,
                   zlim=c(log(zlim_min),log(zlim_max)),
                   col = rgb(1,1,1,0),
                   main=paste("Mixing Ratio of Hydrometeors (Y=",
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


	##### Read QR
	QR <- ncvar_get(ncin,"QR_sd")
	QR_units <- ncatt_get(ncin,"QR_sd","units")
	dimnames(QR)[[4]] <- alltimes

	##### Read QC
	QC <- ncvar_get(ncin,"QC_sd")
	QC_units <- ncatt_get(ncin,"QC_sd","units")
	dimnames(QC)[[4]] <- alltimes

	if(sdm_cold == ".true."){
		##### Read QG
		QG <- ncvar_get(ncin,"QG_sd")
		QG_units <- ncatt_get(ncin,"QG_sd","units")
		dimnames(QG)[[4]] <- alltimes

		##### Read QI
		QI <- ncvar_get(ncin,"QI_sd")
		QI_units <- ncatt_get(ncin,"QI_sd","units")
		dimnames(QI)[[4]] <- alltimes

		##### Read QS
		QS <- ncvar_get(ncin,"QS_sd")
		QS_units <- ncatt_get(ncin,"QS_sd","units")
		dimnames(QS)[[4]] <- alltimes
	}

	##### Close
	nc_close(ncin)

	#### Plot Data
#	for(time in alltimes[which(as.numeric(alltimes)==3600)]){
	for(time in alltimes[seq(i0,length(alltimes),by=skip)]){
		 #### set the device number
		 dev.set(dev_num[as.character(time)])

		 if(sdm_cold == ".true."){
			allhydlist <- c("graupel", "rain", "cloud water", "ice", "snow")
		 }else{
			allhydlist <- c("rain", "cloud water")
		 }

		 for(hyd_category in allhydlist){
		 	if(hyd_category=="graupel"){
				makeActiveBinding("qhyd", function() QG,.GlobalEnv)
			}else if(hyd_category=="ice"){
				makeActiveBinding("qhyd", function() QI,.GlobalEnv)
			}else if(hyd_category=="snow"){
				makeActiveBinding("qhyd", function() QS,.GlobalEnv)
			}else if(hyd_category=="rain"){
				makeActiveBinding("qhyd", function() QR,.GlobalEnv)
			}else if(hyd_category=="cloud water"){
				makeActiveBinding("qhyd", function() QC,.GlobalEnv)
			}

		 	par(pin=c(w_plot_area*2,h_plot_area*2))
		 	par(new=T)
		 	image(CX[3:(IMAX+2)], CZ[3:(KMAX+2)],
		   		log(pmax(pmin(qhyd[1:IMAX,1:1,1:KMAX,as.character(time)],zlim_max),zlim_min)),
                   		xlim=c(xlim_min,xlim_max),
                   		ylim=c(ylim_min,ylim_max),
		   		asp = 1,
				zlim=c(log(zlim_min),log(zlim_max)),
                   		col=pal[[hyd_category]](100),
                   		main=paste("Mixing Ratio of Hydrometeors (Y=",
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
    }

    #### Close all the files
    graphics.off()

}
