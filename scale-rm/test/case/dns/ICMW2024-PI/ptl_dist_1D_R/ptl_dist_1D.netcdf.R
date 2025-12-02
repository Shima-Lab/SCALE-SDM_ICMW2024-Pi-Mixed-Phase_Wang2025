## This is a R script to plot 1D number density distribution 
## using kernel density estimation
# Load required libraries
library(ncdf4)

# Make the list of files
tmp_files = dir("../",pattern="^SD_output_NetCDF_")
tmp = strsplit(tmp_files,"\\SD_output_NetCDF_|\\.000.pe")
alltimes = unique(matrix(unlist(tmp),nrow=3)[2,])
allmpiranks = unique(matrix(unlist(tmp),nrow=3)[3,])

allfiles = matrix(tmp_files,ncol=length(alltimes))
rownames(allfiles) = allmpiranks
colnames(allfiles) = alltimes

# set parameters
VALID2INVALID <- -999.0 

# read parameters
tmp <- read.table(textConnection(gsub("="," ",gsub(",", " ", readLines("../run.conf")))),header=FALSE,fill=TRUE)
DX <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="DX",2]))
DY <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="DY",2]))
DZ <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="DZ",2]))
IMAX <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="IMAX",2]))
JMAX <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="JMAX",2]))
KMAX <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="KMAX",2]))
sdm_zlower <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="sdm_zlower",2]))
sdm_zupper <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="sdm_zupper",2]))
sdm_zupper <- min(sdm_zupper,DZ*KMAX)
KMAX_sdm <- round((sdm_zupper - sdm_zlower)/DZ)
PRC_NUM_X <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="PRC_NUM_X",2])) 
PRC_NUM_Y <- as.numeric(gsub(".[dD]",".e",tmp[tmp$V1=="PRC_NUM_Y",2])) 

tmp <- read.table(textConnection(gsub("="," ",gsub(",", " ", readLines("../run.conf")))),header=FALSE,fill=TRUE,flush=T)
if( is.element("sdm_cold", tmp$V1) ){
    sdm_cold <- tmp[tmp$V1=="sdm_cold",2]
}else{
    sdm_cold <- ".false."
}

sink("numbers.txt")
sink()

# check the number of MPI processes
if(length(allmpiranks)!=(PRC_NUM_X*PRC_NUM_Y)){
	 cat(sprintf("Number of MPI processes is not consistent\n"))
	 cat(sprintf("%s != %s * %s\n",length(allmpiranks),PRC_NUM_X,PRC_NUM_Y))
	 quit()
}

# loop of time
for(time in alltimes){
#for(time in alltimes[which(alltimes=="00000101-005000")]){
	 hh <- substr(time, 10, 11)
	 mm <- substr(time, 12, 13)
	 ss <- substr(time, 14, 15)
	 time_sec <- 3600*as.integer(hh) + 60*as.integer(mm) + as.integer(ss)
	 cat(sprintf("processing the time = %s:%s:%s \n",hh,mm,ss))

	 data <- NULL
	 # loop of MPI rank
	 for(mpirank in allmpiranks){
		# Read data
	 	file <- allfiles[mpirank,time]
		ncin <- nc_open(paste("../",file,sep=""))	

		allvarnames <- attributes(ncin$var)$names
		data_tmp <- NULL
		for(varname in allvarnames){
			data_tmp <- cbind(data_tmp,ncvar_get(ncin,varname))
		}
		colnames(data_tmp) <- allvarnames
		data <- rbind(data,data_tmp)
		nc_close(ncin)
	}
	data <- data.frame(data)
	valid_data <- data[data$sd_z>VALID2INVALID,]

	# number of SDs
	SD_num <- length(valid_data$sd_x)
	SD_num_dens <- SD_num/IMAX/JMAX/KMAX_sdm/length(allmpiranks)

	# number density of non-IN aerosols (soluble, (soluble+IN inactive mineral dust))
	if(sdm_cold == ".true."){
		nonIN_data <- valid_data[valid_data$sdi_tf<(-37.99999),]
	}else{
		nonIN_data <- valid_data
	}
	nonIN_num <- sum(nonIN_data$sd_n)
	nonIN_num_dens <- nonIN_num/IMAX/JMAX/KMAX_sdm/length(allmpiranks)/DX/DY/DZ
	
	# number density of IN active aerosols
	if(sdm_cold == ".true."){
		IN_data <- valid_data[valid_data$sdi_tf>(-37.99999),]
	}else{
		IN_data <- valid_data[valid_data$sd_z<VALID2INVALID,]
	}
	IN_num <- sum(IN_data$sd_n)
	IN_num_dens <- IN_num/IMAX/JMAX/KMAX_sdm/length(allmpiranks)/DX/DY/DZ

	sink("numbers.txt", append=TRUE)
	cat(sprintf("###### SD/RD numbers at time = %s [YYYYMMDD-hhmmss] ######\n",time))
	cat(sprintf("SD number = %d\n",SD_num))
	cat(sprintf("SD number density = %g [/grid]\n",SD_num_dens))
	cat(sprintf("total aerosol number = %g\n",nonIN_num+IN_num))
	cat(sprintf("total aerosol number density = %g [/m^3]\n",nonIN_num_dens+IN_num_dens))
	cat(sprintf("non IN aerosol number = %g\n",nonIN_num))
	cat(sprintf("non IN aerosol number density = %g [/m^3]\n",nonIN_num_dens))
	cat(sprintf("IN aerosol number = %g\n",IN_num))
	cat(sprintf("IN aerosol number density = %g [/m^3]\n\n",IN_num_dens))
	sink()

	# Probability Density of Freezing Temperature
	if(sdm_cold == ".true."){
		prob_dens_tf <- density(IN_data$sdi_tf,bw="nrd",weight=IN_data$sd_n/IN_num)
		#jpeg(paste("prob_dens_freezing_temp.",time,".jpg",sep=""))
		pdf(paste("prob_dens_freezing_temp.",time,".pdf",sep=""))
        	plot(prob_dens_tf,
                   xlim=c(-38.0,-10.0),
		   main=paste("Probability Density of Freezing Temperature ([/degC]) (T=",
		              hh,":",mm,":",ss,")"),
                   xlab="Freezing temperature [degC])",
                   ylab="Probability Density [/degC])")
        	dev.off()
	}

	# Number density distribution of ammonium bisulfate\n (sum of nonIN and IN aerosols))"
	rho_amsul <- 1.78E+6 # density of ammonium bisulfate [g/m3]
	dry_radi_amsul <- (valid_data$sd_asl/(4.0/3.0)/pi/rho_amsul)**(1.0/3.0)
	num_dens_amsul <- density(log(dry_radi_amsul),bw="nrd",weight=valid_data$sd_n/(nonIN_num+IN_num))
	num_dens_amsul$y <- (nonIN_num_dens+IN_num_dens)*num_dens_amsul$y
	num_dens_amsul$x <- exp(num_dens_amsul$x)
	#jpeg(paste("num_dens_amsul.",time,".jpg",sep=""))
	pdf(paste("num_dens_amsul.",time,".pdf",sep=""))
        plot(num_dens_amsul,
		   log="x",
		   xlim=c(1e-8,1e-6),
		   main=paste("Number density distribution of ammonium bisulfate (T=",
		              hh,":",mm,":",ss,")\n (sum of nonIN and IN aerosols))"),
                   xlab="Dry radius of ammonium bisulfate [m])",
                   ylab="Number density distribution dN/dlogr ([/unit logr/m^3]")
        dev.off()

	# Mass density distribution of droplets"
	drop_data <- valid_data[valid_data$sd_liqice==1,]
	drop_radi <- drop_data$sd_r
	drop_num  <- sum(drop_data$sd_n)
	drop_num_dens <- drop_num/IMAX/JMAX/KMAX_sdm/length(allmpiranks)/DX/DY/DZ

	rho_liqw <- 1.0E+6 # density of liquid water [g/m3]
	mass_dens_drop <- density(log(drop_radi),bw="nrd",weight=drop_data$sd_n/drop_num)
	mass_dens_drop$x <- exp(mass_dens_drop$x)
	mass_dens_drop$y <- drop_num_dens*rho_liqw*(4.0/3.0)*pi*mass_dens_drop$x**3*mass_dens_drop$y
	pdf(paste("mass_dens_drop.",time,".pdf",sep=""))
        plot(mass_dens_drop,
		   log="x",
		   xlim=c(1e-9,1e-2),
		   main=paste("Mass density distribution of droplets (T=",
		              hh,":",mm,":",ss,")"),
                   xlab="Radius of droplet [m])",
                   ylab="Mass density distribution of m(r)dN/dlogr ([g/unit logr/m^3]")
        dev.off()

	# terminal velocity of droplets
#	drop_data <- valid_data[valid_data$sd_liqice==1,]
#	drop_radi <- drop_data$sd_r
	drop_termv<- drop_data$sd_vz # terminal velocity of droplets [m/s]
	n_data_skip <- max(round(length(drop_termv)/10000),1)

	pdf(paste("term_vel_drop.",time,".pdf",sep=""))
        plot(2.0*drop_radi[seq(1, length(drop_radi), n_data_skip)],drop_termv[seq(1, length(drop_termv), n_data_skip)],
		   log="x",pch='.',
		   xlim=c(1e-9,1e-1),
		   main=paste("Terminal velocity of droplets (T=",
		              hh,":",mm,":",ss,")"),
                   xlab="Diameter of droplet [m]",
                   ylab="Terminal velocity of droplets [m/s]")
        dev.off()

	# x-z position of SDs
	if(SD_num>0){
		valid_x <- valid_data$sd_x # x coordinate of valid SD [m]
		valid_z <- valid_data$sd_z # z coordinate of valid SD [m]
		n_data_skip <- max(round(length(valid_z)/10000),1)

		pdf(paste("xz_SD.",time,".pdf",sep=""))
        	plot(valid_x[seq(1, length(valid_x), n_data_skip)],valid_z[seq(1, length(valid_z), n_data_skip)],
		   pch='.',
		   xlim=c(0.0,PRC_NUM_X*DX*IMAX),
		   ylim=c(0.0,DZ*KMAX),
		   main=paste("Positions of SDs on x-z plane (T=",
		              hh,":",mm,":",ss,")"),
                   xlab="x [m]",
                   ylab="z [m]")
        	dev.off()
	}
}
