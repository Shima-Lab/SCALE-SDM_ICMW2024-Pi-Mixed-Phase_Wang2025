## This is a R script to plot 1D number density distribution 
## using kernel density estimation
# Load required libraries
library(ncdf4)
library(plotrix)
library(doParallel)

# Set the number of parallelization
cluster_num <- 8

# Make the list of files
tmp_files = dir("../",pattern="^SD_selected_history.")
tmp = strsplit(tmp_files,"\\SD_selected_history.")
#tmp_files = dir("../",pattern="^SD_all_history.")
#tmp = strsplit(tmp_files,"\\SD_all_history.")
allmpiranks = unique(matrix(unlist(tmp),nrow=2)[2,])

allfiles = tmp_files
names(allfiles) = allmpiranks

# Read times
mpirank <- allmpiranks[1]
file <- allfiles[mpirank]
ncin <- nc_open(paste("../",file,sep=""))	
varname <- "sd_dmp_time"
tmp <- ncvar_get(ncin,varname)
alltimes <- tmp
nc_close(ncin)

##### setup plot time
plottimes <- alltimes
#plottimes <- alltimes[which(as.numeric(alltimes)==7200)]

# set parameters
VALID2INVALID <- -999.0 

# read parameters
tmp <- read.table(textConnection(gsub("="," ",gsub(",", " ", readLines("../run.conf")))),header=FALSE,fill=TRUE)
DX <- as.numeric(gsub("[dD]","e",tmp[tmp$V1=="DX",2]))
DY <- as.numeric(gsub("[dD]","e",tmp[tmp$V1=="DY",2]))
DZ <- as.numeric(gsub("[dD]","e",tmp[tmp$V1=="DZ",2]))
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

# check the number of MPI processes
if(length(allmpiranks)!=(PRC_NUM_X*PRC_NUM_Y)){
	 cat(sprintf("Number of MPI processes is not consistent\n"))
	 cat(sprintf("%s != %s * %s\n",length(allmpiranks),PRC_NUM_X,PRC_NUM_Y))
	 quit()
}

cluster <- makeCluster(cluster_num)
clusterEvalQ(cluster, library(ncdf4))
clusterApply(cluster,1:cluster_num, function(x){
	cluster_id <<- x
})
clusterEvalQ(cluster,{
	tmp_filename <- sprintf("numbers_%03d.txt",cluster_id)
	sink(tmp_filename)
	sink()
})

# loop of time
registerDoParallel(cluster)
#for(time in plottimes){
foreach(time = plottimes, .packages="plotrix") %dopar% {
	 time_count <- which(alltimes==time) 
	 cat(sprintf("processing the time = %s [s]\n",time))

	 data <- NULL
	 # loop of MPI rank
	 for(mpirank in allmpiranks){
		# Read data
	 	file <- allfiles[mpirank]
		ncin <- nc_open(paste("../",file,sep=""))	

		allvarnames <- attributes(ncin$var)$names
		allvarnames <- allvarnames[grep(sprintf("_%04d",time_count),allvarnames)]
		data_tmp <- NULL
		varname <- allvarnames[1]
		data_tmp <- cbind(data_tmp,ncvar_get(ncin,varname))
		sd_num <- length(data_tmp)
		if( sd_num == 0 ){
		    dim(data_tmp) <- c(0,length(allvarnames))
		}else{
		    for(varname in allvarnames[2:length(allvarnames)]){
			data_tmp <- cbind(data_tmp,ncvar_get(ncin,varname))
		    }
		}
		colnames(data_tmp) <- strsplit(allvarnames[grep(time_count,allvarnames)],sprintf("_%04d",time_count))
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

	tmp_filename <- sprintf("numbers_%03d.txt",cluster_id)
	sink(tmp_filename, append=TRUE)
	cat(sprintf("###### SD/RD numbers at time = %s [s] ######\n",time))
	cat(sprintf("SD number = %d\n",SD_num))
	cat(sprintf("SD number density = %g [/grid]\n",SD_num_dens))
	cat(sprintf("total aerosol number = %g\n",nonIN_num+IN_num))
	cat(sprintf("total aerosol number density = %g [/m^3]\n",nonIN_num_dens+IN_num_dens))
	cat(sprintf("non IN aerosol number = %g\n",nonIN_num))
	cat(sprintf("non IN aerosol number density = %g [/m^3]\n",nonIN_num_dens))
	cat(sprintf("IN aerosol number = %g\n",IN_num))
	cat(sprintf("IN aerosol number density = %g [/m^3]\n\n",IN_num_dens))
	sink(file=NULL)

	# Probability Density of Freezing Temperature
	if(sdm_cold == ".true."){
		if(length(IN_data$sdi_tf)>1){
	                if(var(IN_data$sdi_tf)==0){
				prob_dens_tf <- density(IN_data$sdi_tf,bw=1.0e-2,weight=IN_data$sd_n/IN_num)
			}else{
				prob_dens_tf <- density(IN_data$sdi_tf,bw="nrd",weight=IN_data$sd_n/IN_num)
			}
		}else{
			prob_dens_tf <- density(c(-1000.0,1000.0),bw="nrd")
			prob_dens_tf$y <- 0.0*prob_dens_tf$y
		}
		#jpeg(paste("prob_dens_freezing_temp.",sprintf("%05d",time),".jpg",sep=""))
		pdf(paste("prob_dens_freezing_temp.",sprintf("%05d",time),".pdf",sep=""))
        	plot(prob_dens_tf,
                   xlim=c(-38.0,-10.0),
		   main=paste("Probability Density of Freezing Temperature ([/degC]) (T=",
		              sprintf("%05d",time)," [s] )"),
                   xlab="Freezing temperature [degC])",
                   ylab="Probability Density [/degC])")
        	dev.off()
	}

	# Number density distribution of ammonium bisulfate\n (sum of nonIN and IN aerosols))"
	rho_amsul <- 1.78E+6 # density of ammonium bisulfate [g/m3]
	dry_radi_amsul <- (valid_data$sd_asl/(4.0/3.0)/pi/rho_amsul)**(1.0/3.0)
	if(length(dry_radi_amsul)>1){
		if(var(dry_radi_amsul)==0){
			num_dens_amsul <- density(log(dry_radi_amsul),width=1.0e-2,weight=valid_data$sd_n/(nonIN_num+IN_num))
		}else{
			num_dens_amsul <- density(log(dry_radi_amsul),bw="nrd",weight=valid_data$sd_n/(nonIN_num+IN_num))
		}
	}else{
		num_dens_amsul <- density(log(c(1.0e-10,1.0e2)),bw="nrd")
		num_dens_amsul$y <- 0.0*num_dens_amsul$y
	}
	num_dens_amsul$y <- (nonIN_num_dens+IN_num_dens)*num_dens_amsul$y
	num_dens_amsul$x <- exp(num_dens_amsul$x)
	#jpeg(paste("num_dens_amsul.",time,".jpg",sep=""))
	pdf(paste("num_dens_amsul.",sprintf("%05d",time),".pdf",sep=""))
        plot(num_dens_amsul,
		   log="x",
		   xlim=c(1e-8,1e-6),
		   main=paste("Number density distribution of ammonium bisulfate (T=",
		              sprintf("%05d",time)," [s])\n (sum of nonIN and IN aerosols))"),
                   xlab="Dry radius of ammonium bisulfate [m])",
                   ylab="Number density distribution dN/dlogr ([/unit logr/m^3]")
        dev.off()

	# Mass density distribution of droplets"
	drop_data <- valid_data[valid_data$sd_liqice==1,]
	drop_radi <- drop_data$sd_r
	drop_num  <- sum(drop_data$sd_n)
	drop_num_dens <- drop_num/IMAX/JMAX/KMAX_sdm/length(allmpiranks)/DX/DY/DZ

	rho_liqw <- 1.0E+6 # density of liquid water [g/m3]
	if(length(drop_radi)>1){
		mass_dens_drop <- density(log(drop_radi),bw="nrd",weight=drop_data$sd_n/drop_num)
	}else{
		mass_dens_drop <- density(log(c(1.0e-10,1.0e2)),bw="nrd")
		mass_dens_drop$y <- 0.0*mass_dens_drop$y
	}
	mass_dens_drop$x <- exp(mass_dens_drop$x)
	mass_dens_drop$y <- drop_num_dens*rho_liqw*(4.0/3.0)*pi*mass_dens_drop$x**3*mass_dens_drop$y
	pdf(paste("mass_dens_drop.",sprintf("%05d",time),".pdf",sep=""))
        plot(mass_dens_drop,
		   log="x",
		   xlim=c(0.5e-6,1e-2),
		   main=paste("Mass density distribution of droplets (T=",
		              sprintf("%05d",time)," [s])"),
                   xlab="Radius of droplet [m])",
                   ylab="Mass density distribution of m(r)dN/dlogr ([g/unit logr/m^3]")
        dev.off()

	# Number density distribution of droplets"
	if(length(drop_radi)>1){
		#num_dens_drop <- density(2.0*1.0e6*drop_radi,bw="nrd",weight=drop_data$sd_n/drop_num)
		num_dens_drop <- density(1.0e6*drop_radi,bw="nrd",weight=drop_data$sd_n/drop_num)
		num_dens_drop$y <- drop_num_dens*num_dens_drop$y * 1e-6 # [/m^3]->[/cm^3]
	}else{
		num_dens_drop <- density(c(1.0e-10,1.0e2),bw="nrd")
		num_dens_drop$y <- 1.0+0.0*num_dens_drop$y
	}
	pdf(paste("num_dens_drop.",sprintf("%05d",time),".pdf",sep=""))
       plot(num_dens_drop,
             log="y",
             main=paste("Number density distribution of droplets (T=",sprintf("%05d",time)," [s])"),
             #xlab="diameter [um]",
             xlab="radius [um]",
             ylab="Concentration density [cm^-3/um]")
        dev.off()

	# Droplet size distribution
    rstep <- 0.2
	pdf(paste("hist_dsd.",sprintf("%05d",time),".pdf",sep=""))
	if(length(drop_radi)>1){
		hist_dsd <- weighted.hist(1.0e6*drop_radi, drop_data$sd_n/drop_num/rstep,
                  breaks=seq(0,50,by=rstep),
                  main=paste("Droplet size distribution (T=", sprintf("%05d",time)," [s])"),
		  ylab="PDF [/um]",
                  xlab="radius [um]")
	}else{
		hist_dsd <- weighted.hist(seq(0,50,by=rstep), rep(1,length(seq(0,50,by=rstep))),
                  breaks=seq(0,50,by=rstep),
                  main=paste("Droplet size distribution (T=", sprintf("%05d",time)," [s])"),
		  ylab="PDF [/um]",
                  xlab="radius [um]")
	}
        dev.off()

    fill <- -9.9999e+30
	rdim <- ncdim_def("r", "um", head(hist_dsd$breaks+rstep/2, n=length(hist_dsd$counts)),
                      longname="radius", unlim=TRUE)
	dsd  <- ncvar_def("DSD", "/um", list(rdim), fill, longname="Droplet size distribution")

	outfile <- nc_create(filename=paste("./LES_Nxxx_DSD_", sprintf("%05d",time), sep=""),
                         force_v4=TRUE, list(dsd))

	ncvar_put(outfile, dsd,  hist_dsd$counts)
	ncatt_put(outfile, 0, "author", "Shima-Lab.")
	nc_close(outfile)

	# Droplet size distribution
    rstep <- 0.5
	pdf(paste("hist2_dsd.",sprintf("%05d",time),".pdf",sep=""))
	if(length(drop_radi)>1){
		hist_dsd <- weighted.hist(1.0e6*drop_radi, drop_data$sd_n/IMAX/JMAX/KMAX_sdm/length(allmpiranks)/DX/DY/DZ,
                  breaks=seq(0,50,by=rstep),
                  main=paste("Droplet size distribution (T=", sprintf("%05d",time)," [s])"),
		  ylab="DSD [/m^3]",
                  xlab="radius [um]")
	}else{
		hist_dsd <- weighted.hist(seq(0,50,by=rstep), rep(1,length(seq(0,50,by=rstep))),
                  breaks=seq(0,50,by=rstep),
                  main=paste("Droplet size distribution (T=", sprintf("%05d",time)," [s])"),
		  ylab="DSD [/m^3]",
                  xlab="radius [um]")
	}
        dev.off()

    fill <- -9.9999e+30
	rdim <- ncdim_def("r", "um", head(hist_dsd$breaks+rstep/2, n=length(hist_dsd$counts)),
                      longname="radius", unlim=TRUE)
	dsd  <- ncvar_def("DSD", "/m^3", list(rdim), fill, longname="Droplet size distribution")

	outfile <- nc_create(filename=paste("./LES2_Nxxx_DSD_", sprintf("%05d",time), sep=""),
                         force_v4=TRUE, list(dsd))

	ncvar_put(outfile, dsd,  hist_dsd$counts)
	ncatt_put(outfile, 0, "author", "Shima-Lab.")
	nc_close(outfile)

	# ice size distribution
	ice_data <- valid_data[valid_data$sd_liqice==10,]
	ice_radi <- pmax(ice_data$sdi_re,ice_data$sdi_rp)
	ice_num  <- sum(ice_data$sd_n)
	ice_num_dens <- ice_num/IMAX/JMAX/KMAX_sdm/length(allmpiranks)/DX/DY/DZ

    rstep <- 0.5
	pdf(paste("hist_isd.",sprintf("%05d",time),".pdf",sep=""))
	if(length(ice_radi)>1){
		hist_isd <- weighted.hist(1.0e6*ice_radi, ice_data$sd_n/IMAX/JMAX/KMAX_sdm/length(allmpiranks)/DX/DY/DZ,
                  breaks=seq(0,50,by=rstep),
                  main=paste("Ice size distribution (T=", sprintf("%05d",time)," [s])"),
		  ylab="ISD [/m^3]",
                  xlab="radius [um]")
	}else{
		hist_isd <- weighted.hist(seq(0,50,by=rstep), rep(1,length(seq(0,50,by=rstep))),
                  breaks=seq(0,50,by=rstep),
                  main=paste("Ice size distribution (T=", sprintf("%05d",time)," [s])"),
		  ylab="ISD [/m^3]",
                  xlab="radius [um]")
	}
        dev.off()

    fill <- -9.9999e+30
	rdim <- ncdim_def("r", "um", head(hist_isd$breaks+rstep/2, n=length(hist_isd$counts)),
                      longname="radius", unlim=TRUE)
	isd  <- ncvar_def("ISD", "/m^3", list(rdim), fill, longname="Ice size distribution")

	outfile <- nc_create(filename=paste("./LES2_Nxxx_ISD_", sprintf("%05d",time), sep=""),
                         force_v4=TRUE, list(isd))

	ncvar_put(outfile, isd,  hist_isd$counts)
	ncatt_put(outfile, 0, "author", "Shima-Lab.")
	nc_close(outfile)

	# terminal velocity of droplets
#	drop_data <- valid_data[valid_data$sd_liqice==1,]
	drop_radi <- drop_data$sd_r
	drop_termv<- drop_data$sd_vz # terminal velocity of droplets [m/s]
	n_data_skip <- max(round(length(drop_termv)/10000),1)

	if(length(drop_termv)==0){ #dummy data
		drop_radi  <- 1.0e2
		drop_termv <- 0.0e0
	}
	pdf(paste("term_vel_drop.",sprintf("%05d",time),".pdf",sep=""))
        plot(drop_radi[seq(1, length(drop_radi), n_data_skip)],drop_termv[seq(1, length(drop_termv), n_data_skip)],
		   log="x",pch='.',
		   xlim=c(0.5e-6,1e-2),
		   main=paste("Terminal velocity of droplets (T=",
		              sprintf('%05d',time)," [s])"),
                   xlab="Radius of droplet [m]",
                   ylab="Terminal velocity of droplets [m/s]")
        dev.off()

	# x-z position of SDs
	valid_x <- valid_data$sd_x # x coordinate of valid SD [m]
	valid_z <- valid_data$sd_z # z coordinate of valid SD [m]
	n_data_skip <- max(round(length(valid_z)/10000),1)

	if(length(valid_x)==0){ #dummy data
		valid_x  <- -1.0e2
		valid_z  <- -1.0e2
	}
	pdf(paste("xz_SD.",sprintf("%05d",time),".pdf",sep=""))
        plot(valid_x[seq(1, length(valid_x), n_data_skip)],valid_z[seq(1, length(valid_z), n_data_skip)],
		   pch='.',
		   xlim=c(0.0,PRC_NUM_X*DX*IMAX),
		   ylim=c(0.0,DZ*KMAX),
		   main=paste("Positions of SDs on x-z plane (T=",
		              sprintf('%05d',time)," [s])"),
                   xlab="x [m]",
                   ylab="z [m]")
        dev.off()

        # x-y position of SDs
        valid_x <- valid_data$sd_x # x coordinate of valid SD [m]
        valid_y <- valid_data$sd_y # y coordinate of valid SD [m]
        n_data_skip <- max(round(length(valid_y)/10000),1)

        if(length(valid_x)==0){ #dummy data
                valid_x  <- -1.0e2
                valid_y  <- -1.0e2
        }
        pdf(paste("xy_SD.",sprintf("%05d",time),".pdf",sep=""))
        plot(valid_x[seq(1, length(valid_x), n_data_skip)],valid_y[seq(1, length(valid_y), n_data_skip)],
                   pch='.',
                   xlim=c(0.0,PRC_NUM_X*DX*IMAX),
                   ylim=c(0.0,PRC_NUM_Y*DY*JMAX),
                   main=paste("Positions of SDs on x-y plane (T=",
                              sprintf('%05d',time)," [s])"),
                   xlab="x [m]",
                   ylab="y [m]")
        dev.off()

}

stopCluster(cluster)

# Merge LES_Nxxx_dsd(by plottime) to single file
DSD_ALL <- NULL
isFirst <- TRUE
for(time in plottimes){
    ncfile <- paste("./LES_Nxxx_DSD_", sprintf("%05d",time), sep="")
    ncin <- nc_open(ncfile)
    if(isFirst){
        r <- ncvar_get(ncin,"r")
        ifFirst <- FALSE
    }
    DSD  <- ncvar_get(ncin,"DSD")
    DSD_ALL <- append(DSD_ALL, DSD)
    nc_close(ncin)
    file.remove(ncfile)
}
DSD_ALL <- array(DSD_ALL,dim=c(length(DSD), length(plottimes)))

rdim <- ncdim_def("r", "um", r)
tdim <- ncdim_def("t", "sec.", plottimes, unlim=TRUE)
fillvalue <- -9.9999e+30
DSDvar    <- ncvar_def("DSD", "cm^-3/um", list(rdim, tdim), fillvalue,
            longname="Droplet size distribution")
outfile <- nc_create(filename="./LES_Nxxx_DSD", force_v4=TRUE, list(DSDvar))
ncvar_put(outfile, DSDvar, DSD_ALL)
ncatt_put(outfile, 0, "author", "Shima-Lab.")
nc_close(outfile)
