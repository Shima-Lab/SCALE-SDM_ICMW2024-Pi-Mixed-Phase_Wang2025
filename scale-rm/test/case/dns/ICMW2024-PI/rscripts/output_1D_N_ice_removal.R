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
RESULT_FILE <- "./LES_Nxxx_1D_iremoval" # result file name
############################################################################
########## Read Data from txts
library(ncdf4)


############################################################################
########## Read Data from txts
cat("Start plotting N_ice_removal", "\n")
directory <- "../"
file_list <- list.files(path = directory, pattern = "^virtual_snow_disdro.pe000[0-9][0-9][0-9]$", full.names = TRUE)
#file_list <- list.files(path = directory, pattern = "^virtual_snow_disdro.pe000000$", full.names = TRUE)

cat("file", "\n")

# Setting of bins
nbins <- 180
min_time <- 0.0    # [s]
max_time <- 5400   # [s]
PI <- acos(-1)     #[]

time.bin <- seq(min_time, max_time, length = nbins + 1)
d_time.bin <- time.bin[2] - time.bin[1]

#initialize histogram
hist_ice <- array(rep(0,nbins))
hist_mass <- array(rep(0,nbins))
#hist_drop <- numeric(nbins)

  for (file in file_list) {
    #file_data <- read.table(file, header = TRUE, sep = ",", stringsAsFactors = FALSE)  # 適宜変更
    file_data <- read.table(file, header = TRUE, sep = "", stringsAsFactors = FALSE, fill = TRUE)
    #time <- read.table(file, header = TRUE, sep = ",","time")
    colnames(file_data) <- c("time", "sd_x", "sd_y", "sd_re", "sd_rp", "sd_rho", "sd_n", "n")

    time <- file_data$time
    sd_x <- file_data$sd_x
    sd_y <- file_data$sd_y
    sd_re <- file_data$sd_re
    sd_rp <- file_data$sd_rp
    sd_rho <- file_data$sd_rho
    sd_n <- file_data$sd_n
    n    <- file_data$n

    mass <- sd_rho * ( 4.0/3.0  * PI * sd_rp * (sd_re ** 2.0) )

    #findInterval(x, vec, rightmost.closed = FALSE, all.inside = FALSE)

    time.bin_index <- findInterval(time, time.bin, all.inside = TRUE)

    # Add the sd_n values to the corresponding bins
    for (i in seq_along(time.bin_index)) {
      hist_ice[time.bin_index[i]] <- hist_ice[time.bin_index[i]] + sd_n[i]
      hist_mass[time.bin_index[i]] <- hist_mass[time.bin_index[i]] + mass[i] * sd_n[i]
    }
    cat(file, "\n")
  }

  # Define dimensions
  time_dim <- ncdim_def("time", "seconds", time.bin[-1],  unlim=TRUE)

  fillvalue <- -9.9999e+30
  N_ice_removalvar <- ncvar_def("N_ice_removal", "m^-3S^-1", list(time_dim), fillvalue,
                        longname="ice removal rate")
  M_ice_removalvar <- ncvar_def("M_ice_removal", "kg m^-3S^-1", list(time_dim), fillvalue,
                        longname="ice mass removal rate")
  ncdflist <- list(N_ice_removalvar, M_ice_removalvar)
  outfile <- nc_create(filename=RESULT_FILE, force_v4=TRUE, ncdflist)

  # Write data to variables
  ncvar_put(outfile, N_ice_removalvar, hist_ice * 0.25 / d_time.bin, 1)
  ncvar_put(outfile, M_ice_removalvar, hist_mass * 0.25 / d_time.bin, 1)

  # Close the NetCDF file
  nc_close(outfile)

  cat("NetCDF file created: N_drop_removal.nc\n")
cat("done.\n")