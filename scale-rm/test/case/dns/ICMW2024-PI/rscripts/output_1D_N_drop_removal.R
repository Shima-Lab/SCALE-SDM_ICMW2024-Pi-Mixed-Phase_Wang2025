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
RESULT_FILE <- "./LES_Nxxx_1D_dremoval" # result file name
############################################################################
########## Read Data from txts
library(ncdf4)


############################################################################
########## Read Data from txts
cat("Start plotting N_drop_removal", "\n")
directory <- "../"
file_list <- list.files(path = directory, pattern = "^virtual_disdro.pe000[0-9][0-9][0-9]$", full.names = TRUE)

cat("file", "\n")

# Setting of bins
nbins <- 180
min_time <- 0.0    # [s]
max_time <- 5400   # [s]
PI <- acos(-1)     #[]
rho_d <- 1.0e+03 #density of water [kgm^-3]

time.bin <- seq(min_time, max_time, length = nbins + 1)
d_time.bin <- time.bin[2] - time.bin[1]

#initialize histogram
hist_drop <- array(rep(0,nbins))
hist_mass <- array(rep(0,nbins))
#hist_drop <- numeric(nbins)

for (file in file_list) {
  file_data <- read.table(file, header = TRUE, sep = "", stringsAsFactors = FALSE, fill = TRUE)
  colnames(file_data) <- c("time", "sd_x", "sd_y", "sd_r", "sd_n", "n")
  file_data <- file_data[file_data$sd_r > 3.5e-6, ]

  time <- file_data$time
  sd_x <- file_data$sd_x
  sd_y <- file_data$sd_y
  sd_r <- file_data$sd_r
  sd_n <- file_data$sd_n
  n    <- file_data$n

  mass <- rho_d * (4.0/3.0 * PI * (sd_r ** 3.0))

  time.bin_index <- findInterval(time, time.bin, all.inside = TRUE)

  # Add the sd_n values to the corresponding bins
  for (i in seq_along(time.bin_index)) {
    hist_drop[time.bin_index[i]] <- hist_drop[time.bin_index[i]] + sd_n[i]
    hist_mass[time.bin_index[i]] <- hist_mass[time.bin_index[i]] + mass[i] * sd_n[i]
  }
  cat(file, "\n")
}

# Define dimensions
time_dim <- ncdim_def("time", "seconds", time.bin[-1],  unlim=TRUE)

fillvalue <- -9.9999e+30
N_drop_removalvar <- ncvar_def("N_drop_removal", "m^-3S^-1", list(time_dim), fillvalue,
                      longname="droplet removal rate")
M_drop_removalvar <- ncvar_def("M_drop_removal", "kg m^-3S^-1", list(time_dim), fillvalue,
                      longname="liquid mass removal rate")
ncdflist <- list(N_drop_removalvar, M_drop_removalvar)
outfile <- nc_create(filename=RESULT_FILE, force_v4=TRUE, ncdflist)

# Write data to variables
ncvar_put(outfile, N_drop_removalvar, hist_drop * 0.25 / d_time.bin, 1)
ncvar_put(outfile, M_drop_removalvar, hist_mass * 0.25 / d_time.bin, 1)

# Close the NetCDF file
nc_close(outfile)

cat("NetCDF file created: N_drop_removal.nc\n")
cat("done.\n")
