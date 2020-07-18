#******************************************************************#
#* this code analyzes the April 2017 vintage of the PMNM prior
#* with Nr = 10, Np = 3, level, rolling
#* where some of the draws were extreme (abs(x) >> 1000)
#* A copy of the mat file is stored in 
#* C:/Users/Philipp/Documents/Dissertation/sparse nowcasting/estim
#* and I am currently re-estimating that particular spec
#* both locally and on the HPC server (26.5.2020)
#******************************************************************#

rm(list = ls())

library(rmatio)

dirname <- "C:/Users/Philipp/Documents/Dissertation/sparse nowcasting/estim/extreme draws PMNM"
filename <- "PH_US_v209_prior3_Nr10_Np3_rolling_level_extremedraws"

tmp <- read.mat(paste0(dirname, "/", filename, ".mat"))

draws <- unlist(tmp$draws$nowcast)
flag_phi <- unlist(tmp$draws$flag_phi_prev)

sum(flag_phi)
draws[flag_phi == 1]
plot(1:length(draws), draws)

draws_NA <- draws
ind_outlier <- abs(draws_NA) > 1000 
draws_NA[ind_outlier] <- NA
draws[ind_outlier]
sum(is.na(draws_NA))

plot(1:length(draws_NA), draws_NA)
