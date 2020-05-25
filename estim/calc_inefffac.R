rm(list = ls())
###################################################

# libraries
library(coda)
library(rmatio)
library(ggplot2)
library(dplyr)
library(doParallel)

# parallel
registerDoParallel(cores=8)

# number of draws
Ndraws <- 1000
Nthin <- 1 # thin argument in mcmc does NOT have an affect on the results from effectivesize!

# spec details 
Nrs <- seq(1, 10)
Npriors <- c(5, 1, 2, 3, 4)
Nsurveys <- c("level", "diff")
Nsamples <- c("rec", "rolling")
Nmod <- seq(1, length(Nsurveys) + length(Nsamples))
Ncountries <- c("GER", "US")
Nps <- c(1, 3) 

# prior names
priors <- c("NIG", "MG", "PMNM", "HS+", "Nd")

# initialize dataframe
ineff_facs <- data.frame(surveysample = vector(), 
                          Np = vector(),
                          prior = vector(),
                          Nr = vector(),
                          country = vector(),
                          vintage = vector(),
                          value = vector(),
                          count_flagged_draws = vector(),
                          filename = vector()
                          )

# loop over files
start_time <- Sys.time()
for (country in Ncountries)
{
  
  # number of vintages and path to matfiles depends on country
  if (country == "GER")
  {
    Nvintages <- seq(1, 157)
    #dirname <- "C:/Users/Philipp/Documents/Dissertation/sparse nowcasting/estim/PH_GER/"
    dirname <- "../../PH_GER/"
  } 
  else if (country == "US")
  {
    Nvintages <- seq(1, 229) 
    #dirname <- "C:/Users/Philipp/Documents/Dissertation/sparse nowcasting/estim/PH_US/"
    dirname <- "../../PH_US/"
  }
  
  # merge spec details into one big loop
  Nrs_loop = rep(Nrs, length(Npriors) * length(Nvintages) * length(Nps) * length(Nmod))
  Npriors_loop = rep(Npriors %x% rep(1, length(Nrs)), length(Nvintages) * length(Nps) * length(Nmod))
  Nvintages_loop = rep(Nvintages %x% rep(1, length(Nrs) * length(Npriors)), length(Nps) * length(Nmod))
  Nps_loop = rep(Nps %x% rep(1, length(Nrs) * length(Npriors) * length(Nvintages)),length(Nmod))
  Nmod_loop = Nmod %x% rep(1, length(Nrs) * length(Npriors) * length(Nvintages) * length(Nps))
  
  tmp_df <- foreach (ind = seq(1, length(Nrs_loop)), .combine = rbind)%dopar%
  {
      # back out survey and sample from Nmod_loop!
      if (Nmod_loop[ind] == 1)
      {
        sample <- Nsamples[1]
        survey <- Nsurveys[1]
      }
      else if (Nmod_loop[ind] == 2)
      {
        sample <- Nsamples[1]
        survey <- Nsurveys[2]
      }
      else if (Nmod_loop[ind] == 3)
      {
        sample <- Nsamples[2]
        survey <- Nsurveys[1]
      }
      else
      {
        sample <- Nsamples[2]
        survey <- Nsurveys[2]
      }

      # paste together filename
      filename <- paste0("PH_", country, 
                         "_v", Nvintages_loop[ind],
                         "_prior", Npriors_loop[ind] ,
                         "_Nr", Nrs_loop[ind],
                         "_Np", Nps_loop[ind],
                         '_', sample, "_", survey,
                         ".mat")
    
      # read in mat file
      x <- rmatio::read.mat(paste0(dirname, filename))
      
      # extract draws of nowcast
      temp <- cbind(x$draws$nowcast[[1]])
      
      # check if forecasts exist
      if ("forecast_mean" %in% names(x$draws))
        temp <- cbind( temp, x$draws$forecast[[1]])

      # multiply US draws with 100
      if (country == "US")
        temp <- temp * 100 
      
      # convert to matrix & remove temp
      draws_mat <- as.matrix(temp)
      rm(temp)
      
      # compute inefficiency factors
      ineff <- Ndraws / coda::effectiveSize(coda::as.mcmc(draws_mat, thin = Nthin)) 
      
      # append to df
      df_out <-   data.frame(surveysample = rep(paste0(survey, ", ", sample), ncol(draws_mat)),
                              Np = rep(paste("P =", Nps_loop[ind]), ncol(draws_mat)),
                              Nr = rep(Nrs_loop[ind], ncol(draws_mat)),
                              prior = rep(priors[Npriors_loop[ind]], ncol(draws_mat)),
                              country = rep(country, ncol(draws_mat)),
                              vintage = rep(Nvintages_loop[ind], ncol(draws_mat)),
                              value = ineff,
                              count_flagged_draws = sum(unlist(x$draws$flag_phi_prevx)),
                              filename = filename,
                              stringsAsFactors = FALSE
                              )
  }
  
  # merge df into master df
  ineff_facs <- rbind(ineff_facs, tmp_df)
}
end_time <- Sys.time()
print(end_time - start_time)

# convert Np, prior, country and surveysample to factors
ineff_facs$surveysample <- as.factor(ineff_facs$surveysample)
ineff_facs$Np <- as.factor(ineff_facs$Np)
ineff_facs$country <- as.factor(ineff_facs$country)
ineff_facs$prior <- as.factor(ineff_facs$prior)

# save to file
#save(ineff_facs, file = "C:/Users/Philipp/Documents/Dissertation/sparse nowcasting/estim/ineff_facs.Rda")
save(ineff_facs, file = "ineff_facs.Rda")
