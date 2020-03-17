rm(list = ls(all = TRUE))
###################################################

# libraries
library(coda)
library(rmatio)
library(ggplot2)
library(dplyr)

# number of draws
Ndraws <- 1000
Nthin <- 1 # thin argument in mcmc does NOT have an affect on the results from effectivesize!

# spec details 
Nrs <- c(1:10)
Npriors <- c(5, 1, 2, 3, 4)
Nsurveys <- c("level")
Nsamples <- c("rec")
Ncountries <- c("GER", "US")
Nps <- 1 

# prior names
priors <- c("NIG", "MG", "PMNM", "HS+", "Nd")

# initialize dataframe
effsampsize <- data.frame(surveys = vector(), 
                          sample = vector(),
                          Np = vector(),
                          prior = vector(),
                          country = vector(),
                          vintage = vector(),
                          value = vector()
                          )

# loop over files
for (country in Ncountries){
  if (country == "GER")
  {
    Nvintages <- 157
    dirname <- "C:/Users/Philipp/Documents/Dissertation/sparse nowcasting/estim/PH_GER/"
  } 
  else if (country == "US")
  {
    Nvintages <- 229 
    dirname <- "C:/Users/Philipp/Documents/Dissertation/sparse nowcasting/estim/PH_US/"
  }
  for (survey in Nsurveys){
    for (sample in Nsamples){
      for (Np in Nps){
        for (r in 1 : length(Nrs)) {
          Nr <- Nrs[ r ] 
          for (prior in Npriors){
            for (Nvintage in 1 : Nvintages){
              # paste together filename
              filename <- paste0("PH_", country, 
                                 "_v", Nvintage,
                                 "_prior", prior ,
                                 "_Nr", Nr,
                                 "_Np", Np,
                                 '_', sample, "_", survey,
                                 ".mat")
            
              # read in mat file
              x <- read.mat(paste0(dirname, filename))
              
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
              ineff <- Ndraws / effectiveSize( as.mcmc( draws_mat , thin = Nthin ) ) 
              
              # append to df
              effsampsize <- rbind(effsampsize, 
                                   data.frame(surveys = rep(survey, ncol(draws_mat)),
                                              sample = rep(sample, ncol(draws_mat)),
                                              Np = rep(Np, ncol(draws_mat)),
                                              prior = rep(priors[prior], ncol(draws_mat)),
                                              country = rep(country, ncol(draws_mat)),
                                              vintage = rep(Nvintage, ncol(draws_mat)),
                                              value = ineff
                                              )
                                  )
            }
          }
        }
      }
    }    
  }      
}    

# plots
library(forcats)
ggplot(effsampsize) + 
  geom_boxplot(aes(x=fct_rev(prior),y=value))+
  facet_wrap( ~ country, scales = "free_y")+
  labs(x = "", y = "inefficiency factor")+
  coord_flip()