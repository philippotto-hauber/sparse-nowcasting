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
Nrs <- seq(1, 3)
Npriors <- c(5, 1, 2, 3, 4)
Nsurveys <- c("level", "diff")
Nsamples <- c("rec", "rolling")
Ncountries <- c("GER", "US")
Nps <- c(1, 3) 

# prior names
priors <- c("NIG", "MG", "PMNM", "HS+", "Nd")

# initialize dataframe
ineff_facs <- data.frame(surveysample = vector(), 
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
              ineff_facs <- rbind(ineff_facs, 
                                  data.frame(surveysample = rep(paste0(survey, ", ", sample), ncol(draws_mat)),
                                              Np = rep(paste("P =", Np), ncol(draws_mat)),
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
ggplot(ineff_facs) + 
  geom_boxplot(aes(x=fct_rev(prior),y=value))+
  facet_wrap( ~ country, scales = "free_y")+
  labs(x = "", y = "inefficiency factor")+
  coord_flip()

# alternative plots
ggplot(ineff_facs) + 
  geom_boxplot(aes(x=prior,y=value, color = country), outlier.size=1)+
  labs(x = "", y = "", caption = "Blablabla")+
  facet_grid(Np ~ surveysample, scales = "free_y")+
  theme_minimal()+
  theme(strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        legend.title = element_blank())

# save file as full size pdf
#ggsave(file="a4_output.pdf", width = 0.5*210, height = 0.5*297, units = "mm")
ggsave(file="plot_ineff_facs.pdf", width = 150, height = 150, units = "mm")
