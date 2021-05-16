setwd("C:/Users/Philipp/Documents/GitHub/sparse-nowcasting/eval/test crps log score")

dat <- read.csv('data_test_crps_logS.csv', header = FALSE)


truegdp <- dat[1, 1]
dens <- dat[1, 2:length(dat)]

library(scoringRules)

crps_sample(y = truegdp, dat = as.matrix(dens))

logs_sample(y = truegdp, dat = as.matrix(dens))
