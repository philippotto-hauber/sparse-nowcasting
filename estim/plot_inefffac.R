rm(list=ls())
library(forcats)
library(ggplot2)
library(dplyr)
dirname <- "C:/Users/Philipp/Documents/Dissertation/sparse nowcasting/estim/"
filename_in <- "ineff_facs.Rda"
filename_out <- "plot_ineff_facs.pdf"

# load data
load(paste0(dirname, filename_in))

# filter Np == 1
ineff_facs %>% filter(Np == "P = 1") -> ineff_facs

# plots
ggplot(ineff_facs) + 
  geom_boxplot(aes(x=prior,y=value, color = country), outlier.size = 2, outlier.shape = 16)+
  #labs(x = "", y = "", caption = "rec/rolling: recursive or rolling estimation window, level/diff: survey variables in levels or first differences\nPriors: HS+ = Horseshoe plus, MG = multiplicative Gamma, NIG = Normal Inverse Gamma, \nNd = Normal diffuse, PMNM = point mass normal mixture. For details, see main text.")+
  labs(x = "", y = "", caption = "") +
  facet_grid(sample ~ survey)+
  scale_y_continuous(breaks = c(50, 100), minor_breaks = c(50, 100))+
  theme_minimal()+
  theme(strip.text.x = element_text(size = 10),
        strip.text.y = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = "top",
        panel.grid.major.x = element_blank(),
        plot.caption = element_text(size = 8))

# save file as full size pdf
#ggsave(file="a4_output.pdf", width = 0.5*210, height = 0.5*297, units = "mm")
ggsave(file = paste0(dirname, filename_out), width = 150, height = 150, units = "mm")
