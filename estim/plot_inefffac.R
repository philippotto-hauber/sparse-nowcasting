rm(list=ls())
library(forcats)
library(ggplot2)
dirname <- "C:/Users/Philipp/Documents/Dissertation/sparse nowcasting/estim/"
filename_in <- "ineff_facs.Rda"
filename_out <- "plot_ineff_facs.pdf"

# load data
load(paste0(dirname, filename_in))

# plots
ggplot(ineff_facs) + 
  geom_boxplot(aes(x=fct_rev(prior),y=value))+
  facet_wrap( ~ country, scales = "free_y")+
  labs(x = "", y = "inefficiency factor")+
  coord_flip()

# alternative plots
ggplot(ineff_facs) + 
  geom_boxplot(aes(x=prior,y=value, color = country), outlier.size = 1, outlier.shape = 16)+
  labs(x = "", y = "", caption = "Blablabla")+
  facet_grid(surveysample ~ Np)+
  scale_y_continuous(breaks = c(50, 100), minor_breaks = c(50, 100))+
  theme_minimal()+
  theme(strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        legend.title = element_blank(),
        panel.grid.major.x = element_blank())

# save file as full size pdf
#ggsave(file="a4_output.pdf", width = 0.5*210, height = 0.5*297, units = "mm")
ggsave(file = paste0(dirname, filename_out), width = 150, height = 150, units = "mm")
