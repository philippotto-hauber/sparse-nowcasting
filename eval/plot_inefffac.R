rm(list=ls())
load("C:/Users/Philipp/Documents/Dissertation/sparse nowcasting/estim/ineff_facs.Rda")

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
  facet_grid(surveysample ~ Np, scales = "free_y")+
  theme_minimal()+
  theme(strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8),
        legend.title = element_blank())

# save file as full size pdf
#ggsave(file="a4_output.pdf", width = 0.5*210, height = 0.5*297, units = "mm")
ggsave(file="plot_ineff_facs.pdf", width = 150, height = 150, units = "mm")