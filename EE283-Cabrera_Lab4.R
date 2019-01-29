library(ggplot2)
library(datasets)
library(tidyverse)
library(wesanderson)
library(colorblindr)

mtcars$vs[mtcars$vs == 0] <- "V-Shaped Engine"
mtcars$vs[mtcars$vs == 1] <- "Straight Engine"
mtcars$am[mtcars$am == 0] <- "Automatic Transmission"
mtcars$am[mtcars$am == 1] <- "Manual Transmission"
mtcars$gear <- sapply(as.character(mtcars$gear), paste, "gears", sep=" ")
mg <- ggplot(mtcars, aes(x = mpg, y = wt, color = gear, shape = gear)) + geom_point(size=3) + 
  xlab("Miles/(US) gallon") + 
  ylab("Weight (1000 lbs)") +
  scale_color_manual(values=wes_palette(n=3, name="FantasticFox1"))
mg + facet_grid(am ~ vs, labeller = labeller(vs = label_value, am = label_value)) +
  theme_bw(base_size = 18) +
  theme(legend.position = c(0.87,0.85), axis.text.x = element_text(angle=45),
        legend.background = element_rect(colour = NA, fill = NA))
cvd_grid(mg)


ds <- ggplot(mpg, aes(displ, hwy, colour = class)) +
  geom_point(aes(shape=class)) +
  geom_smooth(aes(color=class,fill=class), alpha = 0.2) +
  theme_bw(base_size = 18) +
  ylab("Highway Miles per Gallon") + 
  xlab("Engine Displacement (Liters)") + 
  #scale_color_manual(values = wes_palette(n=7, "FantasticFox1")) +
  #scale_fill_manual(values = wes_palette(n=7, "FantasticFox1")) +
  labs(color = "Car class", fill="Car class", shape = "Car class") 
ds
ggsave(mpg-to-ed, plot = last_plot())
