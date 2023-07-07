library(ggplot2)
library(dplyr)
library(latex2exp)

datapoints = read.csv("~/Data/StrongSelection.dat", sep="", as.is = T )
datapoints$TIME = as.numeric(datapoints$TIME)
datapoints$M = as.numeric(datapoints$M)
datapoints$H = as.numeric(datapoints$H)
datapoints$S = as.numeric(datapoints$S)


dp = datapoints %>%
      filter(RES == "Fix") %>%
      group_by(EM, Ndeme, Nind, S, H, M, MU) %>%
      summarise(MEAN.TIME = mean(TIME),
                SD.TIME = sd(TIME), 
                N = n(),
                SE.TIME = SD.TIME/sqrt(N) )

dp %>%
  ggplot() + 
  aes(x = M, y = MEAN.TIME, col = as.factor(H)) + 
  geom_point(size = 1, alpha = 0.8) +
  geom_line(size = 1, alpha = 0.7) + 
  geom_errorbar(aes(ymin = MEAN.TIME-SD.TIME, ymax = MEAN.TIME+SD.TIME), width = 0.1, alpha = 0.6) + 
  facet_grid(cols = vars(S), scales = "free_y") + 
  labs(y = TeX("$T_\\mathrm{fix}$"),
       x = TeX("Dispersal rate, \\textit{d}"),
       col = TeX("Genetic dominance, \\textit{h}:")) + 
  my.theme + x.dispersal + y.meantime + dominance.color
