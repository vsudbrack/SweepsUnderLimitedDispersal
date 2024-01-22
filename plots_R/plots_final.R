################################################
############ Figures from the paper ############
################################################

###### Both simulations and analytical comes from .csv files
###### See the folder simulations or analytical to codes to generate these files

library("RColorBrewer")
library(ggplot2)
library(dplyr)
library(scales)
library(latex2exp)
library(patchwork)
setwd("~/sync/LimitedDispersal/")
source("RozeFrancois.R")

# Define theme of ggplot2 ----

my.theme = theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        legend.position="top",
        text = element_text(size = 16))

y.meantime = scale_y_continuous(trans = log10_trans(), 
                                breaks = trans_breaks("log10", function(x) 10^x), 
                                labels = trans_format("log10", math_format(10^.x)),
                                #limits = c(10**3.7, 10**5),
)

x.dominance = scale_x_continuous(breaks =c(0, 0.25, 0.5, 0.75, 1),
                                 labels = c("0\nRecessive", "0.25", 
                                            "0.5\nAdditive","0.75", 
                                            "1\nDominant")          )

y.dispersal = scale_y_continuous(trans = log10_trans(), 
                                 breaks = c(1, 0.1, 0.01, 0.001, 1E-4), 
                                 labels = c(TeX("$1$"), 
                                            TeX("$10^{-1}$"), 
                                            TeX("$10^{-2}$"), 
                                            TeX("$10^{-3}$"),
                                            TeX("$10^{-4}$")))

x.dispersal = scale_x_continuous(trans = log10_trans(), 
                                 breaks = c(1, 0.1, 0.01, 0.001, 1E-4), 
                                 labels = c(TeX("$1$"), 
                                            TeX("$10^{-1}$"), 
                                            TeX("$10^{-2}$"), 
                                            TeX("$10^{-3}$"),
                                            TeX("$10^{-4}$")))

dominance.color = scale_color_viridis_d(end = 0.8)
dominance.fill =  scale_fill_viridis_d(end = 0.8)

x.FST = scale_x_continuous(trans = identity_trans(),
                           limits = c(0, 1),
                           breaks = 0.01*c(0, 2.5, 20, 25, 50, 70, 75, 96, 100), 
                           labels = c("Well mixed", "\nd = 0.1", 
                                      "\nd = 0.01", "25", "50",
                                      "\nd = 0.001", "75", "\nd = 0.0001", 
                                      "100"))


h.labs <- c("Recessive", "Additive", "Dominant")
names(h.labs) <- c("0", "0.5", "1")

phi.labs <- c("All colonizers from\ndifferent demes", "All colonizers from\nthe same deme")
names(phi.labs) <- c("0", "1")

case.labs <- c("Soft sweep\nwithout dominance shift", "Soft sweep\nwith dominance shift")
names(case.labs) <- c("noshift", "shift")

mu.labs <- c(TeX("10^{-8}"), TeX("10^{-6}"), TeX("10^{-4}"))
names(mu.labs) <- c("1e-8", "1e-06", "1e-04")

sel.labs <- c("Soft selection", "Hard selection")
names(sel.labs) <- c("soft", "hard")


ratio_log = function(a, b){
  # Relates a ratio to a value in logscale
  10 ^ ( (a / (a + b) ) * log10( a + b) )
}

source(../calculations_R/Model.R) 

####################
#### Fig1 #########
###################

#### A.
data.raw = read.csv("xxxx", 
                     sep="")
data = data.raw %>%
  group_by(Ndeme, Nind, S, H, M, MU) %>% 
  summarise(TF.mean = mean(TF), 
            TF.var  = var(TF), 
            TF.sd  = sd(TF) )

data.theory = read.csv("xxxx", 
                    sep=",")

data %>%
  ggplot(aes(x = M*(0.95+0.1*H), y = TF.mean, color=factor(H) )) + 
  geom_point(alpha=0.9, size=2.0 ) + 
  geom_errorbar(aes(ymin = TF.mean - TF.sd, ymax = TF.mean + TF.sd), 
                width = 0.075, alpha = 0.7, size=1.0) + 
  geom_line(data = data.theory, aes(x = m, y = fixtime, color=factor(h)), 
            alpha=0.9, size=1.3) + 
  xlab(TeX("Dispersal rate, \\textit{d}")) + ylab("Generations") +
  labs(col=TeX("Dominance, \\textit{h}:"), title = TeX("$Tfix$")) +  
  my.theme + y.meantime + x.dispersal + dominance.color + 
  theme(legend.position = "none")

#### B.
M = 2^seq(0, -14, -0.25) #c(1, 1E-1, 1E-2, 1E-3) # #c(1, 1E-1, 1E-2, 1E-3) #
datapoints = expand.grid(m = M)  
aux_Neff = function(ma){
  a = Roze(s=0.01, h=0.5, m=ma, mu=0, n=200, N=100)
  a$Neff
}
datapoints$Neff = mapply(aux_Neff,ma = datapoints$m)
datapoints %>%
  ggplot(aes(x = m, y = Neff/(200 * 100) )) + 
  geom_line(alpha = 0.7, size = 1.45) + 
  xlab(TeX("Dispersal rate, \\textit{d}")) + 
  ylab(TeX("$\\textit{N}_{\\mathrm{e}}\ / \ nN$")) + 
  labs(title = "Effective population size") +
  my.theme + x.dispersal + 
  scale_y_continuous(trans = log10_trans(), 
                     breaks = c(1, 2, 5, 10, 50), 
                     labels = c("1", "2", "5", "10", "50") )

#### C. 
data.raw = read.csv("xxxx", 
                    sep="")

aux = function(disp) {
  data.raw.fixed = data.raw %>%
    filter( M == disp ) %>%
    filter( ((H == 0.0)&(ID %in% unique(data.raw %>% filter(H == 0.0, M==disp, FREQ>0.95) %>% select(ID))$ID)) |
            ((H == 0.5)&(ID %in% unique(data.raw %>% filter(H == 0.5, M==disp, FREQ>0.95) %>% select(ID))$ID)) |
            ((H == 1.0)&(ID %in% unique(data.raw %>% filter(H == 1.0, M==disp, FREQ>0.95) %>% select(ID))$ID))  ) 
  
  sample.recessive = sample(unique(data.raw.fixed%>%filter(H == 0.0)%>%select(ID))$ID, 5)
  sample.additive  = sample(unique(data.raw.fixed%>%filter(H == 0.5)%>%select(ID))$ID, 5)
  sample.dominant  = sample(unique(data.raw.fixed%>%filter(H == 1.0)%>%select(ID))$ID, 5)
  
  data.raw.fixed %>%
    filter( ((H == 0.0)&(ID %in% sample.recessive)) |
            ((H == 0.5)&(ID %in% sample.additive )) |
            ((H == 1.0)&(ID %in% sample.dominant ))  )  %>% 
    ggplot(aes(x = TIME, y = FREQ, color=factor(H), group = factor(ID))) + 
    geom_line(alpha=0.7, size=0.8) + 
    xlab("Generations") + ylab("Allele frequency") +
    labs(col=TeX("Dominance, \\textit{h}:")) +  
    my.theme + dominance.color + 
    scale_x_continuous(trans = log10_trans(),
                       limits = c(500, 10000))
}; aux(1) / aux(0.001)


####################
#### Fig2 #########
###################

#### A.

data.raw1 = read.csv("xxxx", 
                     sep="")
data.raw2 = read.csv("xxx", 
                     sep="")
data.raw3 =read.csv("xxx", 
                    sep="")
data.raw4 = read.csv("xxx", 
                     sep="")
data.raw5 = read.csv("xxx", 
                     sep="")

data.raw1$ID = 1
data.raw2$ID = 1
data.raw3$E = 0
data.raw4$E = 0
data.raw5$E = 0

data.raw = rbind(data.raw1, data.raw2, data.raw3, 
                 data.raw4, data.raw5)

data = data.raw %>%
  group_by(Ndeme, Nind, S, H, M, MU, E) %>% 
  summarise(TF.mean = mean(TF), 
            TF.var  = var(TF), 
            TF.sd  = sd(TF) )

datapointsHS = read.csv("xxxx")
data_th = datapointsHS %>% filter(h %in% c(0, 0.5, 1))
data_th$MU = data_th$mu

data %>%
  filter( MU == 1e-6) %>%
  ggplot(aes(x = M, y = TF.mean, color=factor(H) )) + 
  geom_point(alpha=0.9, size=1.8 ) + 
  #geom_point(aes(y = TF.mean - TF.sd), size = 0.45, alpha = 0.75) + 
  #geom_point(aes(y = TF.mean + TF.sd), size = 0.45, alpha = 0.75) + 
  geom_line(data = data_th %>% filter(MU == 1e-6), aes(x = m, y = meantime, 
                                                       color=factor(h)),
            size = 1.2, alpha = 0.8) +
  xlab(TeX("Dispersal rate, \\textit{d}")) + ylab("Generations") +
  labs(col=TeX("Dominance, \\textit{h}:"), title = TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +  
  my.theme + y.meantime + x.dispersal + dominance.color

#### B.
datapointsHS = read.csv("xxxx")
aux = function(ha, ea, phia){
  aa = datapointsHS %>% filter(m == 1) %>% filter(h == ha) %>% filter(e == ea) %>% filter(phi == phia)
  aa$meantime
}
datapointsHS$meantimeWM = mapply(aux, ha=datapointsHS$h, ea=datapointsHS$e, phia=datapointsHS$phi)
datapointsHS$ratio = datapointsHS$meantime / datapointsHS$meantimeWM
datapointsHS$ratio = datapointsHS$ratio * (datapointsHS$ratio<5.4) + 5.4*(datapointsHS$ratio > 5.4)
datapointsHS %>% filter(e == 0, phi == 0) %>%
  ggplot(aes(x = h, y = m, fill = ratio)) +
  geom_raster( interpolate = T ) +
  xlab(TeX("Genetic dominance, \\textit{h}")) + ylab(TeX("Dispersal rate, \\textit{d}")) +
  labs(fill= NULL,
       title = TeX("$\\tau_{HS} / \\tau_{HS}^{\\mathrm{WM}}}$")) + 
  geom_contour(aes(z = ratio), breaks = c(0.25, 0.5, 0.8, 2, 4), 
               colour = "black", linetype = "dashed",
               size = 0.6, alpha = 1) + 
  geom_contour(aes(z = ratio), breaks = c(1.000001), 
               colour = "black",
               size = 1) + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 16),
        legend.position = "top",
        plot.margin = margin(0.1,1.5,0.1,0.1, "cm")) +
  guides(fill = guide_colorbar(barwidth = 20, barheight = 0.5)) + 
  scale_x_continuous(breaks =c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0\nRecessive", "0.25", "0.5\nAdditive","0.75", "1\nDominant")) +
  scale_y_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001, 1E-4), 
                     labels = c("1", TeX("$10^{-1}$"), 
                                TeX("$10^{-2}$"), TeX("$10^{-3}$"),
                                TeX("$10^{-4}$"))) + 
  scale_fill_gradient2(low = "#7A00FF", #"#E65164", 
                       mid = "#FAFAFA", #"#FCFDBF", 
                       high= "#008D72", #"#646D7E", # "#1D1147", 
                       limits = c(0.18, 5.4),
                       trans = log2_trans(), 
                       breaks = c(0.25, 0.5, 1, 2, 4), 
                       labels = c("4x\nquicker", "2x", " \nequal", "2x", "4x\nslower") ) +
  coord_cartesian(xlim=c(0,1), ylim = c(10^{-4.25}, 1), expand = F)



####################
#### FigF1 ########
####################
aux_Charlesworth = function(sa, ha, ma, mua, sela, ea, ka, phia){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, n=200, N=100, selection=sela, e=ea, extc_param=c(ka, phia, phia, 0))
  tau = tryCatch({adaptation.meantime.Charlesworth2020(a)}, error=function(e) {return(NA)})
  return(tau)
}

data_th = read.csv("xxxx")
data_th$meantimeC = mapply(aux_Charlesworth,  sa = 0.01, 
                           ha = data_th$h, ma = data_th$m, 
                           mua = data_th$mu,  sela = "soft", 
                           ea=data_th$e, ka=100 , phia=data_th$phi)
data_th %>% 
  filter(h %in% c(0, 0.1, 0.5, 0.9, 1))  %>% filter(e == 0) %>% 
  ggplot(aes(x = m, y = meantimeC, color=factor(h))) + 
  geom_line(alpha=0.7, size=1.2) + 
  geom_line(aes(y = meantime), alpha=0.7, size=1.0, lty = "dashed" ) +
  # facet_grid(cols = vars(m), scales = "free_y",
  #            labeller = labeller(h = h.labs)) +
  xlab(TeX("Dispersal rate, d")) + 
  ylab(TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +
  labs(col="Genetic dominance, h:") +  
  my.theme + y.meantime + x.dispersal + dominance.color

##################
#### Fig3 ########
##################

#### A.
datapointsSS = read.csv("xxx")
datapoints.extra = expand.grid(s = 0.01, h = c(0, 0.5, 1), m = c(5E-4), mu = 1E-4, sel = "soft")
aux_time = function(sa, ha, ma, mua){
  adaptation.meantime.softsweep(s = c(-1E-3, 0.01), h=c(ha, ha), 
                                m=ma, n=200, N=100, mu=mua)
}
datapoints.extra$meantime = mapply(aux_time,  sa = datapoints.extra$s, 
                                   ha = datapoints.extra$h, ma = datapoints.extra$m, 
                                   mua = datapoints.extra$mu)
datapointsSS = rbind(datapointsSS, datapoints.extra)

data.raw1 =read.csv("xxxx", sep="")
data.raw2 =read.csv("xxxx", sep="")
data.raw3 =read.csv("xxx", sep="")
data.raw = rbind(data.raw2, data.raw3, data.raw1)
data.TFs = data.raw  %>%
  group_by(H, M, ID) %>% 
  summarise( TF = max(T) ) %>%
  group_by(H, M) %>% 
  summarise( TF.mean = mean(TF - 1E4),
             TF.sd = sd(TF - 1E4) )

datapointsSS %>% 
  filter(h %in% c(0, 0.5, 1) , case == "noshift" ) %>%
  filter(mu > 1e-05) %>%
  ggplot(aes(x = m, y = meantime)) + 
  geom_line(aes(color=factor(h)), alpha=0.7, size=1.1 ) + 
  geom_point(data = data.TFs, 
             aes(x = M, y = TF.mean, color=factor(H)),
             size = 2) + 
  # geom_errorbar(data = data.TFs, 
  #             aes(x = M, y = TF.mean, ymin = TF.mean - TF.sd, 
  #                 ymax = TF.mean + TF.sd, color=factor(H)),
  #             width = 0.1, alpha = 0.5) + 
  xlab(TeX("Dispersal rate, $\\textit{d}$")) + 
  ylab("Generations") +
  labs(col=TeX("Dominance, $\\textit{h}$:"), 
       title = TeX("Mean time of a soft sweep, $\\tau_{SS}$")) +  
  my.theme + theme(legend.position = "none") + 
  scale_y_continuous(trans = log10_trans(), 
                     breaks = trans_breaks("log10", function(x) 10^x), 
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(10**3, 10**4),
  ) + 
  dominance.color + x.dispersal


#### B.
M = c(1.0, 0.1, 0.01, 0.001) 
H = c(0, 0.5, 1.0) 
p0 = seq(0, 1, by = 0.01)
params = expand.grid(H = H, M = M)  

aux = function(haa, maa){
  aa = Roze(s = -1E-3, h=haa, m=maa, n = 200, N=100, mu = 1E-4, selection = "soft")
  stationary.FPEq(aa, res = 0.01, plot=F)
}

theoryDist = data.frame()
for (row in 1:nrow(params)) {
  ha = params[row, "H"]
  ma  = params[row, "M"]
  pts = aux(ha, ma)
  theoryDist = rbind(theoryDist, data.frame(H=ha, M=ma, P=pts[,1], FREQ=pts[,2]))
}

aux = function(mm){ 
  data.raw %>%
    filter(T < 1E4, T > 5E3, M == mm) %>%
    ggplot(aes(x = FREQ,fill = factor(H), group = factor(H))) + 
    geom_histogram(aes(y=stat(count/sum(count))), position = "identity", alpha = 0.5, size = 0.0, bins = 20) + 
    geom_line(data = filter(theoryDist, M == mm), aes(x = P, y = FREQ/100, color = factor(H)), size = 1.2) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "none",
          text = element_text(size = 16)) +
    xlab(NULL) +  ylab(NULL) + dominance.fill + dominance.color + 
    scale_y_continuous(limits = c(0, 0.12)) + 
    scale_x_continuous(
      limits = c(0, 0.5), 
      breaks = c(0, 0.2, 0.4), 
      labels = c("0", TeX("0.2"), TeX("0.4"))
    ) 
}; (aux(1) + aux(0.1))/(aux(0.01)+aux(0.001))

#### C. 
data.TFs = data.raw  %>%
  group_by(H, M, ID) %>% 
  summarise( TF = max(T) ) %>%
  group_by(H, M) %>% 
  summarise( TF.mean = mean(TF) )
aux = function(mm){ 
  data.raw %>%
    filter(M == mm, ID < 5) %>%
    ggplot(aes(x = T, y = FREQ, color = factor(H), group = interaction(H, ID))) +
    geom_line(alpha=0.7, size=0.8) + 
    geom_vline(xintercept = 1E4, lty = "dashed", color = "black") + 
    xlab("Generations") + ylab("Allele frequency") +
    labs(col=TeX("Dominance, \\textit{h}:")) +  
    my.theme + dominance.color + 
    scale_x_continuous(trans = log10_trans(), 
                       breaks = 1E4 + c(0, 10**3.2, 10**3.4, 10**3.6, 10**3.8), 
                       labels = c("   Env.\nchange", TeX("$10^{3.2}"), TeX("$10^{3.4}"), TeX("$10^{3.6}"), TeX("$10^{3.8}")),
                       limits = c(8.5E3, 1.8E4) )
}; aux(1) / aux(0.001)


##################
#### Fig4 ########
##################

#### A.
data.raw2 =read.csv("xxxx",
                    sep="")
data.raw3 =read.csv("xxxx",
                    sep="")
data.raw = rbind(data.raw2, data.raw3) 

aux = function(time, h){
  if(time <= 10000) return(h)
  else return(1 - h)
}
data.raw$Hnew = mapply(aux, time = data.raw$TIME, h = data.raw$H)

data.TFs = data.raw  %>%
  group_by(Hnew, M, ID) %>% 
  summarise( TF = max(T) ) %>%
  group_by(Hnew, M) %>% 
  summarise( TF.mean = mean(TF) )
aux = function(mm){ 
  data.raw %>%
    filter(M == mm, ID < 5) %>%
    ggplot(aes(x = T, y = FREQ, color = factor(Hnew), group = interaction(Hnew, ID))) + 
    geom_line(alpha=0.7, size=0.8) + 
    geom_vline(xintercept = 1E4, lty = "dashed", color = "black") + 
    xlab("Generations") + ylab("Allele frequency") +
    labs(col=TeX("Dominance, \\textit{h}:")) +  
    my.theme +  
    scale_x_continuous(trans = log10_trans(), 
                       breaks = 1E4 + c(0, 10**3.2, 10**3.4, 10**3.6, 10**3.8), 
                       labels = c("   Env.\nchange", TeX("$10^{3.2}"), TeX("$10^{3.4}"), TeX("$10^{3.6}"), TeX("$10^{3.8}")),
                       limits = c(8.5E3, 1.8E4) ) + 
    scale_color_manual(values = c("#FA7F5E", "#57157E") )
}; aux(1)/(aux(0.001) + xlab("Generations"))


#### B.
aux = function(mm){ 
  data.raw %>%
    filter(M == mm, ID < 10) %>%
    ggplot(aes(x = T, y = FREQ, color = factor(H), group = interaction(H, ID))) + 
    geom_line(alpha = 0.4, size = .8 ) + 
    geom_vline(aes(xintercept = 1E4), color = "white", 
               alpha = 1, size = 2.5 ) +
    geom_vline(aes(xintercept = 1E4), color = "black", 
               alpha = 0.75, size = 1) +
    geom_vline(data = filter(data.TFs, M == mm ), 
               aes(xintercept = TF.mean, color = factor(H)), 
               alpha = 1, size = 2, linetype = "dashed" ) + 
    geom_point(data = filter(data.FREQ0s, M == mm ), 
               aes(x = T, y = FREQ0.mean, color = factor(H), group = factor(H)), 
               alpha = 0.85, size = 2, pch = 16) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
          legend.position = "none",
          text = element_text(size = 16)) +
    xlab(NULL) +  ylab(NULL) + 
    scale_y_continuous(#trans = logit_trans(),
      limits = c(0, 1), 
      breaks = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c("Rare", TeX("$25\\%"), TeX("$50\\%"), TeX("75\\%"), " Fixed")
    ) + 
    scale_x_continuous(trans = log10_trans(), 
                       breaks = 1E4 + c(0, 10**3.2, 10**3.4, 10**3.6, 10**3.8), 
                       labels = c(TeX("Env. change"), TeX("$10^{3.2}"), TeX("$10^{3.4}"), TeX("$10^{3.6}"), TeX("$10^{3.8}")),
                       limits = c(8.5E3, 1.8E4) )
}; aux(1)/(aux(0.001) + xlab("Generations"))


##################
#### Fig5 ########
##################

#### A.
datapointsHS_ext$fixtime = datapointsHS_ext$timehs - 0.5/200/100/mu/datapointsHS_ext$fixprob

datapointsHS_ext %>% 
  filter(h==0 | h==0.5 | h==1) %>% 
  filter(e < 0.15) %>%
  filter(m %in% c(0.1, 0.001)) %>%
  filter(phi == 1) %>%
  ggplot(aes(x = e)) + 
  geom_ribbon(aes(ymin=1, 
                  ymax= ratio_log(0.5/200/100/mu/fixprob, fixtime)), 
              fill="gray", color = "transparent",
              alpha=0.5 ) +
  geom_ribbon(aes(ymin= ratio_log(0.5/200/100/mu/fixprob, fixtime), 
                  ymax = 0.5/200/100/mu/fixprob + fixtime), 
              fill="black", alpha=0.5, color = "transparent" ) + 
  geom_line(aes(y = 0.5/200/100/mu/fixprob + fixtime),
            color="black") + 
  facet_grid(m~h) +
  xlab(TeX("Extinction rate, \\textit{e}")) + 
  ylab("Generations") +
  labs(col=TeX("Dispersal rate, \\textit{d}:")) +  
  my.theme + x.dispersal +
  theme(legend.position = "none") + 
  scale_y_continuous(trans = log10_trans(),
                     limits = c(1, 3E6),
                     breaks = c(10, 1000, 1E5),
                     labels = c(TeX("$10$"), TeX("$10^3$"), TeX("$10^5$")) ) + 
  coord_cartesian(expand = F)

#### C.
S = 0.01
H = c(0, 0.5, 1)
M = c(1E-1, 1E-3) 
PHI = c(0, 1)
E = 2^-seq(1,12, by = 0.25)
datapoints = expand.grid(m = M, phi = PHI, e = E)  
aux_Neff = function(ma, ea, phia){
  a = Roze(s=0.01, h=0.5, m=ma, mu=0, n=200, N=100, e=ea, extc_param=c(100, phia, phia, 0))
  a$Neff
}
datapoints$Neff = mapply(aux_Neff, 
                         ma = datapoints$m,
                         phia = datapoints$phi,
                         ea = datapoints$e)
datapointsHS_ext %>%
  filter(e < 1) %>%
  ggplot(aes(x = e, y = Neff/(200 * 100), lty = factor(m), col = factor(m))) + 
  geom_hline(yintercept = 1, color = "black") + 
  geom_line(alpha = 0.7, size = 1.45) + 
  xlab(TeX("Extinction rate, \\textit{e}")) + 
  ylab(TeX("$\\textit{N}_{\\mathrm{e}}\ / \ nN$")) + 
  facet_grid(cols = vars(phi)) + 
  labs(title = "Effective population size",
       col = "Dispersal rate",
       lty = "Dispersal rate") +
  my.theme + x.dispersal + 
  scale_y_continuous(trans = log10_trans(), 
                     limits = c(0.04, 2.4),
                     breaks = c(0.05, 0.1, 0.5, 1, 2), 
                     labels = c("1/20", "1/10", "1/2", "1", "2") ) +
  scale_colour_brewer(palette = "Greys")


#### B.
M = c(1E-1, 1E-3) 
PHI = c(0, 1)
E = 2^-seq(1,12, by = 0.25)
datapointsHS_ext = expand.grid(s = S, h = H, m = M, phi = PHI, e = E)  
aux_prob = function(ha, ma, ea, phia){
  a = Roze(s=0.01, h=ha, m=ma, mu=0, n=200, N=100, e=ea, extc_param=c(100, phia, phia, 0))
  res = fixation.probability(a, res = 0.5/100/200)
  res[2,2]
}
datapointsHS_ext$probfix = mapply(aux_prob, 
                         ha = datapointsHS_ext$h,
                         ma = datapointsHS_ext$m,
                         phia = datapointsHS_ext$phi,
                         ea = datapointsHS_ext$e)

datapointsHS_ext %>%
  filter(e < 0.2) %>%
  ggplot(aes(x = e, y = prob*(200 * 200), lty = factor(m), col = factor(h), alpha = factor(m))) + 
  geom_line(size = 1.45) + 
  geom_hline(yintercept = 1, color = "black") +
  xlab(TeX("Extinction rate, \\textit{e}")) + 
  ylab(TeX("$\\textit{P}_{\\mathrm{fix}}$")) + 
  facet_grid(cols = vars(phi)) + 
  labs(title = "Probability of fixation",
       col = "h",
       lty = "d") +
  my.theme + x.dispersal + 
  scale_y_continuous(trans = log10_trans(), 
                     breaks = c(1, 10, 100, 500), 
                     labels = c("1", "10", "100", "500") ) +
  scale_colour_brewer(palette = "Dark2")


####################
###### Fig 1E ######
####################
datapoints_hardselection = read.csv("xxxxx")

datapoints_hardselection1 = read.csv("xxx", sep="")
datapoints_hardselection2 = read.csv("xxx", sep="")
datapoints_hardselection3 = read.csv("xxx", sep="")
datapoints.sim = rbind(datapoints_hardselection1,
                       datapoints_hardselection2,
                       datapoints_hardselection3)

dp.sim = datapoints.sim %>% 
  filter(FREQ > 0.9999) %>%
  group_by(Ndeme, Nind, S, H, M) %>%
  summarise(MEAN.TIME = mean(T),
            SD.TIME = sd(T) )

datapoints_hardselection %>%
  filter( mu == 1e-6, e == 0) %>%
  filter(h %in% c(0, 0.5, 1.0) ) %>% 
  ggplot(aes(x = m, y = meantime, color=factor(h), lty = factor(sel))) + 
  geom_line(size = 1.2, alpha = 0.8) +
  #geom_point(data = dp.sim, aes(x = M, y = MEAN.TIME, col = factor(H), lty = "soft") ) + 
  xlab(TeX("Dispersal rate, \\textit{d}")) + ylab("Generations") +
  labs(col=TeX("Dominance, \\textit{h}:"), title = TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +  
  my.theme + y.meantime + x.dispersal + dominance.color
  

datapoints_hardselection %>%
  filter( mu == 1e-6, e %in% c(0.1, 0.001)) %>%
  filter(h %in% c(0, 0.5, 1.0) ) %>% 
  ggplot(aes(x = m, y = meantime, color=factor(h), lty = factor(sel) )) + 
  geom_line(size = 1.2, alpha = 0.8) +
  facet_grid(e~phi, scales = "free_y", labeller = labeller(phi = phi.labs)) +
  xlab(TeX("Dispersal rate, \\textit{d}")) + ylab("Generations") +
  labs(col=TeX("Dominance, \\textit{h}:"), title = TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +  
  my.theme + y.meantime + x.dispersal + dominance.color



######################################
######  SM 2 -- Boxplots ####
######################################

#### A.
data.raw = read.csv("XXXXX", 
                    sep="")
data = data.raw %>%
  filter(FREQ > 0.99) %>%
  group_by(Ndeme, Nind, S, H, M, MU, ID) %>% 
  summarise(TF = max(TIME)) %>% 
  group_by(Ndeme, Nind, S, H, M, MU) %>% 
  summarise(TF.mean = mean(TF), 
            TF.median = median(TF),
            TF.var  = var(TF), 
            TF.sd  = sd(TF),
            TF.CV =  TF.sd / TF.mean,
            n = n() )

data.raw %>%
  filter(FREQ > 0.99) %>%
  group_by(Ndeme, Nind, S, H, M, MU, ID) %>% 
  summarise(TF = max(TIME)) %>% 
  group_by(Ndeme, Nind, S, H, M, MU) %>% 
  ggplot(aes(x=interaction(H, log10(M)), y=TF, fill=factor(H))) +
  geom_boxplot(alpha=0.7) +
  stat_summary(fun=n, geom = "text", hjust = 0.5) +
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="white") + my.theme +
  #theme(legend.position="none") +
  dominance.color + dominance.fill + y.meantime


###########################################
##### Comparing ours to Whitlock ##########
###########################################

aux_Roze = function(sa, ha, ma){
  a = Roze(s=sa, h=ha, m=ma)
  tau = tryCatch({fixation.meantime.singlemutant(a)}, error=function(e) {return(NA)})
  return(tau)
}

aux_Whitlock = function(sa, ha, ma){
  a = Whitlock(s=sa, h=ha, m=ma)
  tau = tryCatch({fixation.meantime.singlemutant(a)}, error=function(e) {return(NA)})
  return(tau)
}



data_th = expand.grid( h= c(0, 0.5, 1.0), m = 2^seq(-10, -1) )
data_th$meantimeR = mapply(aux_Roze,  sa = 0.01, 
                           ha = data_th$h, ma = data_th$m)
data_th$meantimeW = mapply(aux_Whitlock,  sa = 0.01, 
                           ha = data_th$h, ma = data_th$m)

data_th %>% 
  ggplot(aes(x = m, y = meantimeR, color=factor(h))) + 
  geom_line(alpha=0.7, size=1.2) + 
  geom_line(aes(y = meantimeW), alpha=0.7, size=1.0, lty = "dashed" ) +
  xlab(TeX("Dispersal rate, d")) + 
  ylab(TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +
  labs(col="Genetic dominance, h:") + 
  my.theme + y.meantime + x.dispersal + dominance.color

