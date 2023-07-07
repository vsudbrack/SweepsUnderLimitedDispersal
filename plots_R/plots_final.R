################################################
############ Figures from the paper ############
################################################

###### Both simulations and analytical comes from .csv files
###### See the folder simulations or analytical to codes to generate these files

ratio_log = function(a, b){
  # Relates a ratio to a value in logscale
  10 ^ ( (a / (a + b) ) * log10( a + b) )
}


####################
#### Fig1 #########
###################

#### A.
data.raw = read.csv("~/Data/tau1.dat", 
                     sep="")
data = data.raw %>%
  group_by(Ndeme, Nind, S, H, M, MU) %>% 
  summarise(TF.mean = mean(TF), 
            TF.var  = var(TF), 
            TF.sd  = sd(TF) )
data.theory = read.csv("~/sync/PhDyear1/times_fixation.csv", 
                    sep=",")

data %>%
  ggplot(aes(x = M*(0.95+0.1*H), y = TF.mean, color=factor(H) )) + 
  geom_point(alpha=0.9, size=2.0 ) + 
  geom_errorbar(aes(ymin = TF.mean - TF.sd, ymax = TF.mean + TF.sd), 
                width = 0.075, alpha = 0.7, size=1.0) + 
  geom_line(data = data.theory, aes(x = m, y = fixtime, color=factor(h)), 
            alpha=0.9, size=1.3) + 
  xlab(TeX("Dispersal rate, \\textit{d}")) + ylab("Generations") +
  labs(col=TeX("Dominance, \\textit{h}:"), title = TeX("Mean time to fixation, $\\tau_{1}$")) +  
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
data.raw = read.csv("~/Data/trajectory.dat", 
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

data.raw1 = read.csv("~/Data/MyDataTest_Recessive.dat", 
                     sep="")
data.raw2 = read.csv("~/Data/DataTest_Recessive2.dat", 
                     sep="")
data.raw3 =read.csv("~/Data/time_fixation_lowmigration_res.dat", 
                    sep="")
data.raw4 = read.csv("~/Data/time_fixation_lowmigration_add.dat", 
                     sep="")
data.raw5 = read.csv("~/Data/time_fixation_lowmigration_dom.dat", 
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

datapointsHS = read.csv("~/sync/LimitedDispersal/datapointsHS_extended_new.csv")
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
datapointsHS = read.csv("~/sync/PhDyear1/datapointsHS_extended_new.csv")
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
#### FigSM5 ########
####################
aux_Charlesworth = function(sa, ha, ma, mua, sela, ea, ka, phia){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, n=200, N=100, selection=sela, e=ea, extc_param=c(ka, phia, phia, 0))
  tau = tryCatch({adaptation.meantime.Charlesworth2020(a)}, error=function(e) {return(NA)})
  return(tau)
}

data_th = read.csv("~/sync/PhDyear1/datapointsHS_extended_new.csv")
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
datapointsSS = read.csv("~/sync/LimitedDispersal/datapointsSS_backup.csv")
datapoints.extra = expand.grid(s = 0.01, h = c(0, 0.5, 1), m = c(5E-4), mu = 1E-4, sel = "soft")
aux_time = function(sa, ha, ma, mua){
  adaptation.meantime.softsweep(s = c(-1E-3, 0.01), h=c(ha, ha), 
                                m=ma, n=200, N=100, mu=mua)
}
datapoints.extra$meantime = mapply(aux_time,  sa = datapoints.extra$s, 
                                   ha = datapoints.extra$h, ma = datapoints.extra$m, 
                                   mua = datapoints.extra$mu)
datapointsSS = rbind(datapointsSS, datapoints.extra)

data.raw1 =read.csv("~/Data/TrajectoriesSS_Recessive.dat", sep="")
data.raw2 =read.csv("~/Data/TrajectoriesSS_Additive.dat", sep="")
data.raw3 =read.csv("~/Data/TrajectoriesSS_Dominant.dat", sep="")
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
data.raw2 =read.csv("~/Data/TrajectoriesSSShift_Recessive.dat",
                    sep="")
data.raw3 =read.csv("~/Data/TrajectoriesSSShift_Dominant.dat",
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

#### B.
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

datapointsHS_ext %>%
  filter(e < 1) %>%
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









####################################
####################################
####################################




#### Fig1C --------

data.raw1 =read.csv("~/Data/Trajectories_Additive.dat", 
                    sep="")
data.raw2 =read.csv("~/Data/Trajectories_Recessive.dat", 
                    sep="") 
data.raw3 =read.csv("~/Data/Trajectories_Dominant.dat", 
                    sep="")
data.raw = rbind(data.raw1, data.raw2, data.raw3)

data.TFs = data.raw  %>%
  group_by(H, M, ID) %>% 
  summarise( TF = max(T) ) %>%
  group_by(H, M) %>% 
  summarise( TF.mean = mean(TF) )

aux = function(mm){ 
  data.raw %>%
    filter(M == mm, ID < 10) %>%
    ggplot(aes(x = T, y = FREQ, color = factor(H), group = interaction(H, ID))) + 
    geom_line(alpha = 0.6, size = .7 ) + 
    geom_point(data = filter(data.TFs, M == mm ), 
               aes(x = TF.mean, y = -0.0, group = factor(H), color = factor(H)), 
               alpha = 1, size = 3.00) + 
    my.theme + 
    theme(legend.position = "none") +
    xlab(NULL) +  ylab(NULL) + 
    scale_y_continuous(#trans = logit_trans(),
      limits = c(0, 1), 
      breaks = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00")
    ) + 
    scale_x_continuous(trans = log10_trans(), 
                       breaks = c(1E2, 1E3, 1E4, 1E5), 
                       labels = c(TeX("$10^2$"), TeX("$10^3$"), TeX("$10^4$"), TeX("$10^5$")),
                       limits = c(1E2, 5E5) ) + dominance.color
}; aux(1)/(aux(0.001) + xlab("Generations"))


#### Fig2 --------

data.raw1 = read.csv("~/Data/MyDataTest_Recessive.dat", 
                     sep="")
data.raw2 = read.csv("~/Data/DataTest_Recessive2.dat", 
                     sep="")
data.raw3 =read.csv("~/Data/time_fixation_lowmigration_res.dat", 
                    sep="")
data.raw4 = read.csv("~/Data/time_fixation_lowmigration_add.dat", 
                     sep="")
data.raw5 = read.csv("~/Data/time_fixation_lowmigration_dom.dat", 
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

datapointsHS = read.csv("~/sync/PhDyear1/datapointsHS_backup.csv")
datapoints.extra = expand.grid(s = 0.01, h = c(0, 0.5, 1), m = c(5E-4, 1E-4), mu = 1E-4, sel="soft")
aux_time = function(sa, ha, ma, mua){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, n=200, N=100, selection="soft")
  adaptation.meantime.singlemutant(a)
}
datapoints.extra$meantime = mapply(aux_time,  sa = datapoints.extra$s, 
                                   ha = datapoints.extra$h, ma = datapoints.extra$m, 
                                   mua = datapoints.extra$mu)
data_th = datapointsHS %>% filter(h %in% c(0, 0.5, 1))
data_th = rbind(data_th, datapoints.extra)
data_th$MU = data_th$mu

data %>%
  filter( MU == 1e-4) %>%
  ggplot(aes(x = M, y = TF.mean, color=factor(H) )) + 
  geom_point(alpha=0.9, size=1.8 ) + 
  geom_errorbar(aes(ymin = TF.mean - TF.sd, ymax =  TF.mean + TF.sd), 
                width = 0.05, size = 0.45, alpha = 0.75) + 
  geom_line(data = data_th %>% filter(MU == 1e-4), aes(x = m, y = meantime, 
                                                       color=factor(h)),
            size = 1.2, alpha = 0.8) +
  xlab(TeX("Dispersal rate, \\textit{d}")) + ylab("Generations") +
  labs(col=TeX("Dominance, \\textit{h}:"), title = TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +  
  my.theme + y.meantime + x.dispersal + dominance.color




### Fig3B shift ----

data.raw2 =read.csv("~/Data/TrajectoriesSSShift_Recessive.dat",
                    sep="")
data.raw3 =read.csv("~/Data/TrajectoriesSSShift_Dominant.dat",
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

##### Fig3A no shift ------
datapointsSS = read.csv("~/sync/PhDyear1/datapointsSS_backup.csv")
datapoints.extra = expand.grid(s = 0.01, h = c(0, 0.5, 1), m = c(5E-4), mu = 1E-4, sel = "soft")
aux_time = function(sa, ha, ma, mua){
  adaptation.meantime.softsweep(s = c(-1E-3, 0.01), h=c(ha, ha), 
                                m=ma, n=200, N=100, mu=mua)
}
datapoints.extra$meantime = mapply(aux_time,  sa = datapoints.extra$s, 
                                   ha = datapoints.extra$h, ma = datapoints.extra$m, 
                                   mua = datapoints.extra$mu)
datapointsSS = rbind(datapointsSS, datapoints.extra)

data.raw1 =read.csv("~/Data/TrajectoriesSS_Recessive.dat", sep="")
data.raw2 =read.csv("~/Data/TrajectoriesSS_Additive.dat", sep="")
data.raw3 =read.csv("~/Data/TrajectoriesSS_Dominant.dat", sep="")
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

#### Fig3C no shift ---------
data.raw1 =read.csv("~/Data/TrajectoriesSS_Recessive.dat", sep="")
data.raw2 =read.csv("~/Data/TrajectoriesSS_Additive.dat", sep="")
data.raw3 =read.csv("~/Data/TrajectoriesSS_Dominant.dat", sep="")
data.raw = rbind(data.raw2, data.raw3, data.raw1)

aux_dist = function(sa, ha, ma){
  res0 = 0.005
  model = Roze(s=-1e-3, m=ma, h=ha, n=200, N=100, mu=1E-4)
  aaa = stationary.FPEq(model, res0)
  data.frame(s = 0.01, H = ha, M = ma, mu = 1E-4, sel="soft", freq=aaa[,1], prob=aaa[,2])
}
datapoints.dist = data.frame()
for (ha in c(0, 0.5, 1)) 
  for (ma in c(1, 0.001))
    datapoints.dist = rbind(datapoints.dist, aux_dist(0, ha, ma))

data.raw  %>%
  filter(M %in% c(1,  0.001) ) %>%
  filter(TIME > 7500, TIME < 1E4) %>% 
  ggplot(aes(x = FREQ, color = NULL, fill=factor(H), group = interaction(H, M))) + 
  geom_histogram(binwidth = 0.02, 
                 alpha=0.5, position = "identity") +
  geom_line(data = datapoints.dist, 
            aes(x = freq, y = 2160*prob/275, color=factor(H)),
            size = 1.0) + 
  facet_grid(rows = vars(M)) + 
  xlab(TeX("Allele frequency at\nenvironmental change")) + 
  ylab(NULL) +
  labs(col=TeX("Dominance, $\\textit{h}$:"), 
       title = TeX("Histogram of\ninitial frequencies, $\\p_{0}$")) +  
  my.theme + theme(legend.position = "none") + 
  dominance.fill + dominance.color + 
  scale_x_continuous(limits = c(0,0.5))



##### Fig3A shift ------
data.raw1 =read.csv("~/Data/TrajectoriesSSShift_Recessive.dat", sep="")
data.raw3 =read.csv("~/Data/TrajectoriesSSShift_Dominant.dat", sep="")
data.raw = rbind(data.raw3, data.raw1)

datapointsSS = read.csv("~/sync/PhDyear1/datapointsSS_backup.csv")
datapoints.extra = expand.grid(s = 0.01, h = c(0, 1), m = c(5E-4), mu = 1E-4, sel = "soft")
aux_time = function(sa, ha, ma, mua){
  adaptation.meantime.softsweep(s = c(-1E-3, 0.01), h=c(1.0-ha, ha), 
                                m=ma, n=200, N=100, mu=mua)
}
datapoints.extra$meantime = mapply(aux_time,  sa = datapoints.extra$s, 
                                   ha = datapoints.extra$h, ma = datapoints.extra$m, 
                                   mua = datapoints.extra$mu)
datapoints.extra$case = "shift"
datapoints.extra$sel = "soft"
datapoints.extra$meantimeWM = NA
datapoints.extra$ratio = NA
datapoints.extra$X = NA
datapointsSS = rbind(datapointsSS, datapoints.extra)

data.TFs = data.raw  %>%
  group_by(H, M, ID) %>% 
  summarise( TF = max(T) ) %>%
  group_by(H, M) %>% 
  summarise( TF.mean = mean(TF - 1E4),
             TF.sd = sd(TF - 1E4) )
datapointsSS %>% 
  filter(h %in% c(0, 1) , case == "shift" ) %>%
  filter(mu > 1e-05) %>%
  ggplot(aes(x = m, y = meantime)) + 
  geom_line(aes(color=factor(h)), alpha=0.7, size=1.1 ) + 
  geom_point(data = data.TFs, 
             aes(x = M, y = TF.mean, color=factor(H)),
             size = 2) + 
  # geom_errorbar(data = data.TFs, 
  #               aes(x = M, y = TF.mean, ymin = TF.mean - TF.sd, 
  #               ymax = TF.mean + TF.sd, color=factor(H) ),
  #               width = 0.1, alpha = 0.5) +   
  xlab(TeX("Dispersal rate, $\\textit{d}$")) + 
  ylab("Generations") +
  labs(col=TeX("Dominance, $h_b$:"), title = TeX("Mean time of a soft sweep, $\\tau_{SS}$")) +  
  my.theme + theme(legend.position = "none") + 
  scale_y_continuous(trans = log10_trans(), 
                     breaks = trans_breaks("log10", function(x) 10^x), 
                     labels = trans_format("log10", math_format(10^.x)),
                     limits = c(10**3, 10**4),
  ) + scale_color_manual(values = c("#57157E", "#FA7F5E") ) + x.dispersal


##########################
##### Fig5 --- Extincitons
##########################



datapointsHS_ext = read.csv("times_extinction.csv")

S = c(1e-2)
M = c(0.1, 0.001) #2^seq(0, -10, -0.5) #
H = c(0, 0.5, 1.0) #seq(0, 1, by=0.1) #
MU = c(1e-6)
SEL = c("soft")
E = 2^seq(3, -10, -0.5) #c(1E-1, 1E-2, 1E-3) #
K = c(100)
PHI = c(0, 1)

datapointsHS_ext = expand.grid(s = S, h = H, m = M, mu = MU, sel = SEL, e = E, k = K, phi = PHI)  

## Function to calculate mean time to adaptation (tauHS)  -------

aux_prob = function(sa, ha, ma, mua, sela, ea, ka, phia){
  a = Roze(s=sa, h=ha, m=ma, mu=0, n=200, N=100, selection=sela, e=ea, extc_param=c(ka, phia, phia, 0))
  fixation.probability(a, 0.5/100/200)[2,2]
}

aux_time = function(sa, ha, ma, mua, sela, ea, ka, phia){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, n=200, N=100, selection=sela, e=ea, extc_param=c(ka, phia, phia, 0))
  tau = tryCatch({adaptation.meantime.singlemutant(a)}, error=function(e) {return(NA)})
  return(tau)
}

aux_Ne = function(sa, ha, ma, mua, sela, ea, ka, phia){
  a = Roze(s=sa, h=ha, m=ma, mu=0, n=200, N=100, selection=sela, e=ea, extc_param=c(ka, phia, phia, 0))
  a$Neff
}

datapointsHS_ext$Neff = mapply(aux_Ne,  sa = datapointsHS_ext$s, ha = datapointsHS_ext$h, ma = datapointsHS_ext$m, 
                               mua = datapointsHS_ext$mu, sela = datapointsHS_ext$sel, ea=datapointsHS_ext$e, 
                               ka=datapointsHS_ext$k, phia=datapointsHS_ext$phi)
datapointsHS_ext$fixprob = mapply(aux_prob,  sa = datapointsHS_ext$s, ha = datapointsHS_ext$h, ma = datapointsHS_ext$m, 
                           mua = datapointsHS_ext$mu,  sela = datapointsHS_ext$sel, ea=datapointsHS_ext$e, 
                           ka=datapointsHS_ext$k , phia=datapointsHS_ext$phi)
datapointsHS_ext$timehs = mapply(aux_time,  sa = datapointsHS_ext$s, ha = datapointsHS_ext$h, ma = datapointsHS_ext$m, 
                               mua = datapointsHS_ext$mu, sela = datapointsHS_ext$sel, ea=datapointsHS_ext$e, 
                               ka=datapointsHS_ext$k, phia=datapointsHS_ext$phi)


datapointsHS_ext = read.csv("datapointsHS_ext.csv")

mu = 1E-6
#datapointsHS_ext$timehs = datapointsHS_ext$fixtime
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
  

aa = datapointsHS %>% 
  filter(h==0 | h==0.5 | h==1) %>% 
  filter(e<0.5) %>%
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
  facet_grid(m~h, #scales = "free_y",
             labeller = labeller(h = h.labs)) +
  xlab(TeX("Extinction rate, \\textit{e}")) + 
  ylab("Generations") +
  labs(col=TeX("Dispersal rate, \\textit{d}:")) +  
  my.theme + x.dispersal +
  theme(legend.position = "none") + 
  scale_y_continuous(trans = log10_trans(),
                     limits = c(1, 9E5),
                     breaks = c(10, 1000, 1E5),
                     labels = c(TeX("$10$"), TeX("$10^3$"), TeX("$10^5$")) ) + 
  coord_cartesian(expand = F)

a / aa



#### Fig5B --------
# Effective population size in the presence of extinctions 
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

datapointsHS_ext %>%
  filter(e < 1) %>%
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

#### Fig SM1 --------
# Mean time to fixation with extinctions
datapointsSS = read.csv("datapointsSS_extinction.csv")
mu = 1E-4

datapointsSS %>% 
  ggplot(aes(x = e, y = meantime)) + 
  geom_line(aes(color=factor(h)), alpha=0.8, size=1.2) + 
  facet_grid(m~phi, scales = "free_y", labeller = labeller(phi = phi.labs)) +
  xlab(TeX("Extinction rate, \\textit{e}")) + 
  ylab("Generations") +
  labs(col=TeX("Genetic dominance, \\textit{h}:")) +  
  my.theme + x.dispersal + dominance.color + 
  theme() + y.meantime

datapointsSS %>% 
  ggplot(aes(x = e, y = meanp0)) + 
  geom_line(aes(color=factor(h)), alpha=0.8, size=1.2) + 
  facet_grid(m~phi, scales = "free_y", labeller = labeller(phi = phi.labs)) +
  xlab(TeX("Extinction rate, \\textit{e}")) + 
  ylab("Mean initial frequency") +
  labs(col=TeX("Genetic dominance, \\textit{h}:")) +  
  my.theme + x.dispersal + dominance.color + 
  theme() + scale_y_continuous(limits = c(0,1))


#####################
### Hard selection 
#####################
datapoints_hardselection = read.csv("~/sync/PhDyear1/datapoints_hardselection.csv")

## Simulations have Ndemes == 200, Nind == 100 (!)
datapoints_hardselection1 = read.csv("~/Data/SExtinction_Additive_hardselection.dat", sep="")
datapoints_hardselection2 = read.csv("~/Data/SExtinction_Recessive_hardselection.dat", sep="")
datapoints_hardselection3 = read.csv("~/Data/SExtinction_Dominant_hardselection.dat", sep="")
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









data.raw1 =read.csv("~/Data/TrajectoriesSS_Recessive.dat", sep="")
data.raw2 =read.csv("~/Data/TrajectoriesSS_Additive.dat", sep="")
data.raw3 =read.csv("~/Data/TrajectoriesSS_Dominant.dat", sep="")
data.raw = rbind(data.raw2, data.raw3, data.raw1)
data.TFs = data.raw  %>%
  group_by(H, M, ID) %>% 
  summarise( TF = max(T) ) %>%
  group_by(H, M) %>% 
  summarise( TF.mean = mean(TF - 1E4),
             TF.sd = sd(TF - 1E4) )



datapointsSS = read.csv("~/sync/LimitedDispersal/datapointsSS_backup.csv")
data.raw2 =read.csv("~/Data/TrajectoriesSSShift_Recessive.dat",
                    sep="")
data.raw3 =read.csv("~/Data/TrajectoriesSSShift_Dominant.dat",
                    sep="")
data.raw = rbind(data.raw2, data.raw3) 

data.raw$TIME = data.raw$T
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
    filter(M == mm, ID < 6) %>%
    ggplot(aes(x = T, y = FREQ, color = factor(Hnew), group = interaction(Hnew, ID))) + 
    geom_line(alpha=0.7, size=0.8) + 
    geom_vline(xintercept = 1E4, lty = "dashed", color = "black") + 
    xlab("Generations") + ylab("Allele frequency") +
    labs(col=TeX("Dominance, \\textit{h}:")) +  
    my.theme +  
    scale_x_continuous(trans = log10_trans(), 
                       breaks = 1E4 + c(0, 1500, 2500, 4000, 6500), 
                       labels = c("0", TeX("1500"), TeX("2500"), TeX("4000"), TeX("6500")),
                       limits = c(8.5E3, 1.8E4) ) + 
    scale_color_manual(values = c("#FA7F5E", "#6aa1b0") )
}; aux(1)/(aux(0.001) + xlab("Generations"))


dataSS.shift   = datapointsSS %>% filter(case == "shift") %>% filter(h == 1)
dataSS.shift$meantimeSHIFT   = dataSS.shift$meantime
dataSS.addtive = datapointsSS %>% filter(case == "shift") %>% filter(h == 0.5)
dataSS.addtive$meantimeNOSHIFT   = dataSS.addtive$meantime

dataALL = merge(dataSS.shift, dataSS.addtive, by.x = c("m", "mu"), 
                                    by.y = c("m", "mu"), all.x = FALSE, all.y = FALSE)

data.rawNS =read.csv("~/Data/TrajectoriesSS_Additive_new.dat", sep="")
data.rawNS$case = "NS"
data.rawSS =read.csv("~/Data/TrajectoriesSSShift_Dominant.dat", sep="")
data.rawSS$case = "SS"
data.raw = rbind(data.rawSS, data.rawNS)

data.TFs = data.raw  %>%
  group_by(H, M, ID, case) %>% 
  summarise( TF = max(T) ) %>%
  group_by(H, M, case) %>% 
  summarise( TF.mean = mean(TF - 1E4),
             TF.sd = sd(TF - 1E4), n = n() ) %>% filter(H != 0)

data.TFs.shift = data.TFs %>% filter(case == "SS") %>% filter(H == 1)
data.TFs.shift$TF.mean.SHIFT   = data.TFs.shift$TF.mean
data.TFs.addtive = data.TFs %>% filter(case == "NS") %>% filter(H == 0.5)
data.TFs.addtive$TF.mean.NOSHIFT   = data.TFs.addtive$TF.mean

dataALL.sim = merge(data.TFs.shift, data.TFs.addtive, by.x = c("M"), 
                by.y = c("M"), all.x = FALSE, all.y = FALSE)


dataALL %>% filter(mu == 1E-4) %>%
  ggplot(aes(x = m, y = meantimeSHIFT/meantimeNOSHIFT)) + 
  geom_line(alpha=0.7, size=1.0) + 
  geom_point(data = dataALL.sim, aes(x = M, y =TF.mean.SHIFT/TF.mean.NOSHIFT)) +  
  geom_hline(yintercept = 1, lty="dashed", col="darkgray") + 
  x.dispersal + 
  labs(x=TeX("Dispersal rate, \\textit{d}:")) +  
  my.theme

dataALL %>% filter(mu == 1E-4) %>%
  ggplot(aes(x = m)) + 
  geom_line(aes(y = meantimeSHIFT), size=1.2, col="#FA7F5E") + 
  geom_line(aes(y = meantimeNOSHIFT), alpha=0.7, size=1.2, col="#6aa1b0ff") +
  geom_point(data = filter(data.TFs, case =="SS", M>=1E-3), aes(x= M, y = TF.mean), alpha=0.7, size=2, col="#FA7F5E") + 
  geom_point(data = filter(data.TFs, case =="NS", M>=1E-3), aes(x= M, y = TF.mean),  alpha=0.7, size=2, col="#6aa1b0ff") +
  x.dispersal + 
  scale_y_continuous(trans = log10_trans(), 
                     limits = c(1100, 9E3) , 
                     breaks = trans_breaks("log10", function(x) 10^x), 
                     labels = trans_format("log10", math_format(10^.x)) )+ 
  labs(x=TeX("Dispersal rate, \\textit{d}")) +  
  my.theme


######################################
######  Boxplots of fixation data ####
######################################

#### A.
data.raw = read.csv("~/Data/trajectory.dat", 
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
