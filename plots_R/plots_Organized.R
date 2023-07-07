## Organized version for the graphs of  
## Sudbrack PhD | Year 1 | Draft

############################################
########### Importing Libraries ############
############################################

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

##############################
########### Plots ############
##############################

# Hard sweeps ----

## Setting values -------
S = c(1e-2)
M = 2^-seq(0, 15, 0.25) #c(1, 0.1, 0.01, 0.001)
H = c(0, 0.5, 1) #seq(0, 1, by=0.05)
MU = c(1e-6)
SEL = c("soft") #, "hard")

datapointsHS = expand.grid(s = S, h = H, m = M, mu = MU, sel = SEL)  

## Function to calculate mean time to adaptation (tauHS)  -------

aux_time = function(sa, ha, ma, mua, sela){
  a = Roze(s = sa, h=ha, m=ma, mu=mua,n=200, N=100, selection = sela)
  adaptation.meantime.singlemutant(a)
}

aux_time_fix = function(sa, ha, ma, mua, sela){
  a = Roze(s = sa, h=ha, m=ma, mu = 0, n=200, N=100, selection = sela)
  fixation.meantime.singlemutant(a)
}

aux_neutral = function(ha, ma, mua, sela){
  a = Roze(s = 0, h=ha, m=ma, mu=mua,n=200, N=100, selection = sela)
  adaptation.meantime.singlemutant(a)
}

aux_neutral_fix = function(ha, ma, mua, sela){
  a = datapointsHS %>% filter(m == ma) %>% filter(h == ha) %>% filter(mu == mua) %>% filter(sel == sela)
  a$meantimeN - 1/mua
}

datapointsHS$meantime = mapply(aux_time,  sa = datapointsHS$s, ha = datapointsHS$h, ma = datapointsHS$m, mua = datapointsHS$mu,  sela = datapointsHS$sel)
datapointsHS$meantimeN = mapply(aux_neutral, ha = datapointsHS$h, ma = datapointsHS$m, mua = datapointsHS$mu, sela = datapointsHS$sel)
datapointsHS$meantimeNfix = mapply(aux_neutral_fix, ha = datapointsHS$h, ma = datapointsHS$m, mua = datapointsHS$mu, sela = datapointsHS$sel)
datapointsHS$meantimeC = mapply(aux_CW, ha = datapointsHS$h, ma = datapointsHS$m, mua = datapointsHS$mu, sela = datapointsHS$sel)

datapointsHS_cut$meantimeFix = mapply(aux_time_fix,  sa = datapointsHS_cut$s, ha = datapointsHS_cut$h, 
                                   ma = datapointsHS_cut$m, mua = datapointsHS_cut$mu,  sela = datapointsHS_cut$sel)

## Plotting -------

#### Plotting mean time to adaptation (tauHS) vs dominance -------
ggplot() +  
  geom_line(data = filter(datapointsHS, m < 1), 
            aes(x = h, y = meantime, color=factor(m)), 
            alpha=0.8, size=1.2) + 
  geom_line(data = filter(datapointsHS, m == 1), 
            aes(x = h, y = meantime), 
            linetype="dotted", alpha=0.6, size=0.8, color="black") + 
  facet_grid(cols = vars(mu), rows = vars(sel), scales = "free_y") +
  xlab(TeX("Dominance, $h$")) + ylab(TeX("Mean time of adaptation, $\\tau_{HS}$")) +
  labs(col="Dispersal, m") +  my.theme + log.y #+
  scale_color_brewer(palette="YlOrRd", direction = -1)


#### Plotting relative mean time to adaptation (tau^0HS) vs dominance -------
datapointsHS$m.round = round(datapointsHS$m, 5)
datapointsHS_cut = datapointsHS %>% 
    filter(m.round %in% c(1.0, 0.10511, 0.01105, 0.00116) )
datapointsHS_cut$m = signif(datapointsHS_cut$m, 1)
  
datapointsHS_cut  %>%
ggplot(aes(x = h, y = meantime/meantimeN, color=factor(m))) + 
geom_line(aes(group = m), alpha=0.8, size=1.5) + 
  #geom_line(data = filter(datapointsHS, m == 1), 
  #          aes(x = h, y = meantime/meantime0, group=sel), 
  #          linetype="dotted", alpha=0.6, size=0.8, color="black") + 
    facet_wrap(vars(mu), labeller = labeller(mu = mu.labs), scales = "free") +
  xlab(TeX("Genetic dominance, $h$")) + ylab(TeX("Relative mean time of adaptation, $\\tau^0_{HS}$")) +
  labs(col=TeX("Dispersal rate, $d$")) +  my.theme +
    scale_x_continuous(breaks =c(0, 0.25, 0.5, 0.75, 1),
                       labels = c("0\nRecessive", "0.25", "0.5\nAdditive","0.75", "1\nDominant"),
                       limits = c(0,1))
    

#### Plotting  mean time to adaptation (tauHS) vs migration ----------
datapoints %>% #filter(e == 0.01) %>%
  filter(h==0 | h==0.5 | h==1, mu == 1E-6) %>%
  ggplot(aes(x = rR, y = meantime)) + 
  geom_line( aes(color=factor(h)), alpha=0.6, size=1.6 ) + 
  #geom_line( aes(y = meantimeC, color=factor(h)), alpha=0.6, size=0.85, lty="dashed") + 
  #facet_grid(cols = vars(e), rows = vars(mu), 
  #           labeller = labeller(mu = mu.labs, sel=sel.labs))  + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5)) + 
  xlab(TeX("$F_{ST}$ (%)")) + # TeX("Dispersal rate, $d$")
  ylab(TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +
  labs(col=TeX("Dominance, $h$:"))  +
    theme(legend.position="top") +
  scale_y_continuous(trans = log10_trans(),
                     limits = c(10**3.5, 10**5),
                     breaks = trans_breaks("log10", function(x) 10^x), 
                     labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_continuous(trans = identity_trans(),
                       limits = c(0, 1),
                       breaks = 0.01*c(0, 2.5, 20, 25, 50, 70, 75, 96, 100), 
                       labels = c("Well mixed", "\nd = 0.1", "\nd = 0.01", "25", "50","\nd = 0.001", "75", "\nd = 0.0001", "100"))
  


#### Plotting relative mean time to adaptation (tau^0HS) vs migration --------
datapointsHS_cut %>% filter(h==0 | h==0.25 | h==0.5 | h==0.75 |h==1) %>%
  ggplot(aes(x = m, y = meantime/meantimeN)) + 
  geom_line( aes(color=factor(h)), alpha=0.8, size=1.1 ) + 
  facet_wrap(facets = vars(mu), labeller = labeller(mu=mu.labs), scales = "free_y") +
  xlab("Dispersal, d") + ylab(TeX("Relative mean time of adaptation, $\\tau^0_{HS}$")) +
  labs(col=TeX("Dominance, $h$"), title = "Regular demes") +  
    my.theme + log.x +
    theme(legend.position="top")
  
  
#### Plotting heatmap of ratio on parameter space dominance vs migration --------
aux = function(ha, ea, phia){
    aa = datapointsHS %>% filter(m == 1) %>% filter(h == ha) %>% filter(e == ea) %>% filter(phi == phia)
    aa$meantime
}
datapointsHS$meantimeWM = mapply(aux, ha=datapointsHS$h, ea=datapointsHS$e, phia=datapointsHS$phi)
datapointsHS$ratio = datapointsHS$meantime / datapointsHS$meantimeWM

# Cut-off
datapointsHS$ratio = datapointsHS$ratio * (datapointsHS$ratio<5.4) + 5.4*(datapointsHS$ratio > 5.4)

datapointsHS %>% #filter(e == 0, phi == 0) %>%
    ggplot(aes(x = h, y = m, fill = ratio)) +
    geom_raster( interpolate = T ) +
    xlab(TeX("Genetic dominance, \\textit{h}")) + ylab(TeX("Dispersal rate, \\textit{d}")) +
    labs(fill= NULL,
         title = TeX("$\\tau_{HS} / \\tau_{HS}^{\\mathrm{WM}}}$")) + 
  geom_contour(aes(z = ratio), breaks = c(0.25, 0.5, 0.8, 2, 4), 
               colour = "white", linetype = "dashed",
               size = 0.6, alpha = 1) + 
  geom_contour(aes(z = ratio), breaks = c(1.000001), 
               colour = "white",
               size = 1) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 16)) +
  scale_x_continuous(breaks =c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0\nRecessive", "0.25", "0.5\nAdditive","0.75", "1\nDominant")) +
  scale_y_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001, 1E-4), 
                     labels = c("1", TeX("$10^{-1}$"), 
                                TeX("$10^{-2}$"), TeX("$10^{-3}$"),
                                TeX("$10^{-4}$"))) + 
  scale_fill_gradient2(low = "#E65164", 
                       mid = "#FCFDBF", 
                       high= "#1D1147", 
                       limits = c(0.18, 5.4),
                       trans = log2_trans(), 
                       breaks = c(0.25, 0.5, 1, 2, 4), 
                       labels = c("4x quicker", "2x", "equal", "2x", "4x slower") ) +
  coord_cartesian(xlim=c(0,1), ylim = c(10^{-4.25}, 1), expand = F)
  

#### Plotting heatmap of meantime on parameter space dominance vs migration --------
p2 = datapointsHS %>% filter(mu == 1E-6) %>%
  ggplot(aes(x = h, y = m, fill = meantime )) + 
  geom_raster( interpolate = T ) +
  xlab(TeX("Genetic dominance, \\textit{h}")) + ylab(TeX("Dispersal rate, \\textit{d}")) +
  labs(fill=TeX("$\\tau_{HS}$"),
       title = "Expected time of a hard sweep") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5)) +
  scale_x_continuous(breaks =c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0\nRecessive", "0.25", "0.5\nAdditive","0.75", "1\nDominant")) +
  scale_y_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001), 
                     labels = c("Well mixed", TeX("$10^{-1}$"), TeX("$10^{-2}$"), TeX("$10^{-3}$"))) + 
  scale_fill_gradient(low = "darkblue", high = "red", 
                      trans = log10_trans(), 
                      breaks = 10000*c(0.5, 1, 2, 4), 
                      labels = c(TeX("5$\\cdot$10^{3}"), 
                                 TeX("1$\\cdot$10^{4}"), 
                                 TeX("2$\\cdot$10^{4}"), 
                                 TeX("4\\cdot$10^{4}")) ) + 
  coord_cartesian(xlim=c(0,1), ylim = c(10^{-3},1), expand = F)

p1 + p2 

## Time to intermediate frequency, 50% ----

aux_time_half = function(sa, ha, ma, mua){
  a = Roze(s = sa, h=ha, m=ma, mu=0, n=200, N=100)
  meantime.singlemutant(a, threshold = 0.5)
}

datapointsHS_cut$meantime50 = mapply(aux_time_half,  sa = datapointsHS_cut$s, ha = datapointsHS_cut$h, 
                               ma = datapointsHS_cut$m, mua = datapointsHS_cut$mu)

p1 = datapointsHS %>% filter(mu == 1E-4, h==0 | h==0.5 | h==1) %>%
  ggplot(aes(x = m, y = meantime50)) + 
  geom_line( aes(color=factor(h)), alpha=0.6, size=1.6 ) + 
  xlab("Dispersal, d") + ylab(TeX("$\\tau_{1/2}$")) +
  labs(col=TeX("Dominance, $h$:"), title = "Regular demes") +  my.theme + log.y +
  theme(legend.position="top", axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  scale_x_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001), 
                     labels = c("Well mixed      ", TeX("$10^{-1}$"), TeX("$10^{-2}$"), TeX("$10^{-3}$")))
         
p2 = datapointsHS_cut %>% filter(mu == 1E-4, h==0 | h==0.5 | h==1) %>%
  ggplot(aes(x = m, y = meantime50/meantimeFix)) + 
  geom_line( aes(color=factor(h)), alpha=0.6, size=1.6 ) + 
  geom_abline(intercept = 0.5, slope = 0, alpha=0.6, color = "black", lty = "dashed") + 
  xlab("Dispersal, d") + ylab(TeX("$\\tau_{1/2}/\\tau_{1}$")) +
  labs(col=NULL) +  my.theme +
  theme(legend.position="none", strip.text.x = element_blank(), axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) + ylim(c(0, 1)) + 
  scale_x_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001), 
                     labels = c("Well mixed   ", TeX("$10^{-1}$"), TeX("$10^{-2}$"), TeX("$10^{-3}$")))
p1 + p2

# Soft sweeps ----

## Setting values -------

S = c(1e-2)
H = seq(0, 1,by = 0.1) #c(0, 0.5, 1.0) #seq(0, 1,by = 0.1)  
M = 2^-seq(10, 15, 0.5) #c(1.0, 0.1, 0.01, 0.001) 
MU = c(1E-4, 1E-6) 
CASE = c("shift", "noshift")
SEL = c("soft")

datapointsSS = expand.grid(s = S, h = H, m = M, mu = MU, case = CASE, sel = SEL)  

## Function to calculate probability (piSS) and mean time to adaptation (tauSS)  -------

aux_time = function(case, sa,ha,ma,mua){
  if (case == "noshift")  adaptation.meantime.softsweep(s=c(-1E-3, sa), h=c(ha, ha),  m=ma, n=200, N=100, mu=mua)
  else adaptation.meantime.softsweep(s=c(-1E-3, sa), h=c(1-ha, ha),  m=ma, n=200, N=100, mu=mua)
}

aux_prob = function(case, sa,ha,ma, mua){
  if (case == "noshift") adaptation.fixation.softsweep(s=c(-1E-3, sa), h=c(ha, ha), m=ma, n=200, N=100, mu=mua)
  else adaptation.fixation.softsweep(s=c(-1E-3, sa), h=c(1-ha, ha), m=ma, n=200, N=100, mu=mua)

} 

datapointsSS$meantime = mapply(aux_time, case = datapointsSS$case, sa = datapointsSS$s, ha = datapointsSS$h, 
                               ma = datapointsSS$m, mua = datapointsSS$mu)
datapointsSS$prob     = mapply(aux_prob,  case = datapointsSS$case, sa = datapointsSS$s, ha = datapointsSS$h, 
                               ma = datapointsSS$m, mua = datapointsSS$mu)

## Simulations -------

data.raw1 =read.csv("~/Data/TrajectoriesSS_Recessive_new.dat", 
                    sep="")
data.raw2 =read.csv("~/Data/TrajectoriesSS_Additive_new.dat", 
                    sep="")
data.raw3 =read.csv("~/Data/TrajectoriesSS_Dominant_new.dat", 
                    sep="")

data.raw = rbind(data.raw2, data.raw3, data.raw1)

data.TFs = data.raw  %>%
  group_by(H, M, ID) %>% 
  summarise( TF = max(T) ) %>%
  group_by(H, M) %>% 
  summarise( TF.mean = mean(TF - 1E4),
             TF.sd = sd(TF - 1E4) ) # Minus the moment of env change

## Plotting -------

#### Plotting mean time to adaptation (pi_SS) vs migration --------
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
  labs(col=TeX("Dominance, $h_b$:"), title = TeX("Mean time of a soft sweep, $\\tau_{SS}$")) +  
  my.theme + y.meantime + dominance.color + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()) + x.dispersal

#### Plotting mean time to adaptation (pi_SS) vs dominance --------
datapointsSS %>% 
  filter(mu < 1e-05) %>%
  ggplot(aes(x = h, y = prob)) + 
  geom_line(aes(group=factor(m), color=m), alpha=0.75, size=0.5) + 
  facet_grid(rows = vars(mu), cols = vars(case), scales = "free_y",
             labeller = labeller(case=case.labs, mu=mu.labs)) +
  xlab(TeX("Dominance, $h_b$")) + ylab(TeX("Probability of a soft sweep, $\\pi_{SS}$")) +
  labs(col=TeX("Dispersal rate, $d$"), title = "Soft sweeps") +  my.theme + ylim(c(0,1)) +
  theme(legend.position="top", axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  scale_x_continuous(breaks =c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0\n       Recessive", "0.25", "0.5\nAdditive","0.75", "1\nDominant      ")) + 
  scale_color_gradient(low = "darkblue", 
                       high= "gray80", 
                       trans = log10_trans(), 
                       limits = c(0.00075, 1.25),
                       breaks = c(1, 0.1, 0.01, 0.001), 
                       labels = c("      Well mixed", TeX("$10^{-1}$"), TeX("$10^{-2}$"), TeX("$10^{-3}$")) )


#### Plotting mean time to adaptation (tau_SS) vs migration --------
datapointsSS %>% filter(h==0 | h==0.5 | h==1) %>%
  ggplot(aes(x = m, y = meantime)) + 
  geom_line(aes(color=factor(h)), alpha=0.7, size=1.1 ) + 
  facet_grid(cols = vars(mu), rows  = vars(case),
             labeller = labeller(case=case.labs, mu=mu.labs)) +  
  xlab(TeX("Dispersal, $d$")) + ylab(TeX("Mean time of a soft sweep, $\\tau_{SS}$")) +
  labs(col=TeX("Dominance, $h_b$"), title = "Soft sweeps") +  my.theme + y.meantime + x.dispersal



#### Plotting heatmap of dominance vs migration --------
aux = function(ha, mua, casea){
  aa = datapointsSS %>% filter(m == 1) %>% filter(h == ha) %>% filter(mu == mua) %>% filter(case == casea)
  aa$meantime
}
datapointsSS$meantimeWM = mapply(aux, ha=datapointsSS$h, mua=datapointsSS$mu, casea=datapointsSS$case)
datapointsSS$ratio = datapointsSS$meantime / datapointsSS$meantimeWM
datapointsSS %>%
  ggplot(aes(x = h, y = m, fill = ratio )) + 
  geom_raster( interpolate = TRUE ) + 
  facet_grid(rows = vars(mu), cols = vars(case), 
             labeller = labeller(mu = mu.labs, case = case.labs)) +
  xlab(TeX("Genetic dominance, $h_B$")) + ylab("Dispersal rate, d") +
  labs(fill=TeX("$\\tau_{SS}/\\tau_{SS}^{\\mathrm{WM}}}$"), 
       title = "Soft Sweeps") + 
  geom_contour(aes(z = ratio), breaks = c(0.25, 0.5, 0.8, 2, 4), 
               colour = "white", linetype = "dashed",
               size = 0.4, alpha = 1) + 
  geom_contour(aes(z = ratio), breaks = c(1.000001), 
               colour = "white",
               size = 0.99) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() +
  scale_x_continuous(breaks =c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0\nRecessive", "0.25", "0.5\nAdditive","0.75", "1\nDominant")) +
  scale_y_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001), 
                     labels = c("Well mixed", TeX("$10^{-1}$"), TeX("$10^{-2}$"), TeX("$10^{-3}$"))) + 
  scale_fill_gradient2(mid = "gray90", 
                       low = "orange", 
                       high= "purple", 
                       limits = c(0.2, 4.5),
                       trans = log2_trans(), 
                       breaks = c(0.25, 0.5, 1, 2, 4), 
                       labels = c("4x quicker", "2x", "= equal", "2x", "4x slower") )

# Hard sweeps with extinction ----

## Setting values -------
S = c(1e-2)
M = 2^seq(0, -10, -0.5) #c(1, 0.1, 0.01, 0.001) #
H = seq(0, 1, by=0.1) #c(0, 0.5, 1.0) #
MU = c(1e-6)
SEL = c("soft")
E = c(0) #, 1E-1, 1E-2, 1E-3) #2^seq(0, -10, -0.5) #
K = c(100)
PHI = c(0, 1)

datapointsHS = expand.grid(s = S, h = H, m = M, mu = MU, sel = SEL, e = E, k = K, phi = PHI)  

## Function to calculate mean time to adaptation (tauHS)  -------

aux_prob = function(sa, ha, ma, mua, sela, ea, ka, phia){
  a = Roze(s=sa, h=ha, m=ma, mu=0, n=200, N=100, selection=sela, e=ea, extc_param=c(ka, phia, phia, 0))
  fixation.probability(a, 0.5/100/200)[2,2]
}

aux_time = function(sa, ha, ma, mua, sela, ea, ka, phia){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, n=200, N=100, selection=sela, e=ea, extc_param=c(ka, phia, phia, 0))
  adaptation.meantime.singlemutant(a)
}

aux_neutral = function(ha, ma, mua, sela, ea, ka, phia){
  a = Roze(s=0, h=ha, m=ma, mu=mua, n=200, N=100, selection=sela, e=ea, extc_param=c(ka, phia, phia, 0))
  adaptation.meantime.singlemutant(a)
}

aux_Charlesworth = function(sa, ha, ma, mua, sela, ea, ka, phia){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, n=200, N=100, selection=sela, e=ea, extc_param=c(ka, phia, phia, 0))
  adaptation.meantime.Charlesworth2020(a)
}

aux_Ne = function(sa, ha, ma, mua){
  a = Roze(s=sa, h=ha, m=ma, mu=0, n=200, N=100, selection=sela, e=ea, extc_param=c(ka, phia, phia, 0))
  a$Neff
}

datapointsHS$prob = mapply(aux_prob,  sa = datapointsHS$s, ha = datapointsHS$h, ma = datapointsHS$m, mua = datapointsHS$mu,  sela = datapointsHS$sel, ea=datapointsHS$e, ka=datapointsHS$k , phia=datapointsHS$phi)
datapointsHS$meantime = mapply(aux_time,  sa = datapointsHS$s, ha = datapointsHS$h, ma = datapointsHS$m, mua = datapointsHS$mu,  sela = datapointsHS$sel, ea=datapointsHS$e, ka=datapointsHS$k , phia=datapointsHS$phi)
datapointsHS$meantime0 = mapply(aux_neutral, ha = datapointsHS$h, ma = datapointsHS$m, mua = datapointsHS$mu, sela = datapointsHS$sel, ea=datapointsHS$e, ka=datapointsHS$k , phia=datapointsHS$phi)

data_th$meantimeC = mapply(aux_Charlesworth,  sa = 0.01, 
                                ha = data_th$h, ma = data_th$m, 
                                mua = data_th$mu,  sela = "soft", 
                                ea=data_th$e, ka=100 , phia=data_th$phi)

datapointsHS$Neff = mapply(aux_Ne,  sa = datapointsHS$s, ha = datapointsHS$h, ma = datapointsHS$m, mua = datapointsHS$mu,  sela = datapointsHS$sel, ea=datapointsHS$e, ka=datapointsHS$k , phia=datapointsHS$phi)
datapointsHS$r1R = mapply(aux_r1R,  sa = datapointsHS$s, ha = datapointsHS$h, ma = datapointsHS$m, mua = datapointsHS$mu)

## Plotting -------
datapointsHS = read.csv("times_extinction.csv")

#### Comparing with Charlesworth approximation
data_th$meantimeC = mapply(aux_Charlesworth,  sa = 0.01, 
                           ha = data_th$h, ma = data_th$m, 
                           mua = data_th$mu,  sela = "soft", 
                           ea=data_th$e, ka=100 , phia=data_th$phi)
data_th %>% 
  filter(h %in% c(0, 0.1, 0.5, 0.9, 1))  %>%
  ggplot(aes(x = m, y = meantimeC, color=factor(h))) + 
  geom_line(alpha=0.7, size=1.2) + 
  geom_line(aes(y = meantime), alpha=0.7, size=1.0, lty = "dashed" ) +
  # facet_grid(cols = vars(m), scales = "free_y",
  #            labeller = labeller(h = h.labs)) +
  xlab(TeX("Dispersal rate, d")) + 
  ylab(TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +
  labs(col="Genetic dominance, h:") +  
  my.theme + y.meantime + x.dispersal + dominance.color



#### Plotting  mean time to adaptation (tauHS) vs migration ----------
datapointsHS %>% 
  filter(h==0 | h==0.5 | h==1) %>% filter( m == 0.001) %>%
  #filter(mu != 1E-4) %>%
  #filter(meantime < 10^6) %>% 
  ggplot(aes(x = e, y = fixtime)) + 
  geom_line( aes(color=factor(h), lty = factor(phi) ), alpha=0.7, size=1.1 ) + 
  # facet_grid(cols = vars(m), scales = "free_y",
  #            labeller = labeller(h = h.labs)) +
  xlab(TeX("Extinction rate,e")) + 
  ylab(TeX("Mean time of adaptation, $\\tau_{HS}$")) +
  labs(col="Dominance, h:",
       lty = "Recolonization:") +  
  my.theme + y.meantime + x.dispersal + dominance.color

  #+scale_color_grey(start = 0, end = 0.8)

#### Plotting relative mean time to adaptation (tau^0HS) vs migration --------
datapointsHS %>% filter(h==0 | h==0.5 | h==1) %>% #filter(mu != 1E-6) %>%
  #filter(meantime < 10^6) %>% 
  ggplot(aes(x = m, y = meantime/meantime0)) + 
  geom_line( aes(color=factor(e)), alpha=0.5, size=1.1 ) + 
  facet_grid(rows = vars(phi), cols = vars(h), scales = "free_y",
             labeller = labeller(h = h.labs, phi = phi.labs)) +
  xlab("Dispersal, d") + ylab(TeX("Relative mean time of adaptation, $\\tau^0_{HS}$")) +
  labs(col=TeX("Extinction rate,e:")) +  my.theme + y.meantime + x.dispersal + 
  theme(legend.position="top") +
  scale_color_grey(start = 0, end = 0.6)


#### Plotting the probability --------
datapointsHS %>% filter(h==0 | h==0.5 | h==1) %>%
  ggplot(aes(x = m, y = prob)) + 
  geom_line(aes(color = factor(e), lty = factor(phi)), size = 1.7) + 
  facet_wrap(vars(h), labeller = labeller(h=h.labs)) + 
  ylab("Fixation probablity") + xlab(TeX("Dispersal rate, d")) +
  labs(col="Extinction rate, e",
       lty="Recolonization model") +
  theme_bw( ) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 0, vjust = -0.3, hjust=0.5),
        legend.position="top",
        legend.box="vertical") +
  scale_x_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001), 
                     labels = c("Well mixed     ", TeX("$10^{-1}$"), 
                                TeX("$10^{-2}$"), TeX("$10^{-3}$"))) + 
  scale_y_continuous(trans = log10_trans(), 
                     breaks = c(1/200/200, 5e-4, 0.01, 0.02), 
                     labels = c("1/2nN", "s/20", "s", "2s")) + 
  coord_cartesian(ylim = c(1/200/200, 0.02)) +
  scale_color_manual(values = c("tomato1", "tomato3", "tomato4"))

#### Plotting the effective pop size --------
datapointsHS %>% 
  ggplot(aes(x = m, y = Neff)) + 
  #geom_hline(yintercept = 200*100, size = 1) + 
  geom_line(aes(color = factor(e), lty = factor(phi)), size = 1.7) + 
  ylab("Effective population size") + xlab(TeX("Dispersal rate, d")) +
  labs(col="Extinction rate, e",
       lty="Recolonization model") + theme_bw( ) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 0, vjust = -0.3, hjust=0.5),
        legend.position="top",
        legend.box="vertical") +
  scale_x_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001), 
                     labels = c("Well mixed     ", TeX("$10^{-1}$"), 
                                TeX("$10^{-2}$"), TeX("$10^{-3}$"))) + 
  scale_y_continuous(trans = log10_trans(), 
                     breaks = 200*100*c(1/3, 0.5, 1, 2, 3), 
                     labels = c("nN/3", "nN/2", "nN", "2nN", "3nN")) +
  coord_cartesian(ylim = (200*100*c(0.25, 4))) +
  scale_color_manual(values = c("black", "tomato1", "tomato3", "tomato4"))
  

#### Plotting the heatmap with extinction --------
aux = function(ha, mua, ea, phia){
  aa = datapointsHS %>% filter(m == 1) %>% filter(h == ha) %>% filter(mu == mua) %>% filter(e == ea) %>% filter(phi == phia)
  aa$meantime
}
datapointsHS$meantimeWM = mapply(aux, ha = datapointsHS$h, 
                                 mua = datapointsHS$mu, ea=datapointsHS$e, phia=datapointsHS$phi)
datapointsHS$ratio = datapointsHS$meantime / datapointsHS$meantimeWM

datapointsHS$ratio = pmin(datapointsHS$ratio, 8)

datapointsHS %>% filter(mu == 1E-6) %>%
  ggplot(aes(x = e, y = m, fill = ratio )) + 
  geom_raster( interpolate = TRUE ) + 
  facet_grid(rows = vars(phi), cols = vars(h), 
             labeller = labeller(phi = phi.labs, h = h.labs)) +
  xlab(TeX("Extinction rate, $e$")) + ylab(TeX("Dispersal rate, $d$")) +
  labs(fill=TeX("$\\tau_{HS} / \\tau_{HS}^{\\mathrm{WM}}$"), 
       #title = "Propagule model of recolonization"
       title = "Low mutation") +  
  geom_contour(aes(z = ratio), breaks = c(0.25, 0.5, 0.8, 2, 4), 
               colour = "white", linetype = "dashed",
               size = 0.4, alpha = 1) + 
  geom_contour(aes(z = ratio), breaks = c(2-1.000001), 
               colour = "white",
               size = 0.99) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() +
  scale_y_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001), 
                     labels = c("Well mixed", TeX("$10^{-1}$"), TeX("$10^{-2}$"), TeX("$10^{-3}$"))) +
  scale_x_continuous(trans = log10_trans(), 
                     breaks = c(1.0, 0.1, 0.01, 0.001), 
                     labels = c(TeX("$10^{0}$"), TeX("$10^{-1}$"), TeX("$10^{-2}$"), TeX("$10^{-3}$"))) + 
  scale_fill_gradient2(mid = "gray80", 
                       low = "orange", 
                       high= "purple", 
                       limits = c(0.12, 8.5),
                       trans = log2_trans(), 
                       breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8), 
                       labels = c("8x", "4x quicker", "2x", "= equal", "2x", "4x slower", "8x") ) + 
  coord_cartesian(xlim=c(1E-3, 1), ylim = c(1E-3, 1), expand = F)

datapointsHS %>% 
  ggplot(aes(x = h, y = m, fill = ratio )) + 
  geom_raster( interpolate = TRUE ) + 
  facet_grid(cols = vars(e), rows = vars(phi), 
             labeller = labeller(phi = phi.labs)) +
  xlab(TeX("Genetic dominance, $h$")) + ylab(TeX("Dispersal rate, $d$")) +
  labs(fill=TeX("$\\tau_{HS} / \\tau_{HS}^{\\mathrm{WM}}$")) +  
  geom_contour(aes(z = ratio), breaks = c(0.125, 0.25, 0.5, 0.8, 2, 4, 8.1), 
               colour = "white", linetype = "dashed",
               size = 0.4, alpha = 1) + 
  geom_contour(aes(z = ratio), breaks = c(1.000001), 
               colour = "white",
               size = 0.99) + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() +
  scale_y_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001), 
                     labels = c("Well mixed", TeX("$10^{-1}$"), TeX("$10^{-2}$"), TeX("$10^{-3}$"))) +
  scale_x_continuous(breaks =c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0\n            Recessive", "0.25", "0.5\nAdditive","0.75", "1\nDominant            ")) +
  scale_fill_gradient2(mid = "gray70", 
                       low = "orange", 
                       high= "purple", 
                       limits = c(0.12, 8.5),
                       trans = log2_trans(), 
                       breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8), 
                       labels = c("8x", "4x quicker", "2x", "= equal", "2x", "4x slower", "8x") ) + 
  coord_cartesian(xlim=c(0, 1), ylim = c(1E-3, 1), expand = F)

####### 
#### End of extinction
#######

#### Comparing mean time of hard sweeps with simulation --------

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

datapointsHS <- read.csv("~/sync/PhDyear1/datapointsHS_backup.csv")
datapoints.extra = expand.grid(s = 0.01, h = c(0, 0.5, 1), m = c(5E-4, 1E-4), mu = 1E-4)
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
  #geom_point(aes(y = TF.mean - TF.sd), size = 0.45, alpha = 0.75) + 
  #geom_point(aes(y = TF.mean + TF.sd), size = 0.45, alpha = 0.75) + 
  geom_line(data = data_th %>% filter(MU == 1e-4), aes(x = m, y = meantime, 
                                                       color=factor(h)),
            size = 1.2, alpha = 0.8) +
  xlab(TeX("Dispersal rate, \\textit{d}")) + ylab("Generations") +
  labs(col=TeX("Dominance, \\textit{h}:"), title = TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +  
  my.theme + y.meantime + x.dispersal + dominance.color


#### Testing the hard sweeps with extinctions and simulations ----

data.raw =read.csv("~/Data/MyDataExtMigration.dat", 
                    sep="")

data = data.raw %>%
  group_by(Model, Ndeme, Nind, S, H, M, MU, E) %>% 
  summarise(TF.mean = mean(TF), 
            TF.var  = var(TF), 
            TF.sd  = sd(TF) )

## Setting values -------
S = unique(data$S)
M = unique(data$M)
H = unique(data$H)
MU = unique(data$MU)
SEL = c("soft")
E = unique(data$E)
K = c(100)
PHI = c(0, 1)

data_theory = expand.grid(s = S, h = H, m = M, mu = MU, sel = SEL, e = E, k = K, phi = PHI)  

aux_time = function(sa, ha, ma, mua, sela, ea, ka, phia){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, n=200, N=100, selection=sela, 
           e=ea, extc_param=c(ka, phia, phia, 0) )
  adaptation.meantime.singlemutant(a)
}

data_theory$meantime = mapply(aux_time,  sa = data_theory$s, 
                               ha = data_theory$h, ma = data_theory$m, 
                               mua = data_theory$mu,  sela = data_theory$sel, 
                               ea=data_theory$e, ka=data_theory$k , phia=data_theory$phi)
data_theory = datapointsHS %>%
              filter(h %in% c(0, 0.5, 1), round(e,6) == 0.011049)

data_theory$H = data_theory$h
data_theory$M = data_theory$m
data_theory$E = data_theory$e
data_theory$TF.mean = data_theory$meantime

aux = function(phia){
  if (phia == 0) return("POOL")
  else return("PROPAGULE")
}
data_theory$Model = mapply(aux, data_theory$phi)



data %>% filter(E == 0.01, MU == 1E-6) %>%
  ggplot() + 
  aes(x = M, y = TF.mean, color=factor(H), pch = factor(Model)) + 
  geom_point(alpha=0.7, size=2.0) + 
  #geom_errorbar(aes(ymin = TF.mean-TF.sd, ymax = TF.mean+TF.sd),  
  #              alpha=0.5, size = 0.5, width = 0.05) +
  geom_line(data = filter(data_theory, e == 0.01, mu == 1E-6), 
            aes(lty = factor(Model)), alpha=0.7, size=1.4)  + 
  #facet_grid(E ~ Model, scales="free_y", labeller = labeller(MU = mu.labs)) +
  xlab(TeX("Dispersal rate, \\textit{d}")) + ylab(TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +
  labs(col=TeX("Dominance, \\textit{h}:"),
       pch = TeX("Recolonization model:"),
       lty = TeX("Recolonization model:")) +  
  my.theme + theme(legend.box="vertical") +
  y.meantime + x.dispersal


data %>% filter(MU == 1E-6) %>%
  filter(H == 0 | H == 0.5 | H == 1) %>%
  ggplot() + 
  aes(x = M, y = TF.mean, color=factor(H), group = interaction(H)) + 
  geom_point(alpha=0.7, size=1.6) + 
  #geom_line() + 
  #geom_errorbar(aes(ymin = TF.mean-TF.sd, ymax = TF.mean+TF.sd),  
  #              alpha=0.5, size = 0.5, width = 0.05) +
  geom_line(data = filter(data_theory, mu == 1E-6, e>0.001), alpha=0.5, size=1)  + 
  facet_grid(E ~ Model, labeller = labeller(PHI= phi.labs)) +
  xlab(TeX("Dispersal rate, \\textit{d}:")) + ylab(TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +
  labs(col=TeX("Dominance, \\textit{h}"), title = "Regular demes with extinction") +  
  my.theme + y.meantime + x.dispersal


datapointsHS %>% filter(MU == 1E-6) %>%
  filter(H == 0 | H == 0.5 | H == 1) %>%
  ggplot() + 
  aes(x = M, y = TF.mean, color=factor(H), group = interaction(H)) + 
  geom_point(alpha=0.7, size=1.6) + 
  #geom_line() + 
  #geom_errorbar(aes(ymin = TF.mean-TF.sd, ymax = TF.mean+TF.sd),  
  #              alpha=0.5, size = 0.5, width = 0.05) +
  geom_line(data = filter(data_theory, mu == 1E-6, e>0.001), alpha=0.5, size=1)  + 
  facet_grid(E ~ Model, labeller = labeller(PHI= phi.labs)) +
  xlab(TeX("Dispersal rate, \\textit{d}:")) + ylab(TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +
  labs(col=TeX("Dominance, \\textit{h}"), title = "Regular demes with extinction") +  
  my.theme + y.meantime + x.dispersal


aaa= times_extinction %>% filter(e< 0.5) %>%
  ggplot() + 
  aes(x = e, y = fixtime, 
      color=factor(m), lty=factor(phi),
      group = interaction(phi, m)) + 
  geom_line(alpha=0.7, size=1.6) +
  facet_grid(1 ~ h, #scales = "free_y",
             labeller = labeller(h = h.labs)) +
  #xlab(TeX("Extinction rate, \\textit{e}")) + 
  #ylab(TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +
  labs(col=TeX("Dispersal rate, \\textit{d}:")) +  
  my.theme + y.meantime + x.dispersal

bbb= times_extinction %>%filter(e< 0.5) %>%
  ggplot() + 
  aes(x = e, y = 0.5/200/100/1E-6/fixprob, 
      color=factor(m), lty=factor(phi),
      group = interaction(phi, m)) + 
  geom_line(alpha=0.7, size=1.6) +
  facet_grid(1 ~ h, #scales = "free_y",
             labeller = labeller(h = h.labs)) +
  xlab(TeX("Extinction rate, \\textit{e}")) + 
  #ylab(TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +
  labs(col=TeX("Dispersal rate, \\textit{d}:")) +  
  my.theme + y.meantime + x.dispersal +
  theme(legend.position = "none")

aaa/bbb

aux = function(ma){
  mu = 1E-4
  times_extinction %>% filter(e< 0.5, m == ma) %>%
    ggplot() + 
    aes(x = e) +
    geom_ribbon(aes(ymin=1, ymax=0.5/200/100/mu/fixprob), 
                fill="blue", alpha=0.5 ) +
    geom_ribbon(aes(ymin=0.5/200/100/mu/fixprob, ymax = 0.5/200/100/mu/fixprob + fixtime), 
                fill="red", alpha=0.5 ) +
    facet_grid(phi ~ h, #scales = "free_y",
               labeller = labeller(phi = phi.labs, h = h.labs)) +
    xlab(TeX("Extinction rate, \\textit{e}")) + 
    #ylab(TeX("Mean time of hard sweeps, $\\tau_{HS}$")) +
    labs(col=TeX("Dispersal rate, \\textit{d}:")) +  
    my.theme + y.meantime + x.dispersal +
    theme(legend.position = "none")
}

(aux(1) + aux(0.1))/(aux(0.01) + aux(0.001))
  
#### Testing probabilities of coalescence ----

S = c(1e-2)
M = c(1, 0.1, 0.01, 0.001) #2^seq(0, -10, -0.5) #
H = c(0.5) #seq(0, 1, by=0.1)
MU = c(1e-6)
SEL = c("soft")
E = c(0, 2^seq(0, -12, -0.25)) #c(0, 1E-1, 1E-2, 1E-3)
K = c(100)
PHI = c(0, 1)

datapointsHS = expand.grid(s = S, h = H, m = M, mu = MU, sel = SEL, e = E, k = K, phi = PHI)  

aux_aD = function(sa, ha, ma, mua, sela, ea, ka, phia){
  a = Roze(s=sa, h=ha, m=ma, mu=0, n=200, N=100, selection=sela, e=ea, extc_param=c(ka, phia, phia, 0))
  a$a.D
}
aux_cD = function(sa, ha, ma, mua, sela, ea, ka, phia){
  a = Roze(s=sa, h=ha, m=ma, mu=0, n=200, N=100, selection=sela, e=ea, extc_param=c(ka, phia, phia, 0))
  a$c.D
}

datapointsHS$aD = mapply(aux_aD,  sa = datapointsHS$s, ha = datapointsHS$h, ma = datapointsHS$m, mua = datapointsHS$mu,  sela = datapointsHS$sel, ea=datapointsHS$e, ka=datapointsHS$k , phia=datapointsHS$phi)
datapointsHS$cD = mapply(aux_cD,  sa = datapointsHS$s, ha = datapointsHS$h, ma = datapointsHS$m, mua = datapointsHS$mu,  sela = datapointsHS$sel, ea=datapointsHS$e, ka=datapointsHS$k , phia=datapointsHS$phi)

p1 = datapointsHS %>% 
  ggplot(aes(x = e, y = aD, color = factor(m))) + 
  geom_point(alpha = 0.5, size = 1.15) + 
  geom_line(alpha = 0.5, size = 1.15) + 
  facet_grid(rows = vars(phi), 
             labeller = labeller(phi = phi.labs)) +
  xlab(TeX("Extinction rate, $e$")) + ylab(TeX("Threeway prob. of coalescence, $a^\\mathrm{D}$")) +
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        strip.background = element_blank(), strip.text = element_blank(), legend.position="none") +
  theme_bw( ) + guides(color="none") + 
  scale_x_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001), 
                     labels = c(TeX("$10^{0}$"), TeX("$10^{-1}$"), TeX("$10^{-2}$"), TeX("$10^{-3}$"))) + 
  scale_y_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001), 
                     labels = c(TeX("$10^{0}$"), TeX("$10^{-1}$"), TeX("$10^{-2}$"), TeX("$10^{-3}$")))

p2 = datapointsHS %>% 
  ggplot(aes(x = e, y = cD, color = factor(m))) + 
  geom_point(alpha = 0.5, size = 1.15) + 
  geom_line(alpha = 0.5, size = 1.15) +   
  facet_grid(rows = vars(phi), 
             labeller = labeller(phi = phi.labs)) +
  xlab(TeX("Extinction rate, $e$")) + ylab(TeX("Threeway prob. of no coalescence, $c^\\mathrm{D}$")) +
  labs(color=TeX("Dispersal rate, $d$")) +  
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_bw() +
  scale_x_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001), 
                     labels = c(TeX("$10^{0}$"), TeX("$10^{-1}$"), TeX("$10^{-2}$"), TeX("$10^{-3}$"))) + 
  scale_y_continuous(trans = log10_trans(), 
                     breaks = c(1, 0.1, 0.01, 0.001), 
                     labels = c(TeX("$10^{0}$"), TeX("$10^{-1}$"), TeX("$10^{-2}$"), TeX("$10^{-3}$")))


p1+p2

#### Effective population size and selection gradient
S = c(1e-2)
M = 2^seq(0, -15, -0.25) #c(1, 1E-1, 1E-2, 1E-3) # #c(1, 1E-1, 1E-2, 1E-3) #
H = c(0, 0.5, 1)           #seq(0, 1, by=0.1) #
MU = c(0, 1e-4, 1e-6)
SEL = c("soft")
freq = seq(0, 1, by = 0.01)
  
datapoints = expand.grid(s = S, h = H, m = M, mu = MU, sel = SEL, freq=freq)  

aux_Neff = function(sa, ha, ma, mua, sela){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, n=200, N=100, selection=sela)
  a$Neff
}

aux_adv0 = function(sa, ha, ma, mua, sela){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, n=200, N=100, selection=sela)
  a$advection0
}

aux_adv1 = function(sa, ha, ma, mua, sela){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, n=200, N=100, selection=sela)
  a$advection1
}

aux_rD = function(sa, ha, ma, mua, sela){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, n=200, N=100, selection=sela)
  a$r.0.D
}

aux_rR = function(sa, ha, ma, mua, sela){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, n=200, N=100, selection=sela)
  a$r.1.R
}

aux_aR = function(sa, ha, ma, mua, sela){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, n=200, N=100, selection=sela)
  a$a.R
}

datapoints$Neff = mapply(aux_Neff,  sa = datapoints$s, ha = datapoints$h, ma = datapoints$m, 
                       mua = datapoints$mu,  sela = datapoints$sel)

datapoints$adv0 = mapply(aux_adv0,  sa = datapoints$s, ha = datapoints$h, ma = datapoints$m, 
                       mua = datapoints$mu,  sela = datapoints$sel)

datapoints$adv1 = mapply(aux_adv1,  sa = datapoints$s, ha = datapoints$h, ma = datapoints$m, 
                         mua = datapoints$mu,  sela = datapoints$sel)

datapoints$rD = mapply(aux_rD,  sa = datapoints$s, ha = datapoints$h, ma = datapoints$m, 
                         mua = datapoints$mu,  sela = datapoints$sel)

datapoints$rR = mapply(aux_rR,  sa = datapoints$s, ha = datapoints$h, ma = datapoints$m, 
                         mua = datapoints$mu,  sela = datapoints$sel)

datapoints$aR = mapply(aux_aR,  sa = datapoints$s, ha = datapoints$h, ma = datapoints$m, 
                         mua = datapoints$mu,  sela = datapoints$sel)

datapoints %>% filter(mu == 0, h == 0) %>%
  ggplot(aes(x = m, y = Neff/(200 * 100) )) + 
  #geom_point(alpha = 0.5, size = 1.15) + 
  geom_line(alpha = 0.65, size = 1.45) + 
  xlab(TeX("Dispersal rate, \\textit{d}")) + 
  ylab(NULL) + 
  labs(title = TeX("Ratio of effective population size, $\\textit{N}_{\\mathrm{e}}/nN$")) +
  my.theme + x.dispersal + 
  scale_y_continuous(trans = log10_trans(), 
                     breaks = c(1, 2, 5, 10, 50), 
                     labels = c("1", "2", "5", "10", "50") )


aaaaa = datapoints %>% filter(mu == 0) %>%
  ggplot(aes(x = freq, y = (adv0 + adv1 * freq)/0.01, color = factor(m))) + 
  geom_line(aes(group = m), size = 1.45, alpha = 0.65) + 
  facet_grid(cols = vars(h), labeller = labeller(h = h.labs)) + 
  xlab(TeX("Allele frequency, \\textit{p}")) + 
  ylab("Total selection = (a)+(b)-(c)") +
  labs(color=TeX("Dispersal rate, \\textit{d}:")) + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
        strip.text = element_blank(),
        legend.position="none") + 
  scale_y_continuous(breaks = c(0,  1/2, 1), 
                     labels = c("0.0", "0.5", "1.0")) + 
  scale_x_continuous(breaks = c(0, 1/4, 1/2, 3/4, 1), 
                     labels = c("0%", "25%", "50%", "75%", "100%")) #+ 
  scale_color_manual(values = c("green", "blue", "red", "gray20"))


aa = datapoints %>% filter(mu == 0) %>%
  ggplot(aes(x = freq, y = rD + (1-rD) * freq, color = factor(m))) + 
  geom_line(aes(group = m), size = 1.45, alpha = 0.65) + 
  facet_grid(cols = vars(h), labeller = labeller(h = h.labs)) + 
  #xlab(TeX("Allele frequency, \\textit{p}")) + 
  ylab(TeX("(a) Direct effects in AA")) + ylim(-1, 1) + 
  labs(color=TeX("Dispersal rate, \\textit{d}:")) + theme_bw() + 
  theme(plot.title = element_blank(), 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
        legend.position="top") 

 aaa = datapoints %>% filter(mu == 0) %>%
    ggplot(aes(x = freq, y = h*(1-rD)*(1-2*freq), color = factor(m))) + 
    geom_line(aes(group = m), size = 1.45, alpha = 0.65) + 
    facet_grid(cols = vars(h), labeller = labeller(h = h.labs)) + 
    #xlab(TeX("Allele frequency, \\textit{p}")) + 
    ylab(TeX("(b) Direct effects in aA")) + ylim(-1, 1) +
    labs(color=TeX("Dispersal rate, \\textit{d}:")) + theme_bw() + 
    theme(plot.title = element_blank(), 
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.text = element_blank(),
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
          legend.position="none") 
 
 aaaa =  datapoints %>% filter(mu == 0) %>%
   ggplot(aes(x = freq, y = aR + 2*(rR - aR)*(freq + (1-2*freq)*h), color = factor(m))) + 
   geom_line(aes(group = m), size = 1.45, alpha = 0.65) + 
   facet_grid(cols = vars(h), labeller = labeller(h = h.labs)) + 
   #xlab(TeX("Allele frequency, \\textit{p}")) + 
   ylab(TeX("(c) Kin competition")) + ylim(-1, 1) +
   labs(color=TeX("Dispersal rate, \\textit{d}:")) + theme_bw() + 
   theme(plot.title = element_blank(), 
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         strip.text = element_blank(),
         axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
         legend.position="none") 
 
 aa/aaa/aaaa/aaaaa
 
datapoints %>% filter(mu == 0) %>%
   ggplot(aes(x = freq, y = 2*(adv0 + adv1 * freq)*Neff, color = factor(m))) + 
   geom_line(aes(group = m), size = 1.45, alpha = 0.65) + 
   facet_grid(cols = vars(h), labeller = labeller(h = h.labs)) + 
   xlab(TeX("Allele frequency, \\textit{p}")) + 
   ylab("Total selection = (a)+(b)-(c)") +
   labs(color=TeX("Dispersal rate, \\textit{d}:")) + theme_bw() + 
   theme(plot.title = element_text(hjust = 0.5), 
         axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
         axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
         strip.text = element_blank(),
         legend.position="none") + 
   #scale_y_continuous(breaks = c(0,  1/2, 1), 
   #                    labels = c("0.0", "0.5", "1.0")) + 
   scale_x_continuous(breaks = c(0, 1/4, 1/2, 3/4, 1), 
                      labels = c("0%", "25%", "50%", "75%", "100%")) #+ 
 scale_color_manual(values = c("green", "blue", "red", "gray20"))
 
 #### Trajectories of hard sweeps ---- 
 
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

data.TFs.half = data.raw  %>%
  filter(FREQ > 0.45 & FREQ < 0.55) %>%
  group_by(H, M) %>% 
  summarise( TF.mean = mean(T) )

aux = function(mm){ 
  data.raw %>%
     filter(M == mm, ID < 10) %>%
     ggplot(aes(x = T, y = FREQ, color = factor(H), group = interaction(H, ID))) + 
     geom_line(alpha = 0.4, size = .8 ) + 
     geom_vline(data = filter(data.TFs, M == mm ), 
               aes(xintercept = TF.mean, color = factor(H)), 
               alpha = 1, size = 2.00, linetype = "dashed" ) + 
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
 

#### Trajectories of soft sweeps shift---- 

# data.raw1 =read.csv("~/Data/TrajectoriesSS_Additive.dat", 
#                     sep="")
data.raw2 =read.csv("~/Data/TrajectoriesSSShift_Recessive.dat",
                    sep="")
data.raw3 =read.csv("~/Data/TrajectoriesSSShift_Dominant.dat",
                    sep="")
# data.raw1 =read.csv("~/Data/TrajectoriesSS_Additive.dat", 
#                     sep="")
# data.raw2 =read.csv("~/Data/TrajectoriesSS_Recessive.dat", 
#                     sep="")
# data.raw3 =read.csv("~/Data/TrajectoriesSS_Dominant.dat", 
#                     sep="")

phantom = data.frame(
  Ndeme= 200, Nind = 100,      
  S    = -0.001, H    = 0.5,
  M    = c(1, 0.1, 0.01, 0.001),
  MU   = 1e-4, ID   = 1,  
  T    = 1E4, FREQ = 1.1)

data.raw = rbind(data.raw2, data.raw3) # data.raw1
  
data.TFs = data.raw  %>%
  group_by(H, M, ID) %>% 
  summarise( TF = max(T) ) %>%
  group_by(H, M) %>% 
  summarise( TF.mean = mean(TF) )

data.FREQ0s = data.raw  %>%
  filter(T == 1E4) %>% 
  group_by(H, M) %>% 
  summarise( FREQ0.mean = mean(FREQ), T = 1E4)

aux = function(h, t) {
  if (t <= 1E4) {return(1 - h)} else return(h)
}; data.raw$H = mapply(aux, h = data.raw$H, t =  data.raw$T)

### Plot of trajectories without dominance shifts ----- 
aux = function(mm){ 
  data.raw %>%
    filter(M == mm, ID < 8) %>%
    ggplot(aes(x = T, y = FREQ, color = factor(H), group = interaction(H, ID))) + 
         geom_line(alpha = 0.4, size = .8 ) + 
    # geom_vline(aes(xintercept = 1E4), color = "white", 
    #            alpha = 1, size = 2.5 ) +
    geom_vline(aes(xintercept = 1E4), color = "black", 
               alpha = 0.75, size = 1, linetype = "dotted") +
    geom_point(data = filter(data.TFs, M == mm ), 
               aes(x = TF.mean, y = 0, color = factor(H), group = factor(H)), 
               alpha = 1, size = 2 ) + 
    # geom_point(data = filter(data.FREQ0s, M == mm ), 
    #            aes(x = T, y = FREQ0.mean, color = factor(H), group = factor(H)), 
    #            alpha = 0.85, size = 2, pch = 16) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5),
          legend.position = "none",
          text = element_text(size = 16)) +
    xlab(NULL) +  ylab(NULL) + 
    scale_y_continuous(
      limits = c(0, 1), 
      breaks = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00")
    ) + scale_color_manual(values = c("#57157E", "#FA7F5E") ) + 
    scale_x_continuous(trans = log10_trans(), 
                       breaks = 1E4 + c(0, 10**3.2, 10**3.4, 10**3.6, 10**3.8), 
                       labels = c("  Env.\nchange", TeX("$10^{3.2}"), TeX("$10^{3.4}"), TeX("$10^{3.6}"), TeX("$10^{3.8}")),
                       limits = c(8.5E3, 1.8E4) )
}; aux(1)/(aux(0.001) + xlab("Generations"))

### Plot of trajectories with dominance shifts ----- 
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


#### Plot of initial conditions -----
aux = function(mm){ 
  data.raw %>%
    filter(T < 1E4, T > 1E3, M == mm) %>%
    ggplot(aes(x = FREQ, group = factor(H))) + 
    geom_histogram(aes(fill = factor(H)), position = "identity", alpha = 0.5, size = 0.0, bins = 25) + 
    geom_vline(data = filter(data.FREQ0s, M == mm ), 
               aes(xintercept = FREQ0.mean, color = factor(H)), 
               alpha = 0.75, size = 2, lty = "dashed") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.y = element_blank(),
          legend.position = "none",
          text = element_text(size = 16)) +
    xlab(NULL) +  ylab(NULL) + 
    scale_x_continuous(#trans = logit_trans(),
      limits = c(0, 0.5), 
      breaks = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c("Rare", TeX("$25\\%"), TeX("$50\\%"), TeX("75\\%"), " Fixed")
    ) 
}; (aux(1) + aux(0.1))/(aux(0.01)+aux(0.001))

### Plot of trajectiories wth dominance shift -----
aux = function(mm){ 
  data.raw %>%
    filter(M == mm, ID < 11, T>=1E4) %>%
    ggplot(aes(x = T, y = FREQ, color = factor(H), group = interaction(H, ID))) + 
    geom_line(alpha = 0.4, size = .8 ) + 
    geom_line(data = filter(data.raw, M == mm, ID < 11, T<=1E4), 
              aes(color = factor(1-H)),
              alpha = 0.4, size = .8 ) + 
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
}; aux(1)/aux(0.1)/aux(0.01)/(aux(0.001) + xlab("Generations"))

####  ----





#### Histograms of deme frequencies ----

data.raw1 = read.csv("~/Data/DemeFreqsSS_Dominant.dat", 
                    sep="")
data.raw2 = read.csv("~/Data/DemeFreqsSS_Recessive.dat", 
                    sep="")
data.raw3 = read.csv("~/Data/DemeFreqsSS_Additive.dat", 
                     sep="")
data.raw = rbind(data.raw1, data.raw2, data.raw3)

data.FREQ0s = data.raw  %>%
  filter(T > 0) %>% 
  group_by(H, M) %>% 
  summarise( FREQ0.mean = mean(FREQ))

aux = function(mm){ 
  data.raw %>%
    filter(H %in% c(0, 0.5, 1.0), T > 0, M == mm) %>%
    ggplot(aes(x = FREQ, group = factor(H))) + 
    geom_histogram(aes(fill = factor(H), y = ..density../100), 
                   position = "identity", 
                   alpha = 0.4, size = 0.0, 
                   bins = 100
                   ) + 
    geom_vline(data = filter(data.FREQ0s, M == mm ), 
               aes(xintercept = FREQ0.mean, color = factor(H)), 
               alpha = 0.75, size = 1, lty = "dashed") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text.y = element_blank(),
          legend.position = "none",
          text = element_text(size = 16)) +
    xlab(NULL) +  ylab(NULL) + ylim(0, 0.10) + 
    scale_x_continuous(#trans = logit_trans(),
      limits = c(0, 1), 
      breaks = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c("  Rare", TeX("$25\\%"), TeX("$50\\%"), TeX("75\\%"), "Fixed  ")
    )
}; (aux(1) + aux(0.1))/(aux(0.01)+aux(0.001))


#### Different initial conditions on soft sweeps ----

data.raw1 = read.csv("~/Data/SS_cond_res.dat", sep="")
data.raw2 = read.csv("~/Data/SS_cond_add.dat", sep="")
data.raw3 = read.csv("~/Data/SS_cond_dom.dat", sep="")
data.raw = rbind(data.raw1, data.raw2, data.raw3)

data = data.raw  %>%
  group_by(H, M, COND) %>% 
  summarise( TF.mean = mean(TF),
             TF.sd = sd(TF) )

aux_time = function(sa, ha, ma){
  a = Roze( s=sa, h=ha, m=ma, mu=0, n=200, N=100 )
  fixation.meantime(a)[[20,2]]
}

data$meantime = mapply(aux_time,  
                      sa = 0.01, 
                      ha = data$H, 
                      ma = data$M)


data %>%
  ggplot(aes(x = COND, 
             y = TF.mean/1000, 
             fill= factor(COND))) + 
  geom_bar(stat="identity",
           color = "white", orientation = "vertical", 
           position=position_dodge() ) +
  facet_grid(cols = vars(M), rows = vars(H), 
             scales = "free_x",
             labeller = labeller(H = h.labs)) + 
  #geom_errorbar(aes(ymin=TF.mean-TF.sd, ymax=TF.mean+TF.sd), 
  #              width=.2,
  #              position=position_dodge()) + 
  geom_hline( aes(yintercept = meantime/1000), 
              #lty = "dashed", 
              size = 1.2, 
              alpha = 0.6 ) + 
  my.theme +  theme(legend.position="none") + 
  xlab("Initial condition") +  
  ylab("Mean time (thousand generations)") + coord_flip() +
  scale_x_continuous(breaks = c(1,2,3,4), 
                      labels = c("4 demes of\nhomozygotes",
                                 "8 demes of\nheterozygotes",
                                 "200 demes of\n4 heterozygotes",
                                 "200 demes of\n2 homozygotes"))

#### Mean allele frequency at the moment of environmental change ----
M = 2^seq(0, -10, -0.25) #c(1, 1E-1, 1E-2, 1E-3) # #c(1, 1E-1, 1E-2, 1E-3) #
H = c(0, 0.5, 1)           #seq(0, 1, by=0.1) #
MU = c(1e-4)

data = expand.grid(h = H, m = M, mu = MU)

aux_phi_mean = function(h, m){
  a = Roze(s=-0.001, m=m, h=h, n=200, N=100, mu=1E-4, selection="soft")
  res0 = 0.5/200/100
  probdist = stationary.FPEq(a, res0)
  mm = weighted.mean(probdist[,1], probdist[,2]*res0)
  return(mm)
}

aux_phi_sd = function(h, m){
  a = Roze(s=-0.001, m=m, h=h, n=200, N=100, mu=1E-4, selection="soft")
  res0 = 0.5/200/100
  probdist = stationary.FPEq(a, res0)
  mm = weighted.mean(probdist[,1], probdist[,2]*res0)
  mmm = weighted.mean(probdist[,1]^2, probdist[,2]*res0)
  sd = sqrt(mmm - mm^2)
  return(sd) # return(mm) #
}

data$mean = mapply(aux_phi_mean, h = data$h, m = data$m)
data$sd   = mapply(aux_phi_sd  , h = data$h, m = data$m)

data.raw1 = read.csv("~/Data/FreqsSS_Dominant.dat", 
                     sep="")
data.raw2 = read.csv("~/Data/FreqsSS_Recessive.dat", 
                     sep="")
data.raw3 = read.csv("~/Data/FreqsSS_Additive.dat", 
                     sep="")
data.raw = rbind(data.raw1, data.raw2, data.raw3)

data.FREQ0s = data.raw  %>%
  filter(TF > 5000) %>% 
  group_by(H, M) %>% 
  summarise( FREQ0.mean = mean(FREQ))
data.FREQ0s$h = data.FREQ0s$H
data.FREQ0s$m = data.FREQ0s$M
data.FREQ0s$mean = data.FREQ0s$FREQ0.mean

data %>%
  ggplot(aes(x = m, 
             y = mean,
             fill = factor(h),
             color = factor(h)) ) + 
  geom_line(size = 1.5, alpha = 0.8) +
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd), 
             alpha=0.25) +
  geom_point(data = data.FREQ0s) + 
  my.theme +  dominance.color + 
  xlab(TeX("Dispersal rate, \\textit{d}")) + ylab(NULL) + 
  labs(title = "Mean allele frequency at the moment of environmental change",
       fill = TeX("Dominance, \\textit{h}:"),
       color = TeX("Dominance, \\textit{h}:")) + 
  scale_x_continuous(trans = log10_trans(), 
                     limits = c(0.0009, 1),
                     breaks = c(1, 0.1, 0.01, 0.001), 
                     labels = c(TeX("$1$"), 
                                TeX("$10^{-1}$"), 
                                TeX("$10^{-2}$"), 
                                TeX("$10^{-3}$"))) + 
  scale_fill_viridis_d(end = 0.8)

## Soft sweeps with extinciton ----
S = c(1e-2)
H = c(0, 0.5, 1.0) #seq(0, 1,by = 0.1)  
M = c(1.0, 0.1, 0.01, 0.001) 
MU = c(1E-4,) 
CASE = c("shift", "noshift")
SEL = c("soft")
E = c(0, 1E-1, 1E-2, 1E-3)
K = c(100)
PHI = c(0, 1)

datapointsSS = expand.grid(s = S, h = H, m = M, mu = MU, 
                           case = CASE, sel = SEL,
                           e = E, k = K, phi = PHI)  

## Function to calculate mean time to adaptation (tauSS)  -------
aux_time = function(sa, ha, case, ma, mua, sela, ea, ka, phia){
  if (case == "noshift")  adaptation.meantime.softsweep(s=c(-1E-3, sa), h=c(ha, ha),  
                                                        m=ma, n=200, N=100, mu=mua, e=ea, extc_param=c(ka, phia, phia, 0))
  else adaptation.meantime.softsweep(s=c(-1E-3, sa), h=c(1-ha, ha),  
                                     m=ma, n=200, N=100, mu=mua, e=ea, extc_param=c(ka, phia, phia, 0))
}


datapointsSS$meantime = mapply(aux_time, case = datapointsSS$case, sa = datapointsSS$s, 
                               ha = datapointsSS$h,
                               ma = datapointsSS$m, mua = datapointsSS$mu,
                               ea = datapointsSS$e, ka = datapointsSS$k, phia = datapointsSS$phi)

### Soft sweeps with extincito, simulations ----

list_files = list.files(path = "~/Data", pattern = "SSExtinction_*", full.names = T)
aux = function(x){
  temp = read.csv(x, sep="")
  temp$from = x
  temp
}
df2 = do.call(rbind, lapply(list_files, aux))
df2$TF = df2$T
df2$PHI = 1*(1-grepl("pool",  df2$from))

data.simulations = df2  %>%
              filter(TF > 5000, TF < 10000) %>% 
              group_by(H, M, E, PHI) %>% 
              summarise(FREQ0.mean = mean(FREQ))

df2  %>%
  filter(TF > 5000, TF < 10000) %>% 
  group_by(H, M, E, PHI) %>% 
  summarise(FREQ0.mean = mean(FREQ),
            FREQ0.sd = sd(FREQ)      ) %>% 
  filter(M %in% c(0.001, 0.1)) %>%
  ggplot(aes(x = E, y = FREQ0.mean, color = factor(H), fill = factor(H))) + 
  geom_line(size = 1) + 
  geom_ribbon(aes(ymin = FREQ0.mean-FREQ0.sd,
                  ymax = FREQ0.mean+FREQ0.sd), color = "transparent", alpha = 0.4) + 
  facet_grid(cols = vars(M), rows = vars(PHI)) + 
  xlab(TeX("Extinction rate, \\textit{e}")) + 
  ylab("Initial frequency") + dominance.color + dominance.fill + 
  my.theme + x.dispersal + 
  labs(title = "Mean allele frequency at the moment of environmental change",
       color = TeX("Dominance, \\textit{h}:"))

df2  %>%
  group_by(H, M, E, ID, PHI) %>% 
  summarise(TF.max = max(T)) %>%
  group_by(H, M, E, PHI) %>% 
  summarise(TF.mean = mean(TF.max)) %>%
  ggplot(aes(x = M, y = TF.mean, color = factor(H))) +
  geom_line(alpha = 0.65, size = 1.45) + 
  facet_grid(cols = vars(E), rows = vars(PHI))+
  xlab(TeX("Dispersal rate, \\textit{d}")) + 
  ylab("Generations") + dominance.color + 
  my.theme + x.dispersal + #y.meantime + 
  labs(title = "Mean time of soft sweep",
       color = TeX("Dominance, \\textit{h}:"))

df2  %>%
  group_by(H, M, E, ID, PHI) %>% 
  summarise(TF.max = max(T)) %>%
  group_by(H, M, E, PHI) %>% 
  summarise(TF.mean = mean(TF.max)) %>%
  filter(M %in% c(0.001, 0.1)) %>%
  ggplot(aes(x = E, y = TF.mean, color = factor(H))) +
  geom_line(alpha = 0.65, size = 1.45) + 
  facet_grid(cols = vars(M), rows = vars(PHI))+
  xlab(TeX("Extinction rate, \\textit{e}")) + 
  ylab("Generations") + dominance.color + 
  my.theme + x.dispersal + #y.meantime + 
  labs(title = "Mean time of soft sweep TSS",
       color = TeX("Dominance, \\textit{h}:"))


##### Hard selection soft sweeps -----
data.raw1 = read.csv("~/Data/SS_Recessive_hardselection_smalldemes.dat", sep="")
data.raw2 = read.csv("~/Data/SS_Additive_hardselection_smalldemes.dat",  sep="")
data.raw3 = read.csv("~/Data/SS_Dominant_hardselection_smalldemes.dat",  sep="")
data.raw = rbind(data.raw1, data.raw2, data.raw3)


data.raw %>% 
  filter(ID < 10) %>%
  ggplot(aes(x = T, y = FREQ, 
             color = factor(H), group = interaction(H, ID))) +
    geom_line(alpha = 0.5) + 
    facet_wrap(vars(M)) + dominance.color + 
    scale_x_continuous(trans = log10_trans(), 
                       breaks = trans_breaks("log10", function(x) 10^x), 
                       labels = trans_format("log10", math_format(10^.x)),
                       limits = c(100, 10^5)
    )

data.TFs = data.raw  %>%
  filter(FREQ < 1) %>%
  group_by(H, M, E, ID) %>% 
  summarise( TF = max(T) ) %>%
  group_by(H, M) %>% 
  summarise( TF.mean = mean(TF) )


aux_time = function(ha, ma, selectiona){
  tau =  tryCatch({adaptation.meantime.softsweep(s=c(-1E-3, 0.01), 
                                                 h=c(ha, ha), m=ma, n=2000, N=10, mu=1E-4, 
                                                 selection=selectiona)
                  }, error=function(e) {return(NA)})
  return(tau)
}

data.TFs$theoryHard =  mapply(aux_time,
                          ha = data.TFs$H, 
                          ma = data.TFs$M,
                          selectiona="hard")

data.TFs$theorySoft =  mapply(aux_time,
                              ha = data.TFs$H, 
                              ma = data.TFs$M,
                              selectiona="soft")


data.TFs %>%
  ggplot(aes(x = M, y = TF.mean, 
             color = factor(H) )) +
  geom_point(size = 2, alpha = 0.75) + 
  geom_line(size = 1.2, alpha = 0.75) + 
  #geom_line(aes(y = theorySoft), size = 1.2, alpha = 0.75, lty = "dashed") + 
  dominance.color + 
  x.dispersal + 
  y.meantime

  
## Setting values -------
S = c(1e-2)
M = 2^-seq(0, 12, 0.5) #c(1, 0.1, 0.01, 0.001)
H = seq(0, 1, by=0.5)
#MU = c(1e-6)
SEL = c("soft") #, "hard")

datapointsHS = expand.grid(s = S, h = H, m = M)  

## Function to calculate mean time to adaptation (tauHS)  -------
aux_time_fix = function(sa, ha, ma){
  a = Roze(s = sa, h=ha, m=ma, mu = 0, n=200, N=100)
  fixation.meantime.singlemutant(a)
}

aux_prob_fix = function(sa, ha, ma){
  a = Roze(s = sa, h=ha, m=ma, mu = 0, n=200, N=100)
  fixation.probability(a, 0.5/100/200)[2,2]
}

datapointsHS$time.fix = mapply(aux_time_fix, sa = datapointsHS$s, ha = datapointsHS$h, ma = datapointsHS$m)
datapointsHS$prob.fix = mapply(aux_prob_fix, sa = datapointsHS$s, ha = datapointsHS$h, ma = datapointsHS$m)

aux = function(mu){
  datapointsHS %>%
    ggplot(aes(x = M, y = 1/(2*200*100*mu*prob.fix) + time.fix), color = "black") +
    geom_ribbon(aes(ymin=1, ymax=1/(2*200*100*mu*prob.fix)), 
                fill="blue", alpha=0.5 ) +
    geom_ribbon(aes(ymin=1/(2*200*100*mu*prob.fix), ymax = 1/(2*200*100*mu*prob.fix) + time.fix), 
                fill="red", alpha=0.5 ) +
    geom_line(size = 1.5, alpha = 1) +
      x.dispersal + my.theme + y.meantime + 
      labs(y = "Mean time of hard sweep", 
           x = "Dispersal rate",
           title = sprintf("%1.0e", mu))  +
  coord_cartesian(expand = F)
}

aux(1e-8)+aux(1e-7)+aux(1e-6)+aux(1e-5)+aux(1e-4)+aux(1e-3)
aux(1e-7)+aux(1e-5)+aux(1e-3)

#### Probability of fixation as a function of limited dispersal ----
datapointsHS %>%
    ggplot(aes(x = m, y = prob.fix, color = factor(h))) +
    geom_line(size = 1.5, alpha = 1) +
    x.dispersal + my.theme  + dominance.color + 
    theme(legend.position = "none") +
    labs(title = "Probability of fixation", 
         x = TeX("Dispersal rate, $\\textit{d}$"),
         y = NULL) + 
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5, 2)/100,
                     labels = c("0.0%", "0.5%", "1.0%", "1.5%", "2.0%"))

######################

#### Directional mutations ----
data.raw1 = read.csv("~/Data/directionalMutationTeste_add.dat", 
                     sep="")
data.raw2 = read.csv("~/Data/directionalMutationTeste_res.dat", 
                     sep="")
data.raw3 =read.csv("~/Data/directionalMutationTeste_dom.dat", 
                    sep="")

rbind(data.raw1, data.raw2, data.raw3) %>%
  filter(M == 1e-1) %>%
ggplot(aes(x = T, y = FREQ, color = factor(H), group = interaction(H, ID))) + 
  geom_line(alpha = 0.4, size = .8 ) +
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
  scale_x_continuous(trans = log10_trans())


time.to.fix = 
  rbind(data.raw1, data.raw2, data.raw3) %>%
    filter(FREQ < 1) %>%
    group_by(H, M, ID) %>% 
    summarise( TF = max(T) ) %>%
    group_by(H, M) %>% 
    summarise( TF.mean = mean(TF) )
  
  
time.to.fix %>%
  ggplot(aes(x = M, y = TF.mean, color = factor(H))) + 
  geom_point() + 
  geom_line(data = filter(data_th, e==0), 
            aes(x = m, y = meantime, color = factor(h))) + 
  x.dispersal +
  dominance.color + 
  y.meantime



