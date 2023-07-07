#### C. 
data.raw1 =read.csv("~/Data/TrajectoriesSS_Recessive.dat", sep="")
data.raw2 =read.csv("~/Data/TrajectoriesSS_Additive.dat", sep="")
data.raw3 =read.csv("~/Data/TrajectoriesSS_Dominant.dat", sep="")
data.raw = rbind(data.raw2, data.raw3, data.raw1)
data.TFs = data.raw  %>%
  group_by(H, M, ID) %>% 
  summarise( TF = max(T) )

data.raw$TIME = data.raw$T

aux = function(disp) {
  
  sample.recessive = sample(unique(data.raw %>% filter(H == 0.0) %>% select(ID))$ID, 10)
  sample.additive  = sample(unique(data.raw %>% filter(H == 0.5) %>% select(ID))$ID, 10)
  sample.dominant  = sample(unique(data.raw %>% filter(H == 1.0) %>% select(ID))$ID, 10)
  
  data.raw.mean = data.raw %>% 
    filter(M == disp) %>% 
    group_by(Ndeme, Nind, S, H, M, TIME) %>% 
    summarise(h = H, t = TIME,
              Nseg = n(), #length(unique((data.raw%>%filter(H == h, TIME == t)%>%select(ID))$ID)),
              Nrep = length(unique((data.raw%>%filter(H == h, TIME == 0)%>%select(ID))$ID)),
              FREQ.mean = (sum(FREQ) + (Nrep - Nseg)*1.0) / Nrep  )
  
  data.raw %>% 
    filter(M == disp) %>% 
    filter(   ((H == 0.0)&(ID %in% sample.recessive)) |
              ((H == 0.5)&(ID %in% sample.additive))  |
              ((H == 1.0)&(ID %in% sample.dominant))  ) %>% 
    ggplot(aes(x = TIME, y = FREQ, color=factor(H))) + 
    geom_line(aes(group = interaction(ID, H)), alpha=0.4, size=0.25) + 
    geom_line(data = data.raw.mean, aes(y = FREQ.mean), alpha=0.8, size=1.5) + 
    geom_vline(xintercept = 1E4, lty = "dashed", color = "darkgray") + 
    xlab("Generations") + ylab("Allele frequency") +
    labs(col=TeX("Dominance, \\textit{h}:")) +  
    my.theme + dominance.color + theme(legend.position = "none") +      
    scale_x_continuous(trans = log10_trans(), 
                       breaks = 1E4 + c(0, 1000, 2500, 4000, 6000), 
                       labels = c("0", "1000", "2500", "4000", "6000"),
                       limits = c(9500, 16000) )
}; aux(1) / aux(0.001)

