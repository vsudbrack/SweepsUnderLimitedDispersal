#### C. 
data.rawNS =read.csv("~/Data/TrajectoriesSS_Additive_new.dat", sep="")
data.rawNS$case = "NS"
data.rawSS =read.csv("~/Data/TrajectoriesSSShift_Dominant.dat", sep="")
data.rawSS$case = "SS"
data.raw = rbind(data.rawSS, data.rawNS)

data.raw$TIME = data.raw$T

aux = function(time, h){
  if(time <= 1E4) return(h)
  else return(1 - h)
}
data.raw$Hnew = mapply(aux, time = data.raw$TIME, h = data.raw$H)
data.raw$H = data.raw$Hnew 
  
data.TFs = data.raw  %>%
  group_by(H, M, ID) %>% 
  summarise( TF = max(T) )


aux = function(disp) {
  
  sample.additive  = sample(unique(data.raw %>% filter(H == 0.5) %>% select(ID))$ID, 10)
  sample.dominant  = sample(unique(data.raw %>% filter(H == 0.0) %>% select(ID))$ID, 10)
  
  data.raw.mean = data.raw %>% 
    filter(M == disp) %>% 
    group_by(Ndeme, Nind, S, H, M, TIME) %>% 
    summarise(h = H, t = TIME,
              Nseg = n(), #length(unique((data.raw%>%filter(H == h, TIME == t)%>%select(ID))$ID)),
              Nrep = length(unique((data.raw%>%filter(H == h, TIME == 0)%>%select(ID))$ID)),
              FREQ.mean = (sum(FREQ) + (Nrep - Nseg)*1.0) / Nrep  )
  
  data.raw %>% 
    filter(M == disp) %>% 
    filter(     ((H == 0.5)&(ID %in% sample.additive))  |
                ((H == 0.0)&(ID %in% sample.dominant))  ) %>% 
    ggplot(aes(x = TIME, y = FREQ, color=factor(H))) + 
    geom_line(aes(group = interaction(ID, H)), alpha=0.4, size=0.25) + 
    geom_line(data = data.raw.mean, aes(y = FREQ.mean), alpha=0.8, size=1.5) + 
    geom_vline(xintercept = 1E4, lty = "dashed", color = "darkgray") + 
    xlab("Generations") + ylab("Allele frequency") +
    labs(col=TeX("Dominance, \\textit{h}:")) +  
    my.theme + theme(legend.position = "none") +  
    scale_color_manual(values = c("#FA7F5E", "#6aa1b0") ) +    
    scale_x_continuous(trans = log10_trans(), 
                       breaks = 1E4 + c(0, 1000, 2500, 4000, 6000), 
                       labels = c("0", "1000", "2500", "4000", "600"),
                       limits  = c(9500, 16000) )
}; aux(1) / aux(0.001)

