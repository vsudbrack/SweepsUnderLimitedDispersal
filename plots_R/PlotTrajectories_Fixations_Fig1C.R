#### 0C. 
data.raw = read.csv("~/Data/trajectory.dat", 
                    sep="")

aux = function(disp) {
  data.raw.fixed = data.raw %>%
    filter( M == disp ) %>%
    filter( ((H == 0.0)&(ID %in% unique(data.raw %>% filter(H == 0.0, M==disp, FREQ>0.95) %>% select(ID))$ID)) |
            ((H == 0.5)&(ID %in% unique(data.raw %>% filter(H == 0.5, M==disp, FREQ>0.95) %>% select(ID))$ID)) |
            ((H == 1.0)&(ID %in% unique(data.raw %>% filter(H == 1.0, M==disp, FREQ>0.95) %>% select(ID))$ID))  ) 
  
  sample.recessive = sample(unique(data.raw.fixed%>%filter(H == 0.0)%>%select(ID))$ID, 10)
  sample.additive  = sample(unique(data.raw.fixed%>%filter(H == 0.5)%>%select(ID))$ID, 10)
  sample.dominant  = sample(unique(data.raw.fixed%>%filter(H == 1.0)%>%select(ID))$ID, 10)
  
  data.raw.fixed.mean = data.raw.fixed %>% 
                          group_by(Ndeme, Nind, S, H, M, TIME) %>% 
                          summarise(h = H, t = TIME,
                                    Nseg = length(unique((data.raw.fixed%>%filter(H == h, TIME == t)%>%select(ID))$ID)),
                                    Nrep = length(unique((data.raw.fixed%>%filter(H == h, TIME == 0)%>%select(ID))$ID)),
                                    FREQ.mean = (sum(FREQ) + (Nrep - Nseg)*1.0) / Nrep  )
  
  data.raw.fixed %>% 
    filter( ((H == 0.0)&(ID %in% sample.recessive)) |
            ((H == 0.5)&(ID %in% sample.additive)) |
            ((H == 1.0)&(ID %in% sample.dominant))  ) %>% 
    ggplot(aes(x = TIME, y = FREQ, color=factor(H))) + 
    geom_line(aes(group = factor(ID)),alpha=0.4, size=0.25) + 
    geom_line(data = data.raw.fixed.mean, aes(y = FREQ.mean),  alpha=0.8, size=1.5) + 
    xlab("Generations") + ylab("Allele frequency") +
    labs(col=TeX("Dominance, \\textit{h}:")) +  
    my.theme + dominance.color + 
    scale_x_continuous(trans = log10_trans(),
                       limits = c(500, 15000))
}; aux(1) / aux(0.001)

