setwd("./XXXXX")
source("Model.R")

## Setting values -------
S = c(0.01)
H = seq(0, 1, by = 0.5)
MU = c(1e-6)
SEL = c("soft")
M = 2^seq(0, -13, -1) 

datapoints = expand.grid(s = S, h = H, 
                         Ndemes = 200, 
                         Nind = 100, mu = MU,
                         m = M, sel = SEL )

aux_time = function(ha, sa, ma, mua, sela){
  a = Roze(s=sa, h=ha, m=ma, mu=mua, 
           n=200, N=100, selection=sela)
 tau =  tryCatch({adaptation.meantime.singlemutant(a)}, error=function(e) {return(NA)})  
 #print(c(ha, sa, ma, mua, sela, tau))
}

datapoints$meantime = mapply(aux_time,  
                             ha  = datapoints$h, 
                             sa  = datapoints$s, 
                             ma  = datapoints$m, 
                             mua = datapoints$mu,  
                             sela= datapoints$sel)

write.csv(datapoints, file = "XXXXX.csv")


