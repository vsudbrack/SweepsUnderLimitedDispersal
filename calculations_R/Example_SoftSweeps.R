setwd("./XXXXX")
source("Model.R")

## Setting values -------
E = 2^seq(-1, -10, -0.5) 
H = c(0, 0.5, 1.0) 
SEL = c("soft")
M = c(1E-1, 1E-3)
PHI = c(0, 1)

datapointsSS = expand.grid(h = H, m = M, 
                           sel = SEL, e = E, phi = PHI)  

## Function to calculate mean time of soft sweeps -------
aux_time = function(ha, ma, mua, ea, phia){
  extc_param=c(100, phia, phia, 0)
  res0 = 0.5/200/100
  model = Roze(s=-1E-3, m=ma, h=ha, n=200, N=100, mu=1E-4, selection="soft", e = ea, extc_param=extc_param)
  probdist = stationary.FPEq(model, res0)
  
  samples = sample(x = probdist[,1], prob = probdist[,2]*res0, size = 1000, replace=T)
  
  model = Roze(s=1E-2, m=ma, h=ha, n=200, N=100, mu=0, selection="soft", e = ea, extc_param=extc_param)
  time = fixation.meantime.points(model, p0 = samples)
  
  c( mean(time[,2]), mean(samples) )
}

datapointsSS$meantime = NA
datapointsSS$meanp0 = NA

for(i in 1:nrow(datapointsSS)) { # for-loop over rows
  print(paste(datapointsSS[i, "m"], datapointsSS[i, "e"]))
  result = aux_time(ha = datapointsSS[i, "h"], 
                    ma = datapointsSS[i, "m"], 
                    mua = datapointsSS[i, "mu"],
                    ea = datapointsSS[i, "e"], 
                    phia = datapointsSS[i, "phi"] )
  datapointsSS[i, "meantime"] = result[1]
  datapointsSS[i, "meanp0"]   = result[2]
  write.csv(datapointsSS, file = "datapointsSS_extinction.csv", append = F)
}

