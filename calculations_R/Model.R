############################################
########## Analytical model ##############
############################################

### Creates the model ####
Roze = function (s, m, n=200, N=100, h=0.5, selection="soft", alpha=NULL, e=0, mu=0, extc_param=NULL){
  classes = c("s", "h", "m", "n", "N", "alpha", "e",
              "r.0.D", "r.1.D", "r.0.R", "r.1.R", 
              "a.D", "b.D", "c.D", "a.R", "c.R",
              "selection", "dWij.dzij", "dWij.dzik",
              "advection0", "advection1", "diffusion", "mut.rate", 
              "k", "phi", "psi1", "psi2")
  
  model <- vector(mode="list", length=length(classes))
  names(model) = classes
  
  model$s = s
  model$h = h
  model$m = m
  model$n = n
  model$N = N
  model$e = e
  model$mut.rate = mu
  
  if (is.null(alpha) && e == 0 ){
    
    model$alpha = 1/N
    model$r.0.D = (1-m)^2/(2*N - (1-m)^2*(2*N-1))
    model$r.1.D = model$r.0.D
    model$r.0.R = 0.5 + 0.5*model$r.0.D
    model$r.1.R = 0.5/N + 0.5*model$r.0.D/N + (1-1/N)*model$r.1.D
    
    model$Neff = N*n/(1-model$r.1.D)
    
    quo = (2*N - (2*N-1)*(1-m)^2)*(2*N*N-(2*N-1)*(N-1)*(1-m)^3)
    
    model$a.D   = (1-m)^3*(N+(2*N-1)*(1-m)^2)/quo
    model$b.D   = 3*m*(1-m)^2*N*(2*N + (1-m)*(2-m)*(2*N-1))/quo
    model$c.D   = 2*m^2*N^2*(2*N*(3-2*m)+(3-m)*(1-m)^2*(2*N-1))/quo
    
    model$a.R = 1/model$N*model$r.0.D + (1-1/model$N)*model$a.D
    model$c.R = model$c.D*(1-1/N)
    model$b.D = 1 - model$a.D -  model$c.D 
    model$b.R = 1 - model$a.R -  model$c.R
    
  } else if (is.null(alpha)) {
    
    if(is.null(extc_param)||length(extc_param) < 4 ) {print("Error in extc_param"); return(NULL);}
    
    model$k    = extc_param[1]
    model$phi  = extc_param[2]
    model$psi1 = extc_param[3]
    model$psi2 = extc_param[4]
    
    model$alpha = 1/N
    model$r.0.D = ((1-model$e)*(1-m)^2 + model$e*model$phi)/(2*N - ((1-model$e)*(1-m)^2 + model$e*model$phi)*(2*N-1)) #(2*model$k*(1-m)^2*(1-model$e) + 2*N*model$e + (2*model$k-1)*model$e*model$phi) / (2*model$k*(2*N - (2*N-1)*(1-model$e)*(1-m)^2) - (2*model$k-1)*(2*N-1)*model$e*model$phi)
    model$r.1.D = model$r.0.D
    model$r.0.R = 0.5 + 0.5*model$r.0.D
    model$r.1.R = 0.5/N + 0.5*model$r.0.D/N + (1-1/N)*model$r.1.D
    
    model$Neff = (1-e)*n/( 2*model$r.1.R*(1-(1-e)^2*(1-m)^2) )
    
    #zeta = 4*model$k^2*(1-model$e)*(1-m)^3 + model$e*(2*model$k-1)*(2*model$k-2)*model$psi1
    #xi = model$e* ( 4*N^2 + 6*N*(2*model$k-1)*model$psi1*(1+(2*N-1)*model$r.1.D) )
    #model$a.D   = ((1 + 3*(2*N-1)*model$r.1.D)*zeta + xi) / (16*N^2*model$k^2  - zeta*(2*N-1)*(2*N-2) )
    #delta   = 4*N^2*(m^3+3*m^2*(1-m))+3*m*(1-m)^2*2*N*(2*N-1)*(1-model$r.1.D)
    #epsilon = 4*N^2*(1 - model$psi1 - model$psi2) + model$psi2*2*N*(2*N-1)*(1-model$r.1.D)
    #denom = 16*N^2*model$k^2 - 4*model$k^2*(1-model$e)*(1-m)^3*(2*N-1)*(2*N-2) - model$e*(2*model$k-1)*(2*model$k-2)*model$psi1*(2*N-1)*(2*N-2)
    #model$c.D   = (4*model$k^2*(1-model$e)*delta + model$e*(2*model$k-1)*(2*model$k-2)*epsilon) / denom
    
    term1 = (1-model$m)^3*(1-model$e)+model$e*(model$phi)^2
    model$a.D = (1+3*(2*model$N-1)*model$r.1.D)*term1/(4*model$N*(model$N-(1-0.5/model$N)*(model$N-1)*term1))
    
    model$a.R = 1/model$N*model$r.0.D + (1-1/model$N)*model$a.D
    model$c.R =  1 - model$a.R - (model$r.0.D - model$a.R) - 2*(model$r.1.R - model$a.R)
  
    model$c.D = model$c.R/(1-1/N)
    model$b.D   = 1 - model$a.D -  model$c.D 
  
  }

  
  model$dWij.dzij = 2
  model$selection = selection
  
  if (selection == 'soft') model$dWij.dzik = -2/(N-1)
  if (selection == 'hard') model$dWij.dzik = -2*(1-m)^2*(1-e)/(N-1)
  
  model$advection0  = 0.5*s*(model$dWij.dzij*(2*h*model$r.0.R+(1-2*h)*model$r.0.D) + (N-1)*model$dWij.dzik*(2*h*model$r.1.R+(1-2*h)*model$a.R))
  model$advection1  = 0.5*s*(1-2*h)*(model$dWij.dzij*(1-model$r.0.D)+(N-1)*model$dWij.dzik*(1-model$r.0.D-model$c.R))
  model$diffusion   = 0.5/model$Neff
  
  return(model)
}


### Calculates the Meantime to Segregation ####
segregation.meantime = function (model, res=0.001, plot=FALSE){
  p = seq(res,1,by=res)
  time_scale = 1 #/(2*model$Neff)
  
  if( model$advection1 == 0){
    G = exp(-2*model$advection0/model$diffusion*p) 
    C = 2/model$diffusion/(1-exp(-2*model$advection0/model$diffusion))
    u = (0.5*C*(1-G)/model$advection0*model$diffusion - (1/model$advection0)*p)/time_scale
  } else {
    k1 = model$advection0/model$diffusion
    k2 = model$advection1/model$diffusion
    G = res*cumsum(exp(2*k1*p+k2*p^2))
    E = -2/model$diffusion*res*cumsum(G*exp(-2*k1*p-k2*p^2)) 
    B = res*cumsum(exp(-2*k1*p-k2*p^2))
    C = -E[length(E)]/B[length(B)]
    u = (C*B + E)/time_scale
  }
  
  if(plot) plot(c(0,p), c(0,u), type='l', lwd=2, xlab="Initial frequency on population", ylab="Mean time to segregation")
  return(cbind(c(0,p), c(0,u)))
}

### Calculates the Fixation probability ####
fixation.probability = function (model, res=0.001, plot=FALSE){
  p = seq(res,1,by=res)
  G = exp(-2*model$advection0/model$diffusion*p - model$advection1/model$diffusion*(p^2) )
  u = cumsum(G)/sum(G)
  
  if(plot) plot(p, u, type='l', lwd=2, xlab="Initial frequency on population", ylab="Probability of fixation")
  return(cbind(c(0,p), c(0,u)))
}

### Calculates the meantime to Fixation ####
fixation.meantime  = function (model, res=0.001, plot=FALSE){
  p = seq(res,1-res,by=res)

  PSI = function(p) exp(-2*model$advection0/model$diffusion*p - model$advection1/model$diffusion*p^2)
  
  norm = integrate(Vectorize(PSI), 0, 1)$value
  
  t0 = function(q) {
    2.0/(model$diffusion*q*(1-q)*PSI(q))/norm*(integrate(Vectorize(PSI), 0, q)$value)^2
  }
  
  t1 = function(q) {
    2.0/(model$diffusion*q*(1-q)*PSI(q))/norm*integrate(Vectorize(PSI), 0, q)$value*integrate(Vectorize(PSI), q, 1)$value
  }
  
  aux = function(x0){
    t.star0 = integrate(Vectorize(t0), lower=0, upper=x0)$value*integrate(Vectorize(PSI), x0, 1)$value/integrate(Vectorize(PSI), 0, x0)$value
    t.star1 = integrate(Vectorize(t1), lower=x0, upper=1, rel.tol = 1/1000)$value
    return(t.star0 + t.star1)
  }
  
  t.star = sapply(p, aux)
  
  return(cbind(p, t.star))

  if(plot) plot(c(p,1), c(t.star,0), type='l', lwd=2, xlab="Initial frequency on population", ylab="Mean time to fixation")
}

fixation.meantime.points  = function (model, p0){
  PSI = function(p) exp(-2*model$advection0/model$diffusion*p - model$advection1/model$diffusion*p^2)
  
  norm = integrate(Vectorize(PSI), 0, 1)$value
  
  t0 = function(q) {
    2.0/(model$diffusion*q*(1-q)*PSI(q))/norm*(integrate(Vectorize(PSI), 0, q)$value)^2
  }
  
  t1 = function(q) {
    2.0/(model$diffusion*q*(1-q)*PSI(q))/norm*integrate(Vectorize(PSI), 0, q)$value*integrate(Vectorize(PSI), q, 1)$value
  }
  
  aux = function(x0){
    t.star0 = integrate(Vectorize(t0), lower=0, upper=x0)$value*integrate(Vectorize(PSI), x0, 1)$value/integrate(Vectorize(PSI), 0, x0)$value
    t.star1 = integrate(Vectorize(t1), lower=x0, upper=1, rel.tol = 1/1000)$value
    return(t.star0 + t.star1)
  }
  
  t.star = sapply(p0, aux)
  
  return(cbind(p0, t.star))
  
}

fixation.meantime.singlemutant  = function (model){
  x0 = 0.5/model$n/model$N
  epsilon = 0.00001
  
  f = function(p) model$mut.rate*(1-2*p) + p*(1-p)*(model$advection0 + p*model$advection1)
  g = function(p) model$diffusion*p*(1-p)
  integand = function (p){
    if (p==0 && model$mut.rate == 0) return(-2 * model$advection0 / model$diffusion)
    -2*f(p)/g(p)
  }
  integand = Vectorize(integand)
  
  G = function(p) exp(integrate(integand, lower=epsilon, upper=p-epsilon, subdivisions = 100000L)$value)
  G = Vectorize(G)
  norm = integrate(Vectorize(G), epsilon, 1-epsilon)$value

  t = function(q) {
    2.0/(model$diffusion*q*(1-q)*G(q))/norm*integrate(Vectorize(G), epsilon, q, subdivisions = 100000L)$value*integrate(Vectorize(G), q, 1-epsilon, subdivisions = 100000L)$value
  }
  
  t.star = integrate(Vectorize(t), lower=x0+epsilon, upper = 1-epsilon, subdivisions = 100000L)$value
  
  return(t.star)
}

meantime.singlemutant  = function (model, threshold = 1){
  x0 = 0.5/model$n/model$N
  
  G = function(p) exp(-2*(model$advection0/model$diffusion)*p-(model$advection1/model$diffusion)*p^2)
  norm = integrate(G, 0, 1)$value
  
  t = function(q) {
    2.0/(model$diffusion*q*(1-q)*G(q))/norm*integrate(G, 0, q)$value*integrate(G, q, 1)$value
  }
  
  t.star = integrate(Vectorize(t), lower=x0, upper=threshold)$value
  
  return(t.star)
}

adaptation.meantime.singlemutant  = function (model){
  u = fixation.probability(model, res=0.5/(model$N*model$n))[2,2]
  Tadapt = 0.5/(model$N*model$n*model$mut.rate*u*(1-model$e))
  model$mut.rate=0
  Tadapt = Tadapt + fixation.meantime.singlemutant(model)
  return(Tadapt)
}

adaptation.meantime.softsweep = function(s, h, m, n, N, mu, size = 1E3, selection = "soft", e = 0, extc_param=NULL){
  res0 = 0.5/n/N
  model = Roze(s=s[1], m=m, h=h[1], n=n, N=N, mu=mu, selection=selection, e = e, extc_param=extc_param)
  probdist = stationary.FPEq(model, res0)
  
  samples = sample(x = probdist[,1], prob = probdist[,2]*res0, size = size, replace=T)

  model = Roze(s=s[2], m=m, h=h[2], n=n, N=N, mu=0, selection=selection, e = e, extc_param=extc_param)
  time = fixation.meantime.points(model, p0 = samples)
  
  mean(time[,2])
}

adaptation.fixation.softsweep = function(s, h, m, n, N, mu, selection="soft"){
  res0 = 0.5/n/N
  model = Roze(s=s[1], m=m, h=h[1], n=n, N=N, mu=mu, selection="soft")
  probdist = stationary.FPEq(model, res0)

  model = Roze(s=s[2], m=m, h=h[2], n=n, N=N, mu=0, selection="soft")
  probfix = fixation.probability(model, res0)
  
  return( res0*sum(probdist[,2] * probfix[seq(2, length(probfix[,2])-1),2] ) )
}

### Solves the probability distribution in time numerically
solve.FPEq = function(model, ini.cond, tf, res=0.001, plot=FALSE){
  x = seq(0,1,by=res)
  u = rep(0,length(x))
  #p0 = round(ini.cond/res)+1
  #u[p0]=0.5/res;u[p0-1]=0.25/res;u[p0+1]=0.25/res
  u = dnorm(x, mean = ini.cond[1], sd = ini.cond[2])
    
  data = data.frame(time = rep(0, length(u)), 
                    p = x,
                    u = u)
  
  f = function(p) model$mut.rate*(1-2*p) + p*(1-p)*(model$advection0 + p*model$advection1)
  g = function(p) model$diffusion*p*(1-p)
      
  d  = function(x, u) f(x+res)*u[c(2:length(x), length(x))] - f(x-res)*u[c(1, 1:length(x)-1)]
  
  dd = function(x, u) g(x+res)*u[c(2:length(x), length(x))] + g(x-res)*u[c(1, 1:length(x)-1)] - 2*g(x)*u
    
  dt = 0.05*min(c(abs(res/f(0.5)), abs(res^2/g(0.5))))
  ts = seq(0, max(tf)+dt,by=dt)
  
  for (t in ts) {
    k1 =  -0.5*d(x, u)/res + 0.5*dd(x, u)/res^2         
    k2 =  -0.5*d(x, u+k1*dt)/res + 0.5*dd(x, u+k1*dt)/res^2 
    u = u + 0.5*(k1+k2)*dt
    
    # Boundary conditions -- zero flux
    u[1] = (18*u[2]-9*u[3]+2*u[4])/11
    u[length(u)] = (18*u[length(u)-1]-9*u[length(u)-2]+2*u[length(u)-3])/11 
    
    if (any( abs(t - tf) < (dt/2) )){
      rows = data.frame(time = rep(round(t,0), length(u)), 
                        p = x,
                        u = u)
      data = bind_rows(data, rows)  }
  }
  
  if (plot) {
    plot = ggplot(data, aes(x = p, y = u, group=time)) + 
    geom_path(col = "red", alpha = 0.75 )+ 
    guides(colour=TRUE) + xlim(c(0,1)) +
    facet_grid(rows = vars(time)) + 
    labs(title = "Densities") +   theme(plot.title = element_text(hjust = 0.5)) + 
    xlab("Probabilities") +
    ylab("Probability distribution")
  plot
  }
  
  return(data) 
}

### Solves the probability distribution in time numerically
stationary.FPEq = function(model, res=0.001, plot=FALSE){
  epsilon = 1E-5
  x = seq(res,1-res,by=res)

  f = function(p) model$mut.rate*(1-2*p) + p*(1-p)*(model$advection0 + p*model$advection1)
  g = function(p) model$diffusion*p*(1-p)
  
  integand = function (p) 2*f(p)/g(p)
  
  aux = function(x) integrate(integand, lower= epsilon, upper= max(epsilon, min(1-epsilon, x)))$value
  G = sapply(x, aux)
  
  integand2 = function(p) exp(aux(p))/(model$diffusion*p*(1-p))
  norm =  integrate(integand2, lower= epsilon, upper= 1-epsilon, subdivisions = 100000L)$value #sum(exp(G)/(model$diffusion*x*(1-x))) #
  
  u = exp(G)/(model$diffusion*x*(1-x))/norm
  norm2 = res*sum(u)
  
  data = cbind(x, u/norm2)
    
  return(data) 
}

### Calculates the mean time to fixation using the approximation of Charlesworth 2020
fixation.meantime.Charlesworth2020 = function(model){
  q1 = 1 / (4 * model$advection0 * model$Neff)
  p2 = 1 / (4 * (model$advection0+model$advection1) * model$Neff)
  
  if (model$advection0 < 1E-6) { Td = 1/model$advection1*(1/q1 + log((1.0-q1)*(1.0-p2)/q1/p2) - 1.0/(1.0-p2) ) }
  if (model$advection0+model$advection1 < 1E-6) { Td = -1/model$advection1*(1/p2 + log((1.0-q1)*(1.0-p2)/q1/p2) - 1.0/(1.0-q1) ) }
  else Td = (1.0/model$advection0*log(model$advection0/(q1*(model$advection0 + model$advection1))) + 
                             1.0/(model$advection0 + model$advection1)*log((model$advection0 + model$advection1)/(p2*model$advection0))   )
  
  Ts1 = 4*model$Neff*( (1-q1)/q1*log(1-q1) + 1 )
  Ts2 = 2*model$Neff*(1-exp(-1))/(2*(model$advection0 + model$advection1)*model$Neff)
  return(Ts1 + Td + Ts2)
}

adaptation.meantime.Charlesworth2020 = function(model){
  u = fixation.probability(model, res=0.5/(model$N*model$n))[2,2]
  Tadapt = 0.5/(model$N*model$n*model$mut.rate*u*(1-model$e))
  Tadapt = Tadapt + fixation.meantime.Charlesworth2020(model)
  return(Tadapt)
}

### Calculates the mean time to fixation using the approximation of van Herwaarden and van der Wal (2002)
fixation.meantime.Herwaarden2002 = function(model){
  q1 = 1 / (4 * model$advection0 * model$Neff)
  p2 = 1 / (4 * (model$advection0+model$advection1) * model$Neff)
  Td = (1.0/model$advection0*log(model$advection0/(q1*(model$advection0 + model$advection1))) + 
          1.0/(model$advection0 + model$advection1)*log((model$advection0 + model$advection1)/(p2*model$advection0))   )
  
  Ts1 = 4*model$Neff*( (1-q1)/q1*log(1-q1) + 1 )
  Ts2 = 2*model$Neff*(1-exp(-1))*(2*(model$advection0 + model$advection1)* model$Neff)
  return(Ts1 + Td + Ts1)
}
