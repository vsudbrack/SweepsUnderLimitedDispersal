import numpy as np
import numpy.random as rd
import sys

class Subpopulation:
    def __init__(self, nind, inifreq):

        BBind = np.round(inifreq*inifreq*nind)[0]
        bbind = np.round((1.0-inifreq)*(1.0-inifreq)*nind)[0]
        Bbind = nind - BBind - bbind

        self.nind = np.matrix([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, BBind, Bbind], [0, 0, 0, bbind]])
        
        self.gam =  np.array([0.0, 0.0, 2*BBind + Bbind, 2*bbind + Bbind])
        self.migs = np.array([0.0,  0.0, 0.0, 0.0])

class Population:
    def __init__(self, ndemes, nind, s, h, m, mu, e, c):
        self.ndemes = ndemes
        self.nind = nind
        self.selcoeff = s
        self.gencoef = h
        self.migrationrate = m
        self.mutationrate = mu
        self.extincitonrate = e
        self.recombinationfraction = c

        self.totgam = np.array([0.0, 0.0, 0.5, 0.5])
        self.subpops = [Subpopulation(nind, rd.beta(2*nind*m, 2*nind*m, size=1)) for deme in range(ndemes)]

    def gametes(self):
        c = self.recombinationfraction
        s = self.selcoeff
        hs = self.gencoef*self.selcoeff

        for deme in self.subpops:
            ht1 = deme.nind[0,0]
            ht2 = deme.nind[0,1] + deme.nind[1,0]
            ht3 = deme.nind[0,2] + deme.nind[2,0]
            ht4 = deme.nind[0,3] + deme.nind[3,0] 
            ht5 = deme.nind[1,1] 
            ht6 = deme.nind[1,2] + deme.nind[2,1] 
            ht7 = deme.nind[3,1] + deme.nind[1,3] 
            ht8 = deme.nind[2,2]  
            ht9 = deme.nind[2,3] + deme.nind[3,2] 
            ht10 = deme.nind[3,3]
            deme.gam = [(1+s)*ht1 + 0.5*(1+s)*ht2 + 0.5*(1+hs)*ht3 + 0.5*(1-c)*(1+hs)*ht4 + 0.5*c*(1+hs)*ht6    ,\
                        (1+s)*ht5 + 0.5*(1+s)*ht2 + 0.5*(1+hs)*ht7 + 0.5*c*(1+hs)*ht4     + 0.5*(1-c)*(1+hs)*ht6,\
                        1.0*ht8   + 0.5*ht9       + 0.5*(1+hs)*ht3 + 0.5*c*(1+hs)*ht4     + 0.5*(1-c)*(1+hs)*ht6,\
                        1.0*ht10  + 0.5*ht9       + 0.5*(1+hs)*ht7 + 0.5*(1-c)*(1+hs)*ht4 + 0.5*c*(1+hs)*ht6     ]
            # normalization
            deme.gam = deme.gam / np.sum(deme.gam)

    def migrate(self):
        self.totgam[0] = np.sum([deme.gam[0] for deme in self.subpops])
        self.totgam[1] = np.sum([deme.gam[1] for deme in self.subpops])
        self.totgam[2] = np.sum([deme.gam[2] for deme in self.subpops])
        self.totgam[3] = np.sum([deme.gam[3] for deme in self.subpops])

        for deme in self.subpops:
            deme.migs[0] = (self.totgam[0] - deme.gam[0])/(1.0*self.ndemes-1.0)
            deme.migs[1] = (self.totgam[1] - deme.gam[1])/(1.0*self.ndemes-1.0)
            deme.migs[2] = (self.totgam[2] - deme.gam[2])/(1.0*self.ndemes-1.0)
            deme.migs[3] = (self.totgam[3] - deme.gam[3])/(1.0*self.ndemes-1.0)
            
        for deme in self.subpops:
            deme.gam = (1.0-self.migrationrate) * deme.gam + self.migrationrate * deme.migs
            deme.gam =  deme.gam / np.sum(deme.gam)
            if( (deme.gam<0).any() ): print(deme.gam)


    def fuse(self):
        for deme in self.subpops:
            probs = np.outer(deme.gam, deme.gam).flatten()
            deme.nind = rd.multinomial(self.nind, probs, 1).reshape(4, 4, order='F')

    def next_generation(self):
            self.gametes()
            self.migrate()
            self.fuse()

    def mean_freq(self):
        sum1 = np.sum([deme.nind[0,0]+deme.nind[0,1]+deme.nind[1,0]+deme.nind[1,1] for deme in self.subpops])
        sum2 = np.sum([deme.nind[2,0]+deme.nind[0,2]+deme.nind[3,0]+deme.nind[3,0]+deme.nind[2,1]+deme.nind[1,2]+deme.nind[3,1]+deme.nind[3,1] for deme in self.subpops])
        return( (2.0*sum1+sum2)/(2.0 * self.ndemes * self.nind) )

    def ld_global(self):
        tot = self.totgam[0]+self.totgam[1]+self.totgam[2]+self.totgam[3]
        ld = self.totgam[0]/tot - (self.totgam[0]+self.totgam[1])/tot*(self.totgam[0]+self.totgam[2])/tot
        return(ld)
    
    def freqB(self):
        tot = self.totgam[0]+self.totgam[1]+self.totgam[2]+self.totgam[3]
        totB = self.totgam[0] + self.totgam[2]
        return(totB/tot)
    
    def FST(self):
        freqs = [i.gam[0]+ i.gam[2]/(i.gam[0]+ i.gam[1]+i.gam[2]+ i.gam[3]) for i in self.subpops]
        var = np.var(freqs)
        mu = np.mean(freqs)
        if mu == 0 or mu == 1: return(0)
        return( var/mu/(1-mu) )
    
    def fixed(self):
        if ( self.totgam[2]+self.totgam[3] == 0): return(True)

    def lost(self):
        if ( self.totgam[0]+self.totgam[1] == 0): return(True)

    def param_list(self):
        return("%d %d %.2f %.2f %.2e %.2e"%(self.ndemes, self.nind, self.selcoeff, self.gencoef, self.migrationrate, self.recombinationfraction) )

    def __repr__(self):
        return( "Mean freq: %.2f %%"%(100*self.mean_freq()) ) 


if ( len(sys.argv) != 8): 
    print("Error in arguments")
    exit()

# Creates a population
Ndemes  = int(sys.argv[1])
Ninds   = int(sys.argv[2])
selCoef = float(sys.argv[3])
domCoef = float(sys.argv[4])
disRate = float(sys.argv[5])
mutRate = 0.0
recombinationRate = float(sys.argv[6])
idNumber = int(sys.argv[7])

main_pop = Population(Ndemes, Ninds, selCoef, domCoef, disRate, mutRate, 0, recombinationRate)

gen0 = main_pop.subpops[0].nind
gen1 = gen0 + np.matrix([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 1, 0, -1]]) #First mutant: ab/ab -> Ab/ab
gen2 = gen0 + np.matrix([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, -1, 0], [1, 0, 0, 0]]) #First mutant: aB/aB -> AB/aB

if (gen1 >= 0).all() : 
    main_pop.subpops[0].nind = gen1 
elif (gen2 >= 0).all() : 
    main_pop.subpops[0].nind = gen2
else: 
    print("ALL IC HAVE NEGATIVE VALUES")
    exit()

for gen in range(10**6):
    main_pop.next_generation()

    if ( main_pop.fixed() ):
        print("%s %d F %d %f %f"%( main_pop.param_list(), idNumber, gen, main_pop.FST(), main_pop.freqB() ) )
        exit()

    if ( main_pop.lost() ):
        print("%s %d L %d %f %f"%( main_pop.param_list(), idNumber, gen, main_pop.FST(), main_pop.freqB()) )
        exit()

print("%s %d S NA"%( main_pop.param_list(), idNumber, gen) )
