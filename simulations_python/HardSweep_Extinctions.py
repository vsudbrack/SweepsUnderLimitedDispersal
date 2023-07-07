from pickle import TRUE, FALSE 
import numpy as np
import numpy.random as rd

import sys
np.seterr("ignore") # Ignores 0/0 in extinct demes (they will be recolonized later)
import warnings
warnings.filterwarnings('ignore')


class Subpopulation:
    def __init__(self, nind):
        self.nind = np.array([nind, 0, 0])
        self.gam = np.array([nind, 0.0])
        self.migs = np.array([0.0, 0.0])

class Population:
    def __init__(self, ndemes, nind, s, h, m, mu, e, extmodel):
        self.ndemes = ndemes
        self.nind = nind
        self.selcoeff = s
        self.gencoef = h
        self.migrationrate = m
        self.mutationrate = mu
        self.extincitonrate = e
        self.propagule = TRUE if extmodel == "propagule" else FALSE 

        self.totgam = np.array([ndemes*nind, 0.0])
        self.subpops = [Subpopulation(nind) for deme in range(ndemes)]

    def gametes(self):
        for deme in self.subpops:
            deme.gam = [deme.nind[0] + 0.5*(1.0+self.gencoef*self.selcoeff)*deme.nind[1], \
                0.5*(1+self.gencoef*self.selcoeff)*deme.nind[1] + (1.0+self.selcoeff)*deme.nind[2]]
            # normalization
            deme.gam = deme.gam / np.sum(deme.gam)

    def extinct(self):
        self.nro_ext_demes = rd.binomial(self.ndemes, self.extincitonrate)
        if (self.nro_ext_demes == 0): 
            self.ext_demes = []
            return(0)

        self.ext_demes = rd.choice(range(1, self.ndemes), size = self.nro_ext_demes, replace = False)

        for demeExt_nro in self.ext_demes:
            demeExt = self.subpops[demeExt_nro]
            demeExt.gam[0] = np.nan 
            demeExt.gam[1] = np.nan

    def mutate(self):
        for deme in self.subpops: # mutation
            deme.gam = (1.0-self.mutationrate)*deme.gam + self.mutationrate*(1.0 - deme.gam)

    def migrate(self):
        self.totgam[0] = np.nansum([deme.gam[0] for deme in self.subpops])
        self.totgam[1] = np.nansum([deme.gam[1] for deme in self.subpops])

        for deme in self.subpops:
            deme.migs[0] = (self.totgam[0] - deme.gam[0])/(self.ndemes - self.nro_ext_demes - 1.0)
            deme.migs[1] = (self.totgam[1] - deme.gam[1])/(self.ndemes - self.nro_ext_demes - 1.0)

        for deme in self.subpops:
            deme.gam = (1.-self.migrationrate) * deme.gam + self.migrationrate * deme.migs
            deme.gam =  np.divide(deme.gam, np.nansum(deme.gam)) # avoids 0/0 problems in extinct demes

    def fuse(self):
        for deme in self.subpops:
            if( np.isnan(deme.gam[0]) ): continue # skip extinct demes
            probs = [deme.gam[0]**2, 2.0*deme.gam[0]*deme.gam[1], deme.gam[1]**2]
            deme.nind = rd.multinomial(self.nind, probs)

    def recolonize_migrantpool(self):
        mig_pool0 = self.totgam[0] / (self.totgam[0] + self.totgam[1])
        mig_pool1 = self.totgam[1] / (self.totgam[0] + self.totgam[1])

        for deme_nro in self.ext_demes:
            demeExt = self.subpops[deme_nro]
            probs = [mig_pool0**2, 2.0*mig_pool0*mig_pool1, mig_pool1**2]
            demeExt.nind = rd.multinomial(self.nind, probs)

    def recolonize_propagule(self):
        for deme_nro in self.ext_demes:
            demeExt = self.subpops[deme_nro]
            demeCol_nro = rd.choice(range(1, self.ndemes), size = 1)
            while (not np.isnan(self.subpops[demeCol_nro[0]].gam[0])): 
                demeCol_nro = rd.choice(range(1, self.ndemes), size = 1)
            
            demeCol = self.subpops[demeCol_nro[0]]
            probs = [demeCol.gam[0]**2, 2.0*demeCol.gam[0]*demeCol.gam[1], demeCol.gam[1]**2]
            demeExt.nind = rd.multinomial(self.nind, probs)

    def next_generation(self):
            self.gametes()
            self.extinct()
            self.mutate()
            self.migrate()
            self.fuse()
            self.recolonize_propagule() if self.propagule else self.recolonize_migrantpool()

    def mean_freq(self):
            sum1 = np.sum([deme.nind[2] for deme in self.subpops])
            sum2 = np.sum([deme.nind[1] for deme in self.subpops])
            return( (2.0*sum1+sum2)/(2.0 * self.ndemes * self.nind) )

    def fixed(self):
        if ( self.mean_freq() > 1.0 - 0.5/(self.ndemes * self.nind) ): return(TRUE)

    def extincted(self):
        if ( self.mean_freq() < 0.5/(self.ndemes * self.nind) ): return(TRUE)

    def param_list(self):
        return("%s %d %d %.2e %.2e %.2e %.2e %.2e"%("PROP" if self.propagule else "POOL", self.ndemes, self.nind, self.selcoeff, self.gencoef, self.migrationrate, self.mutationrate, self.extincitonrate) )

    def __repr__(self):
        return( "Mean freq: %.2f %%"%(100*self.mean_freq()) ) 


if ( len(sys.argv) != 9): 
    print("Error in arguments")
    exit()

# Creates a population
Ndemes  = int(sys.argv[1])
Ninds   = int(sys.argv[2])
selCoef = float(sys.argv[3])
domCoef = float(sys.argv[4])
disRate = float(sys.argv[5])
mutRate = float(sys.argv[6])
extRate = float(sys.argv[7])
extModel = sys.argv[8]

main_pop = Population(Ndemes, Ninds, selCoef, domCoef, disRate, mutRate, extRate, extModel)
main_pop.subpops[0].nind = np.array([main_pop.nind-1, 1, 0])

for gen in range(10**6):
    main_pop.next_generation()

    if (main_pop.fixed() ) :
        print("%s %s %d"%(main_pop.param_list(), "Fix", gen))
        exit()

    if (main_pop.extincted() ) :
        print("%s %s %d"%(main_pop.param_list(), "Ext", gen))
        exit()

    #if (gen % 10 == 0) :
    #    print("%s %f %d"%(main_pop.param_list(), main_pop.mean_freq(), gen))


