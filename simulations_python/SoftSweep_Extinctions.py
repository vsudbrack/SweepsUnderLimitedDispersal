from pickle import TRUE
import numpy as np
import numpy.random as rd
import sys

class Subpopulation:
    def __init__(self, nind):
        self.nind = np.array([nind, 0, 0])
        self.gam = np.array([nind, 0.0])
        self.migs = np.array([0.0, 0.0])

class Population:
    def __init__(self, ndemes, nind, s, h, m, mu, e, modeext):
        self.ndemes = ndemes
        self.nind = nind
        self.selcoeff = s
        self.gencoef = h
        self.migrationrate = m
        self.mutationrate = mu
        self.extincitonrate = e
        self.extincitonmode = modeext

        self.totgam = np.array([ndemes*nind, 0.0])
        self.subpops = [Subpopulation(nind) for deme in range(ndemes)]

    def gametes(self):
        for deme in self.subpops:
            deme.gam = [deme.nind[0] + 0.5*(1.0+self.gencoef*self.selcoeff)*deme.nind[1], \
                0.5*(1+self.gencoef*self.selcoeff)*deme.nind[1] + (1.0+self.selcoeff)*deme.nind[2]]
            # normalization
            deme.gam = deme.gam / np.sum(deme.gam)

    def mutate(self):
       if np.sign(self.selcoeff) * self.mean_freq() < 0.01: #If it is deleterious, it is alwyas true
            for deme in self.subpops: # mutation
                deme.gam = (1.0-self.mutationrate)*deme.gam + self.mutationrate*(1.0 - deme.gam)

    def migrate(self):
        self.totgam[0] = np.sum([deme.gam[0] for deme in self.subpops])
        self.totgam[1] = np.sum([deme.gam[1] for deme in self.subpops])

        for deme in self.subpops:
            deme.migs[0] = (self.totgam[0] - deme.gam[0])/(1.0*self.ndemes-1.0)
            deme.migs[1] = (self.totgam[1] - deme.gam[1])/(1.0*self.ndemes-1.0)

        for deme in self.subpops:
            deme.gam = (1.-self.migrationrate) * deme.gam + self.migrationrate * deme.migs
            deme.gam =  deme.gam / np.sum(deme.gam)

    def fuse(self):
        for deme in self.subpops:
            probs = [deme.gam[0]**2, 2.0*deme.gam[0]*deme.gam[1], deme.gam[1]**2]
            deme.nind = rd.multinomial(self.nind, probs)

    def extinct(self):
        nro_ext_demes = rd.binomial(self.ndemes, self.extincitonrate)

        if (nro_ext_demes == 0): return(0)

        ext_demes = rd.choice(range(1, self.ndemes), size = nro_ext_demes, replace = False)

        if(self.extincitonmode == "propagule"):
            for deme_nro in ext_demes:
                demeExt = self.subpops[deme_nro]
                demeCol_nro = deme_nro

                while(demeCol_nro in ext_demes):
                    demeCol_nro = rd.choice(range(1, self.ndemes), size = 1)

                demeCol = self.subpops[demeCol_nro[0]]

                probs = [demeCol.gam[0]**2, 2.0*demeCol.gam[0]*demeCol.gam[1], demeCol.gam[1]**2]
                demeExt.nind = rd.multinomial(self.nind, probs)

        elif(self.extincitonmode == "pool"):
            for demeExt_nro in ext_demes:
                demeExt = self.subpops[demeExt_nro]
                self.totgam[0] = self.totgam[0] - demeExt.gam[0] 
                self.totgam[1] = self.totgam[1] - demeExt.gam[1]

            mig_pool0 = self.totgam[0] / (self.totgam[0] + self.totgam[1])
            mig_pool1 = self.totgam[1] / (self.totgam[0] + self.totgam[1])

            for deme_nro in ext_demes:
                demeExt = self.subpops[deme_nro]
                probs = [mig_pool0**2, 2.0*mig_pool0*mig_pool1, mig_pool1**2]
                demeExt.nind = rd.multinomial(self.nind, probs)

        else:
            print("Error in mode of extinction"); exit()

    def next_generation(self):
            self.gametes()
            self.mutate()
            self.migrate()
            self.fuse()
            self.extinct()

    def mean_freq(self):
            sum1 = np.sum([deme.nind[2] for deme in self.subpops])
            sum2 = np.sum([deme.nind[1] for deme in self.subpops])
            return( (2.0*sum1+sum2)/(2.0 * self.ndemes * self.nind) )

    def env_change(self, snew, hnew):
        self.selcoeff = snew
        self.gencoef = hnew

    def fixed(self):
        if ( self.totgam[0] == 0): return(TRUE)

    def param_list(self):
        return("%d %d %.2e %.2e %.2e %.2e %.2e %s"%(self.ndemes, self.nind, self.selcoeff, self.gencoef, \
        self.migrationrate, self.mutationrate, self.extincitonrate, self.extincitonmode) )

    def __repr__(self):
        return( "Mean freq: %.2f %%"%(100*self.mean_freq()) ) 


# Creates a population
Ndemes  = 200
Ninds   = 100
selCoef = -0.001
domCoef = float(sys.argv[1])
disRate = float(sys.argv[2])
mutRate = 1E-4
extRate = float(sys.argv[3])
extMode = sys.argv[4]
idNumber = int(sys.argv[5])

main_pop = Population(Ndemes, Ninds, selCoef , domCoef, disRate, mutRate, extRate, extMode)

for gen in range(10**6):
    main_pop.next_generation()

    if (gen % 200 == 0) :
        print("%s %d %d %.3e"%( main_pop.param_list(), idNumber, gen, main_pop.mean_freq() ) )

    if (gen == 10**4) :
        main_pop.env_change(0.01, domCoef)

    if (main_pop.fixed()) :
        exit()

