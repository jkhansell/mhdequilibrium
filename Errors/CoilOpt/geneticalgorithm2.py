import numpy as np

def probfunc(Imutprob, Nbred): 
    return Imutprob*np.exp(-Nbred)

class GenAlg(object):

    def __init__(self, costfunc, I, ngenes, kfunc, Rred, nmut, mutprobfunc = probfunc, Ipopsize = 30, Fpopsize = 10, dpopsize = 5,
                PXcross = 0.85, Imutprob = 0.9, Rabs=0.0001):

        self.costfunc = costfunc

        self.I = I
        self.Ipopsize = Ipopsize
        self.Fpopsize = Fpopsize
        self.dpopsize = dpopsize
        self.ngenes = ngenes
        self.PXcross = PXcross
        self.kfunc = kfunc
        self.Imutprob = Imutprob
        self.mutprobfunc = mutprobfunc
        self.Rred = Rred
        self.Rabs = Rabs
        self.maxiter = ngenes*5
        self.delta = self.I[1] - self.I[0] 
        self.Rneigh = self.npop*self.ngenes
        self.epsilon = self.delta/self.Rneigh
        self.maxgen = 2*self.ngenes
        self.nmut = nmut


    def normalize_params(self, p):
        phi = self.I[1]
        plo = self.I[0]
        pnorm = (p-plo)/(phi-plo)
        return pnorm

    def denormalize_params(self, pnorm):
        phi = self.I[1]
        plo = self.I[0]
        p = (phi-plo)*pnorm + plo
        return p

    def pop_generation(self):

        npop = self.Ipopsize
        epsilon = self.epsilon
        I = self.I
        ngenes = self.ngenes
        normeps = self.normalize_params(epsilon, I)

        pop = np.zeros((npop, ngenes))
        pop[0] = np.random.rand(ngenes)
        i = 1
        while i < npop: 
            poptemp = np.random.rand(ngenes)
            truthvect = np.zeros((i+1))
            for k in range(i+1):
                ds = np.sqrt(np.sum((poptemp-pop[k])**2))
                if ds >= normeps: 
                    truthvect[k] = True
            
            if truthvect.all(): 
                pop[i] = poptemp 
                i += 1 
        return pop

    def cost_eval(self, pop):
        costfunc = self.costfunc
        I = self.I
        
        npop = len(pop)
        poptemp = self.denormalize_params(pop, I)
        costvect = np.zeros(npop)
        for i in range(npop):
            costvect[i] = costfunc(poptemp[i])
        return costvect

    def selection(self, pop, nparents):
        costvect = self.cost_eval(self.costfunc, pop, self.I)**2
        total = np.sum(costvect)
        pievect = costvect/total
        pievect = np.sort(pievect)
        inds = np.argsort(pievect)
        parents = []
        for j in range(nparents):
            randnum = np.random.rand()
            p = 0
            i = 0
            while p < randnum: 
                p += pievect[i]
                i += 1

            parents.append(pop[inds[i]])
        
        return np.array(parents)

    def recombination(self, pop, parents):

        nchildren = len(pop)-len(parents)
        newpop = parents.copy()

        for i in range(nchildren): 
            ind1 = np.random.randint(0, len(parents))
            ind2 = np.random.randint(0, len(parents))

            child1 = parents[ind1]
            child2 = parents[ind2]

            indchange = np.random.randint(0, len(child1))

            child1[indchange:-1] = child2[indchange:-1]
            child2[indchange:-1] = child1[indchange:-1]

            M = np.random.randint(1,1001)
        
            dx1 = child1[indchange]/M
            dx2 = child2[indchange]/M

            child1[indchange] += dx2 - dx1
            child2[indchange] += dx1 - dx2
    
            if np.random.rand() > 0.5: 
                newpop = np.append(newpop, [child1], axis = 0)
            else: 
                newpop = np.append(newpop, [child2], axis = 0)

        return newpop 

    def mutation(self, ite, pop, s, mutprob):
        nPop = len(pop)
        probs = np.random.rand(nPop)

        for i in range(nPop): 

            if probs[i] > mutprob: 
                randind = np.random.randint(0,len(pop[i]), size=s)
                M = np.random.randint(1,11, size=s)
                dx = self.delta/M
                sign = np.random.randint(-1,1, size=s)
                pop[i, randind] += sign*dx*self.kfunc(ite)
        
        return pop

    def GA(self):
        pop = self.pop_generation()     #Initialization
        iters = 0

        costs = cost_eval()
        while iters < self.maxiter and 




        


