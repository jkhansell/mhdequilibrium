"""
Rutina de código que implementa una versión del algoritmo genético continuo.

Elaborado por: Johansell Villalobos 
Ing. Física, TEC

"""

from geneticalgorithm import denormalize_params
import numpy as np

class GenAlg(object): 

    def __init__(self, costfunc, I, npop, ngenes, 
                    ncross, nmut, epsilon, nparents, Rred, maxgen, 
                    gentol, maxiter = 5000, abstol = 1e-20, 
                    mutprob = 0.2): 

        self.costfunc = costfunc
        self.I = I
        self.staticI = I.copy() 
        self.npop = npop
        self.ngenes = ngenes
        self.ncross = ncross
        self.nmut = nmut
        self.mutprob = mutprob
        self.epsilon = epsilon
        self.nparents = nparents
        self.gentol = gentol
        self.maxiter = maxiter
        self.Rred = Rred
        self.abstol = abstol 
        self.maxgen = maxgen 


    def normalize_params(self, p):
        """
        Función que normaliza los genes respecto al espacio de solución.
        p :: gen/dato
        I :: intervalo que determina el espacio de solución 
        """
        phi = self.I[:,1] 
        plo = self.I[:,0]
        pnorm = (p-plo)/(phi-plo)
        return pnorm


    def denormalize_params(self, pnorm):
        """
        Función que devuelve los genes a sus valores correspondientes en el espacio de solución.
        pnorm :: gen/dato normalizado 
        I :: intervalo que determina el espacio de solución 
        """
        phi = self.I[:,1] 
        plo = self.I[:,0]
        p = (phi-plo)*pnorm + plo
        return p


    def pop_generation(self):
        """
        Función generadora de una población esparcida por el espacio de solución. 

        npop :: cantidad de individuos en la población
        ngenes :: cantidad de variables a optimizar 
        epsilon :: radio de distanciamiento entre individuos para un debido esparcimiento (normalizado).  
        """
        npop = self.npop
        ngenes = self.ngenes 
        epsilon = self.epsilon

        pop = np.zeros((npop, ngenes))
        pop[0] = np.random.rand(ngenes)
        i = 1
        j = 0
        while i < npop: #se comparan los individuos para corrobar su debido esparcimiento
            poptemp = np.random.rand(ngenes)
            truthvect = np.zeros((i+1))
            for k in range(i+1):
                ds = np.sqrt(np.sum((poptemp-pop[k])**2))
                if ds >= epsilon: 
                    truthvect[k] = True
            
            if truthvect.all(): 
                pop[i] = poptemp 
                i += 1 
            j += 1
        
        return pop


    def cost_eval(self, pop):
        """
        Evaluación de la función a optimizar. 

        costfunc :: función a optimizar
        pop :: población de individuos
        I :: intervalo del espacio de solución

        """
        costfunc = self.costfunc
        npop = self.npop
        poptemp = self.denormalize_params(pop)
        costvect = np.zeros(npop)
        for i in range(npop):
            costvect[i] = costfunc(poptemp[i])
        
        return costvect


    def selection(self, pop):

        costs = self.cost_eval(pop)
        sind = np.argsort(costs)
        pop = pop[sind]
        costs = costs[sind]
        parents = pop[0:self.nparents]

        return parents, pop, costs


    def mating(self, parents):
        """
        Función que implementa el cruce entre los 2 individuos escogidos

        pfather, pmother :: individuos escogidos para el cruce.
        n :: cantidad de genes a cruzar  
        """
        npop = self.npop
        ngenes = self.ngenes
        nparents = self.nparents
        ncross = self.ncross

        nextgen = np.zeros((npop, ngenes))        
        nChildren = npop - nparents
        nextgen[0:nparents] = parents

        for i in range(nChildren):
            indParents = np.random.randint(0,nparents, 2)
            pfather = parents[indParents[0]]
            pmother = parents[indParents[1]]
            
            offspring1 = pfather
            offspring2 = pmother

            indGenes = np.random.randint(0,ngenes, ncross)

            for i in range(ncross):
                A = indGenes[i]
                pm = offspring1[A]
                pf = offspring2[A]
                beta = np.random.rand(1)
                pnew1 = pm - beta*(pm-pf)
                pnew2 = pf + beta*(pm-pf)

                offspring1[A] = pnew2
                offspring2[A] = pnew1

            if np.random.rand()>0.5: 
                child = offspring1
            else: 
                child = offspring2
            
            nextgen[nparents+i] = child

        return nextgen


    def mutations(self,pop):
        """
        Función que muta genes de la población

        percent :: porcentaje de genes a mutar
        pop :: población de individuos

        """
        percent = self.mutprob
        a, b = pop.shape
        mutnum = np.round(a*b*percent/100).astype(int)

        for i in range(mutnum):
            ind1 = np.random.randint(0,a)
            ind2 = np.random.randint(0,b)

            pop[ind1, ind2] = np.random.rand()

        return pop 

    def distance_calc(self, pop):

        dpop = self.denormalize_params(pop)
        costs = self.cost_eval(pop)
        inds = np.argsort(costs)
        dpop = dpop[inds]
        dif = dpop[0] - dpop[1:len(dpop)]
        d = np.sqrt(np.sum(dif**2, axis=1))

        return np.max(d)

        #check for 

    def GA(self):
        """
        Función de implementación del algoritmo genético.

        costfunc :: función a optimizar
        epsilon :: radio de distanciamiento entre individuos para un debido esparcimiento. 
        ite :: número de iteraciones
        ngenes :: número de variables del problema 
        percent :: porcentaje de genes a mutar
        I :: intervalo del espacio de solución
        npop :: número de individuos en la población
        nkeep :: número de individuos que se escogen para ser cruzados 
        
        """
        i = 0
        j = 0
        Nred = 0
        pop = self.pop_generation() #definición de población inicial
        coords = []
        bestcost = []
        k = True 
        gens = 1
        while k:
            
            parents, pop, costs = self.selection(pop)#escogencia de los mejores candidatos
            coords.append(parents[0]) 
            bestcost.append(costs[0])
            pop = self.mating(parents)
            pop = self.mutations(pop)
            j += 1
            if j > 10: 
                movavg = np.mean(abs(np.diff(bestcost[-10:len(bestcost)])))

                delta = np.sqrt(np.sum(self.I[:,1] - self.I[:,0])**2)

                if delta > self.abstol: 

                    if movavg <= self.gentol:
                        Nred += 1
                        gens = 0
                        i += 1
                        j = 0
                        #print("Gen update")
                        bestpoint = self.denormalize_params(pop[0])
                        
                        self.I[:,0] = bestpoint - delta/(Nred*self.Rred)
                        self.I[:,1] = bestpoint + delta/(Nred*self.Rred)
                        self.gentol /= 10

                        lowerbool = self.I[:,0] >= self.staticI[:,0] 
                        upperbool = self.I[:,1] <= self.staticI[:,1] 
                        lowerind = np.where(lowerbool == False)
                        upperind = np.where(upperbool == False)

                        self.I[lowerind,0] = self.staticI[lowerind,0] 
                        self.I[upperind,1] = self.staticI[upperind,1]

                        pop = self.pop_generation()           
                        #print(self.I)                 
                else: 
                    #print("D")
                    k = False

            if i > self.maxiter:
                #print("H") 
                k = False    

            gens += 1
            i += 1

        decodedpop = self.denormalize_params(pop) #decodificación de parámetros
        decodedcoords = self.denormalize_params(coords) #decodificación de coordenadas iniciales
        print("Minimum Found after: "+str(i)+" iteraciones")
        return decodedpop[0], costs, decodedcoords, bestcost

