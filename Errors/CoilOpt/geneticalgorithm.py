"""
Rutina de código que implementa una versión del algoritmo genético continuo.

Elaborado por: Johansell Villalobos 
Ing. Física, TEC

"""

import numpy as np

def normalize_params(p, I):
    """
    Función que normaliza los genes respecto al espacio de solución.
    p :: gen/dato
    I :: intervalo que determina el espacio de solución 
    """

    phi = I[1]
    plo = I[0]
    pnorm = (p-plo)/(phi-plo)
    return pnorm

def denormalize_params(pnorm, I):
    """
    Función que devuelve los genes a sus valores correspondientes en el espacio de solución.
    pnorm :: gen/dato normalizado 
    I :: intervalo que determina el espacio de solución 
    """
    phi = I[1]
    plo = I[0]
    p = (phi-plo)*pnorm + plo
    return p

def pop_generation(npop,ngenes,epsilon):
    """
    Función generadora de una población esparcida por el espacio de solución. 

    npop :: cantidad de individuos en la población
    ngenes :: cantidad de variables a optimizar 
    epsilon :: radio de distanciamiento entre individuos para un debido esparcimiento (normalizado).  
    """
    pop = np.zeros((npop, ngenes))
    pop[0] = np.random.rand(ngenes)
    i = 1
    j = 0
    print("hi")
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
    print(j)
    
    return pop

def mating(pfather, pmother, n):
    """
    Función que implementa el cruce entre los 2 individuos escogidos

    pfather, pmother :: individuos escogidos para el cruce.
    n :: cantidad de genes a cruzar  
    """
    a = len(pfather)

    offspring1 = pfather
    offspring2 = pmother

    ind = np.random.randint(0,a, n)

    for i in range(n):
        A = ind[i]
        pm = offspring1[A]
        pf = offspring2[A]
        beta = np.random.rand(1)
        pnew1 = pm - beta*(pm-pf)
        pnew2 = pf + beta*(pm-pf)

        offspring1[A] = pnew2
        offspring2[A] = pnew1

    return offspring1, offspring2

def mutations(percent,pop):
    """
    Función que muta genes de la población

    percent :: porcentaje de genes a mutar
    pop :: población de individuos

    """

    a, b = pop.shape
    mutnum = np.round(a*b*percent/100).astype(int)

    for i in range(mutnum):
        ind1 = np.random.randint(0,a)
        ind2 = np.random.randint(0,b)

        pop[ind1, ind2] = np.random.rand()

    return pop 

def cost_eval(costfunc,pop, I):
    """
    Evaluación de la función a optimizar. 

    costfunc :: función a optimizar
    pop :: población de individuos
    I :: intervalo del espacio de solución

    """

    npop = len(pop)
    poptemp = denormalize_params(pop, I)
    costvect = np.zeros(npop)
    for i in range(npop):
        costvect[i] = costfunc(poptemp[i])
    
    return costvect

def GA(costfunc, epsilon, ite, ngenes, percent, I, npop=8, nkeep=4):
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

    pop = pop_generation(npop,ngenes,epsilon) #definición de población inicial
    coords = np.zeros((ite,2))
    nchildren = npop - nkeep 
    if nchildren % 2 == 0:
        for i in range(ite):
            nextgen = np.zeros(pop.shape)
            costvect = cost_eval(costfunc, pop, I) #evaluación de función a optimizar para escoger los mejores candidatos
            sortedind = np.argsort(costvect, axis=None)[0:nkeep] 
            parents = pop[sortedind] #escogencia de los mejores candidatos
            coords[i] = parents[0] 
            nextgen[0:nkeep] = parents
            children = []
            for i in range(nchildren // 2): #ciclo de cruce de individuos 
                ind1 = np.random.randint(0,len(parents), 2)
                child1, child2 = mating(parents[ind1[0]], parents[ind1[1]], ngenes)
                children.append(child1)
                children.append(child2)

            nextgen[nkeep:len(nextgen)] = np.array(children) #se obtiene la nueva generación
            pop = nextgen #se sobreescribe la nueva generación
            pop = mutations(percent, pop) #ronda de mutación
    
    else: 
        Exception("La resta de npop y nparents es impar, \t el algoritmo funciona partiendo de que esta sea par.")
        
    decodedpop = denormalize_params(pop, I) #decodificación de parámetros
    decodedcoords = denormalize_params(coords, I) #decodificación de coordenadas iniciales

    return decodedpop[sortedind][0], costvect[sortedind][0], decodedcoords




