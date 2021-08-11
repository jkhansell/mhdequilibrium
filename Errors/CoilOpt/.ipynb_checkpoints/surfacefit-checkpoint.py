"""
Rutina de código que obtiene los coeficientes de la serie de Fourier que 
representa la superficie mapeada por la línea de campo calculada. La forma 
matemática es la siguiente, 

∑∑𝑃c*cos(𝑚𝜃−𝑛𝑁𝜙) + 𝑃s*sin(𝑚𝜃−𝑛𝑁𝜙)

Luego de ciertas discretizaciones se puede representar la serie como la 
suma de Einstein, 

𝑆𝑝=𝑃^𝑖𝑗*𝑀_𝑖𝑗. 

𝑆𝑝 -> R(𝜃,𝜙) & Z(𝜃,𝜙)

"""


import numpy as np
from scipy.optimize import minimize


def initialize_data(data, nplanes):
    '''
    Inicializa los datos, en "nplanes" planos seleccionados del intervalo [0, 2pi] 
    para realizar el ajuste de parámetros. 

    data :: datos a analizar. 
    '''
      
    torplanes = np.linspace(0, 2*np.pi, num=nplanes, endpoint=False)
    dphi = 0.01
    phi_up = torplanes + dphi/2
    phi_down = torplanes - dphi/2
    R0 = 0.2477 

    boolarr = [np.logical_and(
                    data[:,1] < phi_up[i],
                    data[:,1] > phi_down[i]) for i in range(len(torplanes))]
                    #comparativo para los nplanes. 
    
    boolarr = np.array(boolarr).any(axis = 0)
    pos = np.where(boolarr)[0]
    data = data[pos]
    
    R = data[:,0]
    phi = data[:,1]
    z = data[:,2]
    r = data[:,3]
    theta = data[:,4]
    B = data[:,5]
    Br = data[:,6]
    Bphi = data[:,7]
    Bz = data[:,8]

    return np.array([R, phi, z, r, theta,B, Br, Bphi, Bz])

def data_setup(dataset, m, n, nfp):
    '''
    Se construyen las matrices de cosenos y senos para obtener las series de fourier. 
    
    dataset:: datos proporcionados por la función initialize_data()
    m, n :: modos de la serie escogida 
    nfp :: periodo de campo 
    '''

    rdata, phi, zdata, _,theta, *_ = dataset
    cosmat = np.zeros((m+1,2*n+1,len(rdata)))
    sinmat = np.zeros((m+1,2*n+1,len(zdata)))
    
    for j in range(m+1): 
        for i in range(-n, n+1):
            cosmat[j,i+n] = np.cos(j*theta-i*nfp*phi)
            sinmat[j,i+n] = np.sin(j*theta-i*nfp*phi)
    
    return [rdata, zdata], np.array([cosmat, sinmat]), [phi, theta]

def matsmesh_setup(m,n,n_fp,num):
    """
    Función que construye matrices de cosenos y senos distribuidas equitativamente
    en el espacio (𝜃,𝜙) para su posterior análisis diferencial.

    m, n :: modos de la serie escogida 
    nfp :: periodo de campo 
    num :: número de datos por dimensión
    """

    phi = np.linspace(0,2*np.pi,num=num)
    tht = np.linspace(0,2*np.pi,num=num)

    p, t = np.meshgrid(phi, tht)

    cosmatpp = np.zeros((m+1,2*n+1,len(phi),len(tht)))
    sinmatpp = np.zeros((m+1,2*n+1,len(phi),len(tht)))

    for j in range(m+1): 
        for i in range(-n, n+1):
            cosmatpp[j,i+n] = np.cos(j*t-i*n_fp*p)
            sinmatpp[j,i+n] = np.sin(j*t-i*n_fp*p)

    return [cosmatpp, sinmatpp], [phi, tht], [p,t]

def errors(data, mats, m, n):
    """
    Función que define las funciones de error respecto a los datos.

    data :: datos proporcionados por la función initialize_data()
    mats :: matrices proporcionadas por la función data_setup()
    m, n :: modos de la serie escogida  
    """
    def square_error_r(par):
        par = np.array(par).reshape(m+1,2*n+1)
        dy = (1/len(data[0]))*(data[0] - np.einsum('ij, ijm->m', par, mats[0],optimize=True))**2
        return np.sum(dy)

    def square_error_s(par):
        par = np.array(par).reshape(m+1,2*n+1)
        dy = (1/len(data[1]))*(data[1] - np.einsum('ij, ijm->m', par,mats[1],optimize=True))**2
        return np.sum(dy)
    return square_error_r, square_error_s

def optim_funcs(funcs, m, n):
    """
    Función que encuentra los coeficientes de la serie de fourier

    funcs :: funciones de error determinadas por errors()
    m, n :: modos de la serie escogida

    """

    p0 = np.random.rand(m+1, 2*n+1)

    res_r = minimize(funcs[0], p0, method='BFGS', tol=1e-5)
    res_z = minimize(funcs[1], p0, method='BFGS', tol=1e-5)
    rcoefs = res_r.x.reshape(m+1, 2*n+1)
    zcoefs = res_z.x.reshape(m+1, 2*n+1)
    


    return rcoefs, zcoefs

