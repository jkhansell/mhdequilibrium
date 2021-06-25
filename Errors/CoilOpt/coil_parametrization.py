"""
Rutina de código para la parametrización en series de fourier de cada filamento.

Se toma un modelo de serie para cada dimensión (x,y,z) de la forma:

x^i = ∑ A_i cos(mt) + B_i sin(mt)

* Se utilizan sumas de Einstein para desarrollar estas series de una 
manera más eficiente.


Elaborado por: Johansell Villalobos
Ing. Física, TEC

"""


import numpy as np
from scipy.optimize import minimize

def cossin(t, m):
    """
    Función que inicializa los tensores de funciones trigonométricas de
    la serie.

    t :: puntos en el espacio theta a evaluar. (Espacio discretizado)
    """
    mvect = np.arange(0,m+1)
    mt = np.einsum('i,j', mvect, t)
    cosm = np.cos(mt)
    sinm = np.sin(mt)
    X = np.array([cosm, sinm])
    
    return X 

def square_error(coefs, X, data):
    """
    Definición de la función de error cuadrado para la minimización.

    coefs :: coeficientes de la serie a minimizar. 
    X :: tensor de funciones trigonométricas. 
    data :: datos discretos a ajustar.

    """
    coefs = coefs.reshape(2,-1)
    err = 1/len(data)*(data - np.einsum('ij, ijk->k', coefs, X))**2

    return np.sum(err)

def obtenerparams(x, y, z, X, m):
    """
    Función que obtiene los coeficientes de la serie respecto a las 
    dimensiones x, y, z de la bobina. 

    x, y, z :: coordenadas de la bobina. 
    X :: tensor de funciones trigonométricas. 
    m :: número de modos de la serie.     
    """

    coefs0 = np.random.rand(2, m+1)
    resx = minimize(square_error, coefs0 ,args = (X, x), method='BFGS', tol=1e-12)
    resy = minimize(square_error, coefs0 ,args = (X, y), method='BFGS', tol=1e-12)
    resz = minimize(square_error, coefs0 ,args = (X, z), method='BFGS', tol=1e-12)
    coefsx = resx.x.reshape(2,m+1)
    coefsy = resy.x.reshape(2,m+1)
    coefsz = resz.x.reshape(2,m+1)

    return coefsx, coefsy, coefsz

def parambobinas(data, m, t):
    """
    Función que implementa los coeficientes definidas anteriormente para devolver 
    los coeficientes óptimos de la serie. 

    data :: datos de las dimensiones de las bobinas. 
    m :: modos de la serie. 
    t :: puntos paramétricos de evaluación. [0, 2𝜋] 
    """
    a, *_ = data.shape
    X = cossin(t, m)
    coefs = np.zeros((a,3,2,m+1))
    for i in range(a): 
        x = data[i, 0,:]
        y = data[i, 1,:]
        z = data[i, 2,:]
        coefsx, coefsy, coefsz = obtenerparams(x,y,z,X,m)
        coefs[i, 0] = coefsx
        coefs[i, 1] = coefsy
        coefs[i, 2] = coefsz
    
    return coefs, X

def coefstodata(coefs, X, cur): 
    """
    Función que calcula la serie de datos proporcionados por la serie de Fourier, y 
    devuelve la estructura de datos con la corriente por bobina. 

    coefs :: coeficientes óptimos de la serie
    X :: Tensor de funciones trigonométricas. 
    cur :: corrientes ordenadas por bobina. 
    """
    _,_, a = X.shape
    b,c,_,_ = coefs.shape 
    data = np.zeros((b,c+1, a))
    for i in range(b):
        xpar = np.einsum('ij, ijk->k', coefs[i,0], X)
        ypar = np.einsum('ij, ijk->k', coefs[i,1], X)
        zpar = np.einsum('ij, ijk->k', coefs[i,2], X)
        data[i,0] = xpar
        data[i,1] = ypar
        data[i,2] = zpar
        data[i,3,0:-1] = cur[i]*np.ones(a-1)
    
    return data