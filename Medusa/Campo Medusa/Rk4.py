"""
Rutina de código para la solución de la ecuación diferenciales de 
línea de campo magnético en coordenadas cartesianas por medio del método
Runge-Kutta 4.

                        dR/|dR| = B/|B|

Elaborado por: Johansell Villalobos
Ing. Física, TEC

"""

from coil_biot import biot

import numpy as np

def RK4(b, N, p0, coildata):
    """
    Función que implementa el algoritmo de Runge-Kutta 4 vectorialmente.

    b :: metros de línea de campo
    N :: número de pasos
    p0 :: punto inicial para el cálculo
    coildata :: datos de bobinas
    """
    p = np.zeros((N+1,7))
    p[0,0:3] = p0
    p[0,4:8] = biot(coildata, p[0,0:3])
    p[0,3] = np.sqrt(np.sum(p[0,4:8]**2))
    ite = np.linspace(0,b,N+1)
    h = b/N
    
    for i in range(1,N+1):
        
        B1 = biot(coildata, p[i-1,0:3])
        k1 = B1/np.sqrt(np.sum(B1**2))
        B2 = biot(coildata, p[i-1,0:3]+h*k1/2)
        k2 = B2/np.sqrt(np.sum(B2**2))
        B3 = biot(coildata, p[i-1,0:3]+h*k2/2)
        k3 = B3/np.sqrt(np.sum(B3**2))
        B4 = biot(coildata, p[i-1,0:3]+h*k3)
        k4 = B4/np.sqrt(np.sum(B4**2))
        
        p[i,0:3] = p[i-1,0:3] + (1/6)*h*(k1+2*k2+2*k3+k4)
        p[i,4:8] = biot(coildata, p[i,0:3])
        p[i,3] = np.sqrt(np.sum(p[i,4:8]**2))

    return p, ite