"""
Rutina de código que implementa la transformación de coordenadas de 
cartesianas a cilíndricas y a cilíndrico-toroidales. 

Elaborado por: Johansell Villalobos
Ing. Física, TEC
"""

import numpy as np

def carttotorcil(data, R0):
    """
    Función que devuelve las estructuras de datos en los sistemas coordenados 
    necesarios para el cálculo de área y de error local/promedio. 

    data :: estructura de datos proporcionada de la solución de la ecuación 
    diferencial. 

    R0: Radio mayor del SCR-1 o dispositivo.
    """

    x, y, z, b, bx, by, bz = data.T
    R = np.sqrt(x**2+y**2)

    phi = np.arctan2(y,x)
    posp = np.where(phi<0)
    phi[posp] = phi[posp]+2*np.pi

    Br = bx*np.cos(phi)+by*np.sin(phi)
    Bphi = -bx*np.sin(phi)+by*np.cos(phi)
    rtor = np.sqrt((R-R0)**2+z**2)
    theta = np.arctan2(z, (R-R0))+np.pi
    post = np.where(phi<0)
    phi[post] = phi[post]+2*np.pi

    datastruct1 = np.array([R, phi, z, rtor, theta, b, Br, Bphi, bz]).T
    datastruct2 = np.array([rtor, phi, theta, b, bx, by, bz]).T
    
    return datastruct1, datastruct2

