"""
Rutina de código que a partir de los coeficientes de salida del script 
'surfacefit.py' realiza cálculos de área y volumen de la superficie. 

Elaborado por: Johansell Villalobos 
Ing. Física, TEC
"""

from surfacefit import matsmesh_setup
import numpy as np 


def normal_vector_area(rc, zc, num, n_fp):
    """
    Función que devuelve el área partiendo de la discretización del cálculo 
    de el área superficial sobre un sistema (𝜃,𝜙) uniformemente espaciado en
    'num' divisiones.
    
    rc :: coeficientes de la serie R(𝜃,𝜙)
    zc :: coeficientes de la serie Z(𝜃,𝜙)
    num :: número de divisiones en cada dimensión
    n_fp :: periodo de campo 

    """
    m, n = rc.shape
    
    m = int(m-1)
    n = int((n-1)/2)

    mats, vects, mesh = matsmesh_setup(m,n, n_fp, num)
    
    dRdTcoefs = np.zeros(rc.shape)
    dRdPcoefs = np.zeros(rc.shape)
    dZdTcoefs = np.zeros(zc.shape)
    dZdPcoefs = np.zeros(zc.shape)

    for j in range(m+1): 
        for i in range(-n, n+1):
            
            dRdTcoefs[j,i+n] = -rc[j,i+n]*j
            dRdPcoefs[j,i+n] = rc[j,i+n]*i*n_fp
            dZdTcoefs[j,i+n] = zc[j,i+n]*j
            dZdPcoefs[j,i+n] = -zc[j,i+n]*i*n_fp
            
    rpp = np.einsum('ij, ijmk->mk', rc, mats[0])

    dRdT = np.einsum('ij, ijmk->mk', dRdTcoefs, mats[1])
    dRdP = np.einsum('ij, ijmk->mk', dRdPcoefs, mats[1])
    
    dZdT = np.einsum('ij, ijmk->mk', dZdTcoefs, mats[0])
    dZdP = np.einsum('ij, ijmk->mk', dZdPcoefs, mats[0])
    
    dPHIdT = np.array([dRdT*np.cos(mesh[0]),
                       dRdT*np.sin(mesh[0]),
                       dZdT])
    
    dPHIdP = np.array([dRdP*np.cos(mesh[0]) - rpp*np.sin(mesh[0]),
                       dRdP*np.sin(mesh[0]) + rpp*np.cos(mesh[0]),
                       dZdP])
            
    normals = np.cross(dPHIdP.T, dPHIdT.T)
    n_norm = np.linalg.norm(normals, axis=2)
    area = np.trapz(np.trapz(n_norm, vects[1]), vects[0])
    
    return area

