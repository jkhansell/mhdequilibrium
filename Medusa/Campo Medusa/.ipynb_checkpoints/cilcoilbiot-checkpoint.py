"""
Rutina de código que utiliza el algoritmo propuesto en el artículo
"Compact expressions for the Biot–Savart fields of a filamentary segment"


Elaborado por: Johansell Villalobos
Ing. Física, TEC

"""

import numpy as np

def biot(coildata, p):
    """
    Función que calcula el campo magnético de una serie de segmentos de bobina 
    definidos por puntos en el espacio.

    coildata :: datos de la bobina de la forma:
        
        Nbobinas
            - x
                -Ndatos
            - y
                -Ndatos
            - z 
                -Ndatos
            - cur
                -Ndatos
    
    el arreglo tiene dimensiones [Nbobinas, Ndimensiones (4), Ndatos].

    Uso: 
    coildata = read_coils('coils.txt')
    point = np.array([x,y,z]) -- [2.660556e-01,1.297578e-02,-6.691423e-04]
    B = biot(coildata, point)
    """
    x = p[0]*np.cos(p[1])
    y = p[0]*np.sin(p[1])
    point = np.array([x,y,p[2]])
    
    mu_0 = 1e-7
    data = coildata.T[0:3].T
    Rf = point - data
    _,b,_ = Rf.shape
    cur = coildata.T[3].T[:,0:b-1]
    coef = mu_0*cur
    coef_n = np.tile(coef.T, (3,1,1)).T
    e_s = np.diff(Rf, axis = 1)
    L = np.sqrt(np.sum(e_s**2, axis=2))
    L_2 = np.tile(L.T, (3,1,1)).T
    e_hat = e_s/L_2
    Rfmag = np.sqrt(np.sum(Rf**2,axis=2))
    num = 2*L*(Rfmag[:,0:b-1]+Rfmag[:,1:b])
    denom = (Rfmag[:,0:b-1]*Rfmag[:,1:b])*((Rfmag[:,0:b-1]+Rfmag[:,1:b])**2-L**2)
    coef2 = num/denom
    coef_m = np.tile(coef2.T, (3,1,1)).T
    e_vects = np.multiply(e_hat,coef_n)
    Ri = np.multiply(Rf[:,0:b-1], coef_m)
    B = np.cross(e_vects, Ri)
    tot = np.sum(B, axis=(1,0))   
    
    
    phi = np.arctan2(point[1],point[0])+np.pi
    Br = tot[0]*np.cos(phi)+tot[1]*np.sin(phi)
    Bphi = -tot[0]*np.sin(phi)+tot[1]*np.cos(phi) 
    
    return np.array([Br, Bphi, tot[2]])
