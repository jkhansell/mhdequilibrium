import numpy as np 
import coil_parametrization as cp
import Rk4 
import coord_transform as ct
import surfacefit as sf
import propcalc as pc
import errtorpath as etp

def cost(par, parshape, X, cur, m, n, nfp, num, nplanes, init_pos, L, N):
    
    """
    Función que devuelve el error promedio como variable a optimizar por 
    el algoritmo genético. Implementa el proceso en conjunto presentado en
    'Optimización de Bobinas.ipynb'. 
    """

    par.reshape(parshape)
    param_data = cp.coefstodata(par, X, cur)
    tddata = np.transpose(param_data, axes=(0,2,1))
    linedata, _ = Rk4.RK4(L, N, init_pos, tddata)
    areadata, errdata = ct.carttotorcil(linedata, 0.2477)
    surfdata = sf.initialize_data(areadata,nplanes)
    fitdata, fitmats, _ = sf.data_setup(surfdata, m, n, nfp)
    errs = sf.errors(fitdata, fitmats, m, n)
    coefs = sf.optim_funcs(errs, m,n)
    area = pc.normal_vector_area(coefs[0], coefs[1], num, nfp)
    orderddata = etp.DatosFinales(errdata)
    erloc = etp.ErrorLocal(orderddata)
    erprom = etp.ErrorSup(erloc, area)

    return erprom

