import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import plotly.graph_objects as go

def initialize_data(filename, nplanes):
    print('\n\tLoading data...\n')
    torplanes = np.linspace(0, 2*np.pi, num=nplanes, endpoint=False)
    dphi = 0.01
    phi_up = torplanes + dphi/2
    phi_down = torplanes - dphi/2
    R0 = 0.2477 
 
    data = np.loadtxt(filename, dtype=float, skiprows=2, delimiter='\t')

    boolarr = [np.logical_and(
                    data[:,1] < phi_up[i],
                    data[:,1] > phi_down[i]) for i in range(len(torplanes))]
    
    boolarr = np.array(boolarr).any(axis = 0)
    pos = np.where(boolarr)[0]
    data = data[pos]
    
    rtor = data[:,0]
    thetator = data[:,2]
    phi = data[:,1]
    R = R0+np.multiply(rtor,np.cos(thetator))
    z = np.multiply(rtor,np.sin(thetator))
    B = data[:,3]
    Br = data[:,4]
    Bphi = data[:,5]
    Bz = data[:,6]
    print('\tFinished!\n')

    return np.array([R, z, thetator, rtor, phi, B, Br, Bphi, Bz])

def data_setup(dataset, m, n, nfp):

    print('\tSetting up matrices...\n')

    rdata, zdata, theta, _, phi, _, _, _, _ = dataset

    cosmat = np.zeros((m+1,2*n+1,len(rdata)))
    sinmat = np.zeros((m+1,2*n+1,len(zdata)))
    
    for j in range(m+1): 
        for i in range(-n, n+1):
            cosmat[j,i+n] = np.cos(j*theta-i*nfp*phi)
            sinmat[j,i+n] = np.sin(j*theta-i*nfp*phi)

    print('\tDone!\n')
    
    return [rdata, zdata], np.array([cosmat, sinmat])

def matsmesh_setup(m,n,n_fp,num):

    phi = np.linspace(0,2*np.pi,num=num)
    tht = np.linspace(0,2*np.pi,num=num)

    p, t = np.meshgrid(phi, tht)

    cosmatpp = np.zeros((m+1,2*n+1,len(phi),len(tht)))
    sinmatpp = np.zeros((m+1,2*n+1,len(phi),len(tht)))

    for j in range(m+1): 
        for i in range(-n, n+1):
            cosmatpp[j,i+n] = np.cos(j*t-i*n_fp*p)
            sinmatpp[j,i+n] = np.sin(j*t-i*n_fp*p)

    return np.array([cosmatpp, sinmatpp]), [phi, tht], [p,t]


def relative_error(rcof, zcof, dataset, m, n, nfp):

    data_o, mats = data_setup(dataset, m, n, nfp)

    rc = rcof.x.reshape(m+1,2*n+1)
    zc = zcof.x.reshape(m+1,2*n+1)
    

    rpp = np.einsum('ij, ijm->m', rc, mats[0])
    zpp = np.einsum('ij, ijm->m', zc, mats[1])

    error_R = abs(data_o[0]-rpp)
    error_Z = abs(data_o[1]-zpp)

    return np.max(error_R), np.max(error_Z)


def errors(data, mats, m, n):

    print('\tDefining error functions...\n')

    def square_error_r(par):
        par = np.array(par).reshape(m+1,2*n+1)
        dy = (1/len(data[0]))*(data[0] - np.einsum('ij, ijm->m', par, mats[0],optimize=True))**2
        return np.sum(dy)

    def square_error_s(par):
        par = np.array(par).reshape(m+1,2*n+1)
        dy = (1/len(data[1]))*(data[1] - np.einsum('ij, ijm->m', par,mats[1],optimize=True))**2
        return np.sum(dy)

    print('\tFinished!\n')
    
    return square_error_r, square_error_s

def optim_funcs(funcs, m, n):

    print('\tBeginning optimization...\n')

    p0 = np.random.rand(m+1, 2*n+1)

    res_r = minimize(funcs[0], p0, method='BFGS', tol=1e-5)
    res_s = minimize(funcs[1], p0, method='BFGS', tol=1e-5)

    print('\tOptimization finished succesfully!\n')

    return res_r, res_s


def normal_vector(rc, zc, num, n_fp):

    m, n = rc.shape
    
    m = int(m-1)
    n = int((n-1)/2)

    mats, _, mesh = matsmesh_setup(m,n, n_fp, num)
    
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
            
    return np.cross(dPHIdP.T, dPHIdT.T)


def fourier_series(X, param): 
       
    return np.einsum('ijk,ijkm', param, X, optimize=True)


def magnetic_error(i, data, m, n, nfp):
    _, _, _, _, _, _, *B = data
    _, X = data_setup(data, m, n, nfp)
    B = np.array(B)

    def vect_error(par):
        param = np.array(par).reshape(2, m+1, 2*n+1)
        dy = (B[i] - fourier_series(X, param))**2
        return (1/(B[i].size))*np.sum(dy)

    return vect_error

def multi_opt(data, m, n, nfp):
    p0 = np.random.rand(2,m+1, 2*n+1)
    ress = []
    for i in range(3): 
        ress.append(magnetic_error(i, data, m, n, nfp))
        
    coef = []
    for i in range(3):    
        coef.append(minimize(ress[i], p0, method='BFGS', tol=1e-5))
        print('done')

    np.savetxt('Br.txt', coef[0])
    np.savetxt('Bphi.txt', coef[1])
    np.savetxt('Bz.txt', coef[2])
     
     
    return coef

def visualize_magfield(files, m, n, nfp, num):
    mats, _, mesh = matsmesh_setup(m, n, nfp, num)

    brcoef = np.loadtxt(files[0]).reshape(2,m+1, 2*n+1)
    bphicoef = np.loadtxt(files[0]).reshape(2,m+1, 2*n+1)
    bzcoef = np.loadtxt(files[0]).reshape(2,m+1, 2*n+1)

    Br = np.einsum('ijk,ijklm->lm', brcoef, mats, optimize=True)
    Bphi = np.einsum('ijk,ijklm->lm', bphicoef, mats, optimize=True)
    Bz = np.einsum('ijk,ijklm->lm', bzcoef, mats, optimize=True)


    layout = go.Layout(scene=dict(aspectmode='data'))

    fig = go.Figure(data=go.Surface(
        x = mesh[0],
        y = mesh[1], 
        z = Br,),layout=layout)

    fig.show()

    fig2 = go.Figure(data=go.Surface(
        x = mesh[0],
        y = mesh[1], 
        z = Bphi,),layout=layout)

    fig2.show()

    fig3 = go.Figure(data=go.Surface(
        x = mesh[0],
        y = mesh[1], 
        z = Bz),layout=layout)

    fig3.show()


def normal_visualization(res_r, res_s, m, n, num, n_fp):

    print('\tPreparing visualization...\n')
                
    rc = np.loadtxt(res_r).reshape(m+1,2*n+1)
    zc = np.loadtxt(res_s).reshape(m+1,2*n+1)

    mats, _, mesh  = matsmesh_setup(m, n, n_fp, num)

    rpp = np.einsum('ij, ijmk->mk', rc, mats[0])
    zpp = np.einsum('ij, ijmk->mk', zc, mats[1])

    xpp = rpp*np.cos(mesh[0])
    ypp = rpp*np.sin(mesh[0])

    layout = go.Layout(
             scene=dict(
                 aspectmode='data'))

    fig2 = go.Figure(data=go.Surface(
        x = xpp,
        y = ypp, 
        z = zpp,),layout=layout)


    normvect = normal_vector(rc, zc, num, n_fp).T

    xflat = xpp.flatten()
    yflat = ypp.flatten()
    zflat = zpp.flatten()

    fig2.add_trace(go.Cone(
        x = xflat, 
        y = yflat, 
        z = zflat,
        u = normvect[0].flatten(),
        v = normvect[1].flatten(),
        w = normvect[2].flatten(),
        sizemode="absolute",
        sizeref=0.2
    ))

    fig2.show()

def main(filename, nplanes, m, n, nfp):

    dataset = initialize_data(filename, nplanes)
    data, mats = data_setup(dataset, m, n, nfp)
    funcs = errors(data, mats, m, n)
    rc, zc = optim_funcs(funcs, m, n)

    #num = int(input('Number of points in one dimension of mesh: '))
    #visualization(rc, zc, num, m, n, nfp)

    np.savetxt('Rc0085.txt',rc.x)
    np.savetxt('Zc0085.txt',zc.x)
    return rc, zc
