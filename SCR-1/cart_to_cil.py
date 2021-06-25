import numpy as np
import os 
import time

def cart_to_cil(filename, save_dir, main_dir):
    a = time.perf_counter()
    data = np.loadtxt(filename, dtype=float, skiprows=1, delimiter=',')
    x, y, z, b, bx, by, bz = data.T
    R = np.sqrt(x**2+y**2)
    phi = np.arctan2(y,x)
    Br = bx*np.cos(phi)+by*np.sin(phi)
    Bphi = -bx*np.sin(phi)+by*np.cos(phi)
    with open('tor'+filename,'w') as torfile:
        torfile.write("r,phi,theta,|B|,BR, Bphi, BZ\n")
        for i in range(len(R)):
            torfile.write('{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\t{:e}\n'.format(R[i], phi[i], z[i], b[i], Br[i], Bphi[i], bz[i]))
        torfile.close()
    end = time.perf_counter()
    print(end-a)


