import numpy as np
import scipy.linalg as SPL


EPS_0 = 8.85418782e-12;  	#C/(V*m), vacuum permittivity
QE = 1.602176565e-19;		#C, electron charge
AMU = 1.660538921e-27;		#kg, atomic mass unit
ME = 9.10938215e-31;		#kg, electron mass
K = 1.380648e-23;			#J/K, Boltzmann constant
C = 3.0e8;					#speed of light
EvToK = QE/K;				#1eV in K ~ 11604
MU_0 = 4*np.pi*1e-7;		#vacuum permisivity

class World(object):
    def __init__(self, nr, nphi, nz):
        self.nr = nr
        self.nphi = nphi
        self.nz = nz
        self.nodeVolumes = np.zeros(nr,nphi,nz)
        self.phi = np.zeros(nr,nphi,nz)
        self.rho = np.zeros(nr,nphi,nz)
        self.E = np.zeros(nr,nphi,nz,3)
    
    def setExtents(self, R0, rm, zm): 
        self.R0 = R0
        self.rm = rm
        self.zm = zm
        self.dh = np.zeros(3)

        self.dh[0] = 2*rm/self.nr
        self.dh[1] = 2*np.pi/self.nphi
        self.dh[2] = 2*zm/self.nz
        self.a = np.sqrt(zm**2+rm**2)

        self.computeNodeVolumes()

    def computeNodeVolumes(self):

        dh = self.dh
        nr = self.nr
        nphi = self.nphi        
        nz = self.nz

        for i in range(nr):
            for j in range(nphi):
                for k in range(nz): 
                    V = i*dh[0]**2*dh[1]*dh[2]
                    
                    if  i == 0 or i ==nr-1:
                        V*=0.5
                    if  j == 0 or j ==nphi-1:
                        V*=0.5
                    if  k == 0 or k ==nz-1:
                        V*=0.5
                    
                    self.nodeVolumes[i,j,k] = V

class Potential_Solver(object):

    def __init__(self, world, max_it, tol): 
        self.world = world
        self.max_it = max_it 
        self.tol = tol

    def solveGSLinear(self, A, x, b)->bool:
        L2 = 0
        converged = False
        P, L, U = SPL.lu(A)
