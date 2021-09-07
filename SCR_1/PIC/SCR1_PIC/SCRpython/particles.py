
import numpy as np

#local files
import coil_biot as cb
import read_coils as rc


class comp_domain(object):
    def __init__(self, nr, nphi, nz, R0, rm, zm):
        self.nr = nr
        self.nphi = nphi
        self.nz = nz
        self.R0 = R0
        self.rm = rm
        self.zm = zm

    def dhs(self):
        self.dh = np.zeros(3) 

        self.dh[0] = 2*self.rm/self.nr
        self.dh[1] = 2*np.pi/self.nphi
        self.dh[2] = 2*self.zm/self.nz
    
def coord_transform(position, vector=None, coord_type="to_cart", entity="position"):

    if coord_type == "to_cart" and entity == "position": 

        R = position[0]
        phi = position[1]
        Z = position[2] 
        x = R*np.cos(phi)
        y = R*np.sin(phi)

        return np.array([x,y,Z])

    elif coord_type == "to_cil" and entity == "position":
        x = position[0]
        y = position[1]
        z = position[2]

        R = np.sqrt(x**2+y**2)
        phi = np.arctan2(y,x)+np.pi

        return np.array([R, phi, z])

    elif coord_type == "to_cart" and entity == "vector":
        if vector.all() != None: 
            R = position[0]
            phi = position[1]
            Z = position[2]

            a_r= vector[0]
            a_phi = vector[1]
            a_z = vector[2] 

            a_x = a_r*np.cos(phi) - a_phi*np.sin(phi)
            a_y = a_r*np.sin(phi) + a_phi*np.cos(phi)
        
        else: 
            raise(ValueError, "Vector parameter is of type: "
                        +str(type(vector)+",\n\tinsert ndarray"))

        return np.array([a_x, a_y, a_z])
        
    elif coord_type == "to_cil" and entity == "vector": 
        if vector.all() != None: 
            x = position[0]
            y = position[1]
            z = position[2]

            phi = np.arctan2(y,x)+np.pi

            a_x = vector[0]
            a_y = vector[1]
            a_z = vector[2] 

            a_r = a_x*np.cos(phi) + a_y*np.sin(phi)
            a_phi = -a_x*np.sin(phi)+a_y*np.cos(phi)

            return np.array([a_r,a_phi,a_z])

        else: 
            raise(ValueError, "Vector parameter is of type: "
                        +str(type(vector)+",\n\tinsert ndarray"))

def scalar_grad(comp_domain, data): 

    """
    Compute the gradient of a scalar field ∇A, if a vector field is passed,
    the algorithm computes ∇|A|. 

    """

    dh = comp_domain.dh
    nr = comp_domain.nr
    nphi = comp_domain.nphi
    nz = comp_domain.nz

    grad = np.zeros((nr,nphi,nz,3))
    Amag = np.linalg.norm(data, axis = 3)

    R_init = comp_domain.R0 - comp_domain.rm
    for i in range(nr):
        for j in range(nphi):
            for k in range(nz):

                if i == 0: 
                    dAdr = (-3*Amag[i,j,k]+4*Amag[i+1,j,k]-Amag[i+2,j,k])/(2*dh[0])                    
                elif i == nr-1: 
                    dAdr = (3*Amag[i,j,k]-4*Amag[i-1,j,k]+Amag[i-2,j,k])/(2*dh[0])
                else: 
                    dAdr = (Amag[i+1,j,k]-Amag[i-1,j,k])/(2*dh[0])
                if j == 0: 
                    dAdphi = (-3*Amag[i,j,k]+4*Amag[i,j+1,k]-Amag[i,j+2,k])/(2*(R_init+i*dh[0])*dh[1])
                elif j == nphi-1: 
                    dAdphi = (3*Amag[i,j,k]-4*Amag[i,j-1,k]+Amag[i,j-2,k])/(2*(R_init+i*dh[0])*dh[1])
                else:
                    dAdphi = (Amag[i,j+1,k]-Amag[i,j-1,k])/(2*(R_init+i*dh[0])*dh[1])

                if k == 0:
                    dAdz = (-3*Amag[i,j,k]+4*Amag[i,j,k+1]-Amag[i,j,k+2])/(2*dh[2]) 
                elif k == nz-1:
                    dAdz = (3*Amag[i,j,k]-4*Amag[i,j,k-1]+Amag[i,j,k-2])/(2*dh[2])
                else:                  
                    dAdz = (Amag[i,j,k+1]-Amag[i,j,k-1])/(2*dh[2])
                
                grad[i,j,k] = [dAdr, dAdphi, dAdz]
    
    return grad

def curl(comp_domain, A):


    """
    
    FIX CURL DEFINITION
    curl vector defined in 3d it's supposed to be in 1d components following dl 
    
    """
    R_init = comp_domain.R0 - comp_domain.rm
    dh = comp_domain.dh
    nr = comp_domain.nr
    nphi = comp_domain.nphi
    nz = comp_domain.nz
    curl = np.zeros((nr-1,nphi-1,nz-1,3))
    for i in range(nr-1):
        for j in range(nphi-1):
            for k in range(nz-1):

                #6 faces averages 

                F_rho1 = 0.25*(A[i,j,k] + A[i,j+1,k] + A[i,j+1,k+1] + A[i,j,k+1])
                F_rho2 = 0.25*(A[i+1,j,k] + A[i+1,j+1,k] + A [i+1,j+1,k+1] + A[i+1,j,k+1])

                F_phi1 = 0.25*(A[i,j,k] + A[i+1,j,k] + A[i+1,j,k+1] + A[i,j,k+1])
                F_phi2 = 0.25*(A[i,j+1,k] + A[i+1,j+1,k] + A[i+1,j+1,k+1] + A[i,j+1,k+1])

                F_z1 = 0.25*(A[i,j,k] + A[i+1,j,k] + A[i+1,j+1,k] + A[i,j+1,k])
                F_z2 = 0.25*(A[i,j,k] + A[i+1,j,k+1] + A[i+1,j+1,k+1] + A[i,j+1,k+1])
                
                #curl components B·dl
                #

                curl_rho = (1/((R_init+(i+0.5)*dh[0])*dh[1]*dh[2]))*(F_z2[1]*(i+0.5)*dh[0]*dh[1]-F_z1[1]*(i+0.5)*dh[0]*dh[1]+F_phi1[2]*dh[2]-F_phi2[2]*dh[2])

                curl_phi = (1/(dh[0]*dh[2])*(F_z2[0]*dh[0]-F_z1[0]*dh[0]+F_rho1[2]*dh[2]-F_rho2[2]*dh[2]))

                curl_z = (1/((R_init+(i+0.5)*dh[0])*dh[1]*dh[0]))*(F_phi1[0]*dh[0]-F_phi2[0]*dh[0]-F_rho1[1]*(i*dh[0]*dh[1])+F_rho2[1]*(i+1)*dh[0]*dh[1])
                 
                curl[i,j,k] = np.array([curl_rho, curl_phi, curl_z])

    return curl


def trilinear_interpolator(pos, data, comp_domain):
        """
        trilinear grid interpolator
        """
        dh = comp_domain.dh

        i = int((pos[0]-(comp_domain.R0-comp_domain.rm))//dh[0])
        j = int(pos[1]//dh[1])
        k = int((pos[2]+comp_domain.zm)//dh[2])

        r_d = (pos[0]-i*dh[0])/dh[0]
        phi_d = (pos[0]*pos[1]-(i*dh[0]*j*dh[1]))/dh[0]*dh[1]
        z_d = (pos[2]-k*dh[2])/dh[2]

        c000 = data[i,j,k]
        c001 = data[i,j,k+1]
        c010 = data[i,j+1,k]
        c011 = data[i,j+1,k+1]

        c100 = data[i+1,j,k]
        c101 = data[i+1,j,k+1]
        c110 = data[i+1,j+1,k]
        c111 = data[i+1,j+1,k+1]

        c00 = c000*(1-r_d)+c100*r_d
        c01 = c001*(1-r_d)+c101*r_d
        c10 = c010*(1-r_d)+c110*r_d
        c11 = c011*(1-r_d)+c111*r_d
        
        c0 = c00*(1-phi_d)+c10*phi_d
        c1 = c01*(1-phi_d)+c11*phi_d

        c = c0*(1-z_d)+c1*z_d

        return c



class magneticfield(object):

    """
    
    """

    def __init__(self, comp_domain, ncomp):
        self.comp_domain = comp_domain
        self.ncomp = ncomp
    
    def initial_condition(self):
        dh = self.comp_domain.dh
        R0 = self.comp_domain.R0
        rm = self.comp_domain.rm
        zm = self.comp_domain.zm
        nr = self.comp_domain.nr
        nphi = self.comp_domain.nphi
        nz = self.comp_domain.nz
        ncomp = self.ncomp
        coildata = rc.read_coil("coils.txt")
        self.data = np.zeros((nr, nphi, nz, ncomp))
        self.magB = np.zeros((nr, nphi, nz, ncomp))
        for i in range(nr):
            for j in range(nphi):
                for k in range(nz):
                    R =  R0 - rm + i*dh[0]
                    phi = j*dh[1]
                    Z = -zm + k*dh[2]
                    cartpos = coord_transform(np.array([R,phi,Z]))
                    bdata = cb.biot(coildata=coildata, point=cartpos)
                    self.data[i,j,k] = coord_transform(position=np.array([R,phi,Z]), vector=bdata, coord_type="to_cil", entity="vector")
                    self.magB[i,j,k] = np.sqrt(np.sum(self.data[i,j,k]**2))*np.ones(3)

        self.gradB = scalar_grad(self.comp_domain, self.data)
        self.curlunitB = curl(self.comp_domain,np.multiply(self.data,1/self.magB))

    def interpolfields(self, pos):
        
        pointB = trilinear_interpolator(pos, self.data, self.comp_domain)
        pointgradB = trilinear_interpolator(pos, self.gradB, self.comp_domain)
        pointcurlunitB = trilinear_interpolator(pos, self.curlunitB, self.comp_domain)

        return [pointB, pointgradB, pointcurlunitB]    

def drift_velocity(x, v, B):

    Bvec, gradB, curlunitB = B.interpolfields(x)
    magB = np.sqrt(np.sum(Bvec**2))
    unitB = Bvec/magB

    vparallelmag = np.dot(v, unitB)
    vparallel = vparallelmag*unitB
    vperp = v - vparallel
    vperpmag = np.sqrt(np.sum(vperp**2))

    # b x ∇b drift
    v_gradB = (vperpmag**2/(2*magB**2))*(np.cross(unitB, gradB))

    # b x (b·∇)b = b x (0.5∇|b|^2 - b x (∇ x b)

    # curv drift

    bdotNb = -np.cross(unitB, curlunitB)
    v_curv = (vparallelmag**2/magB)*(np.cross(unitB, bdotNb))

    v_drift = v_gradB+v_curv

    return v_drift, Bvec


class Particle(object): 
    """
    
    """
    def __init__(self, q, m, v_0, x_0): 
        self.q = q 
        self.m = m
        self.v = v_0
        self.x = x_0

    def update_drift_velocity(self, E, B, dt):
        Bvec, _, _ = B.interpolfields(self.x)
        qm = self.q/self.m
        v_minus = self.v + qm*E*dt/2

        t = qm*Bvec*dt/2
        v_prime = v_minus + np.cross(v_minus,t)
        s = 2*t/(1+np.dot(t,t))
        v_plus = v_minus + np.cross(v_prime,s)

        self.v = v_plus + qm*E*dt 
        self.x += self.v*dt

        
        




        
