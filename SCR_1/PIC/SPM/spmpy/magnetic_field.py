"""
Magnetic field definition
"""

# global libraries

import numpy as np

# local files

import coil_biot as cb 


class comp_domain(object):
    """
    
    """

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













