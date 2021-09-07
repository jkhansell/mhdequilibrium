import numpy as np

class Grid(object):
    def __init__(self, coord, na, nb, nc, Ncomp):
        """
        coord: coordinate system to be used (cartesian :: 0, cylindrical :: 1)

        (a,b,c) -- (x,y,z) -- (R, phi, Z) 

        na :: number of cells in the a direction 
        nb :: number of cells in the b direction
        nc :: number of cells in the c direction

        n :: number of components in the grid  
        """
        self.coord = coord
        self.Ncomp = Ncomp 
        self.grid = np.zeros(na, nb, nc, Ncomp)
        
        if self.coord == 0: 
            self.nx = na
            self.ny = nb 
            self.nz = nc
            
        else:
            self.nR = na
            self.nPhi = nb
            self.nZ = nc
        
        

        def fill_initial_condition(self):