#global libraries
import numpy as np
import xarray as xr

#local files
import coil_biot as cb 
import read_coils as rc 

def pos_mesh(nr, nphi, nz, R0, rm, zm):
    """
    Calculates magnetic field grid to be used in SCR 1 particle calculations.
    """

    dh = np.zeros(3)
    dh[0] = 2*rm/nr
    dh[1] = 2*np.pi/nphi
    dh[2] = 2*zm/nz
    coildata = rc.read_coil("coils.txt")
    magnetic_mesh = np.zeros((nr, nphi, nz, 3))
    position = np.zeros((nr, nphi, nz, 3))

    for i in range(nr):
        for j in range(nphi):
            for k in range(nz):
                R =  R0 - rm + i*dh[0]
                phi = j*dh[1]
                Z = -zm + k*dh[2]
                X = R*np.cos(phi)
                Y = R*np.sin(phi)
                Bcart = cb.biot(coildata, np.array([X,Y,Z]))
                Br = Bcart[0]*np.cos(phi) + Bcart[1]*np.sin(phi)
                Bphi = -Bcart[0]*np.sin(phi) + Bcart[1]*np.cos(phi)
                position[i,j,k] = np.array([R, phi, Z])
                magnetic_mesh[i,j,k] = np.array([Br, Bphi, Bcart[2]])


    return magnetic_mesh, position


def write_magfield(nr, nphi, nz, R0, rm, zm):
    magmesh, position = pos_mesh(nr, nphi, nz, R0, rm, zm)
    Bfield = xr.Dataset(
        data_vars=dict(
            magnetic_field=(["NR", "NPhi", "NZ", "BXYZ"], magmesh)
        ),
        coords=dict(
            position = (["NR", "NPhi", "NZ", "XYZ"], position)
        )
    )

    Bfield.to_netcdf("Bfield.nc", mode="w")






    

    



    

