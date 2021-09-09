"""
trilinear interpolator

"""


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