#include "ParticleMesh.H"

using namespace amrex;

void 
ParticleMeshFuncs::gather_fields
( 
    const Vector<Real>& position, 
    const Array4<const Real>& Bx,
    const Array4<const Real>& By, 
    const Array4<const Real>& Bz,
    const Array4<const Real>& gradBx,
    const Array4<const Real>& gradBy, 
    const Array4<const Real>& gradBz,
    const Array4<const Real>& Jacobian,
    Vector<Real>& Bp, 
    Vector<Real>& gradBp, 
    Vector<Real>& jacobianp,
    const GpuArray<Real, AMREX_SPACEDIM>& plo,
    const GpuArray<Real, AMREX_SPACEDIM>& dxi
)
{   
    Real lx = (position[0]-plo[0])*dxi[0];
    Real ly = (position[1]-plo[1])*dxi[1];     // logical node-centered coordinates
    Real lz = (position[2]-plo[2])*dxi[2];

    int i = floor(lx);
    int j = floor(ly);
    int k = floor(lz);

    Real di = lx - i;
    Real dj = ly - j;
    Real dk = lz - k;

    Real sx[] = {1.- di, di};
    Real sy[] = {1.- dj, dj};
    Real sz[] = {1.- dk, dk};

    // trilinear interpolation algorithm

    constexpr int ixmin = 0;
    constexpr int ixmax = 1;
    constexpr int iymin = 0;
    constexpr int iymax = 1;
    constexpr int izmin = 0;
    constexpr int izmax = 1;

    Bp = {0.0,0.0,0.0};
    gradBp = {0.0,0.0,0.0};
    jacobianp = {0.0,0.0,0.0,
                 0.0,0.0,0.0,
                 0.0,0.0,0.0,};

    for (int kk = izmin; kk <= izmax; ++kk){
        for (int jj = iymin; jj <= iymax; ++jj){
            for (int ii = ixmin; ii <= ixmax; ++ii){
                
                Bp[0] += sx[ii]*sy[jj]*sz[kk]*Bx(i+ii,j+jj,k+kk);
                Bp[1] += sx[ii]*sy[jj]*sz[kk]*By(i+ii,j+jj,k+kk);
                Bp[2] += sx[ii]*sy[jj]*sz[kk]*Bz(i+ii,j+jj,k+kk);

                gradBp[0] += sx[ii]*sy[jj]*sz[kk]*Bx(i+ii,j+jj,k+kk);
                gradBp[1] += sx[ii]*sy[jj]*sz[kk]*By(i+ii,j+jj,k+kk);
                gradBp[2] += sx[ii]*sy[jj]*sz[kk]*Bz(i+ii,j+jj,k+kk);

                jacobianp[0] += sx[ii]*sy[jj]*sz[kk]*Jacobian(i+ii,j+jj,k+kk,0);
                jacobianp[1] += sx[ii]*sy[jj]*sz[kk]*Jacobian(i+ii,j+jj,k+kk,1);
                jacobianp[2] += sx[ii]*sy[jj]*sz[kk]*Jacobian(i+ii,j+jj,k+kk,2);
                jacobianp[3] += sx[ii]*sy[jj]*sz[kk]*Jacobian(i+ii,j+jj,k+kk,3);
                jacobianp[4] += sx[ii]*sy[jj]*sz[kk]*Jacobian(i+ii,j+jj,k+kk,4);
                jacobianp[5] += sx[ii]*sy[jj]*sz[kk]*Jacobian(i+ii,j+jj,k+kk,5);
                jacobianp[6] += sx[ii]*sy[jj]*sz[kk]*Jacobian(i+ii,j+jj,k+kk,6);
                jacobianp[7] += sx[ii]*sy[jj]*sz[kk]*Jacobian(i+ii,j+jj,k+kk,7);
                jacobianp[8] += sx[ii]*sy[jj]*sz[kk]*Jacobian(i+ii,j+jj,k+kk,8);

            }
        }
    }
}
