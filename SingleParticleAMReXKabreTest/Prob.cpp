#include <AMReX_REAL.H>
#include "Prob.h"
#include "Particles.h"

using namespace amrex; 

void init_values(Real3& Bp, Real3& vel, Real3& pos,
                 Real& magvpar, 
                 Real& mu, Real q, Real m)
{
    Real Bmag = mag(Bp);
    Real3 bhat = Bp/Bmag;
    Real vpar_abs = dot(bhat, vel);
    Real3 vpar = vpar_abs*bhat;
    Real3 vperp = vel - vpar;
    Real Omega = q*Bmag/m;
    Real3 rho = -1.*cross(bhat, vperp)/Omega;
    Real mu_part = m*mag(vperp)/(2*Bmag);
    pos += rho;
    magvpar = vpar_abs;
    mu = mu_part;
}


Tuple<Real, Real3>
K_calc
(
    Real3 Bp, 
    Real3 dmagB,
    Vector<Real> jacobianp,
    const Real charge, const Real mass,
    const Real& magvpar, Real& mu
)
{    
    Real magBp = mag(Bp);
    
    // Derivatives
    
    Real dBxdx = jacobianp[0]; Real dBxdy = jacobianp[1]; Real dBxdz = jacobianp[2];
    Real dBydx = jacobianp[3]; Real dBydy = jacobianp[4]; Real dBydz = jacobianp[5];
    Real dBzdx = jacobianp[6]; Real dBzdy = jacobianp[7]; Real dBzdz = jacobianp[8];

    Real dbxdx = (dBxdx*magBp-Bp[0]*dmagB[0])/pow(magBp,2);
    Real dbxdy = (dBxdy*magBp-Bp[0]*dmagB[1])/pow(magBp,2);
    Real dbxdz = (dBxdz*magBp-Bp[0]*dmagB[2])/pow(magBp,2);

    Real dbydx = (dBydx*magBp-Bp[1]*dmagB[0])/pow(magBp,2);
    Real dbydy = (dBydy*magBp-Bp[1]*dmagB[1])/pow(magBp,2);
    Real dbydz = (dBydz*magBp-Bp[1]*dmagB[2])/pow(magBp,2);

    Real dbzdx = (dBzdx*magBp-Bp[2]*dmagB[0])/pow(magBp,2);
    Real dbzdy = (dBzdy*magBp-Bp[2]*dmagB[1])/pow(magBp,2);
    Real dbzdz = (dBzdz*magBp-Bp[2]*dmagB[2])/pow(magBp,2);

    Real3 unitgradB;
    unitgradB[0] = dbxdx;
    unitgradB[1] = dbydy;
    unitgradB[2] = dbzdz;

    Real3 unitB = Bp/magBp;
    
    Real3 curlb;
    curlb[0] = (mass*magvpar/charge)*(dbzdy-dbydz);
    curlb[1] = (mass*magvpar/charge)*(dbxdz-dbzdx);
    curlb[2] = (mass*magvpar/charge)*(dbydx-dbxdy);

    Real3 B_eff = Bp+curlb;
    Real B_effpar = dot(B_eff,unitB);
    
    Real vpardot = -(mu/(mass*B_effpar))*dot(B_eff, dmagB);

    Real3 vpar = magvpar*unitB;

    Real3 Xdot = magvpar*(B_eff)/B_effpar - (mu/charge)*cross(dmagB, (unitB/B_effpar));

    Tuple<Real, Real3> advance = {vpardot, Xdot};

    return advance;
}

void 
derivative_calculation
(
    int i, int j, int k,
    const Array4<Real>& Bx, 
    const Array4<Real>& By, 
    const Array4<Real>& Bz,
    const Array4<Real>& gradBx, 
    const Array4<Real>& gradBy, 
    const Array4<Real>& gradBz,
    const Array4<Real>& jacobian, 
    const GpuArray<Real, AMREX_SPACEDIM> problo,
    const GpuArray<Real, AMREX_SPACEDIM> dx
)
{
    Real magBx_plusi = sqrt(pow(Bx(i+1,j,k),2)+pow(By(i+1,j,k),2)+pow(Bz(i+1,j,k),2));
    Real magBx_minusi = sqrt(pow(Bx(i-1,j,k),2)+pow(By(i-1,j,k),2)+pow(Bz(i-1,j,k),2));
    
    Real magBy_plusi = sqrt(pow(Bx(i,j+1,k),2)+pow(By(i,j+1,k),2)+pow(Bz(i,j+1,k),2));
    Real magBy_minusi = sqrt(pow(Bx(i,j-1,k),2)+pow(By(i,j-1,k),2)+pow(Bz(i,j-1,k),2));

    Real magBz_plusi = sqrt(pow(Bx(i,j,k+1),2)+pow(By(i,j,k+1),2)+pow(Bz(i,j,k+1),2));
    Real magBz_minusi = sqrt(pow(Bx(i,j,k-1),2)+pow(By(i,j,k-1),2)+pow(Bz(i,j,k-1),2));

    gradBx(i,j,k) = (magBx_plusi-magBx_minusi)/(2*dx[0]);
    gradBy(i,j,k) = (magBy_plusi-magBy_minusi)/(2*dx[1]);
    gradBz(i,j,k) = (magBz_plusi-magBz_minusi)/(2*dx[2]);

    jacobian(i,j,k,0) = (Bx(i+1,j,k)-Bx(i-1,j,k))/(2*dx[0]);
    jacobian(i,j,k,1) = (Bx(i,j+1,k)-Bx(i,j-1,k))/(2*dx[1]);
    jacobian(i,j,k,2) = (Bx(i,j,k+1)-Bx(i,j,k-1))/(2*dx[2]);

    jacobian(i,j,k,3) = (By(i+1,j,k)-By(i-1,j,k))/(2*dx[0]);
    jacobian(i,j,k,4) = (By(i,j+1,k)-By(i,j-1,k))/(2*dx[1]);
    jacobian(i,j,k,5) = (By(i,j,k+1)-By(i,j,k-1))/(2*dx[2]);

    jacobian(i,j,k,6) = (Bz(i+1,j,k)-Bz(i-1,j,k))/(2*dx[0]);
    jacobian(i,j,k,7) = (Bz(i,j+1,k)-Bz(i,j-1,k))/(2*dx[1]);
    jacobian(i,j,k,8) = (Bz(i,j,k+1)-Bz(i,j,k-1))/(2*dx[2]);
}