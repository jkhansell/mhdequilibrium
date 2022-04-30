
#include <AMReX_REAL.H>
#include "Prob.H"
#include "AmrParticles.H"

using namespace amrex;



void init_values(Vector<Real>& Bp,
                 Vector<Real>& Xdot,
                 Real& magvpar,
                 Real& mu, Real charge, Real mass)
{
    Real magBp = sqrt(Dot(Bp, Bp));
    Vector<Real> unitB = Bp/magBp;
    magvpar = Dot(Xdot, unitB);
    Vector<Real> vpar = magvpar*unitB;
    Vector<Real> vperp = Xdot - vpar;
    Real magvperp = sqrt(Dot(vperp,vperp));
    mu = (mass/(2*magBp))*pow(magvperp,2);
}

Tuple<Real, Vector<Real>>
advance_calc
(
    Vector<Real> Bp, 
    Vector<Real> dmagB,
    Vector<Real> jacobianp,
    Real charge, Real mass,
    Real& magvpar, Real& mu
)
{    
    Real magBp = sqrt(Dot(Bp, Bp));
    
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

    Vector<Real> unitgradB(3);
    unitgradB[0] = dbxdx;
    unitgradB[1] = dbydy;
    unitgradB[2] = dbzdz;

    Vector<Real> unitB = Bp/magBp;
    
    Vector<Real> curlb(3);
    curlb[0] = (mass*magvpar/charge)*(dbzdy-dbydz);
    curlb[1] = (mass*magvpar/charge)*(dbxdz-dbzdx);
    curlb[2] = (mass*magvpar/charge)*(dbydx-dbxdy);

    Vector<Real> B_eff = Bp+curlb;
    Real B_effpar = Dot(B_eff,unitB);
    
    Real vpardot = Dot(-(mu/(mass*B_effpar))*B_eff, dmagB);

    Vector<Real> vpar = magvpar*unitB;

    Vector<Real> leftvect = (mass/(charge*B_effpar))*unitB;   

    Vector<Real> gradpar = {unitB[0]*dbxdx+unitB[1]*dbxdy+unitB[2]*dbxdz,
                            unitB[0]*dbydx+unitB[1]*dbydy+unitB[2]*dbydz,
                            unitB[0]*dbzdx+unitB[1]*dbzdy+unitB[2]*dbzdz};

    Vector<Real> rightvect = (mu/mass)*dmagB + magvpar*magvpar*gradpar;
    Vector<Real> Xdot = vpar + Cross(leftvect,rightvect);
    
    Tuple<Real, Vector<Real>> advance = {vpardot, Xdot};

    return advance;
}

void 
init_fields
(
    
    const Box& bx,
    Array4<Real> Bx, 
    Array4<Real> By, 
    Array4<Real> Bz,
    const GpuArray<Real, AMREX_SPACEDIM> problo,
    const GpuArray<Real, AMREX_SPACEDIM> dx,
    const GlobalData& Coils, Coil* rmi, Coil* rmf
)  
{
    const auto lo = lbound(bx);
    const auto hi = ubound(bx);

    for (int k = lo.z; k <= hi.z; ++k){
        Real z = problo[2] + k*dx[2];
        for (int j = lo.y; j <= hi.y; ++j){
            Real y = problo[1] + j*dx[1];
            for (int i = lo.x; i <= hi.x; ++i){
                Real x = problo[0] + i*dx[0];
                cartesian p0 = {x, y, z};
                cartesian magfield = magnetic_field(rmi, rmf, Coils, p0);
                
                Bx(i,j,k) = magfield.x;
                By(i,j,k) = magfield.y; 
                Bz(i,j,k) = magfield.z;
            }
        }
    }


}

void 
derivative_calculation
(
    int i, int j, int k,
    const Array4<Real> Bx, 
    const Array4<Real> By, 
    const Array4<Real> Bz,
    Array4<Real> gradBx, 
    Array4<Real> gradBy, 
    Array4<Real> gradBz,
    Array4<Real> jacobian, 
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

void 
init_fields_zero
(
    
    int i, int j, int k,
    Array4<Real> Bx, 
    Array4<Real> By, 
    Array4<Real> Bz
)  
{
    Bx(i,j,k) = 0.0;
    By(i,j,k) = 0.0;
    Bz(i,j,k) = 0.0;

}