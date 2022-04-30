#include <eigen3/Eigen/Dense>

#include "GCM.h"
using namespace Eigen;

void init_GCM_calc
(
    const Matrix<double, 6,1>& IC_0,
    Matrix<double, 5,1>& solvect,
    const double& mass, 
    const double& charge,
    Field Bfield
)
{   
    Vector3d pos = {IC_0(0), IC_0(1), IC_0(2)};
    Vector3d vel = {IC_0(3), IC_0(4), IC_0(5)};
    Vector3d B = {0.0,0.0,0.0};
    Matrix3d jac = Matrix3d::Zero();
    Bfield.trilin_interp(pos, B, jac);
    double Bmag = B.norm();
    Vector3d bhat = B/Bmag;
    double vparmag = vel.dot(bhat);
    Vector3d vpar = vparmag*bhat;
    Vector3d vperp = vel - vpar; 
    Vector3d rho = bhat.cross(vperp);
    double Omega = charge*Bmag/mass; 
    rho *= -1/Omega; 
    Vector3d X = pos + rho; 
    double mu = mass*vperp.squaredNorm()/(2*Bmag);
    solvect(0) = vparmag; 
    solvect(1) = X(0); 
    solvect(2) = X(1); 
    solvect(3) = X(2); 
    solvect(4) = mu;  

}

Matrix<double, 5, 1> GCM_step
(
    Matrix<double, 5,1>& solvect,
    double mass, 
    double charge,
    Field Bfield
)
{
    double mu, udot, B_effpar, Bmag; 
    Vector3d pos, Xdot, B, B_eff, bhat, gradB, curlB;
    Matrix<double, 5, 1> derivs;
    Matrix3d jac = Matrix3d::Zero();

    pos = {solvect(1), solvect(2), solvect(3)};
    //std::cout << pos << "\n";

    Bfield.trilin_interp(pos, B, jac);

    Matrix3d realjac; 
    double R = sqrt(pow(pos(0),2)+pow(pos(1),2)); 
    double phi = atan2(pos(1),pos(0));
    if (phi < 0){phi += 2*PI;}

    realjac(0,0) = jac(0,0); realjac(0,1) = (jac(0,1)/R)-B(1)/R; realjac(0,2) = jac(0,2);
    realjac(1,0) = jac(1,0); realjac(1,1) = (jac(1,1)/R)+B(0)/R; realjac(1,2) = jac(1,2);
    realjac(2,0) = jac(2,0); realjac(2,1) =  jac(2,1)/R;         realjac(2,2) = jac(2,2);

    Matrix3d cart_jac = Bfield.transformjac(realjac, phi);
    Vector3d B_cart = Bfield.transformvect(B, phi);

    Bmag = B_cart.norm();
    //std::cout << Bmag << "\n";

    bhat = B_cart/Bmag;

    curlB = {cart_jac(2,1)-cart_jac(1,2),
             cart_jac(0,2)-cart_jac(2,0),
             cart_jac(1,0)-cart_jac(0,1)};

    gradB = cart_jac.transpose()*B_cart;
    gradB/=Bmag;
    /*
    gradB(0) = (B_cart(0)*cart_jac(0,0)+B_cart(1)*cart_jac(1,0)+B_cart(2)*cart_jac(2,0))/Bmag; 
    gradB(1) = (B_cart(0)*cart_jac(0,1)+B_cart(1)*cart_jac(1,1)+B_cart(2)*cart_jac(2,1))/Bmag; 
    gradB(2) = (B_cart(0)*cart_jac(0,2)+B_cart(1)*cart_jac(1,2)+B_cart(2)*cart_jac(2,2))/Bmag; 
    */
    B_eff = B + (mass/charge)*curlB;
    
    B_effpar = B_eff.dot(bhat);
 
    udot = -B_eff.dot(gradB)*mu/(mass*B_effpar);

    //std::cout << udot << "\n";


    Xdot = solvect[0]*B_eff/B_effpar-(solvect[4]/charge)*gradB.cross(bhat/B_effpar);

    derivs(0) = udot;
    derivs(1) = Xdot(0); 
    derivs(2) = Xdot(1); 
    derivs(3) = Xdot(2); 
    derivs(4) = 0; 

    return derivs;
}

bool check_particle_domain
(
    const Matrix<double, 5, 1>& solvect
)
{
    double R = sqrt(pow(solvect(1),2)+pow(solvect(2),2));
    double r = sqrt(pow(R-R0,2)+pow(solvect(3),2));
    if (r > a-0.02){return false;}
    else{return true;}
}