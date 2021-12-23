#ifndef CONSTANTS_H
#define CONSTANTS_H 

#include <AMReX_REAL.H>
#include <math.h>

using namespace amrex;

namespace PhysConst
{
    static const amrex::Real c   = 299792458.;
    static const amrex::Real ep0 = 8.854187817e-12;
    static const amrex::Real mu0 = 1.2566370614359173e-06;
    static const amrex::Real q_e = 1.6021764620000001e-19;
    static const amrex::Real m_e = 9.10938291e-31;
    static const amrex::Real m_p = 1.6726231000000001e-27;
    static const amrex::Real kBJ = 1.380649e-23;
    static const amrex::Real kBeV = 8.617333262e-5;
    static const Real electron_th_vel0 = sqrt(kBJ*10/(m_e*kBeV));
    static const Real proton_th_vel0 = sqrt(kBJ*1/(m_p*kBeV));

};

namespace MathConstants
{
    static constexpr amrex::Real pi = 3.14159265358979323846;
};

#endif