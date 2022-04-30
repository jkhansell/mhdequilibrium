#ifndef PROB_H_
#define PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

#include "solctra_utils.H"

using namespace amrex;

Tuple<Real, Vector<Real>> advance_calc(Vector<Real> Bp, 
                    Vector<Real> dmagBp,
                    Vector<Real> jacobianp,
                    Real charge, Real mass,
                    Real& magvpar, Real& mu);

void init_values(Vector<Real>& Bp, Vector<Real>& Xdot, 
            Real& magvpar, Real& mu, Real charge, Real mass);

void init_fields
(    
    const Box& bx,
    Array4<Real> Bx, 
    Array4<Real> By, 
    Array4<Real> Bz,
    const GpuArray<Real, AMREX_SPACEDIM> problo,
    const GpuArray<Real, AMREX_SPACEDIM> dx,
    const GlobalData& Coils, Coil* rmi, Coil* rmf
);

void init_fields_zero
(    
    int i, int j, int k,
    Array4<Real> Bx, 
    Array4<Real> By, 
    Array4<Real> Bz
);


void derivative_calculation
(    int i, int j, int k,
    const Array4<Real> Bx, 
    const Array4<Real> By, 
    const Array4<Real> Bz,
    Array4<Real> gradBx, 
    Array4<Real> gradBy, 
    Array4<Real> gradBz,
    Array4<Real> jacobian, 
    const GpuArray<Real, AMREX_SPACEDIM> problo,
    const GpuArray<Real, AMREX_SPACEDIM> dx);


inline Vector<Real> operator + (const Vector<Real>& a, const Vector<Real>& b)
{
    Vector<Real> c(3);
    for (int i = 0; i < c.size(); ++i)
    {
        c[i] = a[i]+b[i];
    }
    return c;
}

inline Vector<Real> operator - (const Vector<Real>& a, const Vector<Real>& b)
{
    Vector<Real> c(3);
    for (int i = 0; i < c.size(); ++i)
    {
        c[i] = a[i]-b[i];
    }
    return c;
}

inline Vector<Real> operator *(const Real& a, const Vector<Real>& b)
{
    Vector<Real> c(3);
    for (int i = 0; i < c.size(); ++i)
    {
        c[i] = a*b[i];
    }
    return c;
}

inline Vector<Real> operator /(const Vector<Real>& a, const Real& b)
{
    Vector<Real> c(3);
    for (int i = 0; i < c.size(); ++i)
    {
        c[i] = a[i]/b;
    }
    return c;
} 

inline Real Dot(const Vector<Real>& a, const Vector<Real>& b)
{
    Real dot; 
    for (int i = 0; i < 3; ++i)
    {
        dot+=a[i]*b[i];
    }
    return dot;
}

inline Vector<Real> Cross(const Vector<Real>& a, const Vector<Real>& b)
{
    Vector<Real> cross(3);
    cross[0] = a[1]*b[2]-a[2]*b[1];
    cross[1] = a[2]*b[0]-a[0]*b[2];
    cross[2] = a[0]*b[1]-a[1]*b[0];
    return cross;
}

#endif


