#ifndef INITBFIELD_H
#define INITBFIELD_H


#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>


using namespace amrex;
#define pi 3.1415926535897932384

void initBfield
(
    int i, int j, int k,
    const Array4<Real>& iBx,
    const Array4<Real>& iBy,
    const Array4<Real>& iBz,
    const GpuArray<Real, AMREX_SPACEDIM> problo,
    const GpuArray<Real, AMREX_SPACEDIM> probhi,
    const GpuArray<Real, AMREX_SPACEDIM> dx
)
{
    Real x = problo[0] + i*dx[0];
    Real y = problo[1] + j*dx[1];
    Real z = problo[2] + k*dx[2];

    Real Lx = probhi[0]-problo[0];
    Real Ly = probhi[1]-problo[1];
    Real Lz = probhi[2]-problo[2];

    iBx(i,j,k) = sin(y*z);
    iBy(i,j,k) = sin(x*y);
    iBz(i,j,k) = 0.0;//x*(z/Lz);

}

#endif