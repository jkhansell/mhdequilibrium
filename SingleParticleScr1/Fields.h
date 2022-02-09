#ifndef FIELD_H
#define FIELD_H

#include <AMReX_Array.H>
#include <AMReX_MultiFab.H>

#include "solctra_utils.h"


using namespace amrex;
#define pi 3.1415926535897932384

void get_initial_fields
(
    int i, int j, int k, 
    const GlobalData& Coils, Coil* rmi, Coil* rmf, 
    const GpuArray<Real, AMREX_SPACEDIM>& plo,
    const GpuArray<Real, AMREX_SPACEDIM>& dx,
    const Array4<Real>& iBx,
    const Array4<Real>& iBy,
    const Array4<Real>& iBz,
    const Array4<int>& imask
);


/*
void initBfield
(
    int i, int j, int k,
    const Array4<Real>& iBx,
    const Array4<Real>& iBy,
    const Array4<Real>& iBz,
    const GpuArray<Real, AMREX_SPACEDIM> problo,
    const GpuArray<Real, AMREX_SPACEDIM> probhi,
    const GpuArray<Real, AMREX_SPACEDIM> dx,
    const Array4<int>& imask
)
{
    if (imask(i,j,k) == 1)
    {
        Real x = problo[0] + i*dx[0];
        Real y = problo[1] + j*dx[1];
        Real z = problo[2] + k*dx[2];

        Real Lx = probhi[0]-problo[0];
        Real Ly = probhi[1]-problo[1];
        Real Lz = probhi[2]-problo[2];

        iBx(i,j,k) = exp(-(z*z+y*y));
        iBy(i,j,k) = exp(-(z*z+x*x));
        iBz(i,j,k) = exp(-(x*x+y*y));
    }
    else
    {
        iBx(i,j,k) = 0.0;
        iBy(i,j,k) = 0.0;
        iBz(i,j,k) = 0.0;
    }

}
*/

std::unique_ptr<iMultiFab> get_masks
(
    const BoxArray& grids, 
    const DistributionMapping& dmap,
    const Geometry& geom,
    const Real R0, const Real a
);


void build_mask
(
    int i, int j, int k,
    Array4<int> mask_arr,
    const Real R0, const Real a,
    const GpuArray<Real, AMREX_SPACEDIM>& plo,
    const GpuArray<Real, AMREX_SPACEDIM>& dx
);


#endif