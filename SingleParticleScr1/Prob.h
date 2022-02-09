#ifndef PROB_H_
#define PROB_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

#include "utils.h"

using namespace amrex;


void init_values(Real3& Bp, const Vector<Real>& init_velocity, Real3& pos,
                 Real& magvpar, 
                 Real& mu, Real q, Real m);

Tuple<Real, Real3>
K_calc
(
    Real3 Bp, 
    Real3 dmagB,
    Vector<Real> jacobianp,
    const Real charge, const Real mass,
    const Real& magvpar, Real& mu
);

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
);

#endif