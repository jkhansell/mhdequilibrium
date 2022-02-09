#ifndef PARTICLEMESH_H
#define PARTICLEMESH_H


#include <AMReX_MultiFabUtil.H>
#include "utils.h"

using namespace amrex;

namespace ParticleMeshFuncs
{
    void gather_fields
    (
        const Real3& position, 
        const Array4<const Real>& Bx,
        const Array4<const Real>& By, 
        const Array4<const Real>& Bz,
        const Array4<const Real>& gradBx,
        const Array4<const Real>& gradBy, 
        const Array4<const Real>& gradBz,
        const Array4<const Real>& Jacobian,
        Real3& Bp, 
        Real3& gradBp, 
        Vector<Real>& jacobianp,
        const GpuArray<Real, AMREX_SPACEDIM>& plo,
        const GpuArray<Real, AMREX_SPACEDIM>& dxi
    );
}

#endif
