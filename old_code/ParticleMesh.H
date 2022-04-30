#include <AMReX_MultiFabUtil.H>

using namespace amrex;

namespace ParticleMeshFuncs
{
    void gather_fields
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
    );
}
