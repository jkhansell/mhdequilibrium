// standard libraries


// AMReX libraries

#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_StructOfArrays.H>
#include <AMReX_Particles.H>


// local modules

#include "solctra_utils.h"
#include "AMR_utils.h"
#include "Fields.h"
#include "Particles.h"
#include "constants.h"


namespace ParticleMeshFuncs
{
    void Scatter(const Array4<Real>& field, const Vector<Real>& l, Real value);

    void GatherFields(const EMParticleContainer::ParticleType& p,
                      const Array4<const Real>& Exarrk, const Array4<const Real>& Exarrk_05, 
                      const Array4<const Real>& Eyarrk, const Array4<const Real>& Eyarrk_05,
                      const Array4<const Real>& Ezarrk, const Array4<const Real>& Ezarrk_05, 
                      const Array4<const Real>& Bxarrk, const Array4<const Real>& Bxarrk_05, 
                      const Array4<const Real>& Byarrk, const Array4<const Real>& Byarrk_05,
                      const Array4<const Real>& Bzarrk, const Array4<const Real>& Bzarrk_05,
                      const Array4<const Real>& gradBxarrk, const Array4<const Real>& gradBxarrk_05, 
                      const Array4<const Real>& gradByarrk, const Array4<const Real>& gradByarrk_05,
                      const Array4<const Real>& gradBzarrk, const Array4<const Real>& gradBzarrk_05, 
                      Vector<Real>& Efieldk, Vector<Real>& Efieldk_05,
                      Vector<Real>& Bfieldk, Vector<Real>& Bfieldk_05,
                      Vector<Real>& gradBfieldk, Vector<Real>& gradBfieldk_05,
                      GpuArray<Real, AMREX_SPACEDIM> plo,
                      GpuArray<Real, AMREX_SPACEDIM> dxi);  
    
}