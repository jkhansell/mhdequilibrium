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


namespace ParticleUpdate
{


    void PushLeapFrog(const EMParticleContainer::ParticleType& p,
                      Real& uxp, Real& uyp, Real& uzp,
                      const Vector<Real>& Efieldk,
                      const Vector<Real>& Bfieldk, 
                      const Vector<Real>& gradBfieldk,
                      const Vector<Real>& Efieldk_05,
                      const Vector<Real>& Bfieldk_05, 
                      const Vector<Real>& gradBfieldk_05,  
                      const Real& dt, const Real& q, const Real& m);

}