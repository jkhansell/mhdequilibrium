// standard libraries


// AMReX libraries

#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_StructOfArrays.H>
#include <AMReX_Particles.H>


// local modules

#include "AMR_utils.h"
#include "Fields.h"
#include "Particles.h"
#include "particle_update.h"
#include "constants.h"


void
ParticleUpdate::PushLeapFrog
(
    const EMParticleContainer::ParticleType& p,
    Real& uxp, Real& uyp, Real& uzp,
    const Vector<Real>& Efieldk,
    const Vector<Real>& Bfieldk, 
    const Vector<Real>& gradBfieldk, 
    const Vector<Real>& Efieldk_05,
    const Vector<Real>& Bfieldk_05, 
    const Vector<Real>& gradBfieldk_05,
    const Real& dt, const Real& q, const Real& m
)
{
    Real magB = sqrt(Bfieldk[0]*Bfieldk[0]+Bfieldk[1]*Bfieldk[1]+Bfieldk[2]*Bfieldk[2]);
    
    Vector<Real> unitB(3);
    Vector<Real> unitgradB(3);
    Vector<Real> vpar(3);
    Vector<Real> vperp(3);

    Real omega = q*magB/m; 
    
    unitB[0] = Bfieldk[0]/magB;
    unitB[1] = Bfieldk[1]/magB;
    unitB[2] = Bfieldk[2]/magB;

    unitgradB[0] = gradBfieldk[0]/magB;
    unitgradB[1] = gradBfieldk[1]/magB;
    unitgradB[2] = gradBfieldk[2]/magB;
    
    Real magvpar = uxp*unitB[0]+uyp*unitB[1]+uzp*unitB[2]; 
    vpar[0] = magvpar*unitB[0];
    vpar[1] = magvpar*unitB[1];
    vpar[2] = magvpar*unitB[2];

    vperp[0] = uxp - vpar[0];
    vperp[1] = uyp - vpar[1];
    vperp[2] = uzp - vpar[2];

       

}

