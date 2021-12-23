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


namespace ParticleMeshFuncs{
    void Scatter(const Array4<Real>& field, const Vector<Real>& l, Real value);
}