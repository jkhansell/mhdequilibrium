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
#include "ParticleMesh.h"
#include "constants.h"


// namespaces
using namespace amrex;
using namespace PhysConst;
// Functions


void ParticleMeshFuncs::Scatter(const Array4<Real>& field, const Vector<Real>& l, Real value)
{
    
    Real lx = l[0];
    Real ly = l[1]; 
    Real lz = l[2]; 

    int i = (int)lx;
    int j = (int)ly;
    int k = (int)lz;

    Real di = lx - (Real)i;
    Real dj = ly - (Real)j;
    Real dk = lz - (Real)k;

    Real sx[] = {Real(1.0)-di, di};
    Real sy[] = {Real(1.0)-dj, dj};
    Real sz[] = {Real(1.0)-dk, dk};

    for (int kk = 0; kk<=1; ++kk){ 
        for (int jj = 0; jj<=1; ++jj){
            for (int ii = 0; ii<=1; ++ii){
                amrex::Gpu::Atomic::AddNoRet(&field(i, j, k), sx[ii]*sy[jj]*sz[kk]*value);
            }
        }
    } 
    // Print() << "scatterPM";
}

