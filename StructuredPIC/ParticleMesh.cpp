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

void 
ParticleMeshFuncs::GatherFields(
    const EMParticleContainer::ParticleType& p,
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
    GpuArray<Real, AMREX_SPACEDIM> dxi)
{
    Real lx = (p.pos(0)-plo[0])*dxi[0];
    Real ly = (p.pos(1)-plo[1])*dxi[1];     // logical node-centered coordinates
    Real lz = (p.pos(2)-plo[2])*dxi[2];

    int i = floor(lx);
    int j = floor(ly);
    int k = floor(lz);

    /*    
                                    (face node node)
    logical coordinates for E-field (node face node)
                                    (node node face)
    */
    
    Real di = lx - i;
    Real dj = ly - j;
    Real dk = lz - k;

    Real sx[] = {1.- di, di};
    Real sy[] = {1.- dj, dj};
    Real sz[] = {1.- dk, dk};

    Real face_di = di - 0.5;
    Real face_dj = dj - 0.5;
    Real face_dk = dk - 0.5;

    Real face_sx[] = {1.- face_di, face_di};
    Real face_sy[] = {1.- face_dj, face_dj};
    Real face_sz[] = {1.- face_dk, face_dk};   


    // trilinear interpolation algorithm

    constexpr int ixmin = 0;
    constexpr int ixmax = 0;
    constexpr int iymin = 0;
    constexpr int iymax = 0;
    constexpr int izmin = 0;
    constexpr int izmax = 0;
/*
    Ex = 0.0;
    for (int kk = izmin; kk <= izmax+1; ++kk){
        for (int jj = iymin; jj <= iymax+1; ++jj){
            for (int ii = ixmin; ii <= ixmax+1; ++ii){
                Ex += face_sx[ii]*sy[jj]*sz[kk]*Exarr(i+ii,j+jj,k+kk);
            }
        }
    } 

    Ey = 0.0;
    for (int kk = izmin; kk <= izmax+1; ++kk){
        for (int jj = iymin; jj <= iymax+1; ++jj){
            for (int ii = ixmin; ii <= ixmax+1; ++ii){
                Ex += sx[ii]*face_sy[jj]*sz[kk]*Eyarr(i+ii,j+jj,k+kk);
            }
        }
    }

    Ez = 0.0;
    for (int kk = izmin; kk <= izmax+1; ++kk){
        for (int jj = iymin; jj <= iymax+1; ++jj){
            for (int ii = ixmin; ii <= ixmax+1; ++ii){
                Ex += sx[ii]*sy[jj]*face_sz[kk]*Ezarr(i+ii,j+jj,k+kk);
            }
        }
    }

    Bx = 0.0;
    for (int kk = izmin; kk <= izmax+1; ++kk){
        for (int jj = iymin; jj <= iymax+1; ++jj){
            for (int ii = ixmin; ii <= ixmax+1; ++ii){
                Bx += sx[ii]*face_sy[jj]*face_sz[kk]*Bxarr(i+ii,j+jj,k+kk);
            }
        }
    }

    By = 0.0;
    for (int kk = izmin; kk <= izmax+1; ++kk){
        for (int jj = iymin; jj <= iymax+1; ++jj){
            for (int ii = ixmin; ii <= ixmax+1; ++ii){
                By += face_sx[ii]*sy[jj]*face_sz[kk]*Byarr(i+ii,j+jj,k+kk);
            }
        }
    }

    Bz = 0.0;
    for (int kk = izmin; kk <= izmax+1; ++kk){
        for (int jj = iymin; jj <= iymax+1; ++jj){
            for (int ii = ixmin; ii <= ixmax+1; ++ii){
                Bz += face_sx[ii]*face_sy[jj]*sz[kk]*Bzarr(i+ii,j+jj,k+kk);
            }
        }
    }
*/ 
}