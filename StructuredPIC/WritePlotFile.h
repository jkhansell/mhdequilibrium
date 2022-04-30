
#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_Vector.H>


#include "Particles.h"
#include "Fields.h"

using namespace amrex;

void WriteScalarData(const Vector<const MultiFab*>& rhs,
                     const Vector<const MultiFab*>& phi,
                     const Vector<Geometry>& geom, Real time, int step);//, Vector<std::unique_ptr<EMParticleContainer>> particles);


void WriteVectorData(const Vector<Array<const MultiFab*, AMREX_SPACEDIM>>& BField,
                     const Vector<Array<const MultiFab*, AMREX_SPACEDIM>>& EField,
                     const Vector<Array<const MultiFab*, AMREX_SPACEDIM>>& jField,
                     const Vector<Geometry>& geom, Real time, int step);

