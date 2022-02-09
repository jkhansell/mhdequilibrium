#include <AMReX_PlotFileUtil.H>

using namespace amrex;

void WritePlotFile
(
    const Array<const MultiFab*, AMREX_SPACEDIM>& B,
    const Array<const MultiFab*, AMREX_SPACEDIM>& gradB,
    const MultiFab& jacobian, const Geometry& geom, 
    Real time, int step
);
