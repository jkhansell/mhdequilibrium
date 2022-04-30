#ifndef TAGGING_H
#define TAGGING_H

#include <AMReX_Array4.H>
#include <AMReX_Print.H>
using namespace amrex;

void
state_error(int i, int j, int k,
            amrex::Array4<char> const& tag,
            char tagval, Real R0, Real dh,
            const GpuArray<Real, AMREX_SPACEDIM>& problo,
            const GpuArray<Real, AMREX_SPACEDIM>& dx)
{
    Real x = problo[0] + i*dx[0];
    Real y = problo[1] + j*dx[1];
    Real R = sqrt(pow(x,2)+pow(y,2));


    if (R > R0-dh && R < R0+dh)
    {
        tag(i,j,k) = tagval;
    }
}

#endif
