#include <AMReX_iMultiFab.H>
#include "Fields.h"

void build_mask
(
    int i, int j, int k,
    Array4<int> mask_arr,
    const Real R0, const Real a,
    const GpuArray<Real, AMREX_SPACEDIM>& plo,
    const GpuArray<Real, AMREX_SPACEDIM>& dx
)
{

    Real x = plo[0] + i*dx[0];
    Real y = plo[1] + j*dx[1];
    Real z = plo[2] + k*dx[2];
    Real R = sqrt(pow(x,2)+pow(y,2));
    Real r = sqrt(pow(R-R0,2)+pow(z,2));

    if (r < a)
    {
        mask_arr(i,j,k) = 1;
    }
    else
    {
        mask_arr(i,j,k) = 0;
    }
}

std::unique_ptr<iMultiFab> get_masks
(
    const BoxArray& grids, 
    const DistributionMapping& dmap,
    const Geometry& geom,
    const Real R0, const Real a
)
{
    std::unique_ptr<iMultiFab> masks;
    BoxArray nba = grids;
    nba.surroundingNodes();
    masks.reset(new iMultiFab(nba, dmap, 1, 0));

    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    
    for (MFIter mfi(*masks); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox(); 
        auto mask_arr = masks->array(mfi); 
        ParallelFor(bx, 
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            build_mask(i,j,k, mask_arr, R0, a, plo, dx);
        });
    }

    return masks;

}

void get_initial_fields
(
    int i, int j, int k, 
    const GlobalData& Coils, Coil* rmi, Coil* rmf, 
    const GpuArray<Real, AMREX_SPACEDIM>& plo,
    const GpuArray<Real, AMREX_SPACEDIM>& dx,
    const Array4<Real>& iBx,
    const Array4<Real>& iBy,
    const Array4<Real>& iBz,
    const Array4<int>& imask
)
{

    if (imask(i,j,k) == 1)
    {
        Real x = plo[0] + i*dx[0];
        Real y = plo[1] + j*dx[1];
        Real z = plo[2] + k*dx[2];

        cartesian p0 = {x, y, z};
        cartesian magfield = magnetic_field(rmi, rmf, Coils, p0);

        iBx(i,j,k) = magfield.x; 
        iBy(i,j,k) = magfield.y; 
        iBz(i,j,k) = magfield.z; 
    }
    else
    {
        iBx(i,j,k) = 0.0; 
        iBy(i,j,k) = 0.0; 
        iBz(i,j,k) = 0.0; 
    }
}