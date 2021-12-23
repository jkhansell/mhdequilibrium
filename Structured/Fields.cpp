// standard libraries


// AMReX libraries

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_Print.H>
#include <AMReX_MLMG.H>
#include <AMReX_BoxArray.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_InterpBndryData.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>


// local modules

#include "solctra_utils.h"
#include "Fields.h"

// namespaces
using namespace amrex;




void get_initial_fields(const GlobalData& Coils, Coil* rmi, Coil* rmf, 
                        const Vector<Geometry>& geom,
                        const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& BField,
                        int n_cell, Real R0, Real dh)
{
    const int num_levels = BField.size();
    for (int lev = 0; lev < num_levels; ++lev)
    {
        
        const Geometry& gm = geom[lev];    
        const auto dx = gm.CellSizeArray();
        
        for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
        {
            Print() << dim; 
            for (MFIter mfi(*BField[lev][dim]); mfi.isValid(); ++mfi)
            {
                const Array4<Real>& Bdim = (*BField[lev][dim]).array(mfi);

                const Box& bx = mfi.validbox();
                const auto lo = lbound(bx); 
                const auto hi = ubound(bx); 

                for (int k = lo.z; k<= hi.z ;++k){
                    for (int j = lo.y; j<= hi.y ;++j){
                        for (int i = lo.x; i<= hi.x ;++i){
                            
                            Real x = R0 - dh + i*dx[0];
                            Real y = -dh + j*dx[1];
                            Real z = -dh + k*dx[2];

                            cartesian p0 = {x, y, z};
                            cartesian magfield = magnetic_field(rmi, rmf, Coils, p0);

                            if (dim == 0)
                            {
                                Bdim(i,j,k) = Real(magfield.x);
                            }                           
                            else if (dim == 1)
                            {
                                Bdim(i,j,k) = Real(magfield.y);
                            }
                            else
                            {
                                Bdim(i,j,k) = Real(magfield.z);
                            }
                        }
                    }
                }
            }
        }
    }    
}


void NodalFieldSolver::init_sigma()
{
    sigma.resize(num_levels);
    for (int lev = 0; lev < num_levels; ++lev)
    {
        sigma[lev].define(NFSgrids[lev], NFSdmap[lev], 1, 0);
        sigma[lev].setVal(1.0);
    }
}

void NodalFieldSolver::solveforPhi(const Vector<const MultiFab *>& rhs, 
                                   const Vector<MultiFab *>& phi)
{   
    int num_levels = rhs.size();

    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setSemicoarsening(semicoarsening);
    info.setMaxCoarseningLevel(max_coarsening_level);
    info.setMaxSemicoarseningLevel(max_semicoarsening_level);


    MLNodeLaplacian linop(NFSgeom, NFSgrids, NFSdmap, info); 

    linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)},
                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)});

    for (int lev = 0; lev <num_levels; ++lev)
    {
        linop.setSigma(lev, sigma[lev]);
    }
    
    MLMG mlmg(linop); 
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);

    mlmg.solve(phi, rhs, reltol, abstol);
    
    for (int lev = 0; lev < num_levels; ++lev) 
    {
        const Geometry& gm = NFSgeom[lev];
        phi[lev]->FillBoundary();
    }
}

void NodalFieldSolver::calculateEField(const Vector<Array<MultiFab *, AMREX_SPACEDIM>>& Efield,
                                   const Vector<const MultiFab *>& phi,
                                   const YeeGrid& EMYeeGrid)
{   


    int num_levels = phi.size();
    int nComp = 1;
    for (int lev = 0; lev < num_levels; ++lev)
    {
        const Geometry& geom = NFSgeom[lev];
        const auto& dxi = geom.InvCellSizeArray(); 

        for (MFIter mfi(*phi[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& Exbox = mfi.tilebox(EMYeeGrid.Efield[0]); 
            const Box& Eybox = mfi.tilebox(EMYeeGrid.Efield[1]);
            const Box& Ezbox = mfi.tilebox(EMYeeGrid.Efield[2]);

            const auto& iPhi = (*phi[lev]).array(mfi);

            const auto& efieldx = Efield[lev][0]->array(mfi); 
            const auto& efieldy = Efield[lev][1]->array(mfi); 
            const auto& efieldz = Efield[lev][2]->array(mfi); 

            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(Exbox, nComp, i, j, k, n,
            {
                efieldx(i,j,k,n) = -dxi[0]*(iPhi(i,j,k,n)-iPhi(i-1,j,k,n));
            });

            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(Eybox, nComp, i, j, k, n,
            {
                efieldy(i,j,k,n) = -dxi[1]*(iPhi(i,j,k,n)-iPhi(i,j-1,k,n));
            });

            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(Exbox, nComp, i, j, k, n,
            {
                efieldz(i,j,k,n) = -dxi[2]*(iPhi(i,j,k,n)-iPhi(i,j,k-1,n));
            });
        }
    }
}

/*
void NodalFieldSolver::evolveBfield(const Vector<Array<const MultiFab*, AMREX_SPACEDIM>>& E,
                                    const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& B,
                                    const YeeGrid& EMYeeGrid, const Vector<Geometry>& geom,
                                    Real dt)
{
    int num_levels = B.size();

    for (int lev = 0; lev < num_levels; ++lev)
    {
        const auto& dxi = geom[lev].InvCellSize();

        const Real dt_dx = dt*dxi[0]; 
        const Real dt_dy = dt*dxi[1]; 
        const Real dt_dz = dt*dxi[2]; 

        for (MFIter mfi(*B[lev][0]); mfi.isValid(); ++mfi)
        {
            const Box& tempbx = mfi.tilebox(EMYeeGrid.Bfield[0]);
            const Box& tempby = mfi.tilebox(EMYeeGrid.Bfield[1]);
            const Box& tempbz = mfi.tilebox(EMYeeGrid.Bfield[2]);

            auto const& iBx = (*B[lev][0]).array(mfi);
            auto const& iBy = (*B[lev][0]).array(mfi);
            auto const& iBz = (*B[lev][0]).array(mfi);
            auto const& iBx = (*E[lev][0]).array(mfi);
            auto const& iBy = (*E[lev][0]).array(mfi);
            auto const& iBz = (*E[lev][0]).array(mfi);

        }
    }
}
*/