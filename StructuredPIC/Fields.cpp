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
#include "constants.h"

// namespaces
using namespace amrex;
using namespace PhysConst;

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


void FieldEvolver::evolve_Bfield(const Vector<Array<const MultiFab*, AMREX_SPACEDIM>>& E,
                                    const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& B,
                                    const YeeGrid& EMYeeGrid, const Vector<Geometry>& geom,
                                    Vector<Real>& dt)
    {
    int num_levels = B.size();

    for (int lev = 0; lev < num_levels; ++lev)
    {
        const auto& dxi = geom[lev].InvCellSize();

        const Real dt_dx = dt[lev]*dxi[0]; 
        const Real dt_dy = dt[lev]*dxi[1]; 
        const Real dt_dz = dt[lev]*dxi[2]; 

        for (MFIter mfi(*B[lev][0]); mfi.isValid(); ++mfi)
        {
            const Box& tempbx = mfi.tilebox(EMYeeGrid.Bfield[0]);
            const Box& tempby = mfi.tilebox(EMYeeGrid.Bfield[1]);
            const Box& tempbz = mfi.tilebox(EMYeeGrid.Bfield[2]);

            auto const& iBx = (*B[lev][0]).array(mfi);
            auto const& iBy = (*B[lev][1]).array(mfi);
            auto const& iBz = (*B[lev][2]).array(mfi);
            auto const& iEx = (*E[lev][0]).array(mfi);
            auto const& iEy = (*E[lev][1]).array(mfi);
            auto const& iEz = (*E[lev][2]).array(mfi);

            AMREX_PARALLEL_FOR_3D(tempbx, i, j, k,
            {
                iBx(i,j,k) = iBx(i,j,k) - dt_dy*(iEz(i,j+1,k)-iEz(i,j,k))
                                        + dt_dz*(iEy(i,j,k+1)-iEy(i,j,k));
            });

            AMREX_PARALLEL_FOR_3D(tempby, i, j, k,
            {
                iBy(i,j,k) = iBy(i,j,k) - dt_dz*(iEx(i,j,k+1)-iEx(i,j,k))
                                        + dt_dx*(iEz(i+1,j,k)-iEz(i,j,k));
            });
            
            AMREX_PARALLEL_FOR_3D(tempbz, i, j, k,
            {
                iBz(i,j,k) = iBz(i,j,k) - dt_dx*(iEy(i+1,j,k)-iEy(i,j,k))
                                        + dt_dy*(iEx(i,j+1,k)-iEx(i,j,k));
            });

        }
    }
}

void FieldEvolver::evolve_Efield(const Vector<Array<const MultiFab*, AMREX_SPACEDIM>>& B_k05,
                                const Vector<Array<const MultiFab*, AMREX_SPACEDIM>>& J_k05,
                                const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& E_k,
                                const YeeGrid& EMYeeGrid, const Vector<Geometry>& geom,
                                Vector<Real>& dt)
{
    int num_levels = E_k.size();

    for (int lev = 0; lev < num_levels; ++lev)
    {
        const auto& dxi = geom[lev].InvCellSize();

        const Real jcoef = 1/ep0;
        const Real CBcoef = c*c;

        for (MFIter mfi(*B_k05[lev][0]); mfi.isValid(); ++mfi)
        {
            const Box& tempEx = mfi.tilebox(EMYeeGrid.Efield[0]);
            const Box& tempEy = mfi.tilebox(EMYeeGrid.Efield[1]);
            const Box& tempEz = mfi.tilebox(EMYeeGrid.Efield[2]);
 
            auto const& iBx = (*B_k05[lev][0]).array(mfi);
            auto const& iBy = (*B_k05[lev][1]).array(mfi);
            auto const& iBz = (*B_k05[lev][2]).array(mfi);
    
            auto const& ijx = (*J_k05[lev][0]).array(mfi);
            auto const& ijy = (*J_k05[lev][1]).array(mfi);
            auto const& ijz = (*J_k05[lev][2]).array(mfi);
            
            auto const& iEx = (*E_k[lev][0]).array(mfi);
            auto const& iEy = (*E_k[lev][1]).array(mfi);
            auto const& iEz = (*E_k[lev][2]).array(mfi);

            AMREX_PARALLEL_FOR_3D(tempEx, i, j, k,
            {
                iEx(i,j,k) = iEx(i,j,k) + dt[lev]*(CBcoef*((iBz(i,j+1,k)-iBz(i,j,k))*dxi[1]
                                        -   (iBy(i,j,k+1)-iBy(i,j,k))*dxi[2]) - jcoef*ijx(i,j,k)); 
            });

            AMREX_PARALLEL_FOR_3D(tempEy, i, j, k,
            {
                iEy(i,j,k) = iEy(i,j,k) + dt[lev]*(CBcoef*((iBx(i,j,k+1)-iBx(i,j,k))*dxi[2]
                                        -   (iBz(i+1,j,k)-iBz(i,j,k))*dxi[0]) - jcoef*ijy(i,j,k));
            });
            
            AMREX_PARALLEL_FOR_3D(tempEz, i, j, k,
            {
                iEz(i,j,k) = iEz(i,j,k) + dt[lev]*(CBcoef*((iBy(i+1,j,k)-iBy(i,j,k))*dxi[0]
                                        -  (iBx(i,j+1,k)-iBx(i,j,k))*dxi[1]) - jcoef*ijz(i,j,k));
            });

        }
    }
}