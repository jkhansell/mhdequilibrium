// standard libraries


// AMReX libraries

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Particles.H>

// local modules

#include "Particles.h"
#include "AMR_utils.h"
#include "constants.h"
#include "ParticleMesh.h"

// namespaces
using namespace std;
using namespace amrex;
using namespace PhysConst;

// Functions

EMParticleContainer::
EMParticleContainer(const Vector<Geometry>            & geom,
                    const Vector<DistributionMapping> & dmap,
                    const Vector<BoxArray>            & ba,
                    const amrex::Vector<int>          & rr,
                    const int                          species_id,
                    const amrex::Real                  charge,
                    const amrex::Real                  mass, 
                    const amrex::Real                  mpw)
    : ParticleContainer<0, 0, PIdx::nattribs, 0> (geom, dmap, ba, rr),
    m_species_id(species_id), m_charge(charge), m_mass(mass), m_mpw(mpw), 
    num_levels(geom.size()), grids(ba), em_dmap(dmap)
{}

void EMParticleContainer::InitParticles(const IntVect& nppc,
                                        const Real th_mom_std,
                                        const RealBox& bounds)
{
    BL_PROFILE("PartContainer::InitParticles");
    const int lev = 0;
    const auto dx = Geom(lev).CellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();
    const int num_ppc = AMREX_D_TERM(nppc[0],*nppc[1],*nppc[2]);

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box = mfi.tilebox(); 
        
        const auto lo = lbound(tile_box); 
        const auto hi = ubound(tile_box);

        auto& particles = GetParticles(lev); 
        auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        
        for (int i = 0; i < num_ppc; ++i)
        {
            //Print() << "20"; 
            Real x = get_gaussian_random_number((bounds.lo(0) + (bounds.hi(0)-bounds.lo(0))/2), (bounds.hi(0)-bounds.lo(0))/2);
            Real y = get_gaussian_random_number((bounds.lo(1) + (bounds.hi(1)-bounds.lo(1))/2), (bounds.hi(1)-bounds.lo(1))/2);
            Real z = get_gaussian_random_number((bounds.lo(2) + (bounds.hi(2)-bounds.lo(2))/2), (bounds.hi(0)-bounds.lo(2))/2);

            if (x >= bounds.hi(0) || x < bounds.lo(0) ||
                y >= bounds.hi(1) || y < bounds.lo(1) ||
                z >= bounds.hi(2) || z < bounds.lo(2))
            {
                i--;
                continue; 
            };

            int procID = ParallelDescriptor::MyProc();
            
            ParticleType p; 

            p.id() = ParticleType::NextID(); 
            p.cpu() = procID; 
            p.pos(0) = x; 
            p.pos(1) = y; 
            p.pos(2) = z; 

            Real ux = get_gaussian_random_number(0, th_mom_std); 
            Real uy = get_gaussian_random_number(0, th_mom_std); 
            Real uz = get_gaussian_random_number(0, th_mom_std); 

            Array<Real, PIdx::nattribs> real_attribs; 

            real_attribs[PIdx::ux] = ux; 
            real_attribs[PIdx::uy] = uy; 
            real_attribs[PIdx::uz] = uz; 

            real_attribs[PIdx::Ex] = 0.0; 
            real_attribs[PIdx::Ey] = 0.0; 
            real_attribs[PIdx::Ez] = 0.0; 

            real_attribs[PIdx::Bx] = 0.0; 
            real_attribs[PIdx::By] = 0.0; 
            real_attribs[PIdx::Bz] = 0.0; 

            real_attribs[PIdx::ginv] = 0.0;
            real_attribs[PIdx::w] = 0.0;

            real_attribs[PIdx::lx] = (p.pos(0)-bounds.lo(0))/dx[0]; 
            real_attribs[PIdx::ly] = (p.pos(1)-bounds.lo(1))/dx[1]; 
            real_attribs[PIdx::lz] = (p.pos(2)-bounds.lo(2))/dx[2]; 

            particle_tile.push_back(p);
            particle_tile.push_back_real(real_attribs); 

        }        
    }
    Redistribute();
}

Geometry EMParticleContainer::get_Geom(int lev){
    return Geom(lev);
}

void EMParticleContainer::Scatter(const Vector<MultiFab*>& V, const Real value)
{   
    Print() << "scatter flag #1";
    const int num_levels = V.size();  
    for (int lev = 0; lev < num_levels; ++lev)
    {
        const Geometry& geom = Geom(lev);

        Print() << "scatter flag #2";
         
        
        for (EMParIter pti(*this, lev); pti.isValid(); ++pti)
        {   
            const Box& box = pti.validbox();

            int np = pti.numParticles();
            auto& particle_attribs = pti.GetStructOfArrays();
            
            RealVector& lx = particle_attribs.GetRealData(PIdx::lx);
            RealVector& ly = particle_attribs.GetRealData(PIdx::ly);
            RealVector& lz = particle_attribs.GetRealData(PIdx::lz);

            const auto& v = (*V[lev]).array(pti);
            Print() << "scatter flag #3";
            ParallelFor(np, [=] AMREX_GPU_DEVICE (int iP) noexcept
            { 
                Vector<Real> l = {lx[iP], lx[iP],lz[iP]};
                ParticleMeshFuncs::Scatter(v, l, value); 
            }); 
        }
    }
}   

void EMParticleContainer::initNumberDensity()
{
    //Print() << "init_ndens #1";
    N.resize(num_levels);
    int Ncomp = 1;
    for (int lev = 0; lev < num_levels; ++lev)
    {
        const DistributionMapping& dmap = em_dmap[lev];
        const BoxArray& nodeba = convert(grids[lev], IntVect::TheNodeVector());   
        N[lev].reset(new MultiFab(nodeba, dmap, Ncomp, 1));
        N[lev]->setVal(0.0); 
    }
}

void EMParticleContainer::computeNumberDensity()
{
    Scatter(GetVecOfPtrs(N), m_mpw);
    int num_levels = N.size();
    // Print() << num_levels;
    for (int lev = 0; lev < num_levels; ++lev)
    {
        //Print() << "CND #1";
        const auto& dxi = Geom(lev).InvCellSizeArray();
        const Real inv_vol = dxi[0]*dxi[1]*dxi[2];
        N[lev]->mult(inv_vol);
    }
}

void EMParticleContainer::computeRHS(const Vector<MultiFab*>& rhs)
{
    int num_levels = rhs.size();
    
    
    for (int lev = 0; lev < num_levels; ++lev)
    {

        for (MFIter mfi(*rhs[lev], TilingIfNotGPU()); mfi.isValid(); ++ mfi)
        {
            const Box& rhs_box = mfi.tilebox(IntVect::TheNodeVector());

            const auto& iRhs = (*rhs[lev]).array(mfi);
            const auto& iN = (*N[lev]).array(mfi);

            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(rhs_box, 1, i, j, k, n,
            {
                iRhs(i,j,k,n) = m_charge*iN(i,j,k,n)/ep0;
            }); 
        }

    }
}

void EMParticleContainer::computejField(const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& jField)
{

}

