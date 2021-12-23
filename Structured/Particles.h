
#ifndef PARTICLES_H
#define PARTICLES_H
// standard libraries

// AMReX libraries

#include <AMReX_Particles.H>
// local modules

#include "constants.h"

// namespaces
using namespace amrex;
using namespace std; 

// classes & structs

struct PIdx
{
    enum {
        ux = 0, 
        uy, uz, w, Ex, Ey, Ez, Bx, By, Bz, ginv, lx, ly, lz,
        nattribs
    };
};


class EMParIter
    : public amrex::ParIter<0,0,PIdx::nattribs,0>
{
public:
    using amrex::ParIter<0,0,PIdx::nattribs,0>::ParIter;

    //    EMParIter (ContainerType& pc, int level);

    const std::array<RealVector, PIdx::nattribs>& GetAttribs () const {
        return GetStructOfArrays().GetRealData();
    }

    std::array<RealVector, PIdx::nattribs>& GetAttribs () {
        return GetStructOfArrays().GetRealData();
    }

    const RealVector& GetAttribs (int comp) const {
        return GetStructOfArrays().GetRealData(comp);
    }

    RealVector& GetAttribs (int comp) {
        return GetStructOfArrays().GetRealData(comp);
    }
};

class EMParticleContainer
    : public amrex::ParticleContainer<0,0,PIdx::nattribs,0>
{
    public:

        EMParticleContainer (const Vector<Geometry>            & geom,
                             const Vector<DistributionMapping> & dmap,
                             const Vector<BoxArray>            & ba,
                             const amrex::Vector<int>          & rr,
                             const int                          a_species_id,
                             const amrex::Real                  a_charge,
                             const amrex::Real                  a_mass, 
                             const amrex::Real                  mpw);

        void InitParticles(const IntVect& nppc,
                            const Real th_mom_std,
                            const RealBox& bounds);

        void initNumberDensity();

        void Scatter(const Vector<MultiFab*>& V, const Real value);

        void computeNumberDensity();

        void computeRHS(const Vector<MultiFab*>& rhs);

        void computejField(const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& jField);


        
        void RedistributeLocal()
        {
            const int lev_min = 0;
            const int lev_max = 0;
            const int nGrow = 0;
            const int local = 1;
            Redistribute(lev_min, lev_max, nGrow, local);
        }

        Geometry get_Geom(int lev); 
        Vector<unique_ptr<MultiFab>> N;
        const Vector<BoxArray>& grids;
        const Vector<DistributionMapping>& em_dmap;
        Real m_charge;

    protected:
        int m_species_id;
        int num_levels;
        amrex::Real m_mass;
        amrex::Real m_mpw;
};

// Functions



#endif