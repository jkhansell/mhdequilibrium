#ifndef PARTICLES_H
#define PARTICLES_H

#include <AMReX_Particles.H>

#include "utils.h"


using namespace amrex; 

struct PIdx
{
    enum{
        ux = 0, 
        uy, uz, w, vpar, mu,
        nattribs
    };
};

class SPPartIter
    : public amrex::ParIter<0,0, PIdx::nattribs, 0>
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

class SPParticleContainer 
    : public amrex::ParticleContainer<0,0,PIdx::nattribs,0>
{
    public:

        SPParticleContainer(const Geometry& geom, 
                            const DistributionMapping& dmap, 
                            const BoxArray& ba, 
                            const Real charge, 
                            const Real mass, 
                            const Real mpw); 

        void InitParticles(const int num_particles,
                           const Vector<Real> initv);


        // this code needs to be compiled in 3D
        void PushParticleTrajectory(const Array<MultiFab*, AMREX_SPACEDIM>& Bfield,
                                    const Array<MultiFab*, AMREX_SPACEDIM>& dmagB, 
                                    const MultiFab& Jacobian,
                                    const Real dt, const int iStep);

        // checks particle's domain in periodic conditions
        void checkparticledomain
        (
            const GpuArray<Real, AMREX_SPACEDIM> plo, 
            const GpuArray<Real, AMREX_SPACEDIM> phi,
            Real3& pos
        );


    private: 

        Real q;
        Real m;

};


#endif

