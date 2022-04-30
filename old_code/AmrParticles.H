#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX.H>

#include "AmrCoreAdv.H"

using namespace amrex;

struct PIdx
{
    enum 
    {
        ux = 0, uy, uz, vpar, mu, nattribs
    };
};

inline Real get_gaussian_random_number(Real v_mean, Real v_std) 
{
    Real vx = amrex::RandomNormal(v_mean, v_std);
    return vx; 
}

class AmrSPPartCont
    : AmrParticleContainer<0,0,PIdx::nattribs,0>
{
    public:
        AmrSPPartCont(AmrCoreAdv* amrcore, Real q, Real m);
        ~AmrSPPartCont(){};

        void InitParticles(const int num_particles,
                                  const Real R0, 
                                  const Real dh);

        void CheckPoint(const std::string& name){
            Checkpoint(name, "particles");
        }
        
        void PushParticleTrajectory(const Vector<MultiFab>& Bx, 
                                    const Vector<MultiFab>& By,
                                    const Vector<MultiFab>& Bz,
                                    const Vector<MultiFab>& gradBx,
                                    const Vector<MultiFab>& gradBy,
                                    const Vector<MultiFab>& gradBz, 
                                    const Vector<MultiFab>& Jacobian,    
                                    Real dt, int iStep);
        
        void CheckParticles(); 
        
        void RedistributeParticles()
        {
            Redistribute();
        }

        private:
            Real charge; 
            Real mass;
};

class AmrSPPartIter
    : public ParIter<0,0,PIdx::nattribs,0>
{
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
