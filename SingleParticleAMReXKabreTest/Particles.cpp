
#include "Particles.h"
#include "utils.h"
#include "ParticleMesh.h"
#include "Prob.h"

using namespace amrex;

SPParticleContainer::SPParticleContainer
(
    const Geometry& geom, 
    const DistributionMapping& dmap, 
    const BoxArray& ba, 
    const Real charge, 
    const Real mass, 
    const Real mpw
)
    : ParticleContainer<0,0,PIdx::nattribs,0>(geom, dmap, ba), q(charge), m(mass)
{}

void SPParticleContainer::InitParticles(const int num_particles, const Vector<Real> initv)
{
    const int lev = 0; 
    const auto dx = Geom(lev).CellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();


#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        auto& particles = GetParticles(lev);
        auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];

        const Box& bx = mfi.tilebox();

        //AMREX_PARALLEL_FOR_1D(num_particles, i, 
        for (int i = 0; i < num_particles; ++i)
        {
            int procID = ParallelDescriptor::MyProc();
            ParticleType p;
            p.id() = ParticleType::NextID(); 
            p.cpu() = procID; 
            p.pos(0) = 0.0;
            p.pos(1) = (phi[1]-plo[1])*Random();
            p.pos(2) = (phi[2]-plo[2])*Random();
            
            //Print() << p.pos(0) << " " << p.pos(1) << " " << p.pos(2) << "\n";  
            Real ux = initv[0]*RandomNormal(0,1);
            Real uy = initv[1]*RandomNormal(0,1);
            Real uz = initv[2]*RandomNormal(0,1);

            Array<Real, PIdx::nattribs> real_attribs; 

            real_attribs[PIdx::ux] = ux; 
            real_attribs[PIdx::uy] = uy; 
            real_attribs[PIdx::uz] = uz; 

            real_attribs[PIdx::vpar] = 0.0; 
            real_attribs[PIdx::mu] = 0.0;
            
            particle_tile.push_back(p); 
            particle_tile.push_back_real(real_attribs);
        }
    }
    Redistribute();
}

void SPParticleContainer::PushParticleTrajectory
(
    const Array<MultiFab*, AMREX_SPACEDIM>& Bfield,
    const Array<MultiFab*, AMREX_SPACEDIM>& dmagB, 
    const MultiFab& Jacobian,
    const Real dt, const int iStep 
)
{
    const int lev = 0; 
    const auto dxi = Geom(lev).InvCellSizeArray(); 
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();

    for (SPPartIter pti (*this, lev); pti.isValid(); ++pti)
    {
        const int np = pti.numParticles(); 
        auto& particles = pti.GetArrayOfStructs();
        auto& partattrbs = pti.GetStructOfArrays(); 

        const auto& iBx = (*Bfield[0]).array(pti); 
        const auto& iBy = (*Bfield[1]).array(pti); 
        const auto& iBz = (*Bfield[2]).array(pti);

        const auto& igradBx = (*dmagB[0]).array(pti); 
        const auto& igradBy = (*dmagB[1]).array(pti); 
        const auto& igradBz = (*dmagB[2]).array(pti);

        const auto& ijacobian = Jacobian.array(pti);

        auto& ux = partattrbs.GetRealData(PIdx::ux);
        auto& uy = partattrbs.GetRealData(PIdx::uy);
        auto& uz = partattrbs.GetRealData(PIdx::uz);

        auto& vpar = partattrbs.GetRealData(PIdx::vpar);
        auto& mu = partattrbs.GetRealData(PIdx::mu);

        for (int i = 0; i < np; ++i)
        {

            ParticleType& p = particles[i];
            Real3 pos = {p.pos(0),
                         p.pos(1),
                         p.pos(2)};
            Real3 temppos; 
            Real tempvpar;
            Real3 gradB;
            Real3 Bp;
            Real3 vel = {ux[i], uy[i], uz[i]};
            Vector<Real> jacobianp(9);

            ParticleMeshFuncs::gather_fields
            (
                pos, iBx, iBy, iBz,
                igradBx, igradBy, igradBz, ijacobian,
                Bp, gradB, jacobianp, plo, dxi
            );

            if (iStep == 0)
            {
                init_values(Bp, vel, pos, vpar[i], mu[i], q, m); 
            }

            Real3 Xdot; 
            Real mu_part;
            mu_part = mu[i]; 
            
            Tuple<Real, Real3> k1;
            Tuple<Real, Real3> k2;
            Tuple<Real, Real3> k3;
            Tuple<Real, Real3> k4;

            checkparticledomain(plo, phi, pos);
            //Print() << pos(0) << " " << pos(1) << " " << pos(2) << "\n";

            k1 = K_calc(Bp, gradB, jacobianp, q, m, vpar[i], mu_part);
            temppos = pos + 0.5*dt*get<1>(k1);
            
            //Print() << temppos(0) << " " << temppos(1) << " " << temppos(2) << "\n";
            checkparticledomain(plo, phi, temppos);

            tempvpar = vpar[i] + 0.5*dt*get<0>(k1);

            ParticleMeshFuncs::gather_fields
            (
                temppos, iBx, iBy, iBz,
                igradBx, igradBy, igradBz, ijacobian,
                Bp, gradB, jacobianp, plo, dxi
            );

            k2 = K_calc(Bp, gradB, jacobianp, q, m, vpar[i], mu_part);
            temppos = pos + 0.5*dt*get<1>(k2);
            
            checkparticledomain(plo, phi, temppos);
            tempvpar = vpar[i] + 0.5*dt*get<0>(k2);
            
            //Print() << temppos(0) << " " << temppos(1) << " " << temppos(2) << "\n";

            ParticleMeshFuncs::gather_fields
            (
                temppos, iBx, iBy, iBz,
                igradBx, igradBy, igradBz, ijacobian,
                Bp, gradB, jacobianp, plo, dxi
            );

            
            k3 = K_calc(Bp, gradB, jacobianp, q, m, vpar[i], mu_part);
            temppos = pos + dt*get<1>(k3);

            checkparticledomain(plo, phi, temppos);
            
            tempvpar = vpar[i] + dt*get<0>(k3);

            Print() << temppos(0) << " " << temppos(1) << " " << temppos(2) << "\n";


            ParticleMeshFuncs::gather_fields
            (
                temppos, iBx, iBy, iBz,
                igradBx, igradBy, igradBz, ijacobian,
                Bp, gradB, jacobianp, plo, dxi
            );

            k4 = K_calc(Bp, gradB, jacobianp, q, m, vpar[i], mu_part);

            Real vparRK = (get<0>(k1)+get<0>(k2)+get<0>(k3)+get<0>(k4))/6.;
            Real3 XdotRK = (get<1>(k1)+get<1>(k2)+get<1>(k3)+get<1>(k4))/6.;
            
            vpar[i] += vparRK*dt;
            
            pos += XdotRK*dt;

            checkparticledomain(plo, phi, pos);

            p.pos(0) = pos(0);
            p.pos(1) = pos(1); 
            p.pos(2) = pos(2); 

            //Print() << pos(0) << " " << pos(1) << " " << pos(2) << "\n";

            

        }
    }
    Redistribute();
}

void SPParticleContainer::checkparticledomain
(  
    const GpuArray<Real, AMREX_SPACEDIM> plo, 
    const GpuArray<Real, AMREX_SPACEDIM> phi,
    Real3& pos
)
{
    Real dx;

    if (plo[0] >= pos(0))
    {
        //Print() << "poof1";
        dx = plo[0] - pos[0];
        pos[0] = phi[0] - dx;
    }
    else if (pos(0) >= phi[0])
    {
        //Print() << "poof2";

        dx = pos[0] - phi[0];
        pos[0] = plo[0] + dx;
    }
    else if (plo[1] >= pos(1))
    {
        //Print() << "poof3";

        dx = plo[1] - pos[1];
        pos[1] = phi[1] - dx;
    }
    else if (pos(1) >= phi[1])
    {
        //Print() << "poof4";

        dx = pos[1] - phi[1];
        pos[1] = plo[1] + dx;
    }
    else if (plo[2] >= pos(2))
    {
        //Print() << "poof5";

        dx = plo[2] - pos[2];
        pos[2] = phi[2] - dx;
    }
    else if (pos(2) >= phi[2])
    {
        //Print() << "poof6";

        dx = pos[2] - phi[2];
        pos[2] = plo[2] + dx;
    }
    else
    {
        //Print() << "oops";
    }
}


/*

            K[0] = K_calc(Bp, gradB, jacobianp, q, m, tempvpar, mu_part); 
            
            tempvpar = vpar[i] + 0.5*dt*get<0>(K[0]);
            temppos = pos + 0.5*dt*get<1>(K[0]);

            ParticleMeshFuncs::gather_fields
            (
                temppos, iBx, iBy, iBz,
                igradBx, igradBy, igradBz, ijacobian,
                Bp, gradB, jacobianp, plo, dxi
            );

            K[1] = K_calc(Bp, gradB, jacobianp, q, m, tempvpar, mu_part); 
            
            tempvpar = vpar[i] + 0.5*dt*get<0>(K[1]);
            temppos = pos + 0.5*dt*get<1>(K[1]);

            ParticleMeshFuncs::gather_fields
            (
                temppos, iBx, iBy, iBz,
                igradBx, igradBy, igradBz, ijacobian,
                Bp, gradB, jacobianp, plo, dxi
            );

            K[2] = K_calc(Bp, gradB, jacobianp, q, m, tempvpar, mu_part); 
            
            tempvpar = vpar[i] + dt*get<0>(K[2]);
            temppos = pos + dt*get<1>(K[2]);

            ParticleMeshFuncs::gather_fields
            (
                temppos, iBx, iBy, iBz,
                igradBx, igradBy, igradBz, ijacobian,
                Bp, gradB, jacobianp, plo, dxi
            );

            K[3] = K_calc(Bp, gradB, jacobianp, q, m, tempvpar, mu_part); 
            
            tempvpar = vpar[i] + dt*get<0>(K[2]);
            temppos = pos + dt*get<1>(K[2]);








*/