// standard c++ includes
#include <iostream>
#include <fstream>



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

void SPParticleContainer::InitParticles(const int num_particles)
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
            p.pos(0) = 0.2477 + 0.02*RandomNormal(0,1);
            p.pos(1) = 0.0;
            p.pos(2) = 0.0;
            
            //Print() << p.pos(0) << " " << p.pos(1) << " " << p.pos(2) << "\n";  

            Array<Real, PIdx::nattribs> real_attribs; 

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
    const MultiFab& Jacobian, const Vector<Real>& init_velocity,
    const Real dt, const int iStep, const Real R0, const Real a,
    const int plotint
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

        auto& vpar = partattrbs.GetRealData(PIdx::vpar);
        auto& mu = partattrbs.GetRealData(PIdx::mu);

        auto& particle_tile = pti.GetParticleTile(); 

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
            Vector<Real> jacobianp(9);

            ParticleMeshFuncs::gather_fields
            (
                pos, iBx, iBy, iBz,
                igradBx, igradBy, igradBz, ijacobian,
                Bp, gradB, jacobianp, plo, dxi
            );

            if (iStep == 0)
            {
                init_values(Bp,  init_velocity, pos, vpar[i], mu[i], q, m); 
            }
            //Print() << "flag!";
            Real3 Xdot; 
            Real mu_part;
            mu_part = mu[i];

            std::vector<bool> domain_checks(4); 
            
            Tuple<Real, Real3> k1;
            Tuple<Real, Real3> k2;
            Tuple<Real, Real3> k3;
            Tuple<Real, Real3> k4;

            domain_checks[0] = checkparticledomain(R0, a, pos);

            k1 = K_calc(Bp, gradB, jacobianp, q, m, vpar[i], mu_part);
            temppos = pos + 0.5*dt*get<1>(k1);
            
            domain_checks[1] = checkparticledomain(R0, a, temppos);
        
            tempvpar = vpar[i] + 0.5*dt*get<0>(k1);

            ParticleMeshFuncs::gather_fields
            (
                temppos, iBx, iBy, iBz,
                igradBx, igradBy, igradBz, ijacobian,
                Bp, gradB, jacobianp, plo, dxi
            );

            k2 = K_calc(Bp, gradB, jacobianp, q, m, vpar[i], mu_part);
            temppos = pos + 0.5*dt*get<1>(k2);
            domain_checks[2] = checkparticledomain(R0, a, temppos);

            tempvpar = vpar[i] + 0.5*dt*get<0>(k2);

            ParticleMeshFuncs::gather_fields
            (
                temppos, iBx, iBy, iBz,
                igradBx, igradBy, igradBz, ijacobian,
                Bp, gradB, jacobianp, plo, dxi
            );

            
            k3 = K_calc(Bp, gradB, jacobianp, q, m, vpar[i], mu_part);
            temppos = pos + dt*get<1>(k3);

            domain_checks[3] = checkparticledomain(R0, a, temppos);
            /*
            if (domain_checks[3]){particles.erase(particles.begin()+i, particles.begin()+i);
                continue;}
            */
            
            tempvpar = vpar[i] + dt*get<0>(k3);

            //Print() << temppos(0) << " " << temppos(1) << " " << temppos(2) << "\n";


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

            bool realcheck = checkparticledomain(R0, a, pos);
            /*
            if (realcheck){particles.erase(particles.begin()+i, particles.begin()+i);
                continue;}
            */

            p.pos(0) = pos(0);
            p.pos(1) = pos(1); 
            p.pos(2) = pos(2); 


            if (iStep % plotint == 0)
            {
                std::string filename = "trajectories/particle"+std::to_string(i)+".txt";
                std::ofstream myfile(filename, std::ios_base::app); 
                if (myfile.is_open())
                {
                    myfile << std::scientific; 
                    myfile <<iStep<<"\t"<<dt*iStep<<"\t"<<p.pos(0)<<"\t"<<p.pos(1)<<"\t"<<p.pos(2)<<"\t"<<vpar[i]<<"\t"<<mu_part<<"\t"<<Bp(0)<<"\t"<<Bp(1)<<"\t"<<Bp(2)<<"\n";
                }
            }
        }
    }
    Redistribute();
}

bool SPParticleContainer::checkparticledomain
(
    const Real R0,
    const Real a,
    Real3& pos
)
{
    Real R = sqrt(pow(pos(0),2)+pow(pos(1),2)); 
    Real r = sqrt(pow(R-R0,2)+pow(pos(2),2));

    if (r > a){return true;}
    else {return false;}
}


/*

Check particle periodical domain

void SPParticleContainer::checkparticledomain
(  
    const GpuArray<Real, AMREX_SPACEDIM> plo, 
    const GpuArray<Real, AMREX_SPACEDIM> phi,
    Real3& pos
)
{
    Real dx;

    if (plo[0] > pos(0)) 
    {
        //Print() << "poof1";

        dx = plo[0] - pos[0];
        pos[0] = phi[0] - dx;
    }
    else if (pos(0) > phi[0])
    {
        //Print() << "poof2";

        dx = pos[0] - phi[0];
        pos[0] = plo[0] + dx;
    }
    else if (plo[1] > pos(1))
    {
        //Print() << "poof3";

        dx = plo[1] - pos[1];
        pos[1] = phi[1] - dx;
    }
    else if (pos(1) > phi[1])
    {
        //Print() << "poof4";

        dx = pos[1] - phi[1];
        pos[1] = plo[1] + dx;
    }
    else if (plo[2] > pos(2))
    {
        //Print() << "poof5";

        dx = plo[2] - pos[2];
        pos[2] = phi[2] - dx;
    }
    else if (pos(2) > phi[2])
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


*/