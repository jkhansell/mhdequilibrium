#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Particles.H>

#include "AmrParticles.H"
#include "ParticleMesh.H"
#include "AmrCoreAdv.H"
#include "Prob.H"

using namespace amrex;


AmrSPPartCont::AmrSPPartCont(AmrCoreAdv* amrcore, Real q, Real m)
    : AmrParticleContainer<0,0,PIdx::nattribs,0>(amrcore), charge(q), mass(m)
{
}

void 
AmrSPPartCont::InitParticles
(
    const int num_particles,
    const Real R0, 
    const Real dh
)
{
    for (int lev = 0; lev < numLevels(); ++lev)
    {
        const auto dx = Geom(lev).CellSizeArray();
        const auto plo = Geom(lev).ProbLoArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto& particles = GetParticles(lev);
            auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];

            for(int np = 0; np < num_particles; ++np)
            {
                Real R = R0 + get_gaussian_random_number(0, 3*dh/5);
                Real Z = get_gaussian_random_number(0, 3*dh/5);

                int procID = ParallelDescriptor::MyProc(); 
                ParticleType p;
                p.id() = ParticleType::NextID();
                p.cpu() = procID;
                p.pos(0) = R;
                p.pos(1) = 0.0; 
                p.pos(2) = Z;

                Real ux = 1e3;
                Real uy = 0.0; 
                Real uz = 0.0;

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
    }
    Redistribute();
}

void 
AmrSPPartCont::PushParticleTrajectory
(
    const Vector<MultiFab>& Bx, 
    const Vector<MultiFab>& By,
    const Vector<MultiFab>& Bz,
    const Vector<MultiFab>& gradBx,
    const Vector<MultiFab>& gradBy,
    const Vector<MultiFab>& gradBz, 
    const Vector<MultiFab>& Jacobian,  
    Real dt, int iStep
)
{
    for (int lev = 1; lev < numLevels(); ++lev)
    {
        const auto dxi = Geom(lev).InvCellSizeArray(); 
        const auto plo = Geom(lev).ProbLoArray(); 

        for (AmrSPPartIter pti(*this, lev); pti.isValid(); ++pti)
        {
            const int np = pti.numParticles();
            
            auto& particles = pti.GetArrayOfStructs();
            auto& partattrbs = pti.GetStructOfArrays();

            const auto& iBx = Bx[lev].array(pti);
            const auto& iBy = By[lev].array(pti);
            const auto& iBz = Bz[lev].array(pti);

            const auto& igradBx = gradBx[lev].array(pti);
            const auto& igradBy = gradBy[lev].array(pti);
            const auto& igradBz = gradBz[lev].array(pti);

            const auto& ijacobian = Jacobian[lev].array(pti);

            auto& ux = partattrbs.GetRealData(PIdx::ux);
            auto& uy = partattrbs.GetRealData(PIdx::uy);
            auto& uz = partattrbs.GetRealData(PIdx::uz);
            auto& vpar = partattrbs.GetRealData(PIdx::vpar);
            auto& mu = partattrbs.GetRealData(PIdx::mu);

            //AMREX_PARALLEL_FOR_1D(np, i,
            for (int i = 0; i < np; ++i)
            {
                Vector<Real> Xdot = {ux[i], uy[i], uz[i]};
                Vector<Real> position(3);
                Vector<Real> Bp(3);
                Vector<Real> gradB(3);
                Vector<Real> jacobianp(9);

                Real vpardot_1; 
                Real vpardot_2;
                Real vpardot_3;
                Real vpardot_4;

                Vector<Real> position1;

                Real Vpar;

                Vector<Real> Xdot_1;
                Vector<Real> Xdot_2;
                Vector<Real> Xdot_3;
                Vector<Real> Xdot_4;

                Tuple<Real, Vector<Real>> sol;

                position[0] = particles[i].pos(0);
                position[1] = particles[i].pos(1);
                position[2] = particles[i].pos(2);

                ParticleMeshFuncs::gather_fields(position,iBx, iBy, iBz,
                                                 igradBx, igradBy, igradBz,
                                                 ijacobian,
                                                 Bp, gradB, jacobianp, plo, dxi);

                if(iStep==0){init_values(Bp, Xdot, vpar[i], mu[i], charge, mass);}

                sol = advance_calc(Bp, gradB, jacobianp, 
                                charge, mass, vpar[i], mu[i]);

                vpardot_1 = std::get<0>(sol);
                Xdot_1 = std::get<1>(sol);

                Vpar = vpar[i] + 0.5*dt*vpardot_1;
                position1 = position + 0.5*dt*Xdot_1;

                ParticleMeshFuncs::gather_fields(position1,iBx, iBy, iBz,
                                                 igradBx, igradBy, igradBz,
                                                 ijacobian,
                                                 Bp, gradB, jacobianp, plo, dxi);

                sol = advance_calc(Bp, gradB, jacobianp, 
                                charge, mass, Vpar, mu[i]);

                vpardot_2 = std::get<0>(sol);
                Xdot_2 = std::get<1>(sol);

                position1 = position + 0.5*dt*Xdot_2;
                Vpar = vpar[i] + 0.5*dt*vpardot_2;

                ParticleMeshFuncs::gather_fields(position1,iBx, iBy, iBz,
                                                 igradBx, igradBy, igradBz,
                                                 ijacobian,
                                                 Bp, gradB, jacobianp, plo, dxi);

                sol = advance_calc(Bp, gradB, jacobianp, 
                                charge, mass, Vpar, mu[i]);

                vpardot_3 = std::get<0>(sol);
                Xdot_3 = std::get<1>(sol);

                position1 = position + 0.5*dt*Xdot_3;
                Vpar = vpar[i] + dt*vpardot_3;
                
                ParticleMeshFuncs::gather_fields(position1,iBx, iBy, iBz,
                                                 igradBx, igradBy, igradBz,
                                                 ijacobian,
                                                 Bp, gradB, jacobianp, plo, dxi);
                
                sol = advance_calc(Bp, gradB, jacobianp, 
                                charge, mass, Vpar, mu[i]);

                vpardot_4 = std::get<0>(sol);
                Xdot_4 = std::get<1>(sol);

                Real vpardotRK = (vpardot_1+2*vpardot_2+2*vpardot_3+vpardot_4)/6;
                Vector<Real> Xdot_RK = (Xdot_1+2*Xdot_2+2*Xdot_3+Xdot_4)/6;
                
                position1 = position + dt*Xdot_RK;
                Vpar += dt*vpardotRK;

                particles[i].pos(0) = position1[0]; 
                particles[i].pos(1) = position1[1];                
                particles[i].pos(2) = position1[2];

                vpar[i] += Vpar; 
            }
        } 
    }
}