#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "Particles.h"
#include "initBfield.h"
#include "Prob.h"
#include "utils.h"
#include "WritePlotFile.h"

using namespace amrex;

void main_loop();

int main(int argc, char* argv[])
{
    Initialize(argc, argv); 
    main_loop();
    Finalize();
    return 0;
}

void main_loop()
{
    int n_cell, max_grid_size, num_particles, max_steps, plot_int; 
    Vector<Real> init_v;
    Real q_e = -1.602e-19;
    Real m_e = 9.10938356e-31;
    Real proton_mass = 1.6726219e-27;
    Real proton_charge = -1*q_e;

    {
        ParmParse pp;
        pp.query("ncell", n_cell);
        pp.query("max_steps", max_steps);
        pp.query("max_grid_size", max_grid_size);
        pp.query("num_particles", num_particles);
        pp.query("plot_int", plot_int);
        pp.getarr("initial_velocity", init_v);
    }

    Print() << "Input File read successfully...\n";

    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};

    RealBox real_box({AMREX_D_DECL(-2,-2,-2)},
                     {AMREX_D_DECL( 2, 2, 2)});
    BoxArray ba; 
    Geometry geom; 
    
    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0)); 
    IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)); 
    Box domain(dom_lo, dom_hi); 

    ba.define(domain);
    ba.maxSize(max_grid_size);

    geom.define(domain, real_box, CoordSys::cartesian, is_periodic); 
    
    int Nghost = 1; 
    int Ncomp = 1; 

    DistributionMapping dm(ba); 

    Print() << "Geometry generation complete. \n";

    Array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> Bfield;
    Array<std::unique_ptr<MultiFab>, AMREX_SPACEDIM> dmagB;
    MultiFab Jacobian;

    Jacobian.define(ba, dm, 9, Nghost);

    for (int dim = 0; dim < AMREX_SPACEDIM; ++dim)
    {
        Bfield[dim].reset(new MultiFab(ba, dm, Ncomp, Nghost));
        dmagB[dim].reset(new MultiFab(ba, dm, Ncomp, Nghost));
    } 

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

    for (MFIter mfi(*Bfield[0]); mfi.isValid(); ++mfi)
    {
        Print() << "\tinitialing fields...\n";
        const Box& bx = mfi.validbox();

        const auto& iBx = (*Bfield[0]).array(mfi);
        const auto& iBy = (*Bfield[1]).array(mfi);
        const auto& iBz = (*Bfield[2]).array(mfi);

        const auto& igradBx = (*dmagB[0]).array(mfi);
        const auto& igradBy = (*dmagB[1]).array(mfi);
        const auto& igradBz = (*dmagB[2]).array(mfi);
        
        const auto& ijacobian = Jacobian.array(mfi);

        const auto& problo = geom.ProbLoArray();
        const auto& probhi = geom.ProbHiArray();

        const auto& dx = geom.CellSizeArray();


        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            initBfield(i, j, k, iBx, iBy, iBz, problo, probhi, dx);
        });

        Print() << "Derivatives!\n";

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            derivative_calculation(i, j, k, 
                                   iBx, iBy, iBz,
                                   igradBx, igradBy, igradBz,
                                   ijacobian, problo, dx);
        });
    }
    
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        Bfield[idim]->FillBoundary(geom.periodicity());
        dmagB[idim]->FillBoundary(geom.periodicity());
    }
    
    WritePlotFile(GetArrOfConstPtrs(Bfield), GetArrOfConstPtrs(dmagB), 
                        Jacobian, geom, 0.0, 0);
    
    std::unique_ptr<SPParticleContainer> particles;
    particles.reset(new SPParticleContainer(geom, dm, ba, proton_charge, proton_mass, 1)); 
    particles->InitParticles(num_particles, init_v);

    Real dt = 7e-10;
    Real time = 0.0;

    for (int iStep = 0; iStep < max_steps; ++iStep)
    {
        //Print() << iStep << "\n";
        particles->PushParticleTrajectory(GetArrOfPtrs(Bfield),
                                          GetArrOfPtrs(dmagB),
                                          Jacobian, dt, iStep);
        if (iStep % plot_int == 0)
        {               
            std::string plotname = Concatenate("./plotfiles/particles_", iStep, 3);
            particles->Checkpoint(plotname, "particles");
        }
        //Print() << "done"; 
        time += dt;
    }
}