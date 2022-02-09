// standard c++ includes
#include <iostream>
#include <fstream> 

//AMReX include files for runtime

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>


//Local include files

#include "Particles.h"
#include "Fields.h"
#include "Prob.h"
#include "utils.h"
#include "WritePlotFile.h"
#include "solctra_utils.h"

using namespace amrex;

void main_loop(GlobalData& Coils, Coil* rmi, Coil* rmf)
{
    int n_cell, max_grid_size, num_particles, max_steps, plot_int; 
    Vector<Real> init_v;
    Real q_e = -1.602e-19;
    Real m_e = 9.10938356e-31;
    Real proton_mass = 1.6726219e-27;
    Real proton_charge = -1*q_e;
    Real R0; Real a;

    {
        ParmParse pp;
        pp.query("ncell", n_cell);
        pp.query("max_steps", max_steps);
        pp.query("max_grid_size", max_grid_size);
        pp.query("num_particles", num_particles);
        pp.query("R0", R0);
        pp.query("a", a);
        pp.query("plot_int", plot_int);
        pp.getarr("initial_velocity", init_v);
    }

    Print() << "Input File read successfully...\n";

    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};

    RealBox real_box({AMREX_D_DECL(-0.3,-0.3,-0.05)},
                     {AMREX_D_DECL( 0.3, 0.3, 0.05)});
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

    std::unique_ptr<iMultiFab> masks = get_masks(ba, dm, geom, R0, a);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

    for (MFIter mfi(*Bfield[0]); mfi.isValid(); ++mfi)
    {
        Print() << "\tinitializing fields...\n";
        const Box& bx = mfi.validbox();

        const auto& imask = (*masks).array(mfi);

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

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            get_initial_fields(i, j, k, Coils, rmi, rmf, problo,
                     dx, iBx, iBy, iBz, imask);
        });

        Print() << "Derivatives!\n";

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {

            derivative_calculation(i, j, k, 
                                   iBx, iBy, iBz,
                                   igradBx, igradBy, igradBz,
                                   ijacobian, problo, dx);
            //Print() << "Derivatives!\n";

        });
        //Print() << "Derivatives!\n";
    }
    
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        Bfield[idim]->FillBoundary(geom.periodicity());
        dmagB[idim]->FillBoundary(geom.periodicity());
    }

    /*
    WritePlotFile(GetArrOfConstPtrs(Bfield), GetArrOfConstPtrs(dmagB), 
                        Jacobian, geom, 0.0, 0);
    */
    std::unique_ptr<SPParticleContainer> particles;
    particles.reset(new SPParticleContainer(geom, dm, ba, proton_charge, proton_mass, 1)); 
    particles->InitParticles(num_particles);

    Real dt = 5e-9;
    Real time = 0.0;

    for (int i = 0; i < num_particles;++i)
    {
        std::string filename = "trajectories/particle"+ std::to_string(i)+".txt";
        std::ofstream outfile(filename);
        outfile << "step\truntime\tx\ty\tz\tvpar\tmu\tvz\tBx\tBy\tBz\n";
        outfile.close();
    }

    for (int iStep = 0; iStep < max_steps; ++iStep)
    {
        

        //Print() << iStep << "\n";
        particles->PushParticleTrajectory(GetArrOfPtrs(Bfield),
                                          GetArrOfPtrs(dmagB),
                                          Jacobian, init_v, dt, iStep, R0, a, plot_int);
        
        //Print() << "Flag!\n";
        /*
        if (iStep % plot_int == 0)
        {               
            std::string plotname = Concatenate("./plotfiles/particles_", iStep, 3);
            particles->Checkpoint(plotname, "particles");
        }
        */
        //Print() << "done"; 

        time += dt;
    }
}

int main(int argc, char* argv[])
{
    Initialize(argc, argv); 
    // BS-Solctra utils initialization 
    GlobalData Coils = initialize_coils();
    Coil rmi[TOTAL_OF_COILS];
    Coil rmf[TOTAL_OF_COILS];
    initializeGlobals(rmi, rmf);
    main_loop(Coils, rmi, rmf);
    free_coils(Coils);
    finishGlobal(rmi, rmf);
    Finalize();
    return 0;
}