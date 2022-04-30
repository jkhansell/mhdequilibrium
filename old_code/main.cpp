
#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include "AmrCoreAdv.H"
#include "solctra_utils.H"
#include "AmrParticles.H"

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {

        int num_particles;
        int max_steps;
        Real R0;
        Real dh;
        Real q_e = 1.60217662e-19;
        Real m_e = 9.10938356e-31; 

        {
            ParmParse pp("SCRparams"); 
            pp.query("R0", R0);
            pp.query("dh", dh);
        }

        {
            ParmParse pp;
            pp.query("num_particles", num_particles);
            pp.query("max_steps", max_steps);

        }


        // timer for profiling
        BL_PROFILE("main()");
        // wallclock time
        const auto strt_total = amrex::second();
        
        GlobalData Coils = initialize_coils();
        Coil rmi[TOTAL_OF_COILS]; 
        Coil rmf[TOTAL_OF_COILS]; 
        initializeGlobals(rmi, rmf);

        AmrCoreAdv amr_core_adv(Coils, rmi, rmf, R0, dh);
        // initialize AMR data
        amr_core_adv.InitData();
        amr_core_adv.WritePlotFile();
        
        free_coils(Coils);
        finishGlobal(rmi, rmf);
        
        /*
        //particle initialization
        std::unique_ptr<AmrSPPartCont> particles;
        
        particles = std::make_unique<AmrSPPartCont>(&amr_core_adv, q_e, m_e);
        particles->InitParticles(num_particles, R0, dh);
        
        Real dt = 7e-10;
        
        for (int iStep = 0; iStep < max_steps; ++iStep)
        {
            particles->PushParticleTrajectory
                (
                    amr_core_adv.Bx, 
                    amr_core_adv.By,
                    amr_core_adv.Bz,
                    amr_core_adv.gradBx,
                    amr_core_adv.gradBy,
                    amr_core_adv.gradBz,
                    amr_core_adv.Jacobian, dt, iStep
                );
        }
        */
        // wallclock time
        auto end_total = amrex::second() - strt_total;

        if (amr_core_adv.Verbose()) {
            // print wallclock time
            ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
            amrex::Print() << "\nTotal Time: " << end_total << '\n';
        }
    }
    
    amrex::Finalize();
}
