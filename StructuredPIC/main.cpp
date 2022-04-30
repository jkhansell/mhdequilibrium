// standard libraries


// AMReX libraries

#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Vector.H>

// local modules

#include "solctra_utils.h"
#include "AMR_utils.h"
#include "constants.h"
#include "Fields.h"
#include "Particles.h"
#include "WritePlotFile.h"


// namespaces
using namespace amrex;
using namespace PhysConst;
// Functions

void main_main(GlobalData& Coils, Coil* rmi, Coil* rmf);

int main(int argc, char* argv[])
{
    
    // BS-Solctra utils initialization 
    GlobalData Coils = initialize_coils();
    Coil rmi[TOTAL_OF_COILS];
    Coil rmf[TOTAL_OF_COILS];
    initializeGlobals(rmi, rmf);

    Initialize(argc, argv); 

    main_main(Coils, rmi, rmf);
    free_coils(Coils);
    Finalize();
    finishGlobal(rmi, rmf);
        
}

void main_main(GlobalData& Coils, Coil* rmi, Coil* rmf)
{
    
    // Parameter definition could be replaced by a Parsing function
    
    Real R0 = 0.2477;
    Real dh = 0.03;
    int n_cell = 50; 
    int max_grid_size = 60;  
    int Nghost = 1; 
    int num_levels = 1;
    int n_buffer = 2;
    int nsteps = 5;

    // logical and Real space definition

    Vector<int> rr(num_levels-1);
    for (int lev = 1; lev < num_levels; lev++)
        rr[lev-1] = 1;

    RealBox real_box({AMREX_D_DECL(R0-dh,-dh,-dh)},
                    {AMREX_D_DECL(R0+dh, dh, dh)});

    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};

    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
    
    Box domain(dom_lo, dom_hi);

    Vector<Geometry> geom(num_levels); 
    geom[0].define(domain, real_box, CoordSys::cartesian, is_periodic);

    Vector<BoxArray> grids(num_levels);
    grids[0].define(domain);

    for (int lev = 1; lev<num_levels; ++lev)
    {
        geom[lev].define(refine(geom[lev-1].Domain(), rr[lev-1]),
                                real_box, CoordSys::cartesian, is_periodic);

        int n_fine = n_cell*rr[lev-1];
        IntVect refined_lo(D_DECL(3*n_fine/8,3*n_fine/8,3*n_fine/8));
        IntVect refined_hi(D_DECL(5*n_fine/8-1,5*n_fine/8-1,5*n_fine/8-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        grids[lev].define(refined_patch);
    } 
    
    for (int lev = 0; lev < num_levels; lev++)
    {
        grids[lev].maxSize(max_grid_size);
    }


    int Ncomp = 1; 
    YeeGrid EMYeeGrid;
    EMYeeGrid.initialize_YeeGrid();

    Vector<DistributionMapping> dmap(num_levels);
    
    Vector<unique_ptr<MultiFab>> phi(num_levels);
    Vector<unique_ptr<MultiFab>> rhs(num_levels);
    Vector<unique_ptr<MultiFab>> Bmag(num_levels);
    
    Vector<Array<unique_ptr<MultiFab>,AMREX_SPACEDIM>> BField(num_levels);
    Vector<Array<unique_ptr<MultiFab>,AMREX_SPACEDIM>> EField(num_levels);
    Vector<Array<unique_ptr<MultiFab>,AMREX_SPACEDIM>> jField(num_levels);
     
    for (int lev = 0; lev < num_levels; ++lev)
    {
        const BoxArray& nodeba = convert(grids[lev], IntVect::TheNodeVector());

        dmap[lev].define(grids[lev]);

        phi[lev].reset(new MultiFab(nodeba, dmap[lev], Ncomp, 1));
        rhs[lev].reset(new MultiFab(nodeba, dmap[lev], Ncomp, 1));
        
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            EField[lev][idim].reset(new MultiFab(convert(grids[lev], EMYeeGrid.Efield[idim]),
                                                    dmap[lev], Ncomp, 1));                                   
            BField[lev][idim].reset(new MultiFab(convert(grids[lev], EMYeeGrid.Bfield[idim]),
                                                    dmap[lev], Ncomp, 1));
            jField[lev][idim].reset(new MultiFab(convert(grids[lev], EMYeeGrid.jfield[idim]),
                                                    dmap[lev], Ncomp, 1));
            EField[lev][idim]->setVal(0.0);             
            jField[lev][idim]->setVal(0.0);             

        }
        
        rhs[lev]->setVal(0.0); 
        phi[lev]->setVal(0.0);  
    }

    get_initial_fields(Coils, rmi, rmf, geom, 
                        GetVecOfArrOfPtrs(BField),
                        n_cell, R0, dh);

    int num_species = 2; 
    IntVect nppc = {30,30,30};

    Vector<std::unique_ptr<EMParticleContainer>> particles(2); 

    EMParticleContainer* electrons;
    electrons = new EMParticleContainer(geom, dmap, grids, rr, 0, -PhysConst::q_e, PhysConst::m_e, 10);
    electrons->InitParticles(nppc, electron_th_vel0, real_box);
    
    EMParticleContainer* H_ions;
    H_ions = new EMParticleContainer(geom, dmap, grids, rr, 1, PhysConst::q_e, PhysConst::m_p, 10);
    H_ions->InitParticles(nppc, proton_th_vel0, real_box);
    
    particles[0].reset(electrons); 
    particles[1].reset(H_ions);
    // Possible error with same real_box on particles!!! watch out
    
    electrons->Checkpoint("init_particles", "electrons");
    H_ions->Checkpoint("init_particles", "H_ions");
    
    Print() << "Finished initialization.\n";

    for (int sp = 0; sp < num_species; ++sp)
    {
        particles[sp]->InitPartArrVars(EMYeeGrid); 
        particles[sp]->computeNumberDensity();
        particles[sp]->computeRHS(GetVecOfPtrs(rhs));
    }

    NodalFieldSolver NFSpic(geom, grids, dmap);
    //Print() << "Flag main.cpp #1\n "; 
    NFSpic.init_sigma();
    //Print() << "Flag main.cpp #2\n "; 
    NFSpic.solveforPhi(GetVecOfConstPtrs(rhs), GetVecOfPtrs(phi));
    NFSpic.calculateEField(GetVecOfArrOfPtrs(EField), GetVecOfConstPtrs(phi), EMYeeGrid);
    
    for (int iStep = 0; iStep < nsteps; ++iStep)
    {

        Vector<Real> dt = compute_dt(geom);
        for (int lev = 0; lev < num_levels; ++lev)
        {
            dt[lev]*=0.5; // initial pushback timestep
        }

        FieldEvolver::evolve_Bfield(GetVecOfArrOfConstPtrs(EField),
                                    GetVecOfArrOfPtrs(BField), 
                                    EMYeeGrid, geom, dt); 
        
        // Magnetic Field Push back to t = -0.5*dt
                
        for(int lev = 0; lev < num_levels; ++lev)
        {
            for (int idim = 0; idim<AMREX_SPACEDIM; ++idim)
            {
                EField[lev][idim]->FillBoundary();
                BField[lev][idim]->FillBoundary();
                jField[lev][idim]->FillBoundary();
            }
        }


        WriteScalarData(GetVecOfConstPtrs(rhs), 
                    GetVecOfConstPtrs(phi), 
                    geom, 0.0, 0);
        
        WriteVectorData(GetVecOfArrOfConstPtrs(BField),
                        GetVecOfArrOfConstPtrs(EField),
                        GetVecOfArrOfConstPtrs(jField),
                        geom, 0.0, 0);
    
    }
};