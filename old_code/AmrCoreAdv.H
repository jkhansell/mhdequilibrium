#ifndef AmrCoreAdv_H_
#define AmrCoreAdv_H_

#include <string>
#include <limits>
#include <memory>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_BCRec.H>
#include <AMReX_Particles.H>

#include "solctra_utils.H" 

using namespace amrex;

class AmrCoreAdv
    : public AmrCore
{
public:

    ////////////////
    // public member functions

    // constructor - reads in parameters from inputs file
    //             - sizes multilevel arrays and data structures
    AmrCoreAdv(GlobalData& coils, Coil* Rmi, Coil* Rmf,  Real r, Real deltah);
    virtual ~AmrCoreAdv();

    // advance solution to final time
    //void Evolve ();

    // initializes multilevel data
    void InitData ();

    // Make a new level using provided BoxArray and DistributionMapping and
    // fill with interpolated coarse level data.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
                                         const DistributionMapping& dm) override;

    // Remake an existing level using provided BoxArray and DistributionMapping and
    // fill with existing fine and coarse data.
    // overrides the pure virtual function in AmrCore
    virtual void RemakeLevel (int lev, Real time, const BoxArray& ba,
                              const DistributionMapping& dm) override;

    // Delete level data
    // overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

    // Make a new level from scratch using provided BoxArray and DistributionMapping.
    // Only used during initialization.
    // overrides the pure virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
                                          const DistributionMapping& dm) override;

    // tag all cells for refinement
    // overrides the pure virtual function in AmrCore
    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) override;

    // write plotfile to disk
    void WritePlotFile () const;

    Vector<MultiFab> Bx;
    Vector<MultiFab> By;
    Vector<MultiFab> Bz;

    Vector<MultiFab> gradBx;
    Vector<MultiFab> gradBy;
    Vector<MultiFab> gradBz;

    Vector<MultiFab> Jacobian;

private:

    ////////////////
    // private member functions

    // read in some parameters from inputs file
    void ReadParameters();

    // set covered coarse cells to be the average of overlying fine cells
    void AverageDown ();

    // more flexible version of AverageDown() that lets you average down across multiple levels
    void AverageDownTo (int crse_lev);

    // compute a new multifab by coping in phi from valid region and filling ghost cells
    // works for single level and 2-level cases (fill fine grid ghost by interpolating from coarse)
    void FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp);

    // fill an entire multifab by interpolating from the coarser level
    // this comes into play when a new level of refinement appears
    void FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp);

    // utility to copy in data from phi_old and/or phi_new into another multifab
    void GetData (int lev, Real time, Vector<MultiFab*>& data,
                  Vector<Real>& datatime);

    // Advance a level by dt - includes a recursive call for finer levels
    void timeStepWithSubcycling (int lev, Real time, int iteration);

    // Advance all levels by the same dt
    void timeStepNoSubcycling (Real time, int iteration);

    // a wrapper for EstTimeStep
    void ComputeDt ();

    // get plotfile name
    std::string PlotFileName (int lev) const;

    // put together an array of multifabs for writing
    Vector<const MultiFab*> PlotFileMF () const;

    // set plotfile variables names
    Vector<std::string> PlotFileVarNames () const;

    // write checkpoint file to disk
    void WriteCheckpointFile () const;

    // read checkpoint file from disk
    void ReadCheckpointFile ();

    ////////////////
    // private data members

    Vector<int> istep;      // which step?
    Vector<int> nsubsteps;  // how many substeps on each level?

    // keep track of old time, new time, and time step at each level
    Vector<Real> t;
    //Vector<Real> t_old;
    Vector<Real> dt;

    // array of multifabs to store the solution at each level of refinement
    // after advancing a level we use "swap".



    // this is essentially a 2*DIM integer array storing the physical boundary
    // condition types at the lo/hi walls in each direction
    Vector<BCRec> bcs;  // 1-component

    // stores fluxes at coarse-fine interface for synchronization
    // this will be sized "nlevs_max+1"
    // NOTE: the flux register associated with flux_reg[lev] is associated
    // with the lev/lev-1 interface (and has grid spacing associated with lev-1)
    // therefore flux_reg[0] and flux_reg[nlevs_max] are never actually
    // used in the reflux operation
    //Vector<std::unique_ptr<FluxRegister> > flux_reg;

    // Velocity on all faces at all levels
    //Vector< Array<MultiFab, AMREX_SPACEDIM> > facevel;
    ////////////////
    // runtime parameters

    // maximum number of steps and stop time
    int max_step = std::numeric_limits<int>::max();
    Real stop_time = std::numeric_limits<Real>::max();

    // if >= 0 we restart from a checkpoint
    //std::string restart_chkfile = "";

    // advective cfl number - dt = cfl*dx/umax
    //Real cfl = 0.7;

    // how often each level regrids the higher levels of refinement
    // (after a level advances that many time steps)
    int regrid_int = 2;

    // hyperbolic refluxing as part of multilevel synchronization
    //int do_reflux = 1;

    // do we subcycle in time?
    //int do_subcycle = 1;

    // plotfile prefix and frequency
    std::string plot_file {"plt"};
    int plot_int = -1;

    //Solctra resources
    GlobalData Coils; 
    Coil* rmi; 
    Coil* rmf;

    Real R0;
    Real dh;
    
};

#endif
