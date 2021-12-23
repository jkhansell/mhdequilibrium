

#ifndef FIELDS_H
#define FIELDS_H
// standard libraries

// AMReX libraries
#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLLinOp.H>
#include <AMReX_MLMG.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BoxArray.H>

// local modules
#include "solctra_utils.h"
#include "AMR_utils.h"

// namespaces
using namespace amrex;
using namespace std; 


//Functions

void get_initial_fields(const GlobalData& Coils, Coil* rmi, Coil* rmf, 
                        const Vector<Geometry>& geom,
                        const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& BField,
                        int n_cell, Real R0, Real dh);


class NodalFieldSolver
{
    public:
        NodalFieldSolver(const Vector<Geometry>& geom,
                         const Vector<BoxArray>& grids,
                         const Vector<DistributionMapping>& dmap): 
                            NFSgeom(geom), NFSgrids(grids), NFSdmap(dmap){}

        ~NodalFieldSolver(){Print() << "Destructed NFS.";}

        void init_sigma();
        
        void solveforPhi(const Vector<const MultiFab *>& rhs, 
                         const Vector<MultiFab *>& phi);
        
        void calculateEField(const Vector<Array<MultiFab *, AMREX_SPACEDIM>>& Efield,
                         const Vector<const MultiFab *>& phi,
                         const YeeGrid& EMYeeGrid);
        
        void evolveBfield(const Vector<Array<const MultiFab*, AMREX_SPACEDIM>>& E,
                          const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& B,
                          const YeeGrid& EMYeeGrid, const Vector<Geometry>& geom,
                          Real dt);

    private:
        
        const Vector<Geometry>& NFSgeom; 
        const Vector<BoxArray>& NFSgrids; 
        const Vector<DistributionMapping>& NFSdmap; 
        int num_levels = NFSgeom.size();
        Vector<MultiFab> sigma;
        
        int verbose = 2;
        int bottom_verbose = 1;
        int max_iter = 150;
        int max_fmg_iter = 0;
        Real reltol = 1e-5;
        Real abstol = 1e-7; 

        int gpu_regtest = 0;

        bool agglomeration = false;
        bool consolidation = false;
        bool semicoarsening = false;
        int max_coarsening_level = 30;
        int max_semicoarsening_level = 0;

};


#endif
