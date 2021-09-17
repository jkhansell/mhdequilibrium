

// standard modules 

#include <iostream>
#include <fstream>
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>

// AMReX includes

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>


// BS-SOLCTRA modules
#include "solctra_sequential.h"
#include "utils.h"

// local files
#include "initial_cond.h"

using namespace amrex; 
using namespace std; 

GlobalData initialize_coils(){
    /*
    Coil data initializer

    Obtains data from the resources folder which contains the coordinates of SCR-1's modular coils.
    */

    std::string resourcePath = "./resources";   //Coil directory

    const size_t sizeToAllocate = sizeof(double) * TOTAL_OF_GRADES_PADDED * TOTAL_OF_COILS;     //size of allocation


    GlobalData data;    //coil data structure

    // memory allocation
    data.coils.x = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.coils.y = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.coils.z = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    

    // coil data loading
    load_coil_data(data.coils.x, data.coils.y, data.coils.z, resourcePath);
    
    // e_roof memory allocation
    data.e_roof.x = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.e_roof.y = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.e_roof.z = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));
    data.leng_segment = static_cast<double*>(_mm_malloc(sizeToAllocate, ALIGNMENT_SIZE));

    // e_roof calculation
    e_roof(data);


    return data;  

}


// free coil/e_roof memory
void free_coils(GlobalData& data){
    _mm_free(data.coils.x);
    _mm_free(data.coils.y);
    _mm_free(data.coils.z);
    _mm_free(data.e_roof.x);
    _mm_free(data.e_roof.y);
    _mm_free(data.e_roof.z);
    _mm_free(data.leng_segment);
}

// Fill space with magnetic field initial condition
void fill_space(){

    // SIMULATION PARAMETERS

    int n_cell; 
    int max_grid_size; 
    int nsteps; 
    int plot_int; 

    // BS-SOLCTRA compatibility

    GlobalData Coils = initialize_coils();

    Coil rmi[TOTAL_OF_COILS];
    Coil rmf[TOTAL_OF_COILS];
    initializeGlobals(rmi, rmf);

    // Domain Definition

    BoxArray ba;
    Geometry geom; 

    IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));

    Box domain(dom_lo, dom_hi);


    ba.define(domain); 

    ba.maxSize(max_grid_size); 

    RealBox real_box({AMREX_D_DECL(0.2477-0.09, -0.09, -0.09)},
                     {AMREX_D_DECL(0.2477+0.09, 0.09, 0.09)});

    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,1,1)};

    geom.define(domain, real_box, CoordSys::cartesian, is_periodic); 

    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    // Nghost = number of ghost cells for each array
    int Nghost = 1;

    // Ncomp = number of components for each array
    int Ncomp = 3;
     // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    Array<MultiFab,AMREX_SPACEDIM> Bfield;
    MultiFab plotfile_mf;
    plotfile_mf.define(ba, dm, AMREX_SPACEDIM, 0);
    
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Bfield[idim].define(convert(ba, IntVect::TheDimensionVector(idim)), dm, 1, 1);
        for (MFIter mfi(Bfield[idim]); mfi.isValid(); ++mfi){

            const Box& bx = mfi.validbox();

            const Array4<Real>& B_i = Bfield[idim].array(mfi); 

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i ,int j, int k){
                
                Real x = (i+0.5)*dx[0];
                Real y = (j+0.5)*dx[1];
                Real z = (k+0.5)*dx[2];

                cartesian p0 = {x, y, z};
                cartesian magfield = magnetic_field(rmi, rmf, Coils, p0);   // Problem is here!!! 

                B_i(i,j,k,0) = Real(magfield.x);
                B_i(i,j,k,1) = Real(magfield.y);
                B_i(i,j,k,2) = Real(magfield.z);

                }
            );
        }
        Print() << "done";
    }
    free_coils(Coils);

}



