#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_WriteEBSurface.H>

#include <iostream>
#include <cmath>

using namespace amrex; 


void domain_definition(); 

int main(int argc, char* argv[]){

    amrex::Initialize(argc, argv);
    
    domain_definition();
    


    amrex::Finalize();
    return 0; 

}

void domain_definition(){

    double R0 = 0.2477; 
    double a = 0.09; 
    int n_cell = 128; 
    int max_grid_size = 32;  
    double z_lim = 0.11;
    
    //Toroidal boundary SCR-1 vacuum chamber
    Array<Real, AMREX_SPACEDIM> center{0.0, 0.0, 0.0}; 
    bool inside = true; 
    EB2::TorusIF torus_bdry(Real(R0), Real(a), center, inside);
    auto shop = EB2::makeShop(torus_bdry); 

    Box domain(IntVect{AMREX_D_DECL(       0,        0,        0)},
            IntVect{AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)});

    RealBox real_box({AMREX_D_DECL(-Real(R0+a),-Real(R0+a),-z_lim)},
                 {AMREX_D_DECL(Real(R0+a), Real(R0+a), z_lim)});

    int coord = 0; 
    Array<int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(0,0,0)};


    Geometry geom(domain, real_box, coord, is_periodic); 
    
    EB2::Build(shop, geom, 0,0);
    
    Print() << "Finished building geometry."; 



}