#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>

#include "solctra_sequential.h"
#include "initial_cond.h"

using namespace amrex; 


int main(int argc, char* argv[]){

    Initialize(argc, argv);

    fill_space();
    Finalize(); 
    return 0; 
}
