
#include <AMReX.H>
#include <AMReX_Print.H>
#include <omp.h>

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        amrex::Print() << "Hello world from AMReX version " << amrex::Version() << "\n";
    }
    amrex::Finalize();
}

