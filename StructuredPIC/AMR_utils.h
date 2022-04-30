#ifndef AMRutils_H
#define AMRutils_H

// standard libraries
#include <vector>

// AMReX libraries
#include <AMReX_MultiFab.H>
#include <AMReX_IntVect.H>
#include <AMReX_REAL.H>
#include <AMReX_Geometry.H>

// local modules
#include "solctra_utils.h"

// namespaces

using namespace amrex;
using namespace std; 

// Functions

struct YGEnum
{
    enum {
        Bx = 0, By, Bz
    };
    enum {
        Ex = 0, Ey, Ez
    };
    enum {
        jx = 0, jy, jz
    };
};

class YeeGrid
{
    public:
    
        YeeGrid(){};
        ~YeeGrid(){};

        Array<IntVect, AMREX_SPACEDIM> Bfield;
        Array<IntVect, AMREX_SPACEDIM> Efield;
        Array<IntVect, AMREX_SPACEDIM> jfield;
        
        void initialize_YeeGrid();
};



void get_gaussian_random_momentum(Real* u, Real u_mean, Real u_std);

void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part);

Real get_gaussian_random_number(Real v_mean, Real v_std);

Vector<Real> compute_dt(const Vector<Geometry>& Geom);

namespace LinAlg{

    Vector<Real> dot(Vector<Real> a, Vector<Real> b);

}





#endif