#include "AMR_utils.h"
#include "constants.h"

using namespace std; 
using namespace PhysConst;

void YeeGrid::initialize_YeeGrid()
{
    Bfield[YGEnum::Bx] = IntVect(1,0,0);
    Bfield[YGEnum::By] = IntVect(0,1,0);
    Bfield[YGEnum::Bz] = IntVect(0,0,1);

    Efield[YGEnum::Ex] = IntVect(0,1,1);
    Efield[YGEnum::Ey] = IntVect(1,0,1);
    Efield[YGEnum::Ez] = IntVect(1,1,0);

    jfield[YGEnum::jx] = IntVect(0,1,1);
    jfield[YGEnum::jy] = IntVect(1,0,1);
    jfield[YGEnum::jz] = IntVect(1,1,0);
};


void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
{
    int nx = nppc[0];
    int ny = nppc[1];
    int nz = nppc[2];

    int ix_part = i_part/(ny * nz);
    int iy_part = (i_part % (ny * nz)) % ny;
    int iz_part = (i_part % (ny * nz)) / ny;

    r[0] = (0.5+ix_part)/nx;
    r[1] = (0.5+iy_part)/ny;
    r[2] = (0.5+iz_part)/nz;
}

void get_gaussian_random_momentum(Real* u, Real u_mean, Real u_std) {
    Real ux_th = amrex::RandomNormal(0.0, u_std);
    Real uy_th = amrex::RandomNormal(0.0, u_std);
    Real uz_th = amrex::RandomNormal(0.0, u_std);

    u[0] = u_mean + ux_th;
    u[1] = u_mean + uy_th;
    u[2] = u_mean + uz_th;
}

Real get_gaussian_random_number(Real v_mean, Real v_std) 
{
    Real vx = amrex::RandomNormal(v_mean, v_std);
    return vx; 
}

Vector<Real> compute_dt(const Vector<Geometry>& Geom)
{
    int num_levels = Geom.size();
    Vector<Real> dt_vec;
    for (int lev = 0; lev < num_levels; ++lev)
    {
        const Real* dx = Geom[lev].CellSize();
        const Real dt = 1./(sqrt(1./(dx[0]*dx[0]+1./(dx[1]*dx[1])+1./(dx[2]*dx[2])))*c);
        dt_vec.push_back(dt);
    }

    return dt_vec;
}