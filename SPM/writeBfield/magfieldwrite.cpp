#include <netcdf.h>
#include <vector>
#include <iostream>
#include <omp.h>

#include "solctra_utils.h"

using namespace std;

#define FILE_NAME1 "SCR1FIELD.nc"
#define FILE_NAME2 "Jacobian.nc"
#define NDIMS 3
#define NCOMP 4
#define NR 64
#define NPHI 360
#define NZ 64
#define R0 0.2477
#define a 0.06
//#define PI 3.14159265358979323846

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

void mag_field
(
    int i, int j, int k, 
    GlobalData& coils, Coil* rmi, Coil* rmf,
    double (&B)[NR][NPHI][NZ][NDIMS], 
    std::vector<double> dh,
    std::vector<double> plo
);

void calculate_jacobian
(
    int i, int j, int k,
    double (&J)[NR][NPHI][NZ][NDIMS][NDIMS],
    const double (&B)[NR][NPHI][NZ][NDIMS],
    std::vector<double> plo,
    std::vector<double> dh
);

int main()
{
    GlobalData Coils = initialize_coils();
    Coil rmi[TOTAL_OF_COILS]; 
    Coil rmf[TOTAL_OF_COILS]; 
    initializeGlobals(rmi, rmf);


    int ncid, r_dimid, phi_dimid, z_dimid, comp_dimid, varid, row_dimid, col_dimid;
    int dimids[NCOMP], Jdimids[NCOMP+1];
    static double magnetic_field[NR][NPHI][NZ][NDIMS];
    static double Jacobian[NR][NPHI][NZ][NDIMS][NDIMS];
    int retval; 
    
    std::vector<double> dh = {2*a/NR, 2*PI/NPHI, 2*a/NZ};
    std::vector<double> plo = {R0-a, 0, -a};
    std::vector<double> phi = {R0+a, 2*PI, a};

#pragma omp parallel for collapse(3)
    for (int i = 0; i < NR; ++i){
        for (int j = 0; j < NPHI; ++j){
            for (int k = 0; k < NZ; ++k){
                //std::cout << k;
                mag_field(i, j, k, Coils, rmi, rmf, 
                        magnetic_field, dh, plo);
            }
        }
    }
    
#pragma omp parallel for collapse(3) 
    for (int i = 0; i < NR; ++i){
        for (int j = 0; j < NPHI; ++j){
            for (int k = 0; k < NZ; ++k){
                //std::cout << k;
                calculate_jacobian(i, j, k, Jacobian, magnetic_field, plo, dh);
            }
        }
    }

    free_coils(Coils);
    finishGlobal(rmi, rmf);

    printf("Coils freed!");

    if ((retval = nc_create(FILE_NAME1, NC_CLOBBER, &ncid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "r", NR, &r_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "phi", NPHI, &phi_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "z", NZ, &z_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "comp", NDIMS, &comp_dimid)))
        ERR(retval);

    dimids[0] = r_dimid;
    dimids[1] = phi_dimid;
    dimids[2] = z_dimid;
    dimids[3] = comp_dimid;

    if ((retval = nc_def_var(ncid, "data", NC_DOUBLE, NCOMP, dimids, &varid)))
        ERR(retval);
    
    if ((retval = nc_enddef(ncid)))
        ERR(retval)

    if ((retval = nc_put_var_double(ncid, varid, &magnetic_field[0][0][0][0])))
        ERR(retval);
    if ((retval = nc_close(ncid)))
        ERR(retval);

    printf("*** SUCCESS!\n");

    if ((retval = nc_create(FILE_NAME2, NC_CLOBBER, &ncid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "r", NR, &r_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "phi", NPHI, &phi_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "z", NZ, &z_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "row", NDIMS, &row_dimid)))
        ERR(retval);
    if ((retval = nc_def_dim(ncid, "col", NDIMS, &col_dimid)))
        ERR(retval);

    Jdimids[0] = r_dimid;
    Jdimids[1] = phi_dimid;
    Jdimids[2] = z_dimid;
    Jdimids[3] = row_dimid;
    Jdimids[4] = col_dimid;

    if ((retval = nc_def_var(ncid, "data", NC_DOUBLE, NCOMP+1, Jdimids, &varid)))
        ERR(retval);
    
    if ((retval = nc_enddef(ncid)))
        ERR(retval)

    if ((retval = nc_put_var_double(ncid, varid, &Jacobian[0][0][0][0][0])))
        ERR(retval);
    if ((retval = nc_close(ncid)))
        ERR(retval);

   printf("*** SUCCESS!\n");
   return 0;
}

void mag_field
(
    int i, int j, int k, 
    GlobalData& coils, Coil* rmi, Coil* rmf,
    double (&B)[NR][NPHI][NZ][NDIMS], 
    std::vector<double> dh,
    std::vector<double> plo
)
{
    
    double R = plo[0] + i*dh[0];
    double phi = plo[1] + j*dh[1];
    double Z = plo[2] + k*dh[2];

    double X = R*cos(phi); 
    double Y = R*sin(phi);

    cartesian p0 = {X,Y,Z}; 
    cartesian B_vect = magnetic_field(rmi, rmf, coils, p0);

    // Stores magnetic field in cylindrical coordinates
    B[i][j][k][0] = B_vect.x*cos(phi)+B_vect.y*sin(phi);
    B[i][j][k][1] =-B_vect.x*sin(phi)+B_vect.y*cos(phi);
    B[i][j][k][2] = B_vect.z;
}

void calculate_jacobian
(
    int i, int j, int k,
    double (&J)[NR][NPHI][NZ][NDIMS][NDIMS],
    const double (&B)[NR][NPHI][NZ][NDIMS],
    std::vector<double> plo,
    std::vector<double> dh
)
{
    double R = plo[0] + i*dh[0];
    double phi = plo[1] + j*dh[1];
    double Z = plo[2] + k*dh[2];

    double dBrdr, dBrdphi, dBrdz,
           dBphidr, dBphidphi, dBphidz,
           dBzdr, dBzdphi, dBzdz;    

    if (i == NR-1)
    {   
        dBrdr   = (3*B[i][j][k][0]-4*B[i-1][j][k][0]+B[i-2][j][k][0])/(2*dh[0]);
        dBphidr = (3*B[i][j][k][1]-4*B[i-1][j][k][1]+B[i-2][j][k][1])/(2*dh[0]);
        dBzdr   = (3*B[i][j][k][2]-4*B[i-1][j][k][2]+B[i-2][j][k][2])/(2*dh[0]);
    }
    else if (i == 0)
    {
        dBrdr   = (-3*B[i][j][k][0]+4*B[i+1][j][k][0]-B[i+2][j][k][0])/(2*dh[0]);
        dBphidr = (-3*B[i][j][k][1]+4*B[i+1][j][k][1]-B[i+2][j][k][1])/(2*dh[0]);
        dBzdr   = (-3*B[i][j][k][2]+4*B[i+1][j][k][2]-B[i+2][j][k][2])/(2*dh[0]);
    
    }
    else if (j == NPHI-1)
    {
        dBrdphi   = (3*B[i][j][k][0]-4*B[i][j-1][k][0]+B[i][j-2][k][0])/(2*dh[1]);
        dBphidphi = (3*B[i][j][k][1]-4*B[i][j-1][k][1]+B[i][j-2][k][1])/(2*dh[1]);
        dBzdphi   = (3*B[i][j][k][2]-4*B[i][j-1][k][2]+B[i][j-2][k][2])/(2*dh[1]);
    }
    else if (j == 0)
    {
        dBrdphi   = (-3*B[i][j][k][0]+4*B[i][j+1][k][0]-B[i][j+2][k][0])/(2*dh[1]);
        dBphidphi = (-3*B[i][j][k][1]+4*B[i][j+1][k][1]-B[i][j+2][k][1])/(2*dh[1]);
        dBzdphi   = (-3*B[i][j][k][2]+4*B[i][j+1][k][2]-B[i][j+2][k][2])/(2*dh[1]);
    }
    else if (k == NZ-1)
    {
        dBrdz   = (3*B[i][j][k][0]-4*B[i][j][k-1][0]+B[i][j][k-2][0])/(2*dh[2]);
        dBphidz = (3*B[i][j][k][1]-4*B[i][j][k-1][1]+B[i][j][k-2][1])/(2*dh[2]);
        dBzdz   = (3*B[i][j][k][2]-4*B[i][j][k-1][2]+B[i][j][k-2][2])/(2*dh[2]);
    }
    else if (k == 0)
    {
        dBrdz   = (-3*B[i][j][k][0]+4*B[i][j][k-1][0]-B[i][j][k-2][0])/(2*dh[2]);
        dBphidz = (-3*B[i][j][k][1]+4*B[i][j][k-1][1]-B[i][j][k-2][1])/(2*dh[2]);
        dBzdz   = (-3*B[i][j][k][2]+4*B[i][j][k-1][2]-B[i][j][k-2][2])/(2*dh[2]);
    }
    else
    {
        dBrdr   = (B[i+1][j][k][0]-B[i-1][j][k][0])/(2*dh[0]);
        dBrdphi = (B[i][j+1][k][0]-B[i][j-1][k][0])/(2*dh[1]);
        dBrdz   = (B[i][j][k+1][0]-B[i][j][k-1][0])/(2*dh[2]);

        dBphidr   = (B[i+1][j][k][1]-B[i-1][j][k][1])/(2*dh[0]);
        dBphidphi = (B[i][j+1][k][1]-B[i][j-1][k][1])/(2*dh[1]);
        dBphidz   = (B[i][j][k+1][1]-B[i][j][k-1][1])/(2*dh[2]);

        dBzdr   = (B[i+1][j][k][2]-B[i-1][j][k][2])/(2*dh[0]);
        dBzdphi = (B[i][j+1][k][2]-B[i][j-1][k][2])/(2*dh[1]);
        dBzdz   = (B[i][j][k+1][2]-B[i][j][k-1][2])/(2*dh[2]);
    }

    J[i][j][k][0][0] = dBrdr;   J[i][j][k][0][1] = dBrdphi;   J[i][j][k][0][2] = dBrdz; 
    J[i][j][k][1][0] = dBphidr; J[i][j][k][1][1] = dBphidphi; J[i][j][k][1][2] = dBphidz;   
    J[i][j][k][2][0] = dBzdr;   J[i][j][k][2][1] = dBzdphi;   J[i][j][k][2][2] = dBzdz;   

}