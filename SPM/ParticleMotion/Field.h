#ifndef FIELD_H
#define FIELD_H


#include <vector>
#include <eigen3/Eigen/Dense>
#include <netcdf.h>
#include <iostream>


#define FILE_NAME1 "SCR1FIELD.nc"
#define FILE_NAME2 "Jacobian.nc"

#define NDIMS 3
#define NCOMP 4
#define NR 64
#define NPHI 360
#define NZ 64
#define R0 0.2477
#define a 0.08
#define PI 3.14159265358979323846

#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}

using namespace Eigen;

class Field
{
    public:
        Field()
        {
            dh(0) = 2*a/NR;
            dh(1) = 2*PI/NPHI;
            dh(2) = 2*a/NZ;
            pHi(0) = R0+a;
            pHi(1) = 2*PI;
            pHi(2) = a;
            pLo(0) = R0-a;
            pLo(1) = 0;
            pLo(2) = -a;
        }
        ~Field(){}

        void read_file(); 

        void trilin_interp
        (
            const Vector3d& pos, 
            Vector3d& B, 
            Matrix3d& jacobian
        );

        Matrix3d transformjac(const Matrix3d& cyl_jac, double phi);

        Vector3d transformvect(const Vector3d& cyl_vect, double phi);
    
    private:
        static double Bfield[NR][NPHI][NZ][NDIMS]; 
        static double Jac[NR][NPHI][NZ][NDIMS][NDIMS];
        Vector3d dh;
        Vector3d pHi;
        Vector3d pLo; 
};

#endif