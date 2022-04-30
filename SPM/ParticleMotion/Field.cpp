#include <vector>

#include "Field.h"

void Field::trilin_interp
(
    const Vector3d& pos, 
    Vector3d& B,
    Matrix3d& jacobian
)
{

    Matrix<double, 3,1> l;
    Matrix<int, 3,1> ind;

    double R = sqrt(pow(pos(0),2)+pow(pos(1),2)); 
    double phi = atan2(pos(1), pos(0));
    if (phi < 0){phi += 2*PI;}
    
    //std::cout << R << "\n";
    //std::cout << phi << "\n";
    //std::cout << pos(2)<< "\n";
    
    l(0) = (R-pLo(0))/dh(0);
    l(1) = (phi-pLo(1))/dh(1);
    l(2) = (pos(2)-pLo(2))/dh(2);

    //std::cout << "-------\n";
    ind(0) = floor(l[0]);
    ind(1) = floor(l[1]);
    ind(2) = floor(l[2]);

    std::cout << ind << "\n"; 

    //std::cout << ind << "\n";

    double dr = l[0]-ind[0];
    double dphi = l[1]-ind[1];
    double dz = l[2]-ind[2];

    double sr[] = {1.-dr, dr};
    double sphi[] = {1.-dphi, dphi};
    double sz[] = {1.-dz, dz};

    //trilinear interpolation algorithm

    constexpr int irmin = 0;
    constexpr int irmax = 1;
    constexpr int iphimin = 0;
    constexpr int iphimax = 1;
    constexpr int izmin = 0;
    constexpr int izmax = 1;

    B = Vector3d::Zero();
    jacobian = Matrix3d::Zero();

    for (int kk = izmin; kk <= izmax; ++kk){
        for (int jj = iphimin; jj <= iphimax; ++jj){
            for (int ii = irmin; ii <= irmax; ++ii){
                
                B(0) += sr[ii]*sphi[jj]*sz[kk]*Bfield[ind(0)+ii][ind(1)+jj][ind(2)+kk][0];
                B(1) += sr[ii]*sphi[jj]*sz[kk]*Bfield[ind(0)+ii][ind(1)+jj][ind(2)+kk][1];
                B(2) += sr[ii]*sphi[jj]*sz[kk]*Bfield[ind(0)+ii][ind(1)+jj][ind(2)+kk][2];

                jacobian(0,0) = sr[ii]*sphi[jj]*sz[kk]*Jac[ind(0)+ii][ind(1)+jj][ind(2)+kk][0][0];
                jacobian(0,1) = sr[ii]*sphi[jj]*sz[kk]*Jac[ind(0)+ii][ind(1)+jj][ind(2)+kk][0][1];
                jacobian(0,2) = sr[ii]*sphi[jj]*sz[kk]*Jac[ind(0)+ii][ind(1)+jj][ind(2)+kk][0][2];

                jacobian(1,0) = sr[ii]*sphi[jj]*sz[kk]*Jac[ind(0)+ii][ind(1)+jj][ind(2)+kk][1][0];
                jacobian(1,1) = sr[ii]*sphi[jj]*sz[kk]*Jac[ind(0)+ii][ind(1)+jj][ind(2)+kk][1][1];
                jacobian(1,2) = sr[ii]*sphi[jj]*sz[kk]*Jac[ind(0)+ii][ind(1)+jj][ind(2)+kk][1][2];

                jacobian(2,0) = sr[ii]*sphi[jj]*sz[kk]*Jac[ind(0)+ii][ind(1)+jj][ind(2)+kk][2][0];
                jacobian(2,1) = sr[ii]*sphi[jj]*sz[kk]*Jac[ind(0)+ii][ind(1)+jj][ind(2)+kk][2][1];
                jacobian(2,2) = sr[ii]*sphi[jj]*sz[kk]*Jac[ind(0)+ii][ind(1)+jj][ind(2)+kk][2][2];

            }
        }
    }

    //std::cout << B << "\n";
    //std::cout << jacobian << "\n";

}

void Field::read_file()
{
    int ncid, varid, retval;

    if ((retval = nc_open(FILE_NAME1, NC_NOWRITE, &ncid)))
        ERR(retval);

    if ((retval = nc_inq_varid(ncid, "data", &varid)))
        ERR(retval);
    
    if ((retval = nc_get_var_double(ncid, varid,(&Bfield)[0][0][0][0])))
        ERR(retval);

    if ((retval = nc_close(ncid)))
        ERR(retval);

    printf("*** SUCCESS reading magnetic field %s!\n", FILE_NAME1);

    
    if ((retval = nc_open(FILE_NAME2, NC_NOWRITE, &ncid)))
        ERR(retval);

    if ((retval = nc_inq_varid(ncid, "data", &varid)))
        ERR(retval);
    
    if ((retval = nc_get_var_double(ncid, varid,(&Jac)[0][0][0][0][0])))
        ERR(retval);

    if ((retval = nc_close(ncid)))
        ERR(retval);

    printf("*** SUCCESS reading Jacobian %s!\n", FILE_NAME1);

    //std::cout << [0][0][0][0][0] << "\n";
}

Matrix3d Field::transformjac(const Matrix3d& cyl_jac, double phi)
{
    Matrix3d Q;
    Q(0,0) = cos(phi); Q(0,1) =-sin(phi); Q(0,2) = 0;
    Q(1,0) = sin(phi); Q(1,1) = cos(phi); Q(1,2) = 0; 
    Q(2,0) = 0;        Q(2,1) = 0;        Q(2,2) = 1;  

    return Q*cyl_jac*Q.transpose(); // check stuff
}

Vector3d Field::transformvect(const Vector3d& cyl_vec, double phi)
{
    Matrix3d Q;
    Vector3d b;
    Q(0,0) = cos(phi); Q(0,1) =-sin(phi); Q(0,2) = 0;
    Q(1,0) = sin(phi); Q(1,1) = cos(phi); Q(1,2) = 0; 
    Q(2,0) = 0;        Q(2,1) = 0;        Q(2,2) = 1; 

    b = Q*cyl_vec; 
    return b;
}

double Field::Bfield[NR][NPHI][NZ][NDIMS];
double Field::Jac[NR][NPHI][NZ][NDIMS][NDIMS];