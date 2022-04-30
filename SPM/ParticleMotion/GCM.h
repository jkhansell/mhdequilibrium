#ifndef GCM_H
#define GCM_H


#include <eigen3/Eigen/Dense>

#include "Field.h"

void init_GCM_calc
(
    const Matrix<double, 6,1>& IC_0,
    Matrix<double, 5,1>& solvect,
    const double& mass, 
    const double& charge,
    Field Bfield
);

Matrix<double, 5, 1> GCM_step
(
    Matrix<double, 5,1>& solvect,
    double mass, 
    double charge,
    Field Bfield
);

bool check_particle_domain
(
    const Matrix<double, 5,1>& solvect
);

#endif