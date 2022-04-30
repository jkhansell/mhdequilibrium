#ifndef RK4_H
#define RK4_H

#include <eigen3/Eigen/Dense>
#include <iostream>
#include <fstream>

#include "GCM.h"

using namespace Eigen;

bool RK4_step
(
    Matrix<double, 5,1>& solvect, 
    const double dt,
    const double mass, 
    const double charge,
    Field Bfield
);

void RK4_loop
(
    Matrix<double, 5,1>& solvect, 
    const double dt,
    const double Nsteps,
    const double mass, 
    const double charge,
    Field Bfield, 
    std::ofstream& File
);

#endif
