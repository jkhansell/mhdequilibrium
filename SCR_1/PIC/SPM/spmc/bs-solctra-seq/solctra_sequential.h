//
// Created by lchavarr on 4/17/16.
//

#ifndef SOLCTRA_SOLCTRA_H
#define SOLCTRA_SOLCTRA_H

#include <string>
#include <cmath>
#include "utils.h"
#include <sstream>


void load_coil_data(double* x, double* y, double* z, 
                    const std::string& path);

void e_roof(GlobalData& data);

void R_vectors(const Coil& coil, 
                const cartesian& point, 
                Coil* Rmi, Coil* Rmf);

cartesian magnetic_field(Coil* rmi, Coil* rmf, 
                        const GlobalData& data, 
                        const cartesian& point);

void initializeGlobals(Coil* rmi, Coil* rmf);

void finishGlobals(Coil* rmi, Coil* rmf);

void getMagneticProfile(const GlobalData& data, const int num_points, 
                        const int phi_angle, const std::string& output, 
                        const int dimension);

void RK4(const GlobalData& data, const std::string& output, 
            const cartesian& start_point, const int steps, 
            const double& step_size, const int particle, 
            const int mode, const int print_type);

inline double norm_of(const cartesian& vec){
    return sqrt(( vec.x * vec.x ) + ( vec.y * vec.y ) + ( vec.z * vec.z ));
}

void runParticles(const GlobalData& data, const std::string& output, 
                    const Coil& particles, const int length, const int steps, 
                    const double& step_size, const int mode, const int print_type);


#endif //SOLCTRA_SOLCTRA_H
