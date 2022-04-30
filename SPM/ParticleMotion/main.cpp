#include <vector>
#include <iostream>
#include <fstream>
#include <eigen3/Eigen/Dense>

using namespace Eigen; 

#include "GCM.h"
#include "RK4.h"

#define q_p 1.60217662e-19
#define m_p 1.6726219e-27
#define Nsteps 2
#define dt 1e-10

int main(int argc, char* argv[])
{
    Field SCRmagfield;
    SCRmagfield.read_file();

    Matrix<double, 6, 1> init_conditions = Matrix<double, 6, 1>::Zero();
    Matrix<double, 5, 1> solvect = Matrix<double, 5, 1>::Zero();

    init_conditions(0) = 0.24; //x-coord
    init_conditions(1) = 0.0; //y-coord
    init_conditions(2) = 0.0; //z-coord
    init_conditions(3) = 0.0; //vx
    init_conditions(4) = 1e5; //vy
    init_conditions(5) = 0.0; //vz

    init_GCM_calc(init_conditions, solvect, m_p, q_p, SCRmagfield);

    //  std::cout << solvect << "\n";    
    std::ofstream outfile; 
    outfile.open("output.txt");
    outfile << "Step\t" <<"time\t" << "x\t" << "y\t" << "z\t" << "vpar\t" << "mu\n";
    outfile.close();
    RK4_loop(solvect, dt, Nsteps, m_p, q_p, SCRmagfield, outfile);
    
    return 0;
}