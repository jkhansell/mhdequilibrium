#include "RK4.h"

bool RK4_step
(
    Matrix<double, 5,1>& solvect, 
    const double dt,
    const double mass, 
    const double charge,
    Field Bfield
)
{
    bool INDOMAIN = true;
    Matrix<double, 5,5> Kmat;
    Matrix<double, 5,1> tempsol;
    Vector4d dtcoefs = {0.5,0.5,1,1};
    tempsol = solvect;

    for (int i = 0; i < 4; i++)
    {
        /*
        INDOMAIN = check_particle_domain(tempsol);
        if(!INDOMAIN)
        {
            std::cout << "Particle left.\n";
            return INDOMAIN;
        }
        */

        //std::cout << tempsol << "\n";

        Kmat.col(i) = GCM_step(tempsol, mass, charge, Bfield);
        //std::cout << Kmat.col(i) << "\n--------\n";
        tempsol = solvect + dtcoefs[i]*dt*Kmat.col(i);   
        //std::cout << Kmat.col(i); 
    }

    Kmat.col(4) = (Kmat.col(0)+2*Kmat.col(1)+2*Kmat.col(2)+Kmat.col(3))/6;

    solvect += dt*Kmat.col(4);
    //std::cout << Kmat.col(4) << "\n";
    //std::cout << "---------\n";

    return INDOMAIN; 
}

void RK4_loop
(
    Matrix<double, 5,1>& solvect, 
    const double dt,
    const double Nsteps,
    const double mass, 
    const double charge,
    Field Bfield, 
    std::ofstream& File
)
{
    bool INDOMAIN = true; 
    int iStep = 0;

    while (INDOMAIN && iStep < Nsteps)
    {
        File.open("output.txt", std::ios_base::app);
        File << iStep << "\t" << dt*iStep << "\t" << solvect(1) << "\t" << solvect(2)  << "\t" << solvect(3)  << "\t" << solvect(0)  << "\t" << solvect(4) << "\n";      
        std::cout << solvect << "\n";
        
        INDOMAIN = RK4_step(solvect, dt, mass, charge, Bfield);

        

        iStep++; 
    }

    File.close();

}