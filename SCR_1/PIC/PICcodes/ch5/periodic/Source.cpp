#include "Source.h"


//samples monoenergetic particles according to a prescribed density
//the source area is reduced to [0.4,0.6]*Lx, [0.4,0.6]*Ly
void ColdBeamSource::sample()
{
	double3 dh = world.getDh();
	double3 x0 = world.getX0();

	//area of the XY plane, A=Lx*Ly
	double Lx = dh[0]*(world.ni-1);
	double Ly = dh[1]*(world.nj-1);
	double A = 0.2*Lx*0.2*Ly;

	//compute number of real particles to generate: (#/s) = n*v*A; # = (#/s)*dt
	double num_real = den*v_drift*A*world.getDt();

	//number of simulation particles
	int num_sim = (int)(num_real/sp.mpw0+rnd());

	//inject particles
	for (int i=0;i<num_sim;i++)
	{
		double3 pos {x0[0]+(0.4+0.2*rnd())*Lx, x0[1]+(0.4+0.2*rnd())*Ly, x0[2]};
		double3 vel {0,0,v_drift};
		sp.addParticle(pos,vel);
	}
}


