#include "Source.h"


//samples monoenergetic particles according to a prescribed density
void ColdBeamSource::sample()
{
	double3 dh = world.getDh();
	double3 x0 = world.getX0();

	//area of the XY plane, A=Lx*Ly
	double Lx = dh[0]*(world.ni-1);
	double Ly = dh[1]*(world.nj-1);
	double A = Lx*Ly;

	//compute number of real particles to generate: (#/s) = n*v*A; # = (#/s)*dt
	double num_real = den*v_drift*A*world.getDt();

	//number of simulation particles
	int num_sim = (int)(num_real/sp.mpw0+rnd());

	//inject particles
	for (int i=0;i<num_sim;i++)
	{
		double3 pos {x0[0]+rnd()*Lx, x0[1]+rnd()*Ly, x0[2]};
		double3 vel {0,0,v_drift};
		sp.addParticle(pos,vel);
	}
}

//samples monoenergetic particles according to a prescribed density
void WarmBeamSource::sample()
{
	double3 dh = world.getDh();
	double3 x0 = world.getX0();

	//area of the XY plane, A=Lx*Ly
	double Lx = dh[0]*(world.ni-1);
	double Ly = dh[1]*(world.nj-1);
	double A = Lx*Ly;

	//compute number of real particles to generate: (#/s) = n*v*A; # = (#/s)*dt
	double num_real = den*v_drift*A*world.getDt();

	//number of simulation particles
	int num_sim = (int)(num_real/sp.mpw0+rnd());

	//inject particles
	for (int i=0;i<num_sim;i++)
	{
		double3 pos {x0[0]+rnd()*Lx, x0[1]+rnd()*Ly, x0[2]};

        double theta = 2*Const::PI*rnd();

        double R = -1.0+2*rnd();    //pick a random direction for n[2]
        double a = sqrt(1-R*R);

        double n[3];
        n[0] = cos(theta)*a;
        n[1] = sin(theta)*a;
        n[2] = R;

        double v_th = sp.sampleVth(T);
        double3 vel;
        vel[0] = v_th*n[0];
        vel[1] = v_th*n[1];
        vel[2] = v_th*n[2] + v_drift;

        //reverse if going in wrong direction
        if (vel[2]<0) vel[2]=-vel[2];

		sp.addParticle(pos,vel);
	}
}


    /*load isotropic thermal distribution,
        http://www.particleincell.com/blog/2012/isotropic-velocity/*/
        /*pick a random angle*/
