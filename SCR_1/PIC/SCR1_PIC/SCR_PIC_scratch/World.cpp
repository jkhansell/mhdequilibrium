/*defines the simulation domain*/
#include <random>
#include <math.h>
#include "World.h"
#include "Field.h"
#include "Species.h"
#include <iostream>
	
//make an instance of the Rnd class
Rnd rnd;

using namespace std;
using namespace Const;


/*constructor*/
World::World(int nr, int nphi, int nz):
	nr{nr}, nphi{nphi}, nz{nz},	nn{nr,nphi,nz},nn1{nr-1,nphi-1,nz-1},
	phi(nn),j{nn}, rho(nn),node_vol(nn),B(nn1),E(nn), phi_m(nn),b_m(nn), H(nn), 
	M(nn), object_id(nn){
		time_start =  chrono::high_resolution_clock::now();	//save starting time point
	}

/*sets domain bounding box and computes mesh spacing*/
void World::setExtents(double _R0, double _rm, double _zm) {
	/*set origin and the opposite corner*/
	R0 = _R0;
	rm = _rm;
	zm = _zm; 

	x0[0] = R0-rm;
	x0[1] = 0; 
	x0[2] = -zm; 

	dh[0] = 2*rm/nr;
	dh[1] = PI/nphi; 
	dh[2] = 2*zm/nz; 

	a = sqrt(pow(zm,2)+pow(rm,2)); 
	
	/*recompute node volumes*/
	computeNodeVolumes();
}

/*returns elapsed wall time in seconds*/
double World::getWallTime() {
  auto time_now = chrono::high_resolution_clock::now();
  chrono::duration<double> time_delta = time_now-time_start;
  return time_delta.count();
}

/*computes charge density from rho = sum(charge*den)*/
void World::computeChargeDensity(vector<Species> &species)
{
	rho = 0;
	for (Species &sp:species)
	{
		if (sp.charge==0) continue;	//don't bother with neutrals
		rho += sp.charge*sp.den;
	}
}	

void World::computeCurrentDensity(vector<Species> &species){
	j = 0; 
	for (Species &sp:species){
		if (sp.charge == 0) continue; 
		j+=sp.j; 

	}
}

/*computes node volumes, dx*dy*dz on internal nodes and fractional
 * values on domain boundary faces*/
void World::computeNodeVolumes() {
	for (int i=0;i<nr;i++)
		for (int j=0;j<nphi;j++)
			for (int k=0;k<nz;k++){
				double V = i*dh[0]*dh[0]*dh[1]*dh[2];	//default volume
				if (i==0 || i==nr-1) V*=0.5;	//reduce by two for each boundary index
				if (j==0 || j==nphi-1) V*=0.5;
				if (k==0 || k==nz-1) V*=0.5;
				node_vol[i][j][k] = V;
			}
}

/* computes total potential energy from 0.5*eps0*sum(E^2)*/
double World::getPE() {
	double pe = 0;
	for (int i=0;i<nr;i++)
		for (int j=0;j<nphi;j++)
			for (int k=0;k<nz;k++)
			{
				double3 efn = E(i,j,k);	//ef at this node
				double ef2 = efn[0]*efn[0]+efn[1]*efn[1]+efn[2]*efn[2];
				pe += ef2*node_vol[i][j][k];
			}
	return 0.5*Const::EPS_0*pe;
}


/*checks for steady state by comparing change in mass, momentum, and energy*/
bool World::steadyState(vector<Species> &species) {
	// do not do anything if already at steady state
	if (steady_state) return true;

	double tot_mass = 0;
	double tot_mom = 0;
	double tot_en = getPE();
	for (Species &sp:species)
	{
		tot_mass += sp.getRealCount();	//number of real molecules
		double3 mom = sp.getMomentum();
		tot_mom += mag(mom);		//z-component of momentum
		tot_en += sp.getKE();		//add kinetic energy
	}

	/*compute new values to last*/
	const double tol = 1e-3;
	if (abs((tot_mass-last_mass)/tot_mass)<tol &&
		abs((tot_mom-last_mom)/tot_mom)<tol &&
		abs((tot_en-last_en)/tot_en)<tol) {
		steady_state = true;
		cout<<"Steady state reached at time step "<<ts<<endl;
	}

	/*update prior values*/
	last_mass = tot_mass;
	last_mom = tot_mom;
	last_en = tot_en;
	return steady_state;
}