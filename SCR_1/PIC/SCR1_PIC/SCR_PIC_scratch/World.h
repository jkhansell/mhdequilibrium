#ifndef _WORLD_H
#define _WORLD_H

#include <vector>
#include <random>
#include <chrono>
#include "Field.h"
#include <cmath>

class Species;

/*define constants*/
namespace Const
{
	const double EPS_0 = 8.85418782e-12;  	// C/(V*m), vacuum permittivity
	const double QE = 1.602176565e-19;		// C, electron charge
	const double AMU = 1.660538921e-27;		// kg, atomic mass unit
	const double ME = 9.10938215e-31;		// kg, electron mass
	const double K = 1.380648e-23;			// J/K, Boltzmann constant
	const double PI = 3.141592653;			// pi
	const double C = 3.0e8;					// speed of light
	const double EvToK = QE/K;				// 1eV in K ~ 11604
	const double MU_0 = 4*PI*1e-7;			// vacuum permisivity
}

/*object for sampling random numbers*/
class Rnd {
	public:
		//constructor: set initial random seed and distribution limits
		Rnd(): mt_gen{std::random_device()()}, rnd_dist{0,1.0} {}
		double operator() () {return rnd_dist(mt_gen);}

	protected:
		std::mt19937 mt_gen;	    //random number generator
		std::uniform_real_distribution<double> rnd_dist;  //uniform distribution
};

extern Rnd rnd;		//tell the compiler that an object of type Rnd called rnd is defined somewhere

/*defines the computational domain*/
class World{
	public:	
		/*constructor, allocates memory*/
		World(int nr, int nphi, int nz);

		/*functions to set mesh origin and spacing*/
		void setExtents(const double R0, const double rm, 
					const double zm);
		
		double getR0() const {return double(R0);}
		double getRm() const {return double(rm);}
		double getZm() const {return double(zm);}
		double3 getDh() const {return double3(dh);}
		

		/*functions for accessing time information*/
		int getTs() const {return ts;}
		double getTime() const {return time;}
		double getWallTime();  /*returns wall time in seconds*/
		double getDt() const {return dt;}
		bool isLastTimeStep() const {return ts==num_ts-1;}

		bool inBounds(double3 pos) {

			double R = sqrt(pow(pos[0],2)+pow(pos[2],2));

			if (R >= a) {return true;}
			return false; 
			
		}

		/*sets time step and number of time steps*/
		void setTime(double dt, int num_ts) {this->dt=dt;this->num_ts=num_ts;}
		
		/*advances to the next time step, returns true as long as more time steps remain*/
		bool advanceTime() {time+=dt;ts++;return ts<=num_ts;}

				/*checks and sets a steady state flag*/
		bool steadyState(std::vector<Species> &species);

		/*returns steady state flag*/
		bool isSteadyState() {return steady_state;}

		/*converts physical position to logical coordinate*/
		double3 XtoL(double3 r) const {
			double3 lc;
				lc[0] = (r[0]-x0(0))/dh(0);
				lc[1] = (r[1]-x0(1))/dh(1);
				lc[2] = (r[2]-x0(2))/dh(2);
				return lc;
		}

		int3 RtoIJK(const double3 &x) const{
			double3 lc = XtoL(x); 
			int3 ijk {(int)lc[0],(int)lc[0],(int)lc[0]};
			return ijk; 
		}

		double3 pos(double3 lc) {
			double3 x = x0 + dh*lc; 
			return x;  
		}

		/*computes charge density from rho = sum(charge*den)*/
		void computeChargeDensity(std::vector<Species> &species);

		/*computes current density from j = sum(j_sp)*/
		void computeCurrentDensity(std::vector<Species> &species);

		/*returns the system potential energy*/
		double getPE();

		int U(int i,int j, int k) {return object_id.U(i,j,k);}

		//mesh geometry
		const int nr,nphi,nz;	//number of nodes
		const int3 nn;	//another way to access node counts
		const int3 nn1; 

		Field phi;			//potential
		Field rho;			//charge density
		Field3 j; 			//total current density
		Field node_vol;		//node volumes
		Field3 E;			//electric field
		Field3 B;			//magnetic field
		FieldI object_id;	//object id flat to flag fixed nodes 			

		/*magnetic solver support*/
		Field phi_m;		//magnetic potential	
		Field3 H; 			//H field
		Field3 M; 			//magnetization
		Field b_m; 			
	protected:
		double R0;			//major radius of 
		double rm;			//width of solution space
		double zm;			//height of solution space
		double a;			//radius of solution space -> a^2 = rm^2 + zm^2 
		
		double3 dh;  		//cylindrical spacing
		double3 xm;			//origin-diagonally opposite corner (max bound)
		double3 xc;			//domain centroid
		double3 x0;
		double dt = 0;		//time step length
		double time = 0;	//physical time
		int ts = -1;		//current time step
		int num_ts = 0;		//number of time steps

		std::chrono::time_point<std::chrono::high_resolution_clock> time_start;	//time at simulation start

		bool steady_state = false; 
		double last_mass = 0;	//mass at the prior time step
		double last_mom = 0;	//momentum at the prior time step
		double last_en = 0;		//energy at the prior time step
		void computeNodeVolumes();
};

#endif
