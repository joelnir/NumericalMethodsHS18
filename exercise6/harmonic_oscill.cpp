#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <vector>
#include <iostream>
#include "writer.hpp"
#include <assert.h>

using namespace std;

/// Uses the explicit Euler method to compute y from time 0 to time T
/// where y is a 2x1 vector solving the linear system of ODEs as in the exercise
///
/// @param[out] yT at the end of the call, this will have vNext
/// @param[in] y0 the initial conditions
/// @param[in] zeta the damping factor (see exercise)
/// @param[in] h the step size
/// @param[in] T the final time at which to compute the solution.
///
/// The first component of y (the position) will be stored to y1, the second component (the velocity) to y2. The i-th entry of y1 (resp. y2) will contain the first (resp. second) component of y at time i*h.
///

//----------------explicitEulerBegin----------------
void explicitEuler(std::vector<double> & y1, std::vector<double> &y2, std::vector<double> & time,
	const Eigen::Vector2d& y0,
		   double zeta, double h, double T) {

    int timesteps = (int)(T/h);
    y1.clear();
    y2.clear();
    time.clear();

    y1.push_back(y0(0));
    y2.push_back(y0(1));
    time.push_back(0);

    for(int i = 1; i <= timesteps; ++i){
        y1.push_back(y1.at(i-1) + h*y2.at(i-1));
        y2.push_back(y2.at(i-1) + h*(-y1.at(i-1) -2*zeta*y2.at(i-1)));

        time.push_back(i*h);
    }
}
//----------------explicitEulerEnd----------------

// Implements the implicit Euler. Analogous to explicit Euler, same input and output parameters
//----------------implicitEulerBegin----------------
void implicitEuler(std::vector<double> & y1, std::vector<double> & y2, std::vector<double> & time,
	const Eigen::Vector2d& y0,
		   double zeta, double h, double T) {

    int timesteps = (int)(T/h);
    y1.clear();
    y2.clear();
    time.clear();

    y1.push_back(y0(0));
    y2.push_back(y0(1));
    time.push_back(0);

    for(int i = 1; i <= timesteps; ++i){
        y2.push_back((1/(1 + 2*h*zeta + h*h)) * (y2.at(i-1) - h*y1.at(i-1)));
        y1.push_back(y1.at(i-1) + h*y2.at(i));

        time.push_back(i*h);
    }
}
//----------------implicitEulerEnd----------------


//----------------energyBegin----------------
// Energy computation given the velocity. Assume the energy vector to be already initialized with the correct size.
void Energy(const std::vector<double> & v, std::vector<double> & energy)
{
  assert(v.size()==energy.size());

  energy.clear();

  for(double vel : v){
    energy.push_back(0.5*vel*vel);
  }
}
//----------------energyEnd----------------

int main() {

	double T = 20.0;
	double h = 0.5; // Change this for explicit / implicit time stepping comparison
	const Eigen::Vector2d y0(1,0);
	double zeta=0.2;
	std::vector<double> y1;
	std::vector<double> y2;
	std::vector<double> time;
	explicitEuler(y1,y2,time,y0,zeta,h,T);
	writeToFile("position_expl.txt", y1);
	writeToFile("velocity_expl.txt",y2);
	writeToFile("time_expl.txt",time);
	std::vector<double> energy(y2.size());
	Energy(y2,energy);
	writeToFile("energy_expl.txt",energy);
	y1.assign(y1.size(),0);
	y2.assign(y2.size(),0);
	time.assign(time.size(),0);
	energy.assign(energy.size(),0);
	implicitEuler(y1,y2,time,y0,zeta,h,T);
	writeToFile("position_impl.txt", y1);
	writeToFile("velocity_impl.txt",y2);
	writeToFile("time_impl.txt",time);
	Energy(y2,energy);
	writeToFile("energy_impl.txt",energy);

	return 0;
}
