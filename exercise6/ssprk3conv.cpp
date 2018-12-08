#include <vector>
#include <iostream>
#include <cmath>
#include "writer.hpp"

using namespace std;

double f(double t, double u) {
    return std::exp(-2*t) - 2*u;
}


/// Uses the SSP RK3 method to compute u from time 0 to time T
/// for the ODE $u'=e^{-2t}-2u$
///
/// @param[out] u at the end of the call (solution at all time steps)
/// @param[out] time contains the time levels
/// @param[in] u0 the initial data
/// @param[in] dt the step size
/// @param[in] T the final time up to which to compute the solution.
///

//----------------SSPRK3Begin----------------
void SSPRK3(std::vector<double> & u, std::vector<double> & time,
          const double & u0, double dt, double T) {
    const unsigned int nsteps = std::round(T/dt);
    u.resize(nsteps+1);
    time.resize(nsteps+1);

    // Initial conditions
    u[0] = u0;
    time[0] = 0;

    double k1, k2, k3;
    for(int i = 1; i <= nsteps; ++i){
        time[i] = dt*i;

        k1 = f(time[i-1], u[i-1]);
        k2 = f(time[i], u[i-1]+dt*k1);
        k3 = f(time[i-1] + 0.5*dt, u[i-1]+0.25*dt*(k1 + k2));

        u[i] = u[i-1] + (dt/6.0)*(k1 + k2 + 4*k3);
    }
}
//----------------SSPRK3End----------------

int main(int argc, char** argv) {

    double T = 10.0;

    const double u0 = 0.;
    std::vector<double> error;
    std::vector<double> steps;

    std::vector<double> u;
    std::vector<double> time;

    double dt;
    for(int k = 1; k <= 8; ++k){
        dt = (1.0/pow(2, k));
        SSPRK3(u,time,u0,dt,T);

        steps.push_back(dt);
        error.push_back(abs(u[round(T/dt)] - T*exp(-2*T)));
    }

    writeToFile("error.txt", error);
    writeToFile("dt.txt",steps);

    return 0;
}
