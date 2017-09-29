#include <iostream>
#include <chrono>

#include <math.h>
#include <cstddef>
#include <vector>

#include <fstream>
#include <iomanip>

#include "hep/mc.hpp"

// hamiltonian input: theta(co2), phi, theta(h2), r
#include "co2h2_hamiltonian.hpp"
// double kinetic_energy(double q1, double q2, double q3, double q4, double p1, double p2, double p3, double p4, double Jx, double Jy, double Jz);

// dipole input: r, theta(co2), theta(h2), phi
#include "co2h2_dipole_lr.hpp"
//double dipx(double R, double Theta1, double Theta2, double Phi);
//double dipy(double R, double Theta1, double Theta2, double Phi);
//double dipz(double R, double Theta1, double Theta2, double Phi);

// potential: xr -- in angstroms
// pararead -- intializing procedure
// output of co2h2pes -- in cm^-1
// input: r, theta(co2), theta(h2), phi
extern "C" void pararead( void );
extern "C" void co2h2pes( double *xr, double *xth1, double *xth2, double *xphi, double *potvalue );

using namespace std;
using namespace std::placeholders;


int main( int argc, char* argv[] )
{
	pararead();
	double r = 5.0;
	double th1 = 0.5;
	double th2 = 0.5;
	double p = 0.2;
	double v;
	co2h2pes( &r, &th1, &th2, &p, &v);
	cout << "v: " << v << endl;

	return 0;
}
