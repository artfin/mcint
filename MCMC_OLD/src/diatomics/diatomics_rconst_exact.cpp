#include <iostream>
#include <random>

using namespace std;

const double BOLTZCONST = 1.38064e-23;
// unified atomic mass units to kg
const double DALTON = 1.660539e-27;
// atomic length unit
const double ALU = 5.29177e-11;

// reduced mass of ar and co2 = m(ar) * m(co2) / (m(ar) + m(co2)) in kg
const double MU = 20.952 * DALTON;

// hbar
const double HBAR = 1.0545718e-34;

const double temperature = 300;

// -----------------------------------------
// distance between atoms in ALU
const double RDIST = 20.0;
// -----------------------------------------

static thread_local std::mt19937 generator;

// normally distributed random vector ( = \dot{\vec{r}} )
double nextGaussian( const double &mean, const double &sigma )
{
    normal_distribution<double> d( mean, sigma );
    return d( generator ); 
} 

int main( int argc, char* argv[] )
{
	if ( argc != 2 )
	{
		cout << "USAGE: ./... (int) n" << endl;
		exit( 1 );
	}

	int n = atoi( argv[1] );
	cerr << "Given value of n: " << n << endl;

	cout << "# Diatomics exact sampler. " << endl;
	cout << "# Temperature: " << temperature << endl;

	// nextGaussian take mu and SIGMA!!! not sigma**2!! 

	for ( int i = 0; i < n; i++ )
	{
		long double jx = nextGaussian( 0, sqrt(BOLTZCONST * temperature * MU * pow(RDIST * ALU, 2)) ) / HBAR;
		long double jy = nextGaussian( 0, sqrt(BOLTZCONST * temperature * MU * pow(RDIST * ALU, 2))  ) / HBAR;
		long double pR = nextGaussian( 0, sqrt(BOLTZCONST * temperature * MU) ) / HBAR * ALU;

		cout << jx << "  " << jy << "  " << pR << endl;
	}

	return 0;
}




