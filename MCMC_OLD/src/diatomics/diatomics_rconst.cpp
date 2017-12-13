#include <iostream>
#include <random>
#include <chrono>
#include <cmath>

#include <iomanip> // std::atoi
#include <algorithm> // std::min

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// boltzmann constant
const double BOLTZCONST = 1.38064e-23;
// unified atomic mass units to kg
const double DALTON = 1.660539e-27;
// atomic length unit
const double ALU = 5.29177e-11;
// atomic mass unit
const long double AMU = 9.1093826e-31;
// hartree to joules
const double HTOJ = 4.35974417e-18;

// reduced mass of ar and co2 = m(ar) * m(co2) / ( m(ar) + m(co2) ) in kg
const double MUKG = ( 40.0 * 44.0 ) / ( 40.0 + 44.0 ) * DALTON;
const long double MUAMU = MUKG / AMU; 

// temperature in K
const double temperature = 300;

static mt19937 generator;

// !--------------------------------
// distance between two atoms in ALU
const double RDIST = 20.0;
// !--------------------------------

static double nextDouble( const double &min = 0.0, const double &max = 1.0 )
{
    uniform_real_distribution<double> distribution( min, max );
    return distribution( generator );
}

static void nextGaussianVec( Vector3d &v, Vector3d mean, const double sigma )
{
    normal_distribution<double> d0( mean(0), sigma );
    normal_distribution<double> d1( mean(1), sigma );
    normal_distribution<double> d2( mean(2), sigma );

    v(0) = d0( generator );
    v(1) = d1( generator );
    v(2) = d2( generator );
}

// x = [jx, jy, pR]
double target( Vector3d x )
{
    double jx = x(0);
    double jy = x(1);
    double pR = x(2);
    double h = pow(pR, 2) / ( 2 * MUAMU ) + ( pow(jx, 2) + pow(jy, 2) ) / (2 * MUAMU * pow(RDIST, 2));

    return exp( - h * HTOJ / ( BOLTZCONST * temperature ));
}

Vector3d metro_step( Vector3d x, double alpha )
{
    Vector3d prop;

    // generate random vector
    nextGaussianVec( prop, x, alpha );

    if ( nextDouble() < min( 1.0, target(prop) / target(x) ))
    {
        x = prop;
    }

    return x;
}

int main( int argc, char* argv[] )
{
    if ( argc != 5 )
    {
        cout << "USAGE: ./... (int) moves-to-be-made (int) burn-in-steps (double) alpha (bool) show_vecs" << endl;
        exit( 1 );
    }

    const int nsteps = atoi( argv[1] );
    const int burnin = atoi( argv[2] );
    const double alpha = atof( argv[3] ); 
    const bool show_vecs = atoi( argv[4] );

    setprecision(3);

    cerr << "-----" << endl;
    cerr << "Input parameters: " << endl;
    cerr << ">> moves-to-be-made: " << nsteps << endl;
    cerr << ">> burn-in: " << burnin << endl;
    cerr << ">> alpha: " << alpha << endl;
    cerr << ">> show_vecs: " << show_vecs << endl;

    Vector3d x ( 10.0, 10.0, 10.0 );
    Vector3d xnew;

    cout << "# Metropolis-Hastings sampler for diatomics" << endl;
    cout << "# moves-to-be-made = " << nsteps << endl;
    cout << "# burn-in = " << burnin << endl;
    cout << "# alpha = " << alpha << endl;
    cout << "# temperature = " << temperature << endl;
    cout << "# R (a.u.) = " << RDIST << endl;  
    
	chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();	

	// burn-in
	for ( size_t i = 0; i < burnin; i++ )
	{
		x = metro_step( x, alpha );
	}
	
	cerr << "Burn-in finished. Time elapsed: " << chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - startTime).count() / 1000.0 << " s" << endl;

	int moves = 0;
	int attempted_moves = 0;

	// recording
	while ( moves < nsteps )
    {
        xnew = metro_step( x, alpha );

        if ( xnew != x )
        {
            moves++;
        }
        
		if ( show_vecs == true )
        {
            cout << x(0) << " " << x(1) << " " << x(2) << endl;
        }

        x = xnew;
		attempted_moves++;
	}

    cerr << "-----------------------------------" << endl;
    cerr << "Attempted steps: " << attempted_moves << "; moves: " << moves << "; percent: " << (double) moves / attempted_moves * 100.0 << "%" << endl;
    cerr << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - startTime).count() / 1000.0 << " s" << endl; 
    cerr << "-----------------------------------" << endl;

    return 0;
}
