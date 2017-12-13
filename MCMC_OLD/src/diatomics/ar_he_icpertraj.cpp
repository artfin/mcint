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
// hartree to joules
const double HTOJ = 4.35974417e-18;

// reduced mass
const double HE_MASS = 4.00260325413;
const double AR_MASS = 39.9623831237; 
const double PROTON_TO_ELECTRON_RATIO = 1836.15267389; 

const double MU = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * PROTON_TO_ELECTRON_RATIO; 

// temperature in K
const double temperature = 300;

static mt19937 generator;

// !--------------------------------
// distance between two atoms in ALU
const double RDIST = 50.0;
// !--------------------------------

const double BBOUND = 1.0;

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

// x = [pR, theta, pT]
double target( Vector3d x )
{
    double pR = x(0);
	double theta = x(1);
    double pT = x(2);

	double h = pow(pR, 2) / ( 2 * MU ) + ( pow(pT, 2) ) / (2 * MU * pow(RDIST, 2));

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

double calculate_gunsight( Vector3d x )
{
	double pR = x(0);
	double theta = x(1);
	double pT = x(2);

	return pT / pR;
}

// wrap x -> [0, max)
double wrapMax( double x, double max )
{
	return fmod( max + fmod(x, max), max );
}

int main( int argc, char* argv[] )
{
    if ( argc != 5 )
    {
        cout << "USAGE: ./... (int) vectors-to-write (int) burn-in-steps (double) alpha (bool) show_vecs" << endl;
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

    Vector3d x ( 1.0, 0.0, 1.0 );
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

	int wrote_vectors = 0;
	double b;

	// recording
	while ( wrote_vectors < nsteps )
    {
        xnew = metro_step( x, alpha );

        if ( xnew != x )
        {
            moves++;
        }
      
		xnew(1) = wrapMax( xnew(1), 2 * M_PI );

		// gunsight parameter
		b = calculate_gunsight( xnew );

		// if pR < 0
		if ( show_vecs == true && xnew != x && 
		  	 xnew(0) < 0 && b > 0 && b < BBOUND )
        {
            cout << wrote_vectors + 1 << " " << RDIST << " " << xnew(0) << " " << xnew(1) << " " << xnew(2) << endl;
        
			wrote_vectors++;
		}
		
		attempted_moves++;

        x = xnew;
	}

    cerr << "-----------------------------------" << endl;
    cerr << "Attempted steps: " << attempted_moves << "; moves: " << moves << "; percent: " << (double) moves / attempted_moves * 100.0 << "%" << endl;
	cerr << "Vectors wrote: " << wrote_vectors << "; percentage (to attempted moves): " << (double) wrote_vectors / attempted_moves * 100 << "%; percentage (to moves made): " << (double) wrote_vectors / moves * 100 << "%" << endl;
    cerr << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - startTime).count() / 1000.0 << " s" << endl; 
    cerr << "-----------------------------------" << endl;

    return 0;
}
