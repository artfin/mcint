#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

#include <iomanip> // std::atoi
#include <algorithm> // std::min

#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

const double mu1 = 14579.0; // ?
const double mu2 = 36440.0; // ?
const double l = 4.398;
const double l2 = pow( l, 2 );

// boltzmann constant
const double BOLTZCONST = 1.38064e-23;
// hartree to joules
const double HTOJ = 4.35974417e-18;

// ! ---------------------------
const double RDIST = 20.0;
const double RDIST2 = pow( RDIST, 2 );
// -----------------------------

const int DIM = 6;

// temperature in K
const double temperature = 300;

static mt19937 generator;

static double nextDouble( const double &min = 0.0, const double &max = 1.0 )
{
    uniform_real_distribution<double> distribution( min, max );
    return distribution( generator );
}

static void nextGaussianVec( VectorXf &v, VectorXf mean, const double sigma )
{
    for ( int i = 0; i < DIM; i++ )
    {
        normal_distribution<double> d( mean(i), sigma );
        v(i) = d( generator );
    }
} 

// wrap x -> [0, max)
double wrapMax( double x, double max )
{
	return fmod( max + fmod(x, max), max );
}

// x = [ theta, rdot, theta_dot, omega_x, omega_y, omega_z ]
double target( VectorXf x )
{
	double Theta = x[0];
	double rdot = x[1];
	double Theta_dot = x[2];
	double Omega_x = x[3];
	double Omega_y = x[4];
	double Omega_z = x[5];

	double rdot2 = pow( rdot, 2 );
	double Theta_dot2 = pow( Theta_dot, 2 );
	double sinTheta = sin( Theta );
	double cosTheta = cos( Theta );
	double sinTheta2 = pow( sinTheta, 2 );
	double cosTheta2 = pow( cosTheta, 2 );

	double kinetic_part = 0.5 * mu2 * rdot2 + 0.5 * mu1 * l2 * Theta_dot2;
	double angular_part = 0.5 * ( mu1 * l2 * cosTheta2 + mu2 * RDIST2) * pow(Omega_x, 2) + 0.5 * ( mu1 * l2 + mu2 * RDIST2 ) * pow(Omega_y, 2) + 0.5 * mu1 * l2 * sinTheta2 * pow(Omega_z, 2) - mu1 * l2 * sinTheta * cosTheta * Omega_x * Omega_z;
	double coriolis_part = mu1 * l2 * Omega_y * Theta_dot;

    double ke = kinetic_part + angular_part + coriolis_part; 

	return exp( -ke * HTOJ / (BOLTZCONST * temperature ));
}

VectorXf metro_step( VectorXf x, double alpha )
{
    VectorXf prop( DIM );

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
        cout << "Usage: ./ ... (int) moves_to_be_made (int) burn-in-steps (double) alpha (bool) show_vec" << endl;
        exit( 1 );
    }

    const int nsteps = atoi( argv[1] );
    const int burnin = atoi( argv[2] );
    const double alpha = atof( argv[3] );
    const bool show_vecs = atoi( argv[4] );

    setprecision( 3 );

    cerr << "-----------" << endl;
    cerr << "Input parameters: " << endl;
    cerr << ">> nsteps: " << nsteps << endl;
    cerr << ">> burnin: " << burnin << endl;
    cerr << ">> alpha: " << alpha << endl;
    cerr << ">> show_vecs: " << show_vecs << endl;
    cerr << "-----------" << endl;

    VectorXf x( DIM );
    x << 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;

	VectorXf xnew;

    cout << "# Metropolis-Hasitngs sample for CO2-Ar" << endl;
    cout << "# nsteps = " << nsteps << endl;
    cout << "# burnin = " << burnin << endl;
    cout << "# alpha = " << alpha << endl;
    cout << "# temperature = " << temperature << endl;
    cout << "# R (a. u.) = " << RDIST << endl;
	
	cout << "# file structure: " << endl;
	cout << "# theta pR pT jx jy jz" << endl;

	// burnin cycle
	for ( size_t i = 0; i < burnin; i++ )
	{
		x = metro_step( x, alpha );
	}
	cerr << "Burn-in cycle is finished." << endl;

	// running cycle
	size_t attempted_steps = 0;
	size_t moves = 0;

	int block_counter = 0;
	int block_size = 100000;

	chrono::milliseconds time_for_block;
	chrono::milliseconds time_for_blocks;

	vector<chrono::high_resolution_clock::time_point> block_times;
	block_times.push_back( chrono::high_resolution_clock::now() );

	while ( moves < nsteps )  
	{
		if ( attempted_steps % block_size == 0 && attempted_steps != 0 )
		{
			block_counter++;

			block_times.push_back( chrono::high_resolution_clock::now() );

			time_for_block = chrono::duration_cast<chrono::milliseconds>( block_times.end()[-1] - block_times.end()[-2] );
			time_for_blocks = chrono::duration_cast<chrono::milliseconds>( block_times.end()[-1] - block_times[0] );

			cerr << endl;
			cerr << "Block " << block_counter << " finished. " << endl;
			cerr << "Moves attempted: " << attempted_steps << "; moves made: " << moves << endl;
			cerr << "Time for current block: " << time_for_block.count() / 1000.0 << " s" << endl;
			cerr << "Total time elapsed: " << time_for_blocks.count() / 1000.0 << " s" << endl;
		}
			
		attempted_steps++;
        xnew = metro_step( x, alpha );
		
		// wrapping theta to [0, Pi)
		xnew(0) = wrapMax( xnew(0), M_PI );

        if ( xnew != x )
        {
            moves++;
        	x = xnew;
        } 
       
	   	if ( show_vecs == true )
        {
         	cout << x(0) << " " << x(1) << " " << x(2) << " " << x(3) << " " << x(4) << " " << x(5) << endl;
		}

    }

	chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();

    cerr << "-----------------------------------" << endl;
	cerr << "Burnin: " << burnin << endl;
    cerr << "Total steps: " << attempted_steps  << "; moves: " << moves << "; percent: " << (double) moves / attempted_steps * 100 << "%" << endl;
    cerr << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds>(endTime - block_times[0]).count() / 1000.0 << " s" << endl; 
    cerr << "-----------------------------------" << endl;

	return 0;
}

