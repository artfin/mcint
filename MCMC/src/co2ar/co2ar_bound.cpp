#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

#include <iomanip> // std::atoi
#include <algorithm> // std::min

#include <Eigen/Dense>

// co2ar hamiltonian and potential
#include "co2ar_hamiltonian.h"
#include "co2ar_potential_julia.h"

using namespace Eigen;
using namespace std;

const double mu1 = 14579.0; // ?
const double mu2 = 36440.0; // ?
const double l = 4.398;

// boltzmann constant
const double BOLTZCONST = 1.38064e-23;
// hartree to joules
const double HTOJ = 4.35974417e-18;
// cm^-1 to hartree
const double CMTOH = 4.55633e-6;

// temperature in K
const double temperature = 300;

const int DIM = 7;

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

// x = [ R, theta, pR, pT, Jx, Jy, Jz ]
double target( VectorXf x )
{
    double ke = kinetic_energy( x[0], x[1], x[2], x[3], x[4], x[5], x[6] );
	double p = ab_initio_pot( x[0], x[1] ) * CMTOH;
	double h = ke + p; // hamiltonian

	return exp( -h * HTOJ / (BOLTZCONST * temperature ));
}

VectorXf metro_step( VectorXf x, double alpha )
{
    VectorXf prop( DIM );

    // generate random vector
    nextGaussianVec( prop, x, alpha );

	double prop_kin_energy = kinetic_energy( prop[0], prop[1], prop[2], prop[3], prop[4], prop[5], prop[6] );
	double prop_pot_energy = ab_initio_pot( prop[0], prop[1] ) * CMTOH;
	double prop_energy = prop_kin_energy + prop_pot_energy;

	if ( prop_energy > 0 )
	{
		return x;
	}

	double target_prop = exp( - prop_energy * HTOJ / ( BOLTZCONST * temperature ));
	double target_x = target(x);

    if ( nextDouble() < min( 1.0, target_prop / target_x ))
    {
        return prop;
    }
	else
	{
		return x;
	}
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
    x << 6.0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;
	VectorXf xnew( DIM );

	cout << "# Distributions of Hamiltonian variables for bound CO2-Ar dimers" << endl;
    cout << "# Metropolis-Hasitngs algorithm, initial parameters:" << endl;
    cout << "# nsteps = " << nsteps << endl;
    cout << "# burnin = " << burnin << endl;
    cout << "# alpha = " << alpha << endl;
    cout << "# temperature = " << temperature << endl;
	
	cout << "# file structure: R, theta, pR, pT, jx, jy, jz" << endl;

	// burnin cycle
	for ( size_t i = 0; i < burnin; i++ )
	{
		x = metro_step( x, alpha );
	}

	cerr << "Burn-in is finished" << endl;

	// running cycle
	size_t attempted_steps = 0;
	size_t moves = 0;

	int block_counter = 0;
	int block_size = 100000;

	chrono::milliseconds time_for_block;
	chrono::milliseconds time_for_blocks;

	vector< chrono::high_resolution_clock::time_point> block_times;
	
	block_times.push_back( chrono::high_resolution_clock::now() );

	while ( moves < nsteps )  
	{
		if ( attempted_steps % block_size == 0 and attempted_steps != 0 )
		{
			block_counter++;
		
			block_times.push_back( chrono::high_resolution_clock::now() );

			time_for_block = chrono::duration_cast<chrono::milliseconds>( block_times.end()[-1] - block_times.end()[-2] );
			time_for_blocks = chrono::duration_cast<chrono::milliseconds>( block_times.end()[-1] - block_times[0] );			

			cerr << endl;
			cerr << "Block " << block_counter << " finished." << endl;
			cerr << "Moves attempted: " << attempted_steps << "; moves made: " << moves << endl;
			cerr << "Time for current block: " << time_for_block.count() / 1000.0 << " s" << endl;
		   	cerr << "Total time elapsed: " << time_for_blocks.count() / 1000.0 << " s" << endl << endl;
		}

        xnew = metro_step( x, alpha );

		// wrapping theta to [0, Pi)
		xnew(1) = wrapMax( xnew(1), M_PI );

        if ( xnew != x )
        {
            moves++;
        	x = xnew;
        } 
        
		if ( show_vecs == true )
        {
         	cout << x(0) << " " << x(1) << " " << x(2) << " " << x(3) << " " << x(4) << " " << x(5) << " " << x(6) << endl;
		}

		attempted_steps++;
    }

	chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();

    cerr << "-----------------------------------" << endl;
	cerr << "Burnin: " << burnin << endl;
    cerr << "Total steps: " << attempted_steps  << "; moves: " << moves << "; percent: " << (double) moves / attempted_steps * 100 << "%" << endl;
    cerr << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds>( endTime - block_times[0] ).count() / 1000.0 << " s" << endl;
    cerr << "-----------------------------------" << endl;

	return 0;
}

