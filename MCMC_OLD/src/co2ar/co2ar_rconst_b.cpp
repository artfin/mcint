#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <fstream>

#include <iomanip> 
#include <algorithm> // std::min

// Eigen
#include <Eigen/Dense>

// co2ar hamiltonian and potential
#include "co2ar_hamiltonian.hpp"

#include <gsl/gsl_histogram.h>

using namespace Eigen;
using namespace std;

const double mu1 = 14579.0; // ?
const double mu2 = 36440.0; // ?
const double l = 4.398;

// boltzmann constant
const double BOLTZCONST = 1.38064e-23;
// hartree to joules
const double HTOJ = 4.35974417e-18;

// ! ---------------------------
const double RDIST = 100.0;
// -----------------------------

const int NBINS = 100;

const double BBOUND = 100.0;

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

double cot( double x )
{
	return 1.0 / tan( x );
}

// wrap x -> [0, max)
double wrapMax( double x, double max )
{
	return fmod( max + fmod(x, max), max );
}

// x = [ Theta, pR, pT, Jx, Jy, Jz ]
double target( VectorXf x )
{
    double ke = kinetic_energy( RDIST, x[0], x[1], x[2], x[3], x[4], x[5] );
	return exp( -ke * HTOJ / (BOLTZCONST * temperature ));
}

void gsl_histogram_normalize( gsl_histogram* h )
{
	double max = gsl_histogram_max( h );
	double min = gsl_histogram_min( h );
	double step = (max - min) / NBINS;

	double sum = gsl_histogram_sum( h ) * step;
	gsl_histogram_scale( h, 1.0 / sum );
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

void save_histogram( gsl_histogram *histogram, string filename, int nsteps, int burnin, double alpha )
{
	ofstream file( filename );

	double lower_bound, higher_bound, bin_content;
	for ( int counter = 0; counter < NBINS; counter++ )
	{
		gsl_histogram_get_range( histogram, counter, &lower_bound, &higher_bound );
		bin_content = gsl_histogram_get( histogram, counter );

		file << lower_bound << " " << higher_bound << " " << bin_content << endl;
	}
	
	file.close();	
}


// x = [ Theta, pR, pT, Jx, Jy, Jz ]
double calc_gunsight( VectorXf x, gsl_histogram* quad1_hist, gsl_histogram* quad2_hist, gsl_histogram* quad3_hist, gsl_histogram* cos_phi_histogram, gsl_histogram* denumerator_histogram )
{
	double Theta = x( 0 );
	double pR = x( 1 );
	double pT = x( 2 );
	double Jx = x( 3 );
	double Jy = x( 4 );
	double Jz = x( 5 );

	if ( pR > 0 )
	{
		return -1;
	}

	double pRmu2 = pR / mu2;

	double quad1 = pow( pRmu2, 2 );
	gsl_histogram_increment( quad1_hist, quad1 );
	
	double quad2 = pow( (Jx + Jz * cot(Theta))/(mu2 * RDIST), 2 );
	gsl_histogram_increment( quad2_hist, quad2 );

	double quad3 = pow( (Jy - pT)/(mu2 * RDIST), 2 );
	gsl_histogram_increment( quad3_hist, quad3 );

	double numerator = - pRmu2; 
	
	double denumerator = sqrt( quad1  + quad2 + quad3 );
	gsl_histogram_increment( denumerator_histogram, denumerator );

	double cos_phi = numerator / denumerator;
	gsl_histogram_increment( cos_phi_histogram, cos_phi );

	double gunsight = RDIST * sqrt(1 - pow(cos_phi, 2));
	
	// for R -> \infinity
	//double gunsight = 1.0 / fabs(pR) * sqrt( quad2 + quad3 ); 
	
	return gunsight;
}

int main( int argc, char* argv[] )
{
    if ( argc != 4 )
    {
        cout << "Usage: ./ ... (int) vectors_to_write (int) burn-in-steps (double) alpha" << endl;
        exit( 1 );
    }

    const int nsteps = atoi( argv[1] );
    const int burnin = atoi( argv[2] );
    const double alpha = atof( argv[3] );
    
    setprecision( 3 );

    cerr << "-----------" << endl;
    cerr << "Input parameters: " << endl;
    cerr << ">> nsteps: " << nsteps << endl;
    cerr << ">> burnin: " << burnin << endl;
    cerr << ">> alpha: " << alpha << endl;
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
	cout << "# theta pR pT alpha beta j" << endl;
	
	double jx, jy, jz;
	double j, jtheta, jphi;

	// burnin cycle
	for ( size_t i = 0; i < burnin; i++ )
	{
		x = metro_step( x, alpha );
	}
	cerr << "Burn-in cycle is finished." << endl;
	
	// running cycle
	size_t attempted_steps = 0;
	size_t moves = 0;
	size_t wrote_vectors = 0;

	int block_counter = 0;
	int block_size = 1000000;

	double b;
	
	gsl_histogram* histogram_b = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( histogram_b, 0.0, BBOUND );	

	gsl_histogram* quad1_hist = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( quad1_hist, 0.0, 2e-7 );
	
	gsl_histogram* quad2_hist = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( quad2_hist, 0.0, 2e-7 );
	
	gsl_histogram* quad3_hist = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( quad3_hist, 0.0, 2e-7 );

	gsl_histogram* cos_phi_histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( cos_phi_histogram, -1.0, 1.0 );

	gsl_histogram* pr_histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( pr_histogram, -50.0, 50.0 );

	gsl_histogram* denumerator_histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( denumerator_histogram, 0.0, 1e-3 );

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
			cerr << "Vectors wrote: " << wrote_vectors << endl;
			cerr << "Time for current block: " << time_for_block.count() / 1000.0 << " s" << endl;
			cerr << "Total time elapsed: " << time_for_blocks.count() / 1000.0 << " s" << endl;
		}
			
		attempted_steps++;
        xnew = metro_step( x, alpha );
		
		// wrapping theta to [0, Pi)
		xnew(0) = wrapMax( xnew(0), M_PI );

		b = calc_gunsight( x, quad1_hist, quad2_hist, quad3_hist, cos_phi_histogram, denumerator_histogram );
		gsl_histogram_increment( histogram_b, b );
		
		gsl_histogram_increment( pr_histogram, xnew(1) );

        if ( xnew != x )
        {
            moves++;
        } 
       
	   	//if ( xnew != x && b > 0 && b < BBOUND )
        //{
			//jx = xnew(3);
			//jy = xnew(4);
			//jz = xnew(5);

			//j = sqrt( pow(jx, 2) + pow(jy, 2) + pow(jz, 2) );
			//jtheta = acos( jz / j );
			//jphi = atan ( jy / jx );

			//cout << wrote_vectors + 1 << " " << RDIST << " " << xnew(0) << " " << xnew(1) << " " << xnew(2) << " " << jphi << " " << jtheta << " " << j << endl;
			//cout << jx << " " << jy << " " << jz << endl;

			//wrote_vectors++;
		//}
        
		x = xnew;
    }
	
	gsl_histogram_normalize( histogram_b );
	save_histogram( histogram_b, "b_true.txt", nsteps, burnin, alpha );
	gsl_histogram_free( histogram_b );

	gsl_histogram_normalize( cos_phi_histogram );
	save_histogram( cos_phi_histogram, "cos_phi_true.txt", nsteps, burnin, alpha );
	gsl_histogram_free( cos_phi_histogram );

	gsl_histogram_normalize( quad1_hist );
	save_histogram( quad1_hist, "quad1_true.txt", nsteps, burnin, alpha );
	gsl_histogram_free( quad1_hist );
	
	gsl_histogram_normalize( quad2_hist );
	save_histogram( quad2_hist, "quad2_true.txt", nsteps, burnin, alpha );
	gsl_histogram_free( quad2_hist );
	
	gsl_histogram_normalize( quad3_hist );
	save_histogram( quad3_hist, "quad3_true.txt", nsteps, burnin, alpha );
	gsl_histogram_free( quad3_hist );

	gsl_histogram_normalize( pr_histogram ); 
	save_histogram( pr_histogram, "pr_true.txt", nsteps, burnin, alpha );
	gsl_histogram_free( pr_histogram );

	gsl_histogram_normalize( denumerator_histogram );
	save_histogram( denumerator_histogram, "denumerator_true.txt", nsteps, burnin, alpha );

	chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();

    cerr << "-----------------------------------" << endl;
	cerr << "Burnin: " << burnin << endl;
    cerr << "Total steps: " << attempted_steps  << "; moves: " << moves << "; percent: " << (double) moves / attempted_steps * 100 << "%" << endl;

	// * 1000 not on 100 because every 10th element of chain is taken
    cerr << "Wrote vectors: " << wrote_vectors << "; percent (to attempted): " << (double) wrote_vectors / attempted_steps * 100.0 << "%" << "; percent( to moves ): " << (double) wrote_vectors / moves * 100.0 << "%" << endl;
	cerr << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds>(endTime - block_times[0]).count() / 1000.0 << " s" << endl; 
    cerr << "-----------------------------------" << endl;

	return 0;
}

