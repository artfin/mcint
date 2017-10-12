#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <fstream>

#include <iomanip> // std::atoi
#include <algorithm> // std::min
#include <gsl/gsl_histogram.h>

// Eigen
#include <Eigen/Dense>

// co2ar hamiltonian and potential
#include "co2ar_hamiltonian.hpp"

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
const double RDIST = 20.0;
// -----------------------------

const int NBINS = 500;
const int DIM = 6;

// boundaries for L2
const double LBOUND = 899.5;
const double UBOUND = 900.5;

// temperature in K
const double temperature = 300;

static mt19937 generator;

double cotan( double x )
{
	return 1.0 / tan( x );
}

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

void gsl_histogram_normalize( gsl_histogram* h )
{
	double max = gsl_histogram_max( h );
   	double min = gsl_histogram_min( h );
	double step = ( max - min ) / NBINS;	

	double sum = gsl_histogram_sum( h ) * step; 
	gsl_histogram_scale( h, 1.0 / sum );
}

void save_histogram( gsl_histogram *histogram, string filename, int nsteps, int burnin, double alpha )
{
	ofstream file( filename );

	file << "# Distributions of Hamiltonian variables for bound CO2-Ar dimers" << endl;
    file << "# Metropolis-Hastings algorithm, initial parameters:" << endl;
    file << "# nsteps = " << nsteps << endl;
    file << "# burnin = " << burnin << endl;
    file << "# alpha = " << alpha << endl;
    file << "# temperature = " << temperature << endl;
	

	double lower_bound, higher_bound, bin_content;
	for ( int counter = 0; counter < NBINS; counter++ )
	{
		gsl_histogram_get_range( histogram, counter, &lower_bound, &higher_bound );
		bin_content = gsl_histogram_get( histogram, counter );

		file << lower_bound << " " << higher_bound << " " << bin_content << endl;
	}
	
	file.close();	
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

// x = [ Theta, pR, pT, Jx, Jy, Jz ]
double linear_molecule_momentum( VectorXf x )
{
	double Theta = x( 0 );
	double Jz = x( 5 );
	double pTheta = x( 2 );

	return pow( Jz, 2 ) * ( 1 + pow( cotan(Theta), 2 )) + pow( pTheta, 2 );
}

int main( int argc, char* argv[] )
{
    if ( argc != 5 )
    {
        cout << "Usage: ./ ... (int) vectors_to_write (int) burn-in-steps (double) alpha (bool) show_vec" << endl;
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
	
	double lbs[] = {0.0, -30.0, -30.0, -150.0, -150.0, -150.0, 19.0 };
	vector<double> lower_boundaries ( lbs, lbs + sizeof(lbs) / sizeof(double) );
	double ubs[] =  {M_PI, 30.0, 30.0, 150.0, 150.0, 150.0, 21.0 }; 
	vector<double> upper_boundaries( ubs, ubs + sizeof(ubs) / sizeof(double) );

	vector< gsl_histogram* > histograms;
	for ( int i = 0; i < DIM + 1; i++ )
	{
		gsl_histogram *h = gsl_histogram_alloc( NBINS );
		gsl_histogram_set_ranges_uniform( h, lower_boundaries[i], upper_boundaries[i] );
	  	histograms.push_back( h );	
	}	

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

	double L2;

	chrono::milliseconds time_for_block;
	chrono::milliseconds time_for_blocks;

	vector<chrono::high_resolution_clock::time_point> block_times;
	block_times.push_back( chrono::high_resolution_clock::now() );

	while ( wrote_vectors < nsteps )  
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

		L2 = linear_molecule_momentum( x );

        if ( xnew != x )
        {
            moves++;
        } 
       
	   	if ( show_vecs == true && xnew != x &&
				L2 > LBOUND && L2 < UBOUND )
        {
			jx = xnew(3);
			jy = xnew(4);
			jz = xnew(5);

			j = sqrt( pow(jx, 2) + pow(jy, 2) + pow(jz, 2) );
			jtheta = acos( jz / j );
			jphi = atan ( jy / jx );

			cout << wrote_vectors + 1 << " " << RDIST << " " << xnew(0) << " " << xnew(1) << " " << xnew(2) << " " << jphi << " " << jtheta << " " << j << endl;

			//cout << x(0) << " " << x(1) << " " << x(2) << " " << x(3) << " " << x(4) << " " << x(5) << endl;
		
			//for ( int i = 0; i < DIM; i++ )
			//{
				//gsl_histogram_increment( histograms[i], x(i) ); 
			//}

			//gsl_histogram_increment( histograms[DIM], L2 );
			
			wrote_vectors++;
		}
        
		x = xnew;
    }
	
	string names[] = {"theta.txt", "pr.txt", "pt.txt", "jx.txt", "jy.txt", "jz.txt", "l2.txt" };
		
	for ( int i = 0; i < DIM + 1; i++ )
	{
		//gsl_histogram_normalize( histograms[i] );
		//save_histogram( histograms[i], names[i], nsteps, burnin, alpha );
		gsl_histogram_free( histograms[i] );
	}

	chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();

    cerr << "-----------------------------------" << endl;
	cerr << "Burnin: " << burnin << endl;
    cerr << "Total steps: " << attempted_steps  << "; moves: " << moves << "; percent: " << (double) moves / attempted_steps * 100 << "%" << endl;
    cerr << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds>(endTime - block_times[0]).count() / 1000.0 << " s" << endl; 
    cerr << "-----------------------------------" << endl;

	return 0;
}

