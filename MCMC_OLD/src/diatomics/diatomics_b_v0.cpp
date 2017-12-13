#include <iostream>
#include <random>
#include <chrono>
#include <cmath>

#include <fstream>
#include <string>

#include <gsl/gsl_histogram.h>

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

const double he_mass = 4.0;
const double ar_mass = 40.0;

const double MUKG = ( 4.0 * 40.0 ) / ( 40.0 + 4.0 ) * DALTON;
const long double MUAMU = MUKG / AMU; 

// temperature in K
const double temperature = 300;

const int NBINS = 500;

static mt19937 generator;

static double nextDouble( const double &min = 0.0, const double &max = 1.0 )
{
    uniform_real_distribution<double> distribution( min, max );
    return distribution( generator );
}

static void nextGaussianVec( Vector2d &v, Vector2d mean, const double sigma )
{
    normal_distribution<double> d0( mean(0), sigma );
    normal_distribution<double> d1( mean(1), sigma );

    v(0) = d0( generator );
    v(1) = d1( generator );
}

void gsl_histogram_normalize( gsl_histogram* h )
{
	double max = gsl_histogram_max( h );
	double min = gsl_histogram_min( h );
	double step = (max - min) / NBINS;

	double sum = gsl_histogram_sum( h ) * step;
	gsl_histogram_scale( h, 1.0 / sum );
}

void save_histogram( gsl_histogram *histogram, string filename )
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

// x = [ v0, phi ]
double target( Vector2d x )
{
    double v0 = x(0);
    double phi = x(1);
    
	double h = MUAMU * pow(v0, 2) / 2.0 * ( 1.0 + sin(phi) * sin(phi) ); 
    return exp( - h * HTOJ / ( BOLTZCONST * temperature ));
}

// wrap x -> [0, max)
double wrapMax( double x, double max )
{
	return fmod( max + fmod(x, max), max );
}

Vector2d metro_step( Vector2d x, double alpha )
{
    Vector2d prop;

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
    if ( argc != 4 )
    {
        cout << "USAGE: ./... (int) moves-to-be-made (int) burn-in-steps (double) alpha" << endl;
        exit( 1 );
    }

    const int nsteps = atoi( argv[1] );
    const int burnin = atoi( argv[2] );
    const double alpha = atof( argv[3] ); 

    setprecision(3);

    cerr << "-----" << endl;
    cerr << "Input parameters: " << endl;
    cerr << ">> moves-to-be-made: " << nsteps << endl;
    cerr << ">> burn-in: " << burnin << endl;
    cerr << ">> alpha: " << alpha << endl;

    Vector2d x ( -0.0005, 0.0 );
    Vector2d xnew;

	chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();	

	// burn-in
	for ( size_t i = 0; i < burnin; i++ )
	{
		x = metro_step( x, alpha );
	}
	
	cerr << "Burn-in finished. Time elapsed: " << chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - startTime).count() / 1000.0 << " s" << endl;

	int moves = 0;
	int attempted_steps = 0;

	gsl_histogram* v0_hist = gsl_histogram_alloc( NBINS );
	gsl_histogram* phi_hist = gsl_histogram_alloc( NBINS );
	gsl_histogram* sin_phi_hist = gsl_histogram_alloc( NBINS );

	gsl_histogram_set_ranges_uniform( v0_hist, -5e-3, 5e-3 );
	gsl_histogram_set_ranges_uniform( phi_hist, 0, 2 * M_PI ); 
	gsl_histogram_set_ranges_uniform( sin_phi_hist, -1.0, 1.0 );

	int block_counter = 0;
	int block_size = 1000000;
	
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

        xnew = metro_step( x, alpha );

        if ( xnew != x )
        {
            moves++;
        }
    
		xnew(1) = wrapMax( xnew(1), 2 * M_PI ); 

		gsl_histogram_increment( v0_hist, xnew(0) );
		gsl_histogram_increment( phi_hist, xnew(1) );
		gsl_histogram_increment( sin_phi_hist, sin( xnew(1) ));

        x = xnew;
		attempted_steps++;
	}

	gsl_histogram_normalize( v0_hist );
	gsl_histogram_normalize( phi_hist );
	gsl_histogram_normalize( sin_phi_hist );

	save_histogram( v0_hist, "v0.txt" );
	save_histogram( phi_hist, "phi.txt" );
	save_histogram( sin_phi_hist, "sin_phi.txt" );

	gsl_histogram_free( v0_hist );
	gsl_histogram_free( phi_hist );
	gsl_histogram_free( sin_phi_hist );

    cerr << "-----------------------------------" << endl;
    cerr << "Attempted steps: " << attempted_steps << "; moves: " << moves << "; percent: " << (double) moves / attempted_steps * 100.0 << "%" << endl;
    cerr << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - startTime).count() / 1000.0 << " s" << endl; 
    cerr << "-----------------------------------" << endl;

    return 0;
}

