#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

#include <iomanip> // std::atoi
#include <algorithm> // std::min
#include <numeric> // std::accumulate

#include <fstream>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_errno.h>

// Eigen
#include <Eigen/Dense>

// co2ar hamiltonian and potential
#include "co2ar_hamiltonian.hpp"

using namespace Eigen;
using namespace std;

const double MIN_DEVIATION = 0.1;

const int NBINS =  500;

const double EXTEND = 1.2;

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

const int DIM = 6;

// temperature in K
const double temperature = 300;

struct dimensionRange
{
	double lower_boundary;
	double upper_boundary;

	// mean value of lbs/ubs vectors
	double lbs_mean;
	double ubs_mean;

	// standard deviations of lbs/ubs vectors
	double lbs_stdev;
	double ubs_stdev;

	// relative standard deviations of lbs/ubs vectors
	double lbs_relative_stdev;
	double ubs_relative_stdev;

	// if relative deviation < MIN_DEVIATION then the boundaries of the current dimension thought to be found 
	bool convergence_status = false;

	vector<double> lbs;
	vector<double> ubs;
};

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

inline void check_boundaries( double x, double* lower_boundary, double* higher_boundary )
{
	if ( x < *lower_boundary )
	{
		*lower_boundary = x;
	}
	if ( x > *higher_boundary )
	{
		*higher_boundary = x;
	}
}

void cut_bounary_vectors( dimensionRange *d )
{
	vector<double> interim_ubs( d->ubs.end() - 5, d->ubs.end() );
	vector<double> interim_lbs( d->lbs.end() - 5, d->lbs.end() );

	d->ubs.clear();
	d->lbs.clear();

	for ( int i = 0; i < 5; i++ )
	{
		d->ubs.push_back( interim_ubs[i] );
		d->lbs.push_back( interim_lbs[i] );
	}
}

void run_statistics( dimensionRange *d )
{
	double lbs_sum = accumulate( d->lbs.begin(), d->lbs.end(), 0.0 );
	d->lbs_mean = lbs_sum / d->lbs.size();

	double ubs_sum = accumulate( d->ubs.begin(), d->ubs.end(), 0.0 );
	d->ubs_mean = ubs_sum / d->ubs.size();

	double accum = 0;
	for_each( begin( d->lbs ), end( d->lbs ), [&]( const double x ) {
		accum += ( x - d->lbs_mean ) * ( x - d->lbs_mean );
	});
	
	d->lbs_stdev = sqrt( accum / ( d->lbs.size() - 1 ));

	accum = 0;
	for_each( begin( d->ubs ), end ( d->ubs ), [&]( const double x ) {
		accum += ( x - d->ubs_mean ) * ( x - d->ubs_mean ); 
	});

	d->ubs_stdev = sqrt( accum / ( d->ubs.size() - 1 ));

	d->lbs_relative_stdev = abs( d->lbs_stdev / d->lbs_mean );
	d->ubs_relative_stdev = abs( d->ubs_stdev / d->ubs_mean );

	//cout << "Lower; mean: " << d->lbs_mean << "; stdev: " << d->lbs_stdev << "; relative stdev: " << d->lbs_relative_stdev << endl;
	//cout << "Upper; mean: " << d->ubs_mean << "; stdev: " << d->ubs_stdev << "; relative stdev: " << d->ubs_relative_stdev << endl;

	if ( d->lbs_relative_stdev < MIN_DEVIATION && d->ubs_relative_stdev < MIN_DEVIATION )
	{
		d->convergence_status = true;
	}
	else
	{
		d->convergence_status = false;
	}
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


int main( int argc, char* argv[] )
{
    if ( argc != 4 )
    {
        cout << "Usage: ./ ... (int) moves_to_be_made (double) alpha (bool) show_vec" << endl;
        exit( 1 );
    }

    const int nsteps = atoi( argv[1] );
    const double alpha = atof( argv[2] );
    const bool show_vecs = atoi( argv[3] );

    cerr << "-----------" << endl;
    cerr << "Input parameters: " << endl;
    cerr << ">> nsteps: " << nsteps << endl;
    cerr << ">> alpha: " << alpha << endl;
    cerr << ">> show_vecs: " << show_vecs << endl;
    cerr << "-----------" << endl;

    VectorXf x( DIM );
    x << 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;

	VectorXf xnew;

    cout << "# Metropolis-Hastings sample for CO2-Ar" << endl;
    cout << "# nsteps = " << nsteps << endl;
    cout << "# alpha = " << alpha << endl;
    cout << "# temperature = " << temperature << endl;
    cout << "# R (a. u.) = " << RDIST << endl;
	cout << "# file structure: " << endl;
	cout << "# theta pR pT jx jy jz" << endl;
	
	int burnin_step = 0;
	
	int block_counter = 0;
	int block_size = 1000000;

	chrono::milliseconds time_for_block;
	chrono::milliseconds time_for_blocks;

	vector<chrono::high_resolution_clock::time_point> block_times;
	block_times.push_back( chrono::high_resolution_clock::now() );

	// creating a vector of structures to describe phase space boundaries
	vector< dimensionRange > ranges;
	for ( int i = 0; i < DIM; i++ )
	{
		dimensionRange d;
		d.lower_boundary = d.upper_boundary = x(i);

		ranges.push_back( d ); 
	}


	// burnin cycle
	bool total_status = false;

	while ( !total_status )
	{
		if ( burnin_step % block_size == 0 && burnin_step != 0 )
		{
			block_counter++;
			
			block_times.push_back( chrono::high_resolution_clock::now() );
			
			time_for_block = chrono::duration_cast<chrono::milliseconds>( block_times.end()[-1] - block_times.end()[-2] );
			time_for_blocks = chrono::duration_cast<chrono::milliseconds>( block_times.end()[-1] - block_times[0] );
			
			cerr << endl;
			cerr << "---------------------------------------" << endl;
			cerr << "Burn-in block " << block_counter << " finished. " << endl;
			cerr << "Time for current burn-in block: " << time_for_block.count() / 1000.0 << " s" << endl;
			cerr << "Total time elapsed in burn-in: " << time_for_blocks.count() / 1000.0 << " s" << endl;
		
			// putting current lower/upper boundary in the corresponding vector 
			for ( int i = 0; i < DIM; i++ )
			{
				ranges[i].lbs.push_back( ranges[i].lower_boundary );
				ranges[i].ubs.push_back( ranges[i].upper_boundary );
			}

			if ( block_counter > 10 )
			{
				// if some dimension did not converge, then total status will be set to false inside this cycle
				total_status = true;				
				for ( int i = 1; i < DIM; i++ )
				{
					cut_bounary_vectors( &ranges[i] );
					run_statistics( &ranges[i] );
					
					cout << "Convergence_status of  " << i << ": " << ranges[i].convergence_status << endl;
					
					if ( ranges[i].convergence_status == false )
					{
						total_status = false;
					}
				}
			}
			else
			{	
				for ( int i = 0; i < DIM; i++ )
				{
					run_statistics( &ranges[i] );
				}
			}	
			

			cerr << "---------------------------------------" << endl;
		}	
		
		// step of Metropolis algorithm
		x = metro_step( x, alpha );
		
		// wrapping theta to [0, Pi)
		x(0) = wrapMax( x(0), M_PI );
		
		// checking for possible boundary update
		for ( int i = 0; i < DIM; i++ )
		{
			check_boundaries( x(i), &ranges[i].lower_boundary, &ranges[i].upper_boundary );
		}

		burnin_step++;
	}

	cerr << "-----------------------------------------------------" << endl;
	cerr << "------------ Burn-in cycle is finished.--------------" << endl;

	cerr << "Final phase space boundaries: " << endl;
	for ( int i = 0; i < DIM; i++ )
	{
		cout << "i: " << i << "; " << ranges[i].lower_boundary << " --- " << ranges[i].upper_boundary << endl; 
	}

	// resetting auxiliary variables
	block_counter = 0;
	block_times.clear();
	block_times.push_back( chrono::high_resolution_clock::now() );

	// histograms
	vector< gsl_histogram* > histograms;
	for ( int i = 0; i < DIM; i++ )
	{
		gsl_histogram *h = gsl_histogram_alloc( NBINS );

		if ( i == 0 )
		{
			gsl_histogram_set_ranges_uniform( h, 0, M_PI );
		}
		else
		{
			gsl_histogram_set_ranges_uniform( h, ranges[i].lower_boundary* EXTEND, ranges[i].upper_boundary * EXTEND );
		}

		histograms.push_back( h );
	}

	// running cycle
	size_t attempted_steps = 0;
	size_t moves = 0;

	// status of histogram incrementing
	int status;

	while ( moves < nsteps )  
	{
		if ( attempted_steps % block_size == 0 && attempted_steps != 0 )
		{
			block_counter++;

			block_times.push_back( chrono::high_resolution_clock::now() );

			time_for_block = chrono::duration_cast<chrono::milliseconds>( block_times.end()[-1] - block_times.end()[-2] );
			time_for_blocks = chrono::duration_cast<chrono::milliseconds>( block_times.end()[-1] - block_times[0] );

			cerr << endl;
			cerr << "---------------------------------------" << endl;
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

		for ( int i = 0; i < DIM; i++ )
		{
			status = gsl_histogram_increment( histograms[i], x(i) );
			if ( status == GSL_EDOM )
			{
				cerr << "ERROR BINNING: i = " << i << "; x(i) = " << x(i) << endl;
			}
		}
    }

	cerr << "-----------------------------------------------" << endl;
	cerr << "---------------MCMC finished ------------------" << endl;
	
	save_histogram( histograms[0], "theta.txt" );
	save_histogram( histograms[1], "pr.txt" );
	save_histogram( histograms[2], "ptheta.txt" );
	save_histogram( histograms[3], "jx.txt" );
	save_histogram( histograms[4], "jy.txt" );
	save_histogram( histograms[5], "jz.txt" );

	chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();

    cerr << "-----------------------------------" << endl;
    cerr << "Total steps: " << attempted_steps  << "; moves: " << moves << "; percent: " << (double) moves / attempted_steps * 100 << "%" << endl;
    cerr << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds>(endTime - block_times[0]).count() / 1000.0 << " s" << endl; 
    cerr << "-----------------------------------" << endl;

	return 0;
}

