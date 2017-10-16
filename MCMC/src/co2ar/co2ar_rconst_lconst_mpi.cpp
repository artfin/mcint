#include <mpi.h>

#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <chrono>
#include <gsl/gsl_histogram.h>

#include <Eigen/Dense>

// co2ar hamiltonian
#include "co2ar_hamiltonian.hpp"

const int WORK_TAG = 0;
const int FINAL_TAG = 42;

// temperature in K
const double temperature = 300;

const int NBINS = 200;
// boundaries for L2
const double LBOUND = 19.5;
const double UBOUND = 20.5;

const int DIM = 6;

// boltzmann constant
const double BOLTZCONST = 1.38064e-23;
// hartree to joules
const double HTOJ = 4.35974417e-18;

// ! ---------------------------
const double RDIST = 20.0;
// -----------------------------

using namespace Eigen;
using namespace std;

double nextDouble( const double &min = 0.0, const double &max = 1.0, int rank = 0 )
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	uniform_real_distribution<double> dist( min, max );
	
	mt19937 rng( seed * rank );
	
	return dist( rng );	
}

void nextGaussianVec( VectorXf &v, VectorXf mean, const double sigma, int rank )
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	mt19937 rng( seed * rank );

	for ( int i = 0; i < DIM; i++ )
	{
		normal_distribution<double> d( mean(i), sigma );
		v( i ) = d( rng );
	}
}

// x = [ Theta, pR, pT, Jx, Jy, Jz ]
double target( VectorXf x )
{
    double ke = kinetic_energy( RDIST, x[0], x[1], x[2], x[3], x[4], x[5] );
    
	return exp( -ke * HTOJ / (BOLTZCONST * temperature ));
}

double cotan( double x )
{
	return 1.0 / tan( x );
}

// x = [ Theta, pR, pT, Jx, Jy, Jz ]
double linear_molecule_momentum( VectorXf x )
{
	double Theta = x( 0 );
	double Jz = x( 5 );
	double pTheta = x( 2 );

	double L2 = pow( Jz, 2 ) * ( 1 + pow( cotan(Theta), 2 )) + pow( pTheta, 2 );
	return sqrt( L2 );
}


VectorXf metro_step( VectorXf x, double alpha, int rank )
{
	VectorXf prop( DIM );

	nextGaussianVec( prop, x, alpha, rank );
	if ( nextDouble( 0.0, 1.0, rank ) < min ( 1.0, target(prop) / target(x) ))
	{
		 x = prop;
	}
	
	return x;
}

// wrap x -> [0, max)
double wrapMax( double x, double max )
{
	return fmod( max + fmod(x, max), max );
}

vector< gsl_histogram* > allocate_histograms( void )
{
	double lbs[] = {0.0, -30.0, -30.0, -150.0, -150.0, -150.0 };
	vector<double> lower_boundaries ( lbs, lbs + sizeof(lbs) / sizeof(double) );
	double ubs[] =  {M_PI, 30.0, 30.0, 150.0, 150.0, 150.0 }; 
	vector<double> upper_boundaries( ubs, ubs + sizeof(ubs) / sizeof(double) );

	vector< gsl_histogram* > histograms;

	for ( int i = 0; i < DIM; i++ )
	{
		gsl_histogram *h = gsl_histogram_alloc( NBINS );
		gsl_histogram_set_ranges_uniform( h, lower_boundaries[i], upper_boundaries[i] );
	  	histograms.push_back( h );	
	}	
	
	return histograms;	
}
void send_histograms_to_root( vector< gsl_histogram* > histograms, int tag )
{
	vector<double> contents( NBINS );
	
	for ( int histogram_counter = 0; histogram_counter < DIM; histogram_counter++ )
	{
		contents.clear();

		for ( int bin_counter = 0; bin_counter < NBINS; bin_counter++ )
		{
			contents.push_back( gsl_histogram_get( histograms[histogram_counter], bin_counter ));
		}

		MPI_Ssend( &contents[0], NBINS, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD );
		// setting all elements of histogram to 0
		gsl_histogram_reset( histograms[histogram_counter] );
	}
}

void normalize( gsl_histogram* h )
{
	double max = gsl_histogram_max( h );
   	double min = gsl_histogram_min( h );
	double step = ( max - min ) / NBINS;	

	double sum = gsl_histogram_sum( h ) * step; 
	gsl_histogram_scale( h, 1.0 / sum );
}

void save_histogram( gsl_histogram *histogram, string filename )
{
	ofstream file( filename );

	// copying histogram
	gsl_histogram* cp = gsl_histogram_clone( histogram );
	
	// normalizing copy
	normalize( cp );

	double lower_bound, higher_bound, bin_content;
	for ( int counter = 0; counter < NBINS; counter++ )
	{
		gsl_histogram_get_range( cp, counter, &lower_bound, &higher_bound );
		bin_content = gsl_histogram_get( cp, counter );

		file << lower_bound << " " << higher_bound << " " << bin_content << endl;
	}
	
	file.close();	
}

void save_histograms( vector< gsl_histogram* > histograms )
{
	string names[] = { "theta.txt", "pr.txt", "ptheta.txt", "jx.txt", "jy.txt", "jz.txt" };

	for ( int histogram_counter = 0; histogram_counter < DIM; histogram_counter++ )
	{
		save_histogram( histograms[histogram_counter], names[histogram_counter] );
	}
}

void slave_code( int rank, int vectors_to_select, int burnin, double alpha )
{
	VectorXf x( DIM );
	x << 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;

	VectorXf xnew( DIM );

	double seed = ( rank + M_PI ) / pow( M_PI, 0.23834) + 4344.23;

	for ( size_t i = 0; i < burnin; i++ )
	{
		x = metro_step( x, alpha, seed );
	}

	cout << "Burn-in cycle is finished." << endl;
	
	vector< gsl_histogram* > histograms = allocate_histograms();

	size_t attempted_steps = 0;
	size_t vectors_selected = 0;
	size_t moves = 0;

	int block_counter = 0;
	int block_size = 1000000;

	double L;

	chrono::milliseconds time_for_block;
	chrono::milliseconds time_for_blocks;

	vector<chrono::high_resolution_clock::time_point> block_times;
	block_times.push_back( chrono::high_resolution_clock::now() );

	while ( vectors_selected < vectors_to_select )
	{
		if ( attempted_steps % block_size == 0 && attempted_steps != 0 )
		{
			block_counter++;
			
			block_times.push_back( chrono::high_resolution_clock::now() );

			time_for_block = chrono::duration_cast<chrono::milliseconds>( block_times.end()[-1] - block_times.end()[-2] );
			time_for_blocks = chrono::duration_cast<chrono::milliseconds>( block_times.end()[-1] - block_times[0] );

			cout << endl;
			cout << "Rank: " << rank << "; block " << block_counter << " finished. " << endl;
			cout << "Moves attempted: " << attempted_steps << "; moves made: " << moves << "; percentage: " << (double) moves / attempted_steps * 100.0 << "%" << endl;
			cout << "Vectors selected: " << vectors_selected << "; percentage (from steps attempted): " << (double) vectors_selected / attempted_steps * 100.0 << "%; percentage (of work done): " << (double) vectors_selected / vectors_to_select * 100.0 << "%" << endl;
			cout << "Time for current block: " << time_for_block.count() / 1000.0 << " s" << endl;
			cout << "Total time elapsed: " << time_for_blocks.count() / 1000.0 << " s" << endl;

			send_histograms_to_root( histograms, WORK_TAG );
		}

		attempted_steps++;
		xnew = metro_step( x, alpha, seed );

		// wrapping theta to [0, Pi)
		xnew(0) = wrapMax( xnew(0), M_PI );

		L = linear_molecule_momentum( x );

        if ( xnew != x )
        {
            moves++;
        	x = xnew;
        } 

		if ( L > LBOUND && L < UBOUND )
		{
			for ( int i = 0; i < DIM; i++ )
			{
				gsl_histogram_increment( histograms[i], x(i) );
			}

			vectors_selected++;
		}
	}
	
	chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();

    cout << "-----------------------------------" << endl;
    cout << "Total steps: " << attempted_steps  << "; moves: " << moves << "; percent: " << (double) moves / attempted_steps * 100 << "%" << endl;
    cout << "Vectors to select: " << vectors_to_select << "; percentage (from attempted steps): " << (double) vectors_to_select / attempted_steps * 100.0 << "%" << endl;
	cout << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds>(endTime - block_times[0]).count() / 1000.0 << " s" << endl; 
	cout << "Process " << rank << " sends final histograms." << endl;
    cout << "-----------------------------------" << endl;
	
	// final sending of histograms
	send_histograms_to_root( histograms, FINAL_TAG );	

	// freeing memory
	for ( int i = 0; i < DIM; i++ )
	{
		gsl_histogram_free( histograms[i] );
	}	
}

void add_data_to_histogram( gsl_histogram* h, vector<double> contents )
{
	double lb;
	double ub;

	double x_in_bin;

	for ( int bin_counter = 0; bin_counter < NBINS; bin_counter++ )
	{
		gsl_histogram_get_range( h, bin_counter, &lb, &ub );
	
		// x certainly lies in current bin
		x_in_bin = 0.5 * (lb + ub );

		gsl_histogram_accumulate( h, x_in_bin, contents[bin_counter] );
	}
}	

int main( int argc, char* argv[] )
{
	if ( argc != 4 )
	{
		cout << "USAGE: ./... ( int ) vectors_to_select (int) burn-in-steps (double) alpha" << endl;
		exit( 1 );
	}

	MPI_Init( &argc, &argv );

	const int vectors_to_select = atoi( argv[1] );
	const int burnin = atoi( argv[2] );
	const double alpha = atof( argv[3] );

	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	int world_size;
	MPI_Comm_size( MPI_COMM_WORLD, &world_size );

	int vectors_to_select_per_process = (int) vectors_to_select / ( world_size - 1 );

	if ( rank == 0 )
	{
		cout << "-----------" << endl;
		cout << "Input parameters:" << endl;
		cout << ">> vectors_to_select_in_process: " << vectors_to_select << endl;
		cout << ">> burnin: " << burnin << endl;
		cout << ">> alpha: " << alpha << endl;
		cout << "-----------" << endl;

		int alive = world_size - 1;

		MPI_Status status;
		int source;
		int tag;

		vector< gsl_histogram* > histograms = allocate_histograms();
		vector<double> contents( NBINS );

		while( true )
		{
			if ( alive == 0 )
			{
				// freeing memory
				for ( int i = 0; i < DIM; i++ )
				{
					gsl_histogram_free( histograms[i] );
				}	

				break;
			}

			MPI_Recv( &contents[0], NBINS, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			source = status.MPI_SOURCE;
			tag = status.MPI_TAG;			
			
			if ( tag == FINAL_TAG )
			{
				alive--;
			}

			add_data_to_histogram( histograms[0], contents );

			for ( int histogram_counter = 1; histogram_counter < DIM; histogram_counter++ )
			{
				MPI_Recv( &contents[0], NBINS, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

				add_data_to_histogram( histograms[histogram_counter], contents );

			}
			
			save_histograms( histograms );
		}
	}
	else
	{
		slave_code( rank, vectors_to_select_per_process, burnin, alpha );
	}

	MPI_Finalize();

	return 0;
}
