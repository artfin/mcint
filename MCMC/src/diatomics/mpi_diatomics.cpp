#include <iostream>

#include "constants.h"

#include <random>
#include <chrono>
#include <cmath>

#include <vector>
#include <fstream>
#include <string>

#include <functional>

#include <gsl/gsl_histogram.h>
#include <iomanip>
#include <algorithm> // std::min

#include <Eigen/Dense>

using std::cout;
using std::endl;
using std::function;
using std::vector;
using std::pair;
using std::string;

using Eigen::VectorXf;

const int NBINS = 500;

const double he_mass = 4.0;
const double ar_mass = 40.0;

const double MU = ( he_mass * ar_mass ) / (he_mass + ar_mass) * constants::DALTON / constants::AMU;

const double temperature = 300.0;

static vector<pair<int,double>> DEFAULT_VECTOR;
static int BLOCK_SIZE = 1e6;

class MCMC
{
	public:
		int DIM;
		int MOVES;
		double alpha;
		vector<gsl_histogram*> histograms;
		vector<string> names;

		bool set_histograms = false;

		VectorXf point_after_burnin;
		function<double(VectorXf)> f;
		std::mt19937 generator;
	
		double nextDouble( const double &min, const double& max );
		void nextGaussianVec( VectorXf &v, VectorXf &mean );
		double wrapMax( const double& x, const double& max );

		void initialize_histograms( vector<gsl_histogram*> histograms,
					   				vector<string> names
								  );
		void save_histogram( gsl_histogram *histogram, string filename );
		void gsl_histogram_normalize( gsl_histogram* h );

		VectorXf metro_step( VectorXf& x );
		void burnin( VectorXf& initial_point, const int& burnin_length );
		void run_chain( 
			vector<pair<int, double>>& to_wrap = DEFAULT_VECTOR, 
			int& block_size = BLOCK_SIZE 
					  );

		MCMC( function<double(VectorXf)> f, const int &MOVES, const int &DIM, const double& alpha );
		~MCMC();
};

MCMC::MCMC ( function<double(VectorXf)> f, const int& MOVES, const int& DIM, const double& alpha )
{
	this->DIM = DIM;
	this->MOVES = MOVES;
	this->alpha = alpha;
	this->f = f;
}

MCMC::~MCMC()
{
	if ( set_histograms )
	{
		for ( int i = 0; i < DIM; i++ )
		{
			gsl_histogram_free( histograms[i] );
		}
	}
}

void MCMC::gsl_histogram_normalize( gsl_histogram* h )
{
	double max = gsl_histogram_max( h );
	double min = gsl_histogram_min( h );
	double step = (max - min) / NBINS;

	double sum = gsl_histogram_sum( h ) * step;
	gsl_histogram_scale( h, 1.0 / sum );
}

void MCMC::save_histogram( gsl_histogram *histogram, string filename )
{
	std::ofstream file( filename );

	double lower_bound, higher_bound, bin_content;
	for ( int counter = 0; counter < NBINS; counter++ )
	{
		gsl_histogram_get_range( histogram, counter, &lower_bound, &higher_bound );
		bin_content = gsl_histogram_get( histogram, counter );

		file << lower_bound << " " << higher_bound << " " << bin_content << endl;
	}
	
	file.close();	
}

void MCMC::initialize_histograms( vector<gsl_histogram*> histograms,
			   					  vector<string> names )
{
	this->set_histograms = true;
	this->histograms = histograms;
	this->names = names;
}

double MCMC::nextDouble( const double& min, const double& max )
{
	std::uniform_real_distribution<double> distribution( min, max );
	return distribution( this->generator );
}

void MCMC::nextGaussianVec( VectorXf &v, VectorXf &mean )
{
	for ( int i = 0; i < this->DIM; i++ )
	{
		std::normal_distribution<double> d( mean(i), this->alpha );
		v(i) = d( this->generator );
	}
}

// wrap x -> [0, max)
double MCMC::wrapMax( const double& x, const double& max )
{
	return fmod( max + fmod(x, max), max );
}

VectorXf MCMC::metro_step( VectorXf& x )
{
	VectorXf prop( this->DIM );
	nextGaussianVec( prop, x );

	if ( nextDouble(0.0, 1.0) < std::min( 1.0, this->f(prop) / this->f(x) ))
	{
		return prop;
	}
	
	return x;	
}

void MCMC::burnin( VectorXf& initial_point, const int& burnin_length )
{
	std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();

	VectorXf x = metro_step( initial_point );

	for ( size_t i = 0; i < burnin_length; i++ )
	{
		x = metro_step( x );
	}

	cout << "Burnin finished. Time elapsed: " << 
		std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - startTime).count() / 1000.0 << " s " << endl;

	this->point_after_burnin = x;
}

void MCMC::run_chain( vector<pair<int, double>>& to_wrap, int& block_size )
{
	int moves = 0;
	int attempted_steps = 0;	

	VectorXf x = this->point_after_burnin;
	VectorXf xnew( this->DIM );

	std::chrono::milliseconds time_for_block;
	std::chrono::milliseconds time_for_blocks;

	vector<std::chrono::high_resolution_clock::time_point> block_times;
	block_times.push_back( std::chrono::high_resolution_clock::now() );

	int block_counter = 0;

	while ( moves < this->MOVES )
	{
		if ( attempted_steps % block_size == 0 && attempted_steps != 0 )
		{
			block_counter++;

			block_times.push_back( std::chrono::high_resolution_clock::now() );

			time_for_block = std::chrono::duration_cast<std::chrono::milliseconds>( block_times.end()[-1] - block_times.end()[-2] );
			time_for_blocks = std::chrono::duration_cast<std::chrono::milliseconds>( block_times.end()[-1] - block_times[0] );

			cout << endl << "Block " << block_counter << " finished." << endl;
			cout << "Attempted steps: " << attempted_steps << "; moves made: " << moves << endl;
			cout << "Time for current block: " << time_for_block.count() / 1000.0 << " s" << endl;
			cout << "Total time elapsed: " << time_for_blocks.count() / 1000.0 << " s" << endl;
		}
		
		attempted_steps++;

		xnew = metro_step( x );

		if ( to_wrap != DEFAULT_VECTOR )
		{
			int curr_var; // number of current variable
			int curr_max; // maximum of current variable

			for ( size_t i = 0; i < to_wrap.size(); i++ )
			{
				curr_var = to_wrap[i].first;
				curr_max = to_wrap[i].second;

				xnew(curr_var) = wrapMax( xnew(curr_var), curr_max );
			}
		}

		if ( xnew != x )
		{
			x = xnew;
			moves++;
		}

		if ( set_histograms )
		{
			for ( size_t i = 0; i < DIM; i++ )
			{
				gsl_histogram_increment( this->histograms[i], xnew(i) );
			}
		}
	}
	
	if ( set_histograms )
	{
		for ( size_t i = 0; i < DIM; i++ )
		{
			gsl_histogram_normalize( this->histograms[i] );
			save_histogram( this->histograms[i], this->names[i] );
		}
	}

    cout << "-----------------------------------" << endl;
   	cout << "Attempted steps: " << attempted_steps << "; moves: " << moves << "; percent: " << (double) moves / attempted_steps * 100.0 << "%" << endl;
   	cout << "Time elapsed: " << std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::high_resolution_clock::now() - block_times[0]).count() / 1000.0 << " s" << endl; 
   	cout << "-----------------------------------" << endl;
}

// x = [v0, phi]
double target( VectorXf x )
{
	double v0 = x(0);
	double phi = x(1);

	double h = MU * pow(v0, 2) / 2.0 * (1.0 + sin(phi) * sin(phi));
	return exp( -h * constants::HTOJ / (constants::BOLTZCONST * temperature) );
}

int main ( int argc, char* argv[] )
{
	const int MOVES = 1000;
	const int DIM = 2;
	const double ALPHA = 0.001;
	MCMC diatomic_b_v0( target, MOVES, DIM, ALPHA );

	Eigen::VectorXf x(DIM);
   	x << -0.0005, 0.0;
	
	// running burn-in cycle
	diatomic_b_v0.burnin( x, 1000 );

	// wrapping second argument: phi to [0, 2\pi]
	pair<int, double> p1(1, 2*M_PI); 
	vector<pair<int, double>> to_wrap;
	to_wrap.push_back( p1 );

	// ############################################################
	gsl_histogram* v0_histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( v0_histogram, -1e-3, 1e-3 );

	gsl_histogram* phi_histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( phi_histogram, 0.0, 2 * M_PI );

	vector<gsl_histogram*> histograms;
	histograms.push_back( v0_histogram );
	histograms.push_back( phi_histogram );

	vector<string> names;
	names.push_back( "v0_histogram.txt" );
	names.push_back( "phi_histogram.txt" );
	// ############################################################

	diatomic_b_v0.initialize_histograms( histograms, names );

	diatomic_b_v0.run_chain( to_wrap ); 
}
