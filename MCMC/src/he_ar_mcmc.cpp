#include "mcmc.hpp"
#include "constants.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include <Eigen/Dense>
#include <gsl/gsl_histogram.h>

const double MU = ( MCONST::HE_MASS * MCONST::AR_MASS ) / ( MCONST::HE_MASS + MCONST::AR_MASS) * WCONST::DALTON / WCONST::AMU;
const double MU_SI = MU * WCONST::AMU;

const double temperature = 300.0;

const double RDIST = 30.0; // a0

// x = [ pR, theta, pT ]
double target( VectorXf x )
{
	double pR = x( 0 );
	double theta = x( 1 );
	double pT = x( 2 );

	double h = pow(pR, 2) / (2 * MU) + pow(pT, 2) / (2 * MU * pow(RDIST, 2));
	return exp( -h * WCONST::HTOJ / (WCONST::BOLTZCONST * temperature) );
}

int main ( int argc, char* argv[] )
{
	if ( argc < 3 || argc > 4 )
	{
		cout << "Usage: ./... (int) moves-to-be-made (double) alpha (int; optional) burnin" << endl;
	   	exit( 1 );	
	}

	int MOVES = atoi( argv[1] );
	double ALPHA = atof( argv[2] );
	int burnin = 1e4;	

	if ( argc == 4 ) burnin = atoi( argv[3] );

	const int DIM = 3;
	MCMC diatomic( target, MOVES, DIM, ALPHA );

	Eigen::VectorXf x(DIM);
   	x << 1.0, 0.1, 0.0;

	// running burn-in cycle
	diatomic.burnin( x, burnin );

	// wrapping second argument (argument 1):
	pair<int, double> p1(1, 2*M_PI); 
	vector<pair<int, double>> to_wrap;
	to_wrap.push_back( p1 );

	// ############################################################
	gsl_histogram* pr_histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( pr_histogram, -10.0, 10.0 );

	gsl_histogram* theta_histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( theta_histogram, 0.0, 2*M_PI ); 

	gsl_histogram* pt_histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( pt_histogram, -300.0, 300.0 );

	vector<gsl_histogram*> histograms;
	histograms.push_back( pr_histogram );
	histograms.push_back( theta_histogram );
	histograms.push_back( pt_histogram );

	char s[64];
	sprintf( s, "%.2e", (double) MOVES );
	string ss{ s };

	sprintf( s, "%d", (int) NBINS );
	string stringbins{ s };

	vector<string> names;
	names.push_back( "pr_mcmc_" + ss + "_" + stringbins + "bins.txt" );
	names.push_back( "theta_mcmc_" + ss + "_" + stringbins + "bins.txt" );
	names.push_back( "pt_mcmc_" + ss + "_" + stringbins + "bins.txt" );
	// ############################################################

	diatomic.initialize_histograms( histograms, names );

	diatomic.run_chain( to_wrap ); 

	return 0;
}

