#include "constants.hpp"

#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <sstream>
#include <fstream>
#include <gsl/gsl_histogram.h>

using namespace std;

const double MU = ( MCONST::HE_MASS * MCONST::AR_MASS ) / ( MCONST::HE_MASS + MCONST::AR_MASS) * WCONST::DALTON / WCONST::AMU;

const double temperature = 300;

const int NBINS = 100;

// -----------------------------------------
// distance between atoms in a0 
const double RDIST = 30.0;
const double RDIST2 = pow(RDIST, 2);
// -----------------------------------------

static thread_local std::mt19937 generator;

void save_histogram( gsl_histogram *histogram, string filename )
{
	ofstream file( filename );

	double lower_bound, higher_bound, bin_content;
	for ( int counter = 0; counter < NBINS; counter++ )
	{
		gsl_histogram_get_range( histogram, counter, &lower_bound, &higher_bound );
		bin_content = gsl_histogram_get( histogram, counter );

		//cout << lower_bound << " " << higher_bound << " " << bin_content << endl;
		file << lower_bound << " " << higher_bound << " " << bin_content << endl;
	}
	
	file.close();	
}

void gsl_histogram_normalize( gsl_histogram* h )
{
	double max = gsl_histogram_max( h );
	double min = gsl_histogram_min( h );
	double step = (max - min) / NBINS;

	double sum = gsl_histogram_sum( h ) * step;
	cout << "sum: " << sum << endl;

	gsl_histogram_scale( h, 1.0 / sum );
}

// normally distributed random vector ( = \dot{\vec{r}} )
double nextGaussian( const double &mean, const double &sigma )
{
    normal_distribution<double> d( mean, sigma );
    return d( generator ); 
} 

int main( int argc, char* argv[] )
{
	if ( argc != 2 )
	{
		cout << "USAGE: ./... (int) n" << endl;
		exit( 1 );
	}

	int n = atoi( argv[1] );
	cerr << "Given value of n: " << n << endl;

	clock_t start = clock();

	gsl_histogram* pr_histogram = gsl_histogram_alloc( NBINS ); 
	gsl_histogram_set_ranges_uniform( pr_histogram, -10.0, 10.0 );

	gsl_histogram* jx_histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( jx_histogram, -300.0, 300.0 );

	double jx, jy, pr;

	double sigma_jx = sqrt(WCONST::BOLTZCONST * temperature / WCONST::HTOJ * MU * RDIST2 );
	double sigma_pr = sqrt(WCONST::BOLTZCONST * temperature / WCONST::HTOJ * MU ); 

	//cout << "sigma_jx: " << sigma_jx << endl;
	//cout << "sigma_pr: " << sigma_pr << endl;

	// nextGaussian take mu and SIGMA!!! not sigma**2!! 
	for ( size_t i = 0; i < n; i++ )
	{
		jx = nextGaussian( 0, sigma_jx );
 		pr = nextGaussian( 0, sigma_pr );

		gsl_histogram_increment( jx_histogram, jx );
	    gsl_histogram_increment( pr_histogram, pr );	   
	}
	
	gsl_histogram_normalize( pr_histogram );
	gsl_histogram_normalize( jx_histogram );
	
	char s[64];
	sprintf( s, "%.2e", (double) n );
	string ss{ s };

	sprintf( s, "%d", (int) NBINS );
	string binstring{ s };

	save_histogram( pr_histogram, "pr_exact_" + ss + "_" + binstring + "bins.txt" );
	save_histogram( jx_histogram, "jx_exact_" + ss + "_" + binstring + "bins.txt" );	

	gsl_histogram_free( pr_histogram );
	gsl_histogram_free( jx_histogram );	

	cout << "Time elapsed: " << (clock() - start) / (double) CLOCKS_PER_SEC << " s" << endl;

	return 0;
}




