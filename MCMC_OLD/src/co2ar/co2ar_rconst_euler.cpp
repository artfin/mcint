#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

#include <iomanip> // std::atoi
#include <algorithm> // std::min

// Eigen
#include <Eigen/Dense>

#include "co2ar_hamiltonian.hpp"

#include <fstream>
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
const double RDIST = 20.0;
// -----------------------------

const int NBINS = 500;

const double BBOUND = 1000.0;
const int DIM = 9;

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

double cot( double x )
{
	return 1.0 / tan( x );
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

void transform_l_to_lab( VectorXf x, vector<double>& l_lab )
{
	// x = [Theta, pR, pT, phi, theta, psi, p_phi, p_theta, p_psi ]
	double Theta = x( 0 );
	double pR = x( 1 );
	double pT = x( 2 );
	double phi = x( 3 );
	double theta = x( 4 );
	double psi = x( 5 );
	double p_phi = x( 6 );
	double p_theta = x( 7 );
	double p_psi = x( 8 );
	
	Matrix<double, 3, 3> S;
	
	double sin_phi = sin( phi );
	double cos_phi = cos( phi );

	double sin_theta = sin( theta );
	double cos_theta = cos( theta );

	double sin_psi = sin( psi );
	double cos_psi = cos( psi );

	S(0, 0) = cos_psi * cos_phi - cos_theta * sin_phi * sin_psi;
	S(0, 1) = - sin_psi * cos_phi - cos_theta * sin_phi * cos_psi;
	S(0, 2) = sin_theta * sin_phi;

	S(1, 0) = cos_psi * sin_phi + cos_theta * cos_phi * sin_psi;
	S(1, 1) = - sin_psi * sin_phi + cos_theta * cos_phi * cos_psi;
	S(1, 2) = - sin_theta * cos_phi;

	S(2, 0) = sin_theta * sin_psi;
	S(2, 1) = sin_theta * cos_psi;
	S(2, 2) = cos_theta;

	Matrix<double, 3, 3> W;
	W_matrix( W, theta, psi );

	Vector3d pe( p_phi, p_theta, p_psi );
	Vector3d j_vector = W * pe;
	
	double Jz = j_vector( 2 ); 

	Vector3d l_mol( -Jz * cos(Theta) / sin(Theta), pT, Jz ); 

	Vector3d res = l_mol;
	
	l_lab[0] = res(0);
	l_lab[1] = res(1);
	l_lab[2] = res(2);
}


// x = [Theta, pR, pT, phi, theta, psi, p_phi, p_theta, p_psi ]
double calc_gunsight( VectorXf x, gsl_histogram* jx_hist, gsl_histogram* jy_hist, gsl_histogram* jz_hist )
{
	double Theta = x( 0 );
	double pR = x( 1 );
	double pT = x( 2 );
	double phi = x( 3 );
	double theta = x( 4 );
	double psi = x( 5 );
	double p_phi = x( 6 );
	double p_theta = x( 7 );
	double p_psi = x( 8 );

	Matrix<double, 3, 3> W;
	W_matrix( W, theta, psi );

	Vector3d pe( p_phi, p_theta, p_psi );
	Vector3d j_vector = W * pe;
	
	double Jx = j_vector(0);
	double Jy = j_vector(1);
	double Jz = j_vector(2);

	gsl_histogram_increment( jx_hist, Jx );
	gsl_histogram_increment( jy_hist, Jy );
	gsl_histogram_increment( jz_hist, Jz );

	//cout << "p_phi: " << p_phi << "; p_theta: " << p_theta << "; p_psi: " << p_psi << endl;
	//cout << "phi: " << phi << "; theta: " << theta << "; psi: " << psi << endl;
	//cout << "jx: " << Jx << "; Jy: " << Jy << "; Jz: " << Jz << endl;

	if ( pR > 0 )
	{
		return -1;
	}

	double pRmu2 = pR / mu2;

	double quad1 = pow( pRmu2, 2 );

	double quad2 = pow( (Jx + Jz * cot(Theta)) / (mu2 * RDIST), 2 );

	double quad3 = pow( (Jy - pT) / (mu2 * RDIST), 2 );
	
	double numerator = - pRmu2; 
	double denumerator = sqrt( quad1  + quad2 + quad3 );

	double cos_phi = numerator / denumerator;

	double gunsight = RDIST * sqrt(1 - pow(cos_phi, 2));
	
	return gunsight;
}


// x = [Theta, pR, pT, phi, theta, psi, p_phi, p_theta, p_psi ]
double target( VectorXf x )
{
    double ke = kinetic_energy_euler( RDIST, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8] );
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
    if ( argc != 4 )
    {
        cout << "Usage: ./ ... (int) moves_to_be_made (int) burn-in-steps (double) alpha " << endl;
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
    x << 0.5, 0.5, 0.5, 0.0, 0.0, 0.0, 0.5, 0.5, 0.5;

	VectorXf xnew( DIM );

    cout << "# Metropolis-Hastings sample for CO2-Ar" << endl;
    cout << "# nsteps = " << nsteps << endl;
    cout << "# burnin = " << burnin << endl;
    cout << "# alpha = " << alpha << endl;
    cout << "# temperature = " << temperature << endl;
    cout << "# R (a. u.) = " << RDIST << endl;
	
	cout << "# file structure: " << endl;
	cout << "# number R Theta pR pT phi theta psi p_phi p_theta p_psi" << endl;

	// burnin cycle
	for ( size_t i = 0; i < burnin; i++ )
	{
		x = metro_step( x, alpha );
	}

	cerr << "Burn-in is finished." << endl;
	
	double lbs[] = {0.0, -50.0, -50.0, 0.0, 0.0, 0.0, -1000.0, -1000.0, -100.0};
	vector<double> lower_boundaries ( lbs, lbs + sizeof(lbs) / sizeof(double) );
	double ubs[] =  {M_PI, 50.0, 50.0, 2 * M_PI, M_PI, 2 * M_PI, 1000.0, 1000.0, 100.0}; 
	vector<double> upper_boundaries( ubs, ubs + sizeof(ubs) / sizeof(double) );

	vector< gsl_histogram* > histograms;
	for ( int i = 0; i < DIM; i++ )
	{
		gsl_histogram *h = gsl_histogram_alloc( NBINS );
		gsl_histogram_set_ranges_uniform( h, lower_boundaries[i], upper_boundaries[i] );
	  	histograms.push_back( h );	
	}

	gsl_histogram* histogram_b = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( histogram_b, 0.0, BBOUND );	
	
	gsl_histogram* jx_histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( jx_histogram, -1000.0, 1000.0 );	
	
	gsl_histogram* jy_histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( jy_histogram, -1000.0, 1000.0 );	
	
	gsl_histogram* jz_histogram = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( jz_histogram, -100.0, 100.0 );	
	
	double LBL = 50.0;
	gsl_histogram* lx_lab = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( lx_lab, -LBL, LBL );

	gsl_histogram* ly_lab = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( ly_lab, -LBL, LBL );
	
	gsl_histogram* lz_lab = gsl_histogram_alloc( NBINS );
	gsl_histogram_set_ranges_uniform( lz_lab, -LBL, LBL );

	// running cycle
	size_t attempted_steps = 0;
	size_t moves = 0;
	size_t wrote_vectors = 0;

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
			cerr << "Vectors written: " << wrote_vectors << "; percentage (to be done): " << (double) wrote_vectors / nsteps * 100 << "%" << endl;
			cerr << "Moves attempted: " << attempted_steps << "; moves made: " << moves << endl;
			cerr << "Time for current block: " << time_for_block.count() / 1000.0 << " s" << endl;
			cerr << "Total time elapsed: " << time_for_blocks.count() / 1000.0 << " s" << endl;
		}

        xnew = metro_step( x, alpha );
		
		// x[0] == Theta
		// x[1] == pR
		// x[2] == pT
		// x[3] == phi
		// x[4] == theta
		// x[5] == psi
		// x[6] == p_phi
		// x[7] == p_theta
		// x[8] == p_psi
		
		// wrapping angles to [0, Pi)
		xnew(0) = wrapMax( xnew(0), M_PI );
		xnew(3) = wrapMax( xnew(3), 2 * M_PI );
		xnew(4) = wrapMax( xnew(4), M_PI );
		xnew(5) = wrapMax( xnew(5), 2 * M_PI );	
		
		double gunsight = calc_gunsight( xnew, jx_histogram, jy_histogram, jz_histogram );

		// ##################
		vector<double> l_lab( 3 );
		transform_l_to_lab( xnew, l_lab );

		gsl_histogram_increment( lx_lab, l_lab[0] );
		gsl_histogram_increment( ly_lab, l_lab[1] );
		gsl_histogram_increment( lz_lab, l_lab[2] );
		// ##################
			
        if ( xnew != x )
        {
            moves++;
        } 
		
		for ( int i = 0; i < DIM; i++ )
		{
			gsl_histogram_increment( histograms[i], xnew(i) ); 
		}
       
		//gsl_histogram_increment( histogram_b, gunsight );

		// pR < 0.0
		//if ( xnew != x && gunsight > 0 && gunsight < BBOUND )
		//{
				//cout << wrote_vectors + 1 << " " 
						//<< RDIST << " " 
						//<< xnew(0) << " " 
						//<< xnew(1) << " " 
						//<< xnew(2) << " " 
						//<< xnew(3) << " " 
						//<< xnew(4) << " " 
						//<< xnew(5) << " "
						//<< xnew(6) << " " 
						//<< xnew(7) << " "
						//<< xnew(8) << " " << endl;

				//wrote_vectors++;
		//}

		attempted_steps++;

        x = xnew;
    }
	
	string names[] = {"Theta.txt", "pr.txt", "pt.txt", "phi.txt", "theta.txt", "psi.txt", "p_phi.txt", "p_theta.txt", "p_psi.txt"};
		
	for ( int i = 0; i < DIM; i++ )
	{
		gsl_histogram_normalize( histograms[i] );
		save_histogram( histograms[i], names[i] );
		gsl_histogram_free( histograms[i] );
	}

	gsl_histogram_normalize( histogram_b );
	save_histogram( histogram_b, "b.txt" );
	gsl_histogram_free( histogram_b );
	
	gsl_histogram_normalize( jx_histogram );
	save_histogram( jx_histogram, "jx.txt" );
	gsl_histogram_free( jx_histogram );
	
	gsl_histogram_normalize( jy_histogram );
	save_histogram( jy_histogram, "jy.txt" );
	gsl_histogram_free( jy_histogram );
	
	gsl_histogram_normalize( jz_histogram );
	save_histogram( jz_histogram, "jz.txt" );
	gsl_histogram_free( jz_histogram );

	gsl_histogram_normalize( lx_lab );
	save_histogram( lx_lab, "lx_lab.txt" );
	gsl_histogram_free( lx_lab );
	
	gsl_histogram_normalize( ly_lab );
	save_histogram( ly_lab, "ly_lab.txt" );
	gsl_histogram_free( ly_lab );
	
	gsl_histogram_normalize( lz_lab );
	save_histogram( lz_lab, "lz_lab.txt" );
	gsl_histogram_free( lz_lab );
	
	chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();

    cerr << "-----------------------------------" << endl;
	cerr << "Burnin: " << burnin << endl;
    cerr << "Total steps: " << attempted_steps  << "; moves: " << moves << "; percent: " << (double) moves / attempted_steps * 100 << "%" << endl;
    cerr << "Percentage of points choosed (to moves): " << (double) wrote_vectors / moves * 100 << "%; (to attempted_steps): " << (double) wrote_vectors / attempted_steps * 100 << "%" << endl;
	cerr << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds>(endTime - block_times[0]).count() / 1000.0 << " s" << endl; 
    cerr << "-----------------------------------" << endl;
	
	return 0;
}

