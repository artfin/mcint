#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

#include <iomanip> // std::atoi
#include <algorithm> // std::min

// Eigen
#include <Eigen/Dense>

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
const double RDIST = 50.0;
// -----------------------------

const double BBOUND = 10.0;
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

void W_matrix( Matrix<double, 3, 3> &W, double &theta, double &psi )
{
	double sin_psi = sin(psi);
	double cos_psi = cos(psi);

	double sin_theta = sin(theta);
	double cos_theta = cos(theta);
	
	W(0, 0) = sin_psi / sin_theta;
	W(0, 1) = cos_psi;
	W(0, 2) = - sin_psi * cos_theta / sin_theta;

	W(1, 0) = cos_psi / sin_theta;
	W(1, 1) = - sin_psi;
	W(1, 2) = - cos_psi * cos_theta / sin_theta;

	W(2, 0) = 0;
	W(2, 1) = 0;
	W(2, 2) = 1;
}

// x = [Theta, pR, pT, phi, theta, psi, p_phi, p_theta, p_psi ]
double calc_gunsight( VectorXf x )
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

	cout << "#####" << endl;
	cout << "Jx: " << Jx << endl;
	cout << "Jy: " << Jy << endl;
	cout << "Jz: " << Jz << endl;
	cout << "#####" << endl;

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

	double gunsight = RDIST * sqrt( 1 - pow( cos_phi, 2 ));
	
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
		
		// wrapping Theta to [0, Pi)
		xnew(0) = wrapMax( xnew(0), M_PI );
		xnew(3) = wrapMax( xnew(3), 2 * M_PI );
		xnew(4) = wrapMax( xnew(4), M_PI );
		xnew(5) = wrapMax( xnew(5), 2 * M_PI );	

		double gunsight = calc_gunsight( xnew );

        if ( xnew != x )
        {
            moves++;
        } 
        
		// pR < 0.0
		if ( show_vecs == true && xnew != x && moves % 10 == 0 &&
		  	 gunsight > 0 && gunsight < BBOUND )
        {
            cout << wrote_vectors + 1 << " " 
				 << RDIST << " " 
				 << xnew(0) << " " 
				 << xnew(1) << " " 
				 << xnew(2) << " " 
				 << xnew(3) << " " 
				 << xnew(4) << " " 
				 << xnew(5) << " "
				 << xnew(6) << " " 
				 << xnew(7) << " "
				 << xnew(8) << " " << endl;

			wrote_vectors++;
        }

		attempted_steps++;

        x = xnew;
    }

	chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();

    cerr << "-----------------------------------" << endl;
	cerr << "Burnin: " << burnin << endl;
    cerr << "Total steps: " << attempted_steps  << "; moves: " << moves << "; percent: " << (double) moves / attempted_steps * 100 << "%" << endl;
    cerr << "Percentage of points choosed (to moves): " << (double) wrote_vectors / moves * 100 << "%; (to attempted_steps): " << (double) wrote_vectors / attempted_steps * 100 << "%" << endl;
	cerr << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds>(endTime - block_times[0]).count() / 1000.0 << " s" << endl; 
    cerr << "-----------------------------------" << endl;
	
	return 0;
}

