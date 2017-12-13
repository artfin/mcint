#include <mpi.h>
#include <iostream>
//#include <stdio.h>
//#include <stdlib.h>
#include <chrono>
#include <fstream>
#include <iomanip>

// string conversion
#include <string>
#include <sstream>

#include <algorithm> // std::min

#include <Eigen/Dense>

#include<assert.h>
#include <boost/random.hpp>

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

// reduced mass of ar and co2 = m(ar) * m(co2) / ( m(ar) + m(co2) ) in kg
const double MUKG = ( 40.0 * 44.0 ) / ( 40.0 + 44.0 ) * DALTON;
const long double MUAMU = MUKG / AMU; 

// temperature in K
const double temperature = 300;

static mt19937 generator;

static double nextDouble( const double &min = 0.0, const double &max = 1.0 )
{
	 boost::random::uniform_real_distribution<double> distribution( min, max );
    return distribution( generator );
}

static void nextGaussianVec( VectorXf &v, VectorXf mean, const double sigma )
{
	for ( int i = 0; i < 4; i++ )
	{
			boost::random::normal_distribution<double> d( mean(i), sigma );
		v(i) = d( generator );
	}
}

// x = [R, pR, jx, jy]
double target( VectorXf x )
{
	double R = x(0);
	double pR = x(1);
	double jx = x(2);
	double jy = x(3);

	double h = pow(pR, 2) /  ( 2 * MUAMU ) + ( pow(jx, 2) + pow(jy, 2) ) / (2 * MUAMU * pow(R, 2));

	return exp( -h * HTOJ / (BOLTZCONST * temperature) );
}

VectorXf metro_step( VectorXf x, double alpha )
{
	VectorXf prop(4);

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
		cerr << "USAGE: ./.. (int) moves-to-be-recorded (int) burn-in steps (double) alpha" << endl;
		exit(1);
	}

	int nsteps = atoi( argv[1] );
	int burnin = atoi( argv[2] );
	double alpha = atof( argv[3] );

	// initializing MPI-environment
	MPI_Init( NULL, NULL );

	// getting id of the current process
	int world_rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );
	//cout << "World rank: " << world_rank << endl;

	// getting number of running processes
	//int world_size;
	//MPI_Comm_size( MPI_COMM_WORLD, &world_size );

	ostringstream strs;
	strs << world_rank;
	string filename = "output/" + strs.str() + "_chain.dat";
	cout << "Process " << world_rank << "; Outputing to: " << filename << endl;

	ofstream file;
	file.open( filename );

	file << "# Metropolis-Hastings sampler for full phase space of diatomics" << endl;
	file << "# moves-to-be-made = " << nsteps << endl;
	file << "# burn-in = " << burnin << endl;
	file << "# alpha = " << alpha << endl;
	file << "# temperature = " << temperature << endl;
	file << "# structure: R, pR, jx, jy" << endl;

	size_t moves = 0;
	size_t total_steps = 0;

	auto startTime = chrono::high_resolution_clock::now();
	
	VectorXf x(4);
	VectorXf xnew(4);
	x << 5.0, 1.0, 1.0, 1.0;

	// burn-in 
	for ( size_t i = 0; i < burnin; i++ )
	{
		x = metro_step( x, alpha );
	}

	auto endTime = chrono::high_resolution_clock::now();

	cerr << ">> Process " << world_rank << "; burn-in stage completed. Time elapsed: " << chrono::duration_cast<chrono::milliseconds>( endTime - startTime ).count() / 1000.0 << "s" << endl;

	startTime = chrono::high_resolution_clock::now();

	// recording
	while ( moves < nsteps )
	{
		xnew = metro_step( x, alpha );

		if ( xnew != x )
		{
			moves++;
			file << x(0) << " " << x(1) << " " << x(2) << " " << x(3) << endl;
		}

		x = xnew;
		total_steps++;
	}

	endTime = chrono::high_resolution_clock::now();

	cerr << "--------------------------------" << endl;
	cerr << "Process " << world_rank << endl;
	cerr << "Total steps: " << total_steps << "; moves: " << moves << "; percent: " << (double) moves / total_steps * 100 << "%" << endl;
	cerr << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds>( endTime - startTime ).count() / 1000.0 << "s" << endl;

	MPI_Finalize();

	return 0;	
}
