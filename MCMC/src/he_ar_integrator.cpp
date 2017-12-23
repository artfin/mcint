#include "hep/mc.hpp"
#include "integrator.hpp"
#include "integrand.hpp"
#include "constants.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cassert>
#include <stdexcept>

using namespace std;

const double MU = ( MCONST::HE_MASS * MCONST::AR_MASS ) / ( MCONST::HE_MASS + MCONST::AR_MASS) * WCONST::DALTON / WCONST::AMU;
const double MU_SI = MU * WCONST::AMU;

const double temperature = 300.0; // K
const double RDIST = 30.0; // a0

// x = [ pR, theta, pT ]
double hamiltonian( hep::mc_point<double> const& x )
{
	double pR = x.point()[0];
	double theta = x.point()[1];
	double pT = x.point()[2];

	double h = pow(pR, 2) / (2 * MU) + pow(pT, 2) / (2 * MU * pow(RDIST, 2));
	return exp( -h * WCONST::HTOJ / (WCONST::BOLTZCONST * temperature) );
}

vector<double> create_distribution_x_array( const double& lb, const double& ub, const int& ndots )
{
	vector<double> v;
	
	double step = (ub - lb) / (ndots - 1);
	for ( double s = lb; s <= ub; s += step )
		v.push_back( s );

	return v;
}

void save_distribution( const string& filename, vector<double>& vectorx, vector<double>& vectory )
{
	assert( vectorx.size() == vectory.size() );

	ofstream file( filename );
	
	for ( size_t i = 0; i < vectorx.size(); i++ )
	{
		file << vectorx[i] << " " << vectory[i] << endl;
	}

	file.close();
}

void normalize_distribution( vector<double>& vectorx, vector<double>& vectory )
{
	assert( vectorx.size() == vectory.size() );

	double step = vectorx[1] - vectorx[0];
	double sum = 0;
	
	for ( size_t i = 0; i < vectorx.size(); i++ )
		sum += vectory[i] * step;

	for ( size_t i = 0; i < vectory.size(); i++ )
		vectory[i] /= sum;
}

int main( int argc, char* argv[] )
{
	if ( argc != 5 )
	{
		cout << "USAGE: ./... (string) distribution-name (double) left-bound (double) right-bound (int) ndots" << endl
			<< "Examples: " << endl
			<< ">> pR -20.0 20.0 50" << endl
			<< ">> theta 0.0 3.14 50" << endl 
			<< ">> pT -100.0 100.0 50" << endl;
		exit( 1 );
	}

	string distribution_name = argv[1];
	double lb = atof( argv[2] );
	double rb = atof( argv[3] );
	int ndots = atoi( argv[4] );

	MPI_Init( &argc, &argv );
	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	Integrand integrand( hamiltonian );
	integrand.set_limits()->add_limit( 0, "-inf", "+inf" )
						  ->add_limit( 1, 0.0, 2 * M_PI )
						  ->add_limit( 2, "-inf", "+inf" );

	if ( distribution_name == "pR" )
		integrand.set_argn( 0 );
	else if ( distribution_name == "theta" )
		integrand.set_argn( 1 );
	else if ( distribution_name == "pT" )
		integrand.set_argn( 2 );
	else
		throw invalid_argument( "Unknown distribution name!" );

	Integrator integrator( integrand, rank );
	//integrator.set_callback();

	vector<double> distribution_x_array = create_distribution_x_array( lb, rb, ndots );
	vector<double> distribution_y_array;

	double v;
	for ( size_t i = 0; i < distribution_x_array.size(); i++ ) 
	{
		integrand.set_default_value(  distribution_x_array[i] );
		
		v = integrator.run_integration( 30, 10000 );
		distribution_y_array.push_back( v );
	}

	if ( rank == 0 )
	{
		normalize_distribution( distribution_x_array, distribution_y_array );

		save_distribution( "pr_integrated_distribution.txt", distribution_x_array, distribution_y_array );
	}

	return 0;
}
