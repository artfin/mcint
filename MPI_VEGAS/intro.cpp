#include "hep/mc-mpi.hpp"

#include <mpi.h>
#include <chrono>

#include <cstddef>
#include <iostream>
#include <vector>

using namespace std;

double square( hep::mc_point<double> const& x )
{
	return 3.0 * x.point()[0] * x.point()[0];
}

int main( int argc, char* argv[] )
{
	// Initialize MPI
	MPI_Init( &argc, &argv ); // why passing?

	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	// print reference result, but only for processor with rank zero ( to avoid printing it more than once )
	if ( rank == 0 )
	{
		cout << ">> computing integral of 3 * x**2 from 0 to 1 which is 1.0" << endl << endl;
	}

	// set the verbose callback function
	hep::mpi_vegas_callback<double>( hep::mpi_vegas_verbose_callback<double> );

	auto startTime = chrono::high_resolution_clock::now();

	// perform 5 iteration with 1000 calls each; this function will also call
	// vegas verbose_callback after each iteration which in turn prints the
	// individual iterations
	auto results = hep::mpi_vegas(
		MPI_COMM_WORLD,
		hep::make_integrand<double>( square, 1 ),
		vector<size_t>( 5, 1e7 )
	);

	auto result = hep::cumulative_result0( results.begin() + 1, results.end());
	double chi_square_dof = hep::chi_square_dof0( results.begin() + 1, results.end() );

	auto endTime = chrono::high_resolution_clock::now();

	if ( rank == 0 )
	{
		// print the cumulative result
		cout << ">> cumulative result: " << endl;
		cout << ">> N = " << result.calls() << "; I = " << result.value() << " +- " << result.error() << "; chi^2/dof = " << chi_square_dof << endl;

		cout << "Time needed: " << chrono::duration_cast<chrono::milliseconds>( endTime - startTime ).count() / 1000.0 << "s" << endl;
	}	

	// clean up
	MPI_Finalize();

	return 0;
}
