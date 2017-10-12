#include <iostream>

#include <mpi.h>
#include <cmath>
#include <random>
#include <chrono>

using namespace std;

int main( int argc, char* argv[] )
{
	MPI_Init( &argc, &argv );

	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	int world_size;
	MPI_Comm_size( MPI_COMM_WORLD, &world_size );

	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	uniform_real_distribution<double> dist( 0, 36 );
	mt19937 rng( seed * rank );

	cout << "Process " << rank << "; random = " << dist(rng) << endl;

	MPI_Finalize();

	return 0;
}
