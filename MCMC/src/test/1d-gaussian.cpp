#include <iostream>
#include <random>
#include <chrono>

using namespace std;

static mt19937 generator;

static double nextDouble( const double &min = 0.0, const double &max = 1.0 )
{
    uniform_real_distribution<double> distribution( min, max );
    return distribution( generator );
}

double pdf_gaussian( double x )
{
    return (1 / sqrt( 2 * M_PI )) * exp( -0.5 * pow(x, 2.0));
}

double target( double x )
{
    double s = sin( x );
    double s2 = sin( 2 * x );
    return pow(s, 2) * pow(s2, 2) * pdf_gaussian( x );
}

double metropolis( double x, double alpha )
{
    double y = nextDouble( x - alpha, x + alpha );
    if ( nextDouble() > target( y ) / target( x ))
    {
        y = x;
    } 

    return y;
}

int main( int argc, char* argv[] )
{
	if ( argc != 4 )
	{
		cerr << "USAGE: ./.. (int) moves-to-be-done (int) burn-in (double) alpha" << endl;
		exit( 1 );
	}

    int nsteps = atoi( argv[1] );
	int burnin = atoi( argv[2] );
	double alpha = atof( argv[3] );

	chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();

    double x = 3.14;

	for ( int i = 0; i < burnin; i++ )
	{
		x = metropolis( x, alpha );
	}

	cerr << "Burn-in finished. Time elapsed: " << chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - startTime).count() / 1000.0 << " s" << endl; 

	double xnew;

	int moves = 0;
	int attempted_moves = 0;
    while ( moves < nsteps )
	{
        xnew = metropolis( x, 0.2 );
    
		if ( xnew != x )
		{
			x = xnew;
			cout << x << endl;

			moves++;
		}

		attempted_moves++;
	}

	cerr << "Total time: " << chrono::duration_cast<chrono::milliseconds>( chrono::high_resolution_clock::now() - startTime).count() / 1000.0 << " s" << endl;

	cerr << "Attempted moves: " << attempted_moves << "; moves: " << moves << "; percentage: " << (double) moves / attempted_moves * 100.0 << "%" << endl;

    return 0;
}
