#include <iostream>
#include <iomanip>
#include <random>
#include <ctime>
#include <cmath>

#include <algorithm>

using namespace std;

static mt19937 generator;

static double nextDouble( const double &min = 0.0, const double &max = 1.0 )
{
	uniform_real_distribution<double> distribution( min, max );
	return distribution( generator );
}

static double nextGaussian( const double &mean, const double &sigma )
{
	normal_distribution<double> d( mean, sigma );
	return d( generator );
}

double target( double x )
{
	return exp(- pow(x, 2) );
}

double metro_step( double x, double alpha = 1.0 )
{
	double y = nextGaussian( x, alpha );

	if ( nextDouble() < min( 1.0, target(y) / target(x)) )
	{
		x = y;
	}

	return x;
}

int main( int argc, char* argv[] )
{
	if ( argc != 4 )
	{
		cout << "USAGE: ./... (int) moves-to-made (int) burn-in (double) alpha" << endl;
		exit( 1 );
	}

	int nsteps = atoi( argv[1] );
	int burnin = atoi( argv[2] );
	double alpha = atof( argv[3] );

	cout << "# moves-to-made = " << nsteps << endl;
	cout << "# burn-in = " << burnin << endl;
	cout << "# alpha = " << alpha << endl;

	double init = 0.5;
	
	double x = metro_step( init, alpha );

	// burnin cycle
	for ( int i = 0; i < burnin; i++ )
	{
		x = metro_step( x, alpha );
	}

	double sum = 0;
	double x_min = x;
	double x_max = x;

	for ( int i = 0; i < nsteps; i++ )
	{
		x = metro_step( x, alpha );

		sum += target( x );

		if ( x < x_min )
		{
			x_min = x;
		}
		
		if ( x > x_max )
		{
			x_max = x;
		}
	}
	
	sum = sum * ( x_max - x_min ) / nsteps;	

	cout << "x_min: " << x_min << endl;
	cout << "x_max: " << x_max << endl;
	cout << "x_max - x_min : " << x_max - x_min << endl;

	double exact = sqrt( 3.1415 );
	cout << "Exact: " << exact << endl;
	cout << "Monte Carlo estimate: " << sum << endl;

	cout << "ratio: " << sum / exact << endl;

	return 0;
}
