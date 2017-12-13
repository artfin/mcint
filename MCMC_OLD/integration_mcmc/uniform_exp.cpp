#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <chrono>

using namespace std;

static mt19937 generator;

static double nextDouble( const double &min = 0.0, const double &max = 1.0 )
{
	uniform_real_distribution<double> distribution( min, max );
	return distribution( generator );
}

double integrand( double x )
{
	return exp( -pow(x, 2) );
}

int main( int argc, char* argv[] )
{
	if ( argc != 2 )
	{
		cout << "USAGE: ./.. (int) nsteps" << endl;
		cout << "Program implements dumb Monte-Carlo integration of exp(-x**2)" << endl;
	}

	auto startTime = chrono::high_resolution_clock::now();

	int nsteps = atoi( argv[1] );

	double x_min = -2.5;
	double x_max =  2.5;

	double summ = 0;
	for ( int i = 0; i < nsteps; i++ )
	{
		summ += integrand( nextDouble( x_min, x_max ));
	}

	summ = summ * (x_max - x_min) / nsteps;

	double exact = sqrt( M_PI );

	double abs_error = fabs( exact - summ );
	double rel_error = abs_error / exact;
	
	cout << "Exact result: " << exact << endl;
	cout << "Monte-Carlo estimate: " << summ << endl;
	cout << "abs error: " << abs_error << endl;
	cout << "rel error: " << rel_error * 100 << "%" << endl;

	auto endTime = chrono::high_resolution_clock::now();

	cout << "Time elapsed: " << chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() << "ms" << endl;

	return 0;
}

