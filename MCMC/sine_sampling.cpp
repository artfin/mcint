#include <iostream>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

default_random_engine generator;
uniform_real_distribution<double> distribution( 0.0, 1.0 );

double custom_distribution( void )
{
	double u = distribution( generator ); 
	return acos( 1 - 2 * u);
}

int main( int argc, char* argv[] )
{
	if ( argc != 2 )
	{
		cout << "USAGE: ./.. (int) number-of-samples-to-be-generated" << endl;
		exit( 1 );
	}

	int n = atoi( argv[1] );

	for ( int i = 0; i < n; i++ )
	{
		cout << custom_distribution() << endl;
	}

	return 0;
}
