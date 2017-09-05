#include <iostream>
#include <random>
#include <ctime>

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

double metropolis( double x, double alpha = 1.0 )
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
    int nsteps = atoi( argv[1] );
    //cout << "Given number of steps: " << nsteps << endl;

    double x = 3.14;
    for ( int i = 0; i < nsteps; i++ )
    {
        x = metropolis( x, 0.2 );
        cout << x << endl;
    }

    return 0;
}
