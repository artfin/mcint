#include <iostream>
#include <random>
#include <ctime>
#include <cmath>

#include <iomanip>   // std::atoi
#include <algorithm> // std::min

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

static mt19937 generator;

static double nextDouble( const double &min = 0.0, const double &max = 1.0 )
{
    uniform_real_distribution<double> distribution( min, max );
    return distribution( generator );
}

static void nextGaussianVec( Vector2d &v, Vector2d mean, const double &sigma )
{
    normal_distribution< double > d0( mean(0), sigma );
    normal_distribution< double > d1( mean(1), sigma );
    v(0) = d0( generator );
    v(1) = d1( generator );
}

double target( Vector2d x )
{
    return exp( - 5 * fabs( pow(x(0), 2) + pow(x(1), 2) - 1 ));
}

Vector2d metro_step( Vector2d x, double alpha )
{
    Vector2d prop;
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
    int nsteps = atoi( argv[1] );
    int burnin = atoi( argv[2] );
    // cout << "nsteps: " << nsteps << endl;

    Vector2d ic ( 0.0, 0.0 );
    Vector2d x = metro_step( ic, 1.0 );

    for ( int i = 0; i < nsteps + burnin; i++ )
    {
        x = metro_step( x, 0.3 );
    
        // outputting the result only after burn-in
        if ( i > burnin )
        {
            cout << x(0) << " " << x(1) << endl;
        } 
    }

    return 0;
}
