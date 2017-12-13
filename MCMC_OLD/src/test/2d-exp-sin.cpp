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
    return exp( - pow( x(0), 2 ) / sin( x(1) ) ); 
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
    if ( argc != 5 )
    {
        cerr << "USAGE: ./... (double) alpha (int) nsteps (int) burnin (bool) show_steps" << endl;
        exit( 1 ); 
    }

    double ALPHA = atof( argv[1] ); 
    int nsteps = atoi( argv[2] );
    int burnin = atoi( argv[3] );
    bool showsteps = atoi( argv[4] );

    Vector2d ic ( 0.1, 0.1 );
    Vector2d x = metro_step( ic, 1.0 );
    Vector2d xnew;

    clock_t start = clock();
    int moves = 0;

    for ( int i = 0; i < nsteps + burnin; i++ )
    {
        xnew = metro_step( x, ALPHA );
    
        if ( xnew != x )
        {
            moves++;

            if ( i > burnin && showsteps == true )
            {
                cout << x(0) << " " << x(1) << endl;
            }
        }

        x = xnew;
    }

    cerr << "---------------------------------------" << endl;
    cerr << "total steps: " << nsteps + burnin << "; moves: " << moves << "; percent: " << (double) moves / ( nsteps + burnin ) * 100 << "%" << endl;
    cerr << "Time elapsed: " << (double) ( clock() - start ) / CLOCKS_PER_SEC << "s" << endl;
    cerr << "---------------------------------------" << endl;

    return 0;
}
