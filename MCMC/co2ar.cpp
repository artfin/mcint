#include <iostream>
#include <cmath>
#include <random>
#include <ctime>

#include <iomanip> // std::atoi
#include <algorithm> // std::min

// Eigen
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

const double mu1 = 14579.0; // ?
const double mu2 = 36440.0; // ?
const double l = 4.398;

// boltzmann constant
const double BOLTZCONST = 1.38064e-23;
// hartree to joules
const double HTOJ = 4.35974417e-18;

// ! ---------------------------
const double RDIST = 20.0;
// -----------------------------

// temperature in K
const double temperature = 300;

static mt19937 generator;

static double nextDouble( const double &min = 0.0, const double &max = 1.0 )
{
    uniform_real_distribution<double> distribution( min, max );
    return distribution( generator );
}

static void nextGaussianVec( VectorXf &v, VectorXf mean, const double sigma )
{
    for ( int i = 0; i < 6; i++ )
    {
        normal_distribution<double> d( mean(i), sigma );
        v(i) = d( generator );
    }
} 

void fill_inertia_tensor(Matrix3d &inertia_tensor, double &theta)
{
    double sin_t = sin(theta);
    double cos_t = cos(theta);
    double l2 = l * l;
    double R2 = RDIST * RDIST;

    inertia_tensor(0, 0) = mu1 * l2 * cos_t * cos_t + mu2 * R2; 
    inertia_tensor(1, 0) = 0;
    inertia_tensor(2, 0) = -mu1 * l2 * sin_t * cos_t;

    inertia_tensor(0, 2) = inertia_tensor(2, 0);
    inertia_tensor(1, 2) = 0;
    inertia_tensor(2, 2) = mu1 * l2 * sin_t * sin_t;

    inertia_tensor(0, 1) = 0;
    inertia_tensor(1, 1) = inertia_tensor(0, 0) + inertia_tensor(2, 2);
    inertia_tensor(2, 1) = 0;
}

double hamiltonian( double theta, double pR, double pT, double Jx, double Jy, double Jz )
{
	Matrix2d a;
	a << mu2, 0, 0, mu1 * l * l;

	Matrix2d a_inv = a.inverse();

	Matrix<double, 3, 2> A;
	A(0, 0) = A(0, 1) = 0;
    A(1, 0) = 0;
    A(1, 1) = mu1 * l * l;
    A(2, 0) = A(2, 1) = 0;

	Matrix<double, 3, 3> prod1 = A * a_inv * A.transpose();
	Matrix<double, 3, 2> prod2 = A * a_inv;

	Vector3d j_vector(Jx, Jy, Jz);
	Vector2d p_vector(pR, pT);

	Matrix<double, 3, 3> I;

	fill_inertia_tensor(I, theta);

	Matrix<double, 3, 3> I_inv = I.inverse();

	Matrix<double, 3, 3> G11;
	Matrix<double, 2, 2> G22;
	Matrix<double, 3, 2> G12;

	Matrix<double, 3, 3> t1;
	Matrix<double, 2, 2> t2;

	t1 = I - prod1;
	G11 = t1.inverse();

	t2 = a;
	t2.noalias() -= A.transpose() * I_inv * A;
	G22 = t2.inverse();

	G12.noalias() = - G11 * prod2;

	double ang_term = 0.5 * j_vector.transpose() * G11 * j_vector;
	double kin_term = 0.5 * p_vector.transpose() * G22 * p_vector;
	double cor_term = j_vector.transpose() * G12 * p_vector;

	return ang_term + kin_term + cor_term;
}

// x = [theta, pR, pT, Jx, Jy, Jz ]
double target( VectorXf x )
{
    double h = hamiltonian( x[0], x[1], x[2], x[3], x[4], x[5] );
    return exp( -h * HTOJ / (BOLTZCONST * temperature ));
}

VectorXf metro_step( VectorXf x, double alpha )
{
    VectorXf prop( 6 );

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
        cout << "Usage: ./ ... (int) record_steps (int) burn-in-steps (double) alpha (bool) show_vec" << endl;
        exit( 1 );
    }

    const int nsteps = atoi( argv[1] );
    const int burnin = atoi( argv[2] );
    const double alpha = atof( argv[3] );
    const bool show_vecs = atoi( argv[4] );

    setprecision( 3 );

    cerr << "-----------" << endl;
    cerr << "Input parameters: " << endl;
    cerr << ">> nsteps: " << nsteps << endl;
    cerr << ">> burnin: " << burnin << endl;
    cerr << ">> alpha: " << alpha << endl;
    cerr << ">> show_vecs: " << show_vecs << endl;
    cerr << "-----------" << endl;

    VectorXf x(6);
    x << 0.5, 0.5, 0.5, 0.5, 0.5, 0.5;
    
    VectorXf xnew;
    int moves = 0;

    clock_t start = clock();

    cout << "# Metropolis-Hasitngs sample for CO2-Ar" << endl;
    cout << "# nsteps = " << nsteps << endl;
    cout << "# burnin = " << burnin << endl;
    cout << "# alpha = " << alpha << endl;
    cout << "# temperature = " << temperature << endl;
    cout << "# R (a. u.) = " << RDIST << endl;

    for ( int i = 0; i < nsteps + burnin; i++ )
    {
        xnew = metro_step( x, alpha );

        if ( xnew != x )
        {
            moves++;

            if ( i > burnin && show_vecs == true )
            {
                cout << x(0) << " " << x(1) << " " << x(2) << " " << x(3) << " " << x(4) << " " << x(5) << endl;
            }
        } 

        x = xnew;
    }

    cerr << "-----------------------------------" << endl;
    cerr << "Total steps: " << nsteps + burnin << "; moves: " << moves << "; percent: " << (double) moves / ( nsteps + burnin ) * 100 << "%" << endl;
    cerr << "Time elapsed: " << (double) ( clock() - start ) / CLOCKS_PER_SEC <<"s" << endl; 
    cerr << "-----------------------------------" << endl;

	return 0;
}

