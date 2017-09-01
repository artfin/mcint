#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include <ctime>
#include <math.h>

#include <vector>

#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

const double DX = 0.5;
const double DX2 = pow(DX, 2);
const double X_MIN = -1.0;
const double X_MAX =  1.0;

const int SIZE = (int) (X_MAX - X_MIN) / DX + 1;

const double EPSILON = 10.0;
const double SIGMA = 1.0;

const double DE = 100.0;
const double A = 1.0;
const double RE = 1.0;

vector<double> calculate_levels ( int n )
{
	vector<double> levels;

	double NU0 = A / (2 * M_PI) * sqrt ( 2 * DE );

	for ( int i = 0; i < n; i++ )
	{
		double interim = 2 * M_PI * NU0 * ( i + 0.5 );
		levels.push_back( interim - pow(interim, 2) / (4 * DE) );
	}

	return levels;
}	

double potential( double x )
{
		//double r6 = pow(x / SIGMA, -6);
		//double r12 = pow(r6, 2);

		//return 4 * EPSILON * (r12 - r6);
	
	return 0.5 * x * x;
	
	//return DE * pow(1 - exp( - A * (x - RE)), 2);
}

void fillSparseHamiltonian( SparseMatrix<double> &h, vector<double> pot )
{
	vector< Triplet<double> > tripletList( SIZE * SIZE );

	for ( int i = 0; i < SIZE; i++ )
	{
		for ( int j = 0; j < SIZE; j++ )
		{
			if ( i == j )
			{
				double el = 1 / DX2 + pot[i];
				tripletList.push_back( Triplet<double>(i, j, el) );
			}
			if (  i == j + 1 || i == j - 1 )
			{
				double el = - 0.5 / DX2;
				tripletList.push_back( Triplet<double>(i, j, el) );
			}
		}
	}

	h.setFromTriplets( tripletList.begin(), tripletList.end() );
}

int main()
{
		//vector<double> levels = calculate_levels(20);
		//for (int i = 0; i < levels.size(); i++ )
		//{
		//cout << "Level " << i << ": " << levels[i] << endl;
		//}
	
	
	cout << "SIZE of Hamiltonian matrix is " << SIZE << "x" << SIZE << endl;

	vector<double> pot(SIZE);
	for ( int i = 0; i < SIZE; i++ )
	{
		double x = X_MIN + i * DX;
		pot[i] = potential( x );
	}
	
	clock_t start = clock();

	SparseMatrix<double> h( SIZE, SIZE );
	fillSparseHamiltonian( h, pot );
	
	cout << h << endl;

	cout << "Time needed to fill Sparse hamiiltonian matrix: " << (double) (clock() - start) / CLOCKS_PER_SEC << "s" << endl;

	start = clock();

	SelfAdjointEigenSolver< SparseMatrix<double> > eigensolver( h );
	if ( eigensolver.info() != Success ) abort();

	cout << "The eigenvalues of H are: " << endl;
	vector<double> eigs;
	for ( int i = 0; i < eigensolver.eigenvalues().size(); i++ )
	{
		eigs.push_back( eigensolver.eigenvalues()[i] );
	}

	for ( int i = 0; i < 5; i++ )
	{
			cout << i << ": " << eigs[i] << endl;
	}

	cout << "Time needed to find eigenvalues: " << (double) (clock() - start) / CLOCKS_PER_SEC << "s" << endl;

	return 0;
}
