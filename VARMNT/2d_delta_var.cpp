#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include <ctime>
#include <math.h>

#include <algorithm>
#include <iomanip>
#include <vector>

#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

const double MIN = -10.0;
const double MAX =  10.0;

double potential( double x, double y )
{
	return 0.5 * x * x + 0.5 * y * y;
}

void fillSparseHamiltonian( SparseMatrix<double> &h, vector<double> pot, int SIZE )
{
	double DX = ( MAX - MIN ) / ( SIZE - 1 );
	double DX2 = DX * DX;
	
	int SIZE2 = pow( SIZE, 2 );
	int SIZE4 = pow( SIZE2, 2 );

	vector< Triplet<double> > tripletList( SIZE4 );
	int counter = 0;

	for ( int i = 0; i < SIZE2; i++ )
	{
		for ( int j = 0; j < SIZE2; j++ )
		{
			if ( i == j )
			{
				double el = 2 / DX2 + pot[counter];
				counter++;
				tripletList.push_back( Triplet<double>( i, j, el )); 
			}
			
			if ( i == j + 1 || i == j - 1 || i == j + SIZE || i == j - SIZE )
			{
					//cout << "Non-diagonal el: i = " << i << "; j = " << j << endl;
				double el = - 0.5 / DX2;
				tripletList.push_back( Triplet<double>( i, j, el ));
			}
		}
	}

	h.setFromTriplets( tripletList.begin(), tripletList.end() );
}

int main()
{
	cout << "Variational approach using delta-functions on constant 2D-grid." << endl;

	cout << "X_MIN = Y_MIN: " << MIN << "; X_MAX = Y_MAX = " << MAX << endl;
	cout << "Enter the number of delta-functions on ONE axis to be used:  " << endl;
	int SIZE;
	cin >> SIZE;
	
	int SIZE2 = pow( SIZE, 2 );
	int SIZE4 = pow( SIZE2, 2 );

	cout << "TOTAL NUMBER OF DELTA FUNCTIONS TO BE USED: " << SIZE2 << endl;

	cout << "SIZE of Hamiltonian matrix is " << SIZE2 << "x" << SIZE2 << endl;

	clock_t start;

	double DX = ( MAX - MIN ) / ( SIZE - 1 );

	cout << "Filling potential matrix..." << endl;
	vector<double> pot;
	for ( int i = 0; i < SIZE; i++ )
	{
		for ( int j = 0; j < SIZE; j++ )
		{
			double x = MIN + i * DX;
			double y = MIN + j * DX;

			pot.push_back( potential( x, y ));
			//cout << "x = " << x << "; y = " << y << "; pot: " << potential(x, y) << endl;
		}
	}
	cout << "Potential vector is filled. " << endl;
	
	cout << "Filling Hamiltonian matrix..." << endl;

	start = clock();
	SparseMatrix<double> h( SIZE2, SIZE2 );
	fillSparseHamiltonian( h, pot, SIZE );

	cout << h << endl;

	cout << "Time needed to fill Hamiltonian matrix: " << (double) (clock() - start ) / CLOCKS_PER_SEC << "s" << endl;
	
	cout << "Starting diagonalization..." << endl;
	start = clock();

	SelfAdjointEigenSolver< SparseMatrix<double> > eigensolver( h );
	if ( eigensolver.info() != Success ) abort();

	cout << "Time needed to find eigenvalues: " << (double) (clock() - start) / CLOCKS_PER_SEC << "s" << endl;

	cout << "Number of eigenvalues to be displayed: "<< endl;
	
	int n_values;
	cin >> n_values;

	cout << "Number \t Eigenvalue" << endl;
	cout << "-----------------------------------" << endl;
	for ( int i = 0; i < n_values; i++ )
	{
		cout << i << setw(10) << eigensolver.eigenvalues()[i] << endl;
	}

	return 0;
}
