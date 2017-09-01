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


// why static ? 
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

float pythag( float a, float b )
{
    float absa = fabs(a);
    float absb = fabs(b);

    if ( absa > absb )
        return absa * sqrt( 1.0 + SQR( absb / absa ));
    else
        return ( absb == 0.0 ? 0.0 : absb * sqrt( 1.0 + SQR( absa / absb )));
}

void tqli( double *d, double *e, int n )
{
   register int   m,l,iter,i,k;
   double         s,r,p,g,f,dd,c,b;

   for(i = 1; i < n; i++) e[i-1] = e[i];
     e[n] = 0.0;
   for(l = 0; l < n; l++) {
      iter = 0;
      do {
         for(m = l; m < n-1; m++) {
            dd = fabs(d[m]) + fabs(d[m+1]);
            if((double)(fabs(e[m])+dd) == dd) break;
         }
         if(m != l) {
            if(iter++ == 50) {
               printf("\n\nToo many iterations in tqli.\n");
               exit(1);
            }
            g = (d[l+1] - d[l])/(2.0 * e[l]);
            r = pythag( g, 1.0);
            g = d[m] - d[l] + e[l] / ( g + SIGN( r, g ) );
            s = c = 1.0;
            p = 0.0;
            for(i = m-1; i >= l; i--) {
               f      = s * e[i];
               b      = c*e[i];
               e[i+1] = (r=pythag(f,g));
               if(r == 0.0) {
                  d[i+1] -= p;
                  e[m]    = 0.0;
                  break;
               }
               s      = f/r;
               c      = g/r;
               g      = d[i+1] - p;
               r      = (d[i] - g) * s + 2.0 * c * b;
               d[i+1] = g + (p = s * r);
               g      = c * r - b;
            
            	/*  
               for(k = 0; k < n; k++) {
                  f         = z[k][i+1];
                  z[k][i+1] = s * z[k][i] + c * f;
                  z[k][i]   = c * z[k][i] - s * f;
               } /* end k-loop 
               */
            
            } /* end i-loop */
            if(r == 0.0 && i >= l) continue;
            d[l] -= p;
            e[l]  = g;
            e[m]  = 0.0;
         } /* end if-loop for m != 1 */
      } while(m != l);
   } /* end l-loop */
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
	cout << "Variational approach using delta-functions on constant grid. " << endl;
	cout << "X_MIN : " << X_MIN << "; X_MAX: " << "; DX: " << DX << endl << endl;

	cout << "SIZE of Hamiltonian matrix is " << SIZE << "x" << SIZE << endl;

	clock_t start;

	cout << "Time needed to fill Sparse hamiiltonian matrix: " << (double) (clock() - start) / CLOCKS_PER_SEC << "s" << endl;

	cout << "1: Eigen diagonalization procedure" << endl;
	cout << "2: Custom diagonalization procedure for tridiagonal matrices" << endl;
	cout << "3: Both ( test ) " << endl << endl;

	cout << "Enter the number depending on the diagonalization procedure: " << endl;
	cout << ">>>  ";
	int v;
	cin >> v;
	cout << endl << endl;

	// pre-filling potential array on the grid
	vector<double> pot(SIZE);
	for ( int i = 0; i < SIZE; i++ )
	{
		double x = X_MIN + i * DX;
		pot[i] = potential( x );
	}

	switch( v ) 
	{
		case 1:
			cout << "Filling Eigen SparseMatrix..." << endl;

			start = clock();
			
			SparseMatrix<double> h( SIZE, SIZE );
			fillSparseHamiltonian( h, pot );

			cout << "Time needed to fill Hamiltonian matrix: " << (double) (clock() - start) / CLOCKS_PER_SEC << "s" << endl;

			cout << "Starting diagonalization..." << endl;
			start = clock();

			SelfAdjointEigenSolver< SparseMatrix<double> > eigensolver( h );
			if ( eigensolver.info() != Success ) abort();

			cout << "Time needed to find eigenvalues: " << (double) (clock() - start) / CLOCKS_PER_SEC << "s" << endl;

		case 2:

			cout << "Filling arrays of diagonal and over-diagonal elements..." << endl;

		default:
			cout << "Invalid number!" << endl;
			exit( 1 );

	} 
	
	/*
	vector<double> eigs;
	for ( int i = 0; i < eigensolver.eigenvalues().size(); i++ )
	{
		eigs.push_back( eigensolver.eigenvalues()[i] );
	}

	for ( int i = 0; i < 5; i++ )
	{
			cout << i << ": " << eigs[i] << endl;
	}
	*/

	return 0;
}
