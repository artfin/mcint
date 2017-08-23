#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include <ctime>
#include <math.h>

#include <vector>

// Warn about use of deprecated functions
#define GNUPLOT_DEPRECATE_WARN
#include "gnuplot-iostream.h"

using namespace std;

const double EPSILON = 10.0;
const double SIGMA = 1.0;

const double DX = 0.2;
const double X_MIN = 0.5;
const double X_MAX = 10.0;
const int SIZE = ( X_MAX - X_MIN ) / DX; 

const double RANDSTEP = 0.01;
const double STEPS = 1e6;

double potential( double x )
{
	double r6 = pow(x / SIGMA, -6);
	double r12 = pow(r6, 2);

	return 4 * EPSILON * (r12 - r6);
	//return x * x;
}

double calc_energy( vector<double> psi )
{
	double energy = 0;
	double sum = 0;

	for ( int i = 0; i < psi.size(); i++ )
	{
		double x = X_MIN + i * DX;

		if ( i == 0 || i == psi.size() )
		{
			energy += 0.5 * DX * potential( x ) * pow(psi[i], 2);
			sum += 0.5 * DX * pow(psi[i], 2);
		}
		else
		{
			energy += DX * potential( x ) * pow(psi[i], 2) - 0.5 / DX * psi[i] * (psi[i + 1] + psi[i - 1] - 2 * psi[i]);
			sum += psi[i] * psi[i] * DX;
		}
	}

	return energy / sum;
}

void normalize( vector<double> &psi )
{
	double sum = 0;

	for ( int i = 0; i < psi.size(); i++ )
	{
		sum += DX * psi[i] * psi[i];
	}

	for ( int i = 0; i < psi.size(); i++ )
	{
		psi[i] /= sqrt( sum );
	}
}

void print_vector( vector<double> psi )
{
	for ( int counter = 0; counter < psi.size(); counter++ )
	{
		cout << "psi[" << counter << "]: " << psi[counter] << endl;
	}
}

int main()
{
	srand( 1 );

	vector<double> psi;

	for ( int counter = 0; counter < SIZE; counter++ )
	{
		psi.push_back(1);
	}

	double initial_energy = calc_energy( psi ); 	
	cout << "Initial energy: " << initial_energy << endl;

	int moves = 0;

	clock_t start = clock();

	Gnuplot gp;
	gp << "set xrange [" << X_MIN << ":" << X_MAX << "]" << endl;
	gp << "set yrange [" << -0.2 << ":" << 1.2 << "]" << endl;
	
	for ( int i = 0; i < STEPS; i++ )
	{
		cout << "i: " << i << endl;

		int pos = rand() % SIZE;
		double step = ( rand() % 1000 / 1000.0 - 0.5 ) * RANDSTEP;

		psi[pos] += step;

		double energy = calc_energy( psi );

		// if energy is lowering then keeping the perturbation 
		if ( energy < initial_energy )
		{
			moves++;
			initial_energy = energy;
		}
		// if energy is not lowering then reversing the perturbation
		else
		{
			psi[pos] -= step;
		}

		normalize( psi );

		int counter = 0;
		vector< pair<double,double> > xy_pts;
		for ( double x = X_MIN; x <= X_MAX; x += DX )
		{
			xy_pts.push_back( make_pair( x, psi[counter] ) );
			counter++;
		}
		
		//gp << "plot '-' binary" << gp.binFmt1d(psi, "array") << "with lines notitle\n";
		//gp.sendBinary1d(pts);
		//gp.flush();

		//auto Y_MIN = min_element( psi.begin(), psi.end() );
		//auto Y_MAX = max_element( psi.begin(), psi.end() );
	
		gp << "plot '-' with lines title 'psi'" << endl;
		gp.send1d(xy_pts);

		gp.flush();
		
		if ( i == 0 )
		{
			sleep( 3 );
		}
		
		//sleep( 100 );
	}




	cout << "------------------" << endl << endl;

	cout << "Final energy: " << initial_energy << endl;
	cout << "Moves made: " << moves << endl;

	cout << "Time elapsed: " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "s" << endl;

	return 0;
}
