#include <iostream>

#include <math.h>
#include <ctime>

using namespace std;

double integrand ( double x )
{
	return exp( - x * x );
}

int main()
{
	cout << "Using Gregory formulae integrating exp(-x^2) on the interval [ 0, 10 ]" << endl;

	double lb = 0;
	double rb = 10;

	int KNOTS = 100;

	double step = ( rb - lb ) / KNOTS;

	double x = lb;
	double* vals = new double [KNOTS + 1];
	for ( int counter = 0; counter < KNOTS + 1; counter++ )
	{	
		vals[counter] = integrand( x );
		x += step;
	}

	double* diffs = new double [KNOTS];
	for (int counter = 0; counter < KNOTS; counter++ ) 
	{
		diffs[counter] = vals[counter + 1] - vals[counter];
	}
	
	double integ = 0;
	for ( int counter = 0; counter < KNOTS + 1; counter++ )
	{
		if ( counter == 0 || counter == KNOTS )
			integ += 0.5 * vals[counter];
		else
			integ += vals[counter];
	}

	double perturb1 = 1.0 / 12.0 * ( diffs[0] - diffs[KNOTS] );
	double perturb2 = - 1.0 / 24.0 * ( pow(diffs[0], 2) + pow(diffs[KNOTS - 1], 2) );
	double perturb3 = 19.0 / 720.0 * ( pow(diffs[0], 3) - pow(diffs[KNOTS - 2], 3) );
	double perturb4 = - 3.0 / 160.0 * ( pow(diffs[0], 4) + pow(diffs[KNOTS - 3], 4) );
	
	cout << "perturb1: " << perturb1 << endl;
	cout << "perturb2: " << perturb2 << endl;
	cout << "perturb3: " << perturb3 << endl;
	cout << "perturb4: " << perturb4 << endl;

	integ = integ + perturb1 + perturb2 + perturb3 + perturb4;

	cout << "integ: " << integ << endl;



	return 0;
}
