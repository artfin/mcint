#include <iostream>

#include <math.h>
#include <ctime>

#include <functional>

using namespace std;
using namespace std::placeholders;

double gammainc_integrand( double x, double a )
{
	return pow(x, a - 1) * exp( -x );
}

double simps_primitive( decltype( bind(gammainc_integrand, _1, 1.0) ) f, double x0, double x1 )
{
	double h = x1 - x0;
	
	double f0 =  f( x0 );
	double f1 =  f( x1 );
	double f12 = 0.5 * ( f0 + f1 );

	return h / 6 * ( f0 + 4 * f12 + f1 );
}

double simps( decltype( bind(gammainc_integrand, _1, 1.0) ) f, double x0, double x1, double step )
{
	double integ = 0;

	double lb = x0;
	double rb = x0 + step;

	while ( rb < x1 ) 
	{
		integ += simps_primitive( f, lb, rb );

		lb += step;
		rb += step;
	}

	return integ;
}

double gammainc ( double a, double b )
{
	auto gi = bind( gammainc_integrand, _1, a);

	int n = 6;
	double* steps = new double[n];
	steps[0] = 1e-1;
	steps[1] = 1e-2;
	steps[2] = 1e-3;
	steps[3] = 1e-4;
	steps[4] = 1e-5;	
	steps[5] = 1e-6;

	double r1, r2;

	for ( int counter = 0; counter < n; counter++ )
	{
		clock_t cycle_time = clock();

		double step = steps[counter];

		r1 = simps( gi, 0.0, b, b * step );
		r2 = simps( gi, 0.0, 20.0, step );

		double res = r1 / r2;
		double rel_error = (1 - res / 0.626155) * 100;

		cout << "step: " << step << "; result: " << res << "; rel error: " << rel_error << "% ; time elapsed: " << ( clock() - cycle_time ) / (double) CLOCKS_PER_SEC << "s" << endl;
	}

	return r1 / r2;
}

int main()
{
	/*
	double VAR = 2.5;
	auto gi = bind(gammainc_integrand, _1, VAR); 
	
	double integ = simps( gi, 0.0, 50.0, 0.01 );
	
	double g = 1.32934038818;

	cout << "integral: " << integ << endl;
	cout << "Time elapsed: " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "s" << endl;

	double abs_error = integ - g;
	double rel_error = abs_error / g;

	cout << "abs_error: " << abs_error << endl;
	cout << "rel_error: " << rel_error << endl;

	cout << endl << endl;
	*/

	double x = gammainc( 1.1, 1.1 );
	cout << "gammainc(1.1, 1.1): " << x << endl;

	return 0;
}
