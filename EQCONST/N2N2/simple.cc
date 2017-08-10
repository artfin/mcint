#include <iostream>

#include "hep/mc.hpp"

#include <math.h>
#include <cstddef>
#include <vector>

#include <ctime>

extern "C" void potinit(void);
extern "C" void potn2n2(double* rr, double* theta1, double* theta2, double* phi, double* res);

using namespace std;
using namespace std::placeholders;

// cm^-1 to hartree
const double CMTOH = 1 / 2.1947 * pow(10, -5);
// hartree to joules
const double HTOJ = 4.35974417 * pow(10, -18);
// boltzmann constant
const long double BOLTZCONST = 1.38064852 * pow(10, -23);
// avogadro number
const double AVOGADRO  = 6.022 * pow(10, 23); 
// atomic length unit
const double ALU = 5.291772 * pow(10, -11);
// m3 to cm3
const double M3TOCM3 = pow(10, 6);

struct stop_after_precision
{
	stop_after_precision( double rel_error )
			: rel_error(rel_error)
	{
	}

	bool operator()(std::vector<hep::vegas_result<double>> const& r)
	{
		// hep::vegas_verbose_callback<double>(r);

		//compute cumulative result
		auto const result = hep::cumulative_result0(r.begin(), r.end());

		// check the relative error
		if ( result.error() < rel_error * result.value() )
		{
		//	cout << ">> relative error "
		//		 << ( result.error() / result.value() )
		//		 << " is smaller than the limit " << rel_error << endl;

			// returning false stops all remaining iterations
			return false;
		}

		return true;
	}

	double rel_error;
};


double gammainc_integrand(hep::mc_point<double> const& x, double a, double b)
{
	return pow(b, a) * pow(x.point()[0], a - 1) * exp(- b * x.point()[0]);
}

double gammainc(double a, double b)
{
		auto integrand = bind(gammainc_integrand, _1, a, b);
		
		// stop if error is lower than 0.05%
		hep::vegas_callback<double>(stop_after_precision(0.0005));
		
		auto results = hep::vegas(
			hep::make_integrand<double>(integrand, 1),
			vector<size_t>(20, 1000)
		);

		auto result = hep::cumulative_result0(results.begin() + 1, results.end());

		return result.value();
}

int main()
{
	clock_t start = clock();

	double x = gammainc(1.0, 0.1);
	cout << "x: " << x << endl;

	cout << "time elapsed: " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "s" << endl;
	
	return 0;
}

/*
double integrand(hep::mc_point<double> &x)
{
	double R_new = x.point()[0];
	double Theta1_new = x.point()[1];
	double Theta2_new = x.point()[2];
	double Phi_new = x.point()[3];

	double R = tan(M_PI / 2 * x.point()[0]);
	double Theta1 = M_PI * Theta1_new;
	double Theta2 = M_PI * Theta2_new;
	double Phi = 2 * M_PI * Phi_new;

	if ( R > 4.4 && R < 45.0 )
	{

	}
}
*/
