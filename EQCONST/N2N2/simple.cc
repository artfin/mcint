#include <iostream>

#include "hep/mc.hpp"

#include <math.h>
#include <cstddef>
#include <vector>
#include <ctime>

#include <boost/math/special_functions/gamma.hpp>

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
// pascals to atmospheres
const double PATOATM = 101325.0;
// m3 to cm3
const double M3TOCM3 = pow(10, 6);
// universal gas consant
const double UGASCONST = 8.314;

struct stop_after_precision
{
	stop_after_precision( double rel_error )
			: rel_error(rel_error)
	{
	}

	bool operator()(std::vector<hep::vegas_result<double>> const& r)
	{
		hep::vegas_verbose_callback<double>(r);

		//compute cumulative result
		auto const result = hep::cumulative_result0(r.begin(), r.end());

		// check the relative error
		if ( result.error() < rel_error * result.value() )
		{

			// returning false stops all remaining iterations
			return false;
		}

		return true;
	}

	double rel_error;
};


double integrand_(hep::mc_point<double> const& x, double Temperature)
{
	double R_new = x.point()[0];
	double Theta1_new = x.point()[1];
	double Theta2_new = x.point()[2];
	double Phi_new = x.point()[3];

	double R = tan( M_PI / 2 * R_new );
	double Theta1 = M_PI * Theta1_new;
	double Theta2 = M_PI * Theta2_new;
	double Phi = 2 * M_PI * Phi_new;

	if ( R < 4.4 )
	{
		return 0;
	}

	double potential_value;
	potn2n2(&R, &Theta1, &Theta2, &Phi, &potential_value);

	potential_value = potential_value * CMTOH * HTOJ;

	if ( potential_value < 0 )
	{
		//jacobians
		double jacR = M_PI / 2 * (1 + pow(R, 2));
		double jacTheta1 = M_PI;
		double jacTheta2 = M_PI;
		double jacPhi = 2 * M_PI;

		double U_KT = potential_value / (BOLTZCONST * Temperature);

		return jacR * jacTheta1 * jacTheta2 * jacPhi * pow(R, 2) * sin(Theta1) * sin(Theta2) * boost::math::gamma_p(3.5, - U_KT ) * exp( -U_KT );
	} 

	else 
	{
		return 0;
	}
}

int main()
{
	potinit();

	cout << "--- Computing Simple EQCONST for N2N2 --- " << endl;

	cout << ">> Enter the lower boundary for temperature interval: " << endl;
	double LTEMP;
	cin >> LTEMP;

	cout << ">> Enter the higher boundary for temperature interval: " << endl;
	double HTEMP;
	cin >> HTEMP;

	cout << ">> Enter the step for temperature: " << endl;
	double STEP;
	cin >> STEP;

	vector<double> temperatures;
	vector<double> constants;

	clock_t full_clock = clock();
	clock_t cycle_clock;
	
	// set the verbose vegas callback function
	// hep::vegas_callback<double>(hep::vegas_verbose_callback<double>);
	hep::vegas_callback<double>(stop_after_precision(0.01));
	
	for ( double TEMP = LTEMP; TEMP <= HTEMP; TEMP += STEP ) 
	{
		cycle_clock = clock();
		
		temperatures.push_back(TEMP);

		auto integrand = bind( integrand_, _1, TEMP );

		auto results = hep::vegas(
			hep::make_integrand<double>(integrand, 4),
			std::vector<std::size_t>(10, 10000)
		);
					
		auto result = hep::cumulative_result0(results.begin() + 1, results.end());
		double chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());

		// cout << "Cumulative result (excluding the first iteration): \nN = " << result.calls() << " I = " << result.value() << " +- " << result.error() << " chi^2/dof = " << chi_square_dof << endl;

		double res = AVOGADRO / ( 4 * UGASCONST * TEMP ) * result.value() * pow(ALU, 3) * PATOATM;

		constants.push_back( res );
					
		cout << endl << "Temperature: " << TEMP << "; EQCONST: " << res << "; time elapsed (per cycle): " << (clock() - cycle_clock ) / (double) CLOCKS_PER_SEC << "s " << endl;
	}

	cout << "Time elapsed (full calculation): " << ( clock() - full_clock ) / (double) CLOCKS_PER_SEC << "s" << endl;

	FILE* const_file = fopen("simple_constants.dat", "w");

	for ( int counter = 0; counter < temperatures.size(); counter++ )
	{
		fprintf(const_file, "%.12lf %.12lf \n", temperatures[counter], constants[counter]);
	}

	fclose(const_file);

	return 0;
}

