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
// current temperature
// const double Temperature = 200;
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
		hep::vegas_verbose_callback<double>(r);

		//compute cumulative result
		auto const result = hep::cumulative_result0(r.begin(), r.end());
		
		/*
		cout << "--------" << endl;
		cout << "current relative error: " << result.error() / result.value() << endl;
		cout << "result.error(): " << result.error() << endl;
		cout << "rel_error * result.value(): " << rel_error * result.value() << endl;
		cout << "--------" << endl;
		*/

		// check the relative error
		if ( result.error() < abs(rel_error * result.value()) )
		{
			cout << ">> relative error "
				 << ( result.error() / result.value() )
				 << " is smaller than the limit " << rel_error << endl;

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

	double R = tan(M_PI / 2 * R_new);
	double Theta1 = M_PI * Theta1_new;
	double Theta2 = M_PI * Theta2_new;
	double Phi = 2 * M_PI * Phi_new;
	
	if ( R > 4.4 )
	{
			// transforming potential value from CM-1 to J
			double potential_value;
			potn2n2(&R, &Theta1, &Theta2, &Phi, &potential_value);
			potential_value = potential_value * CMTOH * HTOJ;
			
			double in_braces = ( 1 - exp ( - potential_value / (BOLTZCONST * Temperature) ));
			
			/*
			cout << "---" << endl;
			cout << "R: " << R << endl;
			cout << "Theta1: " << Theta1 << endl;
			cout << "Theta2: " << Theta2 << endl;
			cout << "Phi: " << Phi << endl;

			cout << "potential value: " << potential_value << endl;
			cout << "exp (- potential_value / k T ): " << exp ( - potential_value / (BOLTZCONST * Temperature )) << endl;
			cout << "in braces: " << in_braces << endl;
			*/

			return pow(M_PI, 4) / 4 * in_braces * pow(R, 2) * sin(Theta1) * sin(Theta2) * (1 + pow(R, 2));
	
	} else {
		return pow(M_PI, 4) / 4 * pow(R, 2) * sin(Theta1) * sin(Theta2) * (1 + pow(R, 2));
	}
	
}

int main()
{
	double r = 5.0;
	double theta1 = 0.5;
	double theta2 = 0.5;
	double phi = 0.5;
	double v;
	
	potinit();
	potn2n2(&r, &theta1, &theta2, &phi, &v);
	cout << "v: " << v << endl;

	cout << "--- Computing SVC for N2N2 --- " << endl;
	
	cout << ">> Enter the lower boundary for temperature interval:" << endl;
	double LTEMP;
	cin >>LTEMP;

	cout << ">> Enter the higher boundary for temperature interval:" << endl;
	double HTEMP;
	cin >> HTEMP;

	cout << ">> Enter the step for temperature:" << endl;
	double STEP;
	cin >> STEP;

	int N = (int) (HTEMP - LTEMP) / STEP + 1;
	cout << ">> Number of values to be calculated is : " << N << endl;

	cout << "----------------------------------" << endl;

	vector<double> temperatures;
	vector<double> vir_coeffs;

	clock_t full_clock;
	clock_t cycle_clock;

	full_clock = clock();
	
	// initializing internal parameters for calculating potential value
	potinit();

	for ( double temp = LTEMP; temp <= HTEMP; temp += STEP )
	{
		cycle_clock = clock();

		temperatures.push_back(temp);
			
		// creating integrand function
		auto integrand = bind(integrand_, _1, temp);
			
		hep::vegas_callback<double>(stop_after_precision(0.005));

		auto results = hep::vegas(
			hep::make_integrand<double>(integrand, 4),
			std::vector<std::size_t>(10, 15000)
		);
			
		auto result = hep::cumulative_result0(results.begin() + 1, results.end());
		double chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());

		cout << "Cumulative result (excluding the first iteration): \nN = " << result.calls() << " I = " << result.value() << " +- " << result.error() << " chi^2/dof = " << chi_square_dof << endl;

		double res = AVOGADRO * result.value() * pow(ALU, 3) * M3TOCM3;
			
		vir_coeffs.push_back(res);
		
		cout << "Temperature: " << temp << "; SVC: " << res << "; time elapsed: " << (clock() - cycle_clock ) / (double) CLOCKS_PER_SEC << endl << endl;

	}

	cout << "-----------------------" << endl;
	cout << "Total elapsed time: " << ( clock() - full_clock ) / (double) CLOCKS_PER_SEC << endl;

	FILE* svc_file = fopen("svc.dat", "a");

	for ( int counter = 0; counter < temperatures.size(); counter++ )
	{
		fprintf(svc_file, "%.2lf %.8lf \n", temperatures[counter], vir_coeffs[counter]);
	}

	fclose(svc_file);

	return 0;
}
