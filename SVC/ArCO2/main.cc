#include "ab_initio_potential.h"
#include "hep/mc.hpp"

#include <math.h> 
#include <cstddef>
#include <iostream>
#include <vector>

#include <ctime>

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

// transforming variables:
// R -> 2/pi arctg(R)
// theta -> 1/pi theta
double integrand_(hep::mc_point<double> const& x, double Temperature)
{
	double R_new = x.point()[0];
	double Theta_new = x.point()[1];

	double R = tan(M_PI / 2 * R_new);
	double Theta = M_PI * Theta_new;

	double potential_value = ab_initio_pot(R, Theta) * CMTOH * HTOJ;
	
	double in_braces = ( 1 - exp( - potential_value / (BOLTZCONST * Temperature) ));	

	cout << "---" << endl;
	cout << "R: " << R << endl;
	cout << "Theta: " << Theta << endl;
	cout << "potential value: " << potential_value << endl;
	cout << "exp(-u/kt): " << exp( - potential_value / (BOLTZCONST * Temperature) ) << endl;
	cout << "---" << endl << endl;
	
	return pow(M_PI, 2) / 2 * in_braces * pow(R, 2) * sin(Theta) * (1 + pow(R, 2)); 
}


int main()
{
	cout << "--- Computing second virial coefficient (SVC) for AB INITIO potential Ar-CO2 (Y. Kalugina) --- " << endl;
	cout << "-------------------------------------------------" << endl;

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

	clock_t full_start;
	full_start = clock();

	double res;
	clock_t start;
	
	for ( int temp = LTEMP; temp <= HTEMP; temp += STEP )
	{
		temperatures.push_back(temp);

		// starting clock
     	start = clock();
     	
		// it is called currying
		// (partial function application)
		auto integrand = bind(integrand_, _1, temp);

		// set the verbose vegas callback function
		// hep::vegas_callback<double>(hep::vegas_verbose_callback<double>);

		// this function also calls vegas_verbose_callback after each iteration which in turn prints the individual iterations
		auto results = hep::vegas(
			hep::make_integrand<double>(integrand, 2),
			std::vector<std::size_t>(15, 50000)
		);

		// results contains the estimations for each iteration. We could take the
		// result from last iteration, but here we instead choose to combine the
		// results of all iterations but the first one in a cumulative
		auto result = hep::cumulative_result0(results.begin() + 1, results.end());
		double chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());

		// print the cumulative result
		cout << ">> Cumulative result (excluding the first iteration): \n>> N = " << result.calls() << " I = " << result.value() << " +- " << result.error() << " chi^2/dof = " << chi_square_dof << endl;

		res = M_PI * AVOGADRO * result.value() * pow(ALU, 3) * M3TOCM3;
		vir_coeffs.push_back(res);

		cout << ">> Temperature: " << temperatures.back() <<"; SVC: " << vir_coeffs.back() << "; time needed: " << (clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << endl;
		cout << "----------------------------------" << endl;
	}

	cout << "Total time elapsed: " << (clock() - full_start) / (double)(CLOCKS_PER_SEC ) << " s" << endl;

	FILE *svc_file = fopen("data.txt", "w");
	for ( int i = 0; i < temperatures.size(); i++ ) 
	{
		fprintf(svc_file, "%.4lf %.4lf\n", temperatures[i], vir_coeffs[i]);
	}

	fclose(svc_file);
	
	return 0;
}


