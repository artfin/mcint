#include "hep/mc.hpp"
#include "hamiltonian.h"

#include <math.h>
#include <cstddef>
#include <iostream>
#include <vector>

#include<ctime>

using namespace std;
using namespace std::placeholders;

// hartree to joules
const double HTOJ = 4.35974417 * pow(10, -18);
// boltzmann constant
const long double BOLTZCONST = 1.38064852 * pow(10, -23);
// avogadro number
const double AVOGADRO  = 6.022 * pow(10, 23); 
// planck constant
const long double PLANKCONST = 6.626070040 * pow(10, -34);
const long double PLANKCONST2 = pow(PLANKCONST, 2);

// atomic units to CI:
// -- length
// -- mass
const double ALU = 5.291772 * pow(10, -11);
const long double AMU = 9.1093826 * pow(10, -31);

// dalton to kg
const long double DA = 1.660539040 * pow(10, -27);
// R 
const double UNIGASCONST = 8.314;
// pascals to atmospheres
const double PATOATM = 101325.0;

// reduced masses
const double MU1 = 14583.109;
const double MU2 = 38193.858;
// 
const double CO2_LENGTH = 4.398;

// particles molar masses: g/mol
const double AR_MOLARMASS = 39.948;
const double CO2_MOLARMASS =  44.01; 
const double COMPLEX_MOLARMASS = AR_MOLARMASS + CO2_MOLARMASS;

double integrand_(hep::mc_point<double> const& x, double Temperature)
{
	double R_new = x.point()[0];
	double Theta_new = x.point()[1];
	double pR_new = x.point()[2];
	double pT_new = x.point()[3];
	double Jx_new = x.point()[4];
	double Jy_new = x.point()[5];
	double Jz_new = x.point()[6];

	double R = tan(M_PI / 2 * R_new);
	double Theta = M_PI * Theta_new;
	double pR = tan(M_PI * (pR_new - 0.5));
	double pT = tan(M_PI * (pT_new - 0.5));
	double Jx = tan(M_PI * (Jx_new - 0.5));
	double Jy = tan(M_PI * (Jy_new - 0.5));
	double Jz = tan(M_PI * (Jz_new - 0.5));

	double h = hamiltonian(R, Theta, pR, pT, Jx, Jy, Jz);

	if ( h < 0 )
	{
		double boltz_exp = exp( -h * HTOJ / ( BOLTZCONST * Temperature ));
		
		// jacobians
		double jacR = M_PI / 2 * (1 + pow(R, 2));
		double jacTheta = M_PI;
		double jacpR = M_PI * (1 + pow(pR, 2));
		double jacpT = M_PI * (1 + pow(pT, 2));
		double jacJx = M_PI * (1 + pow(Jx, 2));
		double jacJy = M_PI * (1 + pow(Jy, 2));
		double jacJz = M_PI * (1 + pow(Jz, 2));

		return boltz_exp * jacR * jacTheta * jacpR * jacpT * jacJx * jacJy * jacJz;
	}
	else
	{
		return 0;
	}
}

double calculate_qtr( double mass, double Temperature )
{
	return pow(2 * M_PI * mass * DA * BOLTZCONST * Temperature / PLANKCONST2, 1.5);
}	

int main()
{
	cout << "--- Computing equilibrium constant (EQCONST) for AR-CO2 ---" << endl;
	cout << "-----------------------------------------------------------" << endl;

	cout << ">> Enter the lower boundary for temperature interval: " << endl;
	double LTEMP;
	cin >> LTEMP;

	cout << ">> Enter the higher boundary for temperature interval: " << endl;
	double HTEMP;
	cin >> HTEMP;

	cout << ">> Enter the step for temperature: " << endl;
	double STEP;
	cin >> STEP;

	int N = (int) (HTEMP - LTEMP) / STEP + 1;
	cout << ">> Number of values to be calculated: " << N << endl;

	cout << ">> Switch for callback function: " << endl;
	cout << ">> --- Hit 0 to switch off." << endl;
	cout << ">> --- Hit 1 to switch on." << endl;

	int CALLSWITCH;
	cin >> CALLSWITCH;

	cout << "---------------------------------------" << endl;

	vector<double> temperatures;
	vector<double> eqconsts;

	clock_t full_start = clock();

	double res;
	clock_t start;

	double Qtr_complex;
	double Q_Ar;
	double Qtr_CO2;
	double Qrot_CO2;
    double Q_CO2;
	
	for ( int temp = LTEMP; temp <= HTEMP; temp += STEP )
	{
		temperatures.push_back(temp);

		// starting clock
		start = clock();

		// creating integrand for current temperature
		auto integrand = bind(integrand_, _1, temp);

		if ( CALLSWITCH == 1 )
		{
			// set the verbose vegas callback function
			hep::vegas_callback<double>(hep::vegas_verbose_callback<double>);
		}
		auto results = hep::vegas(
			hep::make_integrand<double>(integrand, 7),
			std::vector<std::size_t>(10, pow(10, 5))
		);

		auto result = hep::cumulative_result0(results.begin() + 1, results.end());
		double chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());

		cout << ">> Cumulative result: " << endl;
		cout << ">> N = " << result.calls() << "; I = " << result.value() << " +- " << result.error() << endl;
		cout << ">> chi^2/dof = " << chi_square_dof << endl;

		cout << "Time needed: " << (clock() - start) / (double)(CLOCKS_PER_SEC) << "s" << endl;
	
		Qtr_complex = calculate_qtr(COMPLEX_MOLARMASS, temp);
		Q_Ar = calculate_qtr(AR_MOLARMASS, temp); 
		Qtr_CO2 = calculate_qtr(CO2_MOLARMASS, temp);
		Qrot_CO2 = 4 * pow(M_PI, 2) * BOLTZCONST * temp / PLANKCONST2 * MU1 * AMU * pow(CO2_LENGTH * ALU, 2);
    	Q_CO2 = Qtr_CO2 * Qrot_CO2;
	
	// cout << "Qtr complex: " << Qtr_complex << endl;
	// cout << "Q_Ar: " << Q_Ar << endl;
	// cout << "Q_CO2: " << Q_CO2 << endl;

		res = AVOGADRO / (UNIGASCONST * temp) * Qtr_complex / Q_Ar / Q_CO2 * result.value() * PATOATM * 8 * pow(M_PI, 2) / pow(2 * M_PI, 5) / 2;
		eqconsts.push_back(res);

    	cout << "T = " << temp << "; EQCONST: " << res << endl;
		cout << "-------------------------------" << endl;	
	}

	cout << "Total time elapsed: " << (clock() - full_start) / (double)(CLOCKS_PER_SEC) << "s" << endl;

	FILE *output = fopen("data.txt", "w");
	for ( int i = 0; i < temperatures.size(); i++ )
	{
		fprintf(output, "%.2lf %.8lf\n", temperatures[i], eqconsts[i]);
	}
	fclose(output);

	return 0;	
}
