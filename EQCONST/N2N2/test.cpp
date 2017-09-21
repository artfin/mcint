#include "hep/mc.hpp"

#include <cstddef>
#include <iostream>
#include <vector>

using namespace std;

double square( hep::mc_point<double> const& x )
{
	return 3.0 * x.point()[0] * x.point()[0];
}

int main()
{
	cout << "computing integral of 3*x**2 from 0 to 1 which 1" << endl;

	hep::vegas_callback<double>( hep::vegas_verbose_callback<double>);

	auto results = hep::vegas(
		hep::make_integrand<double>(square, 1),
		vector<size_t>( 5, 1000 )
	);

	auto result = hep::cumulative_result0( results.begin() + 1, results.end() );

	double chi_square_dof = hep::chi_square_dof0( results.begin() + 1, results.end() );

	cout << ">> cumulative result: " << endl << ">> N = " << result.calls() << " I = " << result.value() << " +- " << result.error() << " chi^2/dof = " << chi_square_dof << endl;

	return 0;
}
