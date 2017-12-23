#include "hep/mc.hpp"

#include <iostream>
#include <cmath>

using namespace std;

class power
{
	public:
		power( double exponent ) : exponent(exponent) { }

		double operator()(hep::mc_point<double> const& x )
		{
			return std::pow( x.point()[0], exponent );
		}

	private:
		double exponent;
};

int main()
{
	power p(9);

    // set the verbose vegas callback function
    hep::vegas_callback<double>(hep::vegas_verbose_callback<double>);
    
	auto results = hep::vegas(
        hep::make_integrand<double>(p, 1),
        std::vector<std::size_t>(10, 1000)
    );
    auto result = hep::cumulative_result0(results.begin() + 1, results.end());
    double chi_square_dof = hep::chi_square_dof0(results.begin() + 1,
        results.end());
    // print the cumulative result
    std::cout << ">> cumulative result (excluding first iteration):\n>> N="
        << result.calls() << " I=" << result.value() << " +- " << result.error()
        << " chi^2/dof=" << chi_square_dof << "\n";

	return 0;	
}
