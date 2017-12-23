#pragma once

//#include "hep/mc.hpp"
#include "hep/mc-mpi.hpp"
#include "integrand.hpp"

#include <iostream>
#include <ctime>

using std::cout;
using std::endl;

class Integrator
{
public:
	Integrand& integrand;
	int rank;

	Integrator( Integrand& integrand, int rank ) : integrand(integrand), rank(rank)
	{
		//cout << "(integrator) integrand constructor" << endl;
	}
	
	void set_callback( void )
	{
		// set the verbose vegas callback function
		hep::vegas_callback<double>(hep::vegas_verbose_callback<double>);
	}
	
	double run_integration( const int& niter, const int& ndots )
	{
		clock_t start = clock();

		//auto results = hep::vegas
			//(
				//hep::make_integrand<double>( integrand, integrand.dim() ),
				//std::vector<std::size_t>( niter, ndots )
			//);

		auto results = hep::mpi_vegas(
			MPI_COMM_WORLD,
			hep::make_integrand<double>( integrand, integrand.dim() ),
			std::vector<std::size_t>( niter, ndots )
		);

		auto result = hep::cumulative_result0(results.begin() + 1, results.end());
		double chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());

		if ( rank == 0 )
		{
			integrand.show_default_value();
			cout << ">> value: " << result.value() << "; error: " << result.error() << "; relative error (%): " << (double) result.error() / result.value() * 100 << "; time elaped: " << (double) ( clock() - start ) / CLOCKS_PER_SEC << " s" << endl; 
		}

		return result.value();
		//cout << ">> N = " << result.calls() << "; I = " << result.value() << " +- " << result.error() << endl;
		//cout << ">> chi^2/dof = " << chi_square_dof << endl;
	}	
};
