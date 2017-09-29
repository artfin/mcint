#include <iostream>
#include <chrono>

#include <math.h>
#include <cstddef>
#include <vector>

#include <fstream>
#include <iomanip>

//#include "hep/mc.hpp"
#include "hep/mc-mpi.hpp"

// hamiltonian input: theta(co2), phi, theta(h2), r
#include "co2h2_hamiltonian.hpp"
// double kinetic_energy(double q1, double q2, double q3, double q4, double p1, double p2, double p3, double p4, double Jx, double Jy, double Jz);

// dipole input: r, theta(co2), theta(h2), phi
#include "co2h2_dipole_lr.hpp"
//double dipx(double R, double Theta1, double Theta2, double Phi);
//double dipy(double R, double Theta1, double Theta2, double Phi);
//double dipz(double R, double Theta1, double Theta2, double Phi);

// potential: xr -- in angstroms
// pararead -- intializing procedure
// output of co2h2pes -- in cm^-1
// input: r, theta(co2), theta(h2), phi
extern "C" void pararead( void );
extern "C" void co2h2pes( double *xr, double *xth1, double *xth2, double *xphi, double *potvalue );

// bohr to meter and to angstrom
const double BOHRTOM = 5.2917721067 * pow(10, -11);
const double BOHRTOANG = BOHRTOM * pow(10, 10);

// cm^-1 to hartree
const double CMTOH = 1 / 2.1947 * pow(10, -5);
// hartree to joules
const double HTOJ = 4.35974417 * pow(10, -18);
// boltzmann constant
const long double BOLTZCONST = 1.38064852 * pow(10, -23);

const double COEFF = 0.0036148;
const double RMAX = 100.0;
const double VOLUME = 1.0 / 3.0 * pow( RMAX, 3 );

using namespace std;
using namespace std::placeholders;

struct stop_after_precision
{
	stop_after_precision( double rel_error )
			: rel_error( rel_error )
	{
	}

	bool operator()(std::vector<hep::vegas_result<double>> const& r)
	{
		hep::vegas_verbose_callback<double>(r);

		// compute cumulative result
		auto const result = hep::cumulative_result0(r.begin(), r.end());
	
		if ( result.error() < abs(rel_error * result.value()) )
		{
			cout << ">> relative error " << (result.error() / result.value() ) << " is smaller than the limit " << rel_error << endl;
			
			return false;
		}

		return true;
	}

	double rel_error;
};

double numerator_integrand_( hep::mc_point<double> const& x, double Temperature )
{
	double R_new = x.point()[0];
	double pR_new = x.point()[1];
	double Theta_co2_new = x.point()[2];
	double pTheta_co2_new = x.point()[3];
	double Theta_h2_new = x.point()[4];
	double pTheta_h2_new = x.point()[5];
	double Phi_new = x.point()[6];
	double pPhi_new = x.point()[7];
	double Jx_new = x.point()[8];
	double Jy_new = x.point()[9];
	double Jz_new = x.point()[10];

	double R = tan(M_PI / 2 * R_new);
	double Theta_co2 = M_PI * Theta_co2_new;
	double Theta_h2 = M_PI * Theta_h2_new;
	double Phi = 2 * M_PI * Phi_new;

	double pR = tan(M_PI * (pR_new - 0.5));
	double pTheta_co2 = tan(M_PI * (pTheta_co2_new - 0.5));
	double pTheta_h2 = tan(M_PI * (pTheta_h2_new - 0.5));
	double pPhi = tan(M_PI * (pPhi_new - 0.5));

	double Jx = tan(M_PI * (Jx_new - 0.5));
	double Jy = tan(M_PI * (Jy_new - 0.5));
	double Jz = tan(M_PI * (Jz_new - 0.5));

	// integrating only to R = 100.0
	if ( R > 100.0 )
	{
		return 0;
	}

	double potential_value;
	double R_angstrom = R * BOHRTOANG; 
	co2h2pes(&R_angstrom, &Theta_co2, &Theta_h2, &Phi, &potential_value);
	
	potential_value = potential_value * CMTOH;

	// hamiltonian input: theta(co2), phi, theta(h2), r
	double h = kinetic_energy(Theta_co2,  Phi, Theta_h2, R, pTheta_co2, pPhi, pTheta_h2, pR, Jx, Jy, Jz) + potential_value;

	double boltz_exp = exp ( - h * HTOJ / (BOLTZCONST * Temperature ));

	// dipole input: r, theta(co2), theta(h2), phi
	double dipx_value = dipx(R, Theta_co2, Theta_h2, Phi);
	double dipy_value = dipy(R, Theta_co2, Theta_h2, Phi);
	double dipz_value = dipz(R, Theta_co2, Theta_h2, Phi);

	double dip_squared = pow(dipx_value, 2) + pow(dipy_value, 2) + pow(dipz_value, 2);

	// jacobians
	double jacR = M_PI / 2 * (1 + pow(R, 2));
	double jacTheta_h2 = M_PI;
	double jacTheta_co2 = M_PI;
	double jacPhi = 2 * M_PI;
	double jacpR = M_PI * (1 + pow(pR, 2));
	double jacpTheta_h2 = M_PI * (1 + pow(pTheta_h2, 2));
	double jacpTheta_co2 = M_PI * (1 + pow(pTheta_co2, 2));
	double jacpPhi = M_PI * (1 + pow(pPhi, 2));
	double jacJx = M_PI * (1 + pow(Jx, 2));
	double jacJy = M_PI * (1 + pow(Jy, 2));
	double jacJz = M_PI * (1 + pow(Jz, 2));

	return dip_squared * boltz_exp * jacR * jacTheta_h2 * jacTheta_co2 * jacPhi * jacpR * jacpTheta_h2 * jacpTheta_co2 * jacpPhi * jacJx * jacJy * jacJz;
}	

double denumerator_integrand_( hep::mc_point<double> const& x, double Temperature )
{
	double R_new = x.point()[0];
	double pR_new = x.point()[1];
	double Theta_co2_new = x.point()[2];
	double pTheta_co2_new = x.point()[3];
	double Theta_h2_new = x.point()[4];
	double pTheta_h2_new = x.point()[5];
	double Phi_new = x.point()[6];
	double pPhi_new = x.point()[7];
	double Jx_new = x.point()[8];
	double Jy_new = x.point()[9];
	double Jz_new = x.point()[10];

	double R = tan(M_PI / 2 * R_new);
	double Theta_co2 = M_PI * Theta_co2_new;
	double Theta_h2 = M_PI * Theta_h2_new;
	double Phi = 2 * M_PI * Phi_new;

	double pR = tan(M_PI * (pR_new - 0.5));
	double pTheta_co2 = tan(M_PI * (pTheta_co2_new - 0.5));
	double pTheta_h2 = tan(M_PI * (pTheta_h2_new - 0.5));
	double pPhi = tan(M_PI * (pPhi_new - 0.5));

	double Jx = tan(M_PI * (Jx_new - 0.5));
	double Jy = tan(M_PI * (Jy_new - 0.5));
	double Jz = tan(M_PI * (Jz_new - 0.5));
	
	// integrating only to R = RMAX 
	if ( R > RMAX )
	{
		return 0;
	}
	
	double potential_value;
	double R_angstrom = R * BOHRTOANG; 
	co2h2pes(&R_angstrom, &Theta_co2, &Theta_h2, &Phi, &potential_value);
	
	potential_value = potential_value * CMTOH;

	// hamiltonian input: theta(co2), phi, theta(h2), r
	double h = kinetic_energy(Theta_co2,  Phi, Theta_h2, R, pTheta_co2, pPhi, pTheta_h2, pR, Jx, Jy, Jz) + potential_value;

	double boltz_exp = exp ( - h * HTOJ / (BOLTZCONST * Temperature ));

	// jacobians
	double jacR = M_PI / 2 * (1 + pow(R, 2));
	double jacTheta_h2 = M_PI;
	double jacTheta_co2 = M_PI;
	double jacPhi = 2 * M_PI;
	double jacpR = M_PI * (1 + pow(pR, 2));
	double jacpTheta_h2 = M_PI * (1 + pow(pTheta_h2, 2));
	double jacpTheta_co2 = M_PI * (1 + pow(pTheta_co2, 2));
	double jacpPhi = M_PI * (1 + pow(pPhi, 2));
	double jacJx = M_PI * (1 + pow(Jx, 2));
	double jacJy = M_PI * (1 + pow(Jy, 2));
	double jacJz = M_PI * (1 + pow(Jz, 2));

	return boltz_exp * jacR * jacTheta_h2 * jacTheta_co2 * jacPhi * jacpR * jacpTheta_h2 * jacpTheta_co2 * jacpPhi * jacJx * jacJy * jacJz;
}	


int main( int argc, char* argv[] )
{
	// initialize MPI
	MPI_Init( &argc, &argv ); // ??

	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	// initialzing potential parameters
	pararead();

	double LTEMP = 300.0;
	double HTEMP = 300.0;
	double STEP = 5.0;

	if ( rank == 0 )
	{	
		cout << "------------ Zeroth moment calculation for CO2-H2 ----------" << endl;
		cout << "LTEMP = " << LTEMP << "; HTEMP = " << HTEMP << "; STEP = " << STEP << endl;
	}

	chrono::high_resolution_clock::time_point startTime = chrono::high_resolution_clock::now();

	ofstream file;
	file.open( "out.txt" );
	
	//hep::vegas_callback<double>(stop_after_precision(0.01));
	hep::mpi_vegas_callback<double>( hep::mpi_vegas_verbose_callback<double> );

	for ( double TEMP = LTEMP; TEMP <= HTEMP; TEMP += STEP )
	{
		chrono::high_resolution_clock::time_point cycleStartTime = chrono::high_resolution_clock::now();
		
		auto numerator_integrand = bind( numerator_integrand_, _1, TEMP );
		auto denumerator_integrand = bind( denumerator_integrand_, _1, TEMP );

		auto numerator_results = hep::mpi_vegas(
			MPI_COMM_WORLD,
			hep::make_integrand<double>(numerator_integrand, 11),
			std::vector<std::size_t>(10, 1e5)
		);

		auto denumerator_results = hep::mpi_vegas(
			MPI_COMM_WORLD,
			hep::make_integrand<double>(denumerator_integrand, 11),
			std::vector<std::size_t>(10, 1e5)
		);

		auto numerator_result = hep::cumulative_result0( numerator_results.begin() + 1, numerator_results.end() );
		auto denumerator_result = hep::cumulative_result0( denumerator_results.begin() + 1, denumerator_results.end() );

		double numerator_chi_square_dof = hep::chi_square_dof0( numerator_results.begin() + 1, numerator_results.end() );
		double denumerator_chi_square_dof = hep::chi_square_dof0( denumerator_results.begin() + 1, denumerator_results.end() );

		double zero_moment = numerator_result.value() / denumerator_result.value() * COEFF * VOLUME;

		if ( rank == 0 )
		{
			cout << ">> Numerator; cumulative result: " << endl;
			cout << ">> N = " << numerator_result.calls() << "; I = " << numerator_result.value() << " +- " << numerator_result.error() << endl;
			cout << ">> chi^2/dof = " << numerator_chi_square_dof << endl << endl;

			cout << ">> Denumerator; cumulative_result: " << endl;
			cout << ">> N = " << denumerator_result.calls() << "; I = " << denumerator_result.value() << " +- " << denumerator_result.error() << endl;
			cout << ">> chi^2/dof = " << denumerator_chi_square_dof << endl << endl;

			file << setprecision(5) << TEMP << " " << setprecision(10) << numerator_result.value() << " " << denumerator_result.value() << " " << setprecision(10) << zero_moment << endl;

			chrono::high_resolution_clock::time_point cycleEndTime = chrono::high_resolution_clock::now();

			cout << "Cycle time: " << chrono::duration_cast<chrono::milliseconds>(cycleEndTime - cycleStartTime).count() / 1000.0 << " s" << endl; 
		}
	}

	file.close();

	chrono::high_resolution_clock::time_point endTime = chrono::high_resolution_clock::now();

	if ( rank == 0 )
	{
		cout << endl << "Total time elapsed: " << chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count() / 1000.0 << " s" << endl;
	}
	
	MPI_Finalize();

	return 0;
}
