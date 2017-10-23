#include "hep/mc-mpi.hpp"

#include <math.h>
#include <cstddef>
#include <vector>

#include <ctime>

#include <fstream>
#include <iomanip> // setprecision

#include "n2n2_hamiltonian.hpp"

extern "C" void potinit(void);
extern "C" void potn2n2(double* rr, double* theta1, double* theta2, double* phi, double* res);

using namespace std;
using namespace std::placeholders;

// cm^-1 to hartree
const double CMTOH = 1.0 / 2.1947 * pow(10, -5);
// hartree to joules
const double HTOJ = 4.35974417 * pow(10, -18);
// boltzmann constant
const long double BOLTZCONST = 1.38064852 * pow(10, -23);
// avogadro number
const double AVOGADRO  = 6.022140857 * pow(10, 23); 
// pascals to atmospheres
const double PATOATM = 101325.0;
// universal gas consant
const double UGASCONST = 8.314;

// planck constant
const long double PLANKCONST = 6.626070040 * pow(10, -34);
const long double PLANKCONST2 = pow(PLANKCONST, 2);

// dalton to kg
const long double DA = 1.660539040 * pow(10, -27);

const double BOHRTOM = 0.529177 * pow( 10, -10 );
const double N2_LENGTH = 2.0744 * BOHRTOM;

const double MASS_COEFF = 10000.0;

// particles molar masses, g/mol
const double N2_MOLARMASS = 28.0 * MASS_COEFF;
const double N_MOLARMASS = 0.5 * N2_MOLARMASS; 
const double COMPLEX_MOLARMASS = 56.0 * MASS_COEFF;

double integrand_( hep::mc_point<double> const& x, double Temperature )
{
	double R_new = x.point()[0];
	double pR_new = x.point()[1];
	double Theta1_new = x.point()[2];
	double pTheta1_new = x.point()[3];
	double Theta2_new = x.point()[4];
	double pTheta2_new = x.point()[5];
	double Phi_new = x.point()[6];
	double pPhi_new = x.point()[7];
	double Jx_new = x.point()[8];
	double Jy_new = x.point()[9];
	double Jz_new = x.point()[10];

	double R = tan(M_PI / 2 * R_new);
	double Theta1 = M_PI * Theta1_new;
	double Theta2 = M_PI * Theta2_new;
	double Phi = 2 * M_PI * Phi_new;

	double pR = tan(M_PI * (pR_new - 0.5));
	double pTheta1 = tan(M_PI * (pTheta1_new - 0.5));
	double pTheta2 = tan(M_PI * (pTheta2_new - 0.5));
	double pPhi = tan(M_PI * (pPhi_new - 0.5));

	double Jx = tan(M_PI * (Jx_new - 0.5));
	double Jy = tan(M_PI * (Jy_new - 0.5));
	double Jz = tan(M_PI * (Jz_new - 0.5));

	if ( R < 4.4 ) 
	{
		return 0;
	}
	
	double potential_value;
	potn2n2(&R, &Theta1, &Theta2, &Phi, &potential_value);

	potential_value = potential_value * CMTOH;

	double h = kinetic_energy(R, Theta1, Theta2, Phi, pR, pTheta1, pTheta2, pPhi, Jx, Jy, Jz) + potential_value;

	if ( h < 0 )
	{
		double boltz_exp = exp ( - h * HTOJ / (BOLTZCONST * Temperature ));
	
		// jacobians
		double jacR = M_PI / 2 * (1 + pow(R, 2));
		double jacTheta1 = M_PI;
		double jacTheta2 = M_PI;
		double jacPhi = 2 * M_PI;
		double jacpR = M_PI * (1 + pow(pR, 2));
		double jacpTheta1 = M_PI * (1 + pow(pTheta1, 2));
		double jacpTheta2 = M_PI * (1 + pow(pTheta2, 2));
		double jacpPhi = M_PI * (1 + pow(pPhi, 2));
		double jacJx = M_PI * (1 + pow(Jx, 2));
		double jacJy = M_PI * (1 + pow(Jy, 2));
		double jacJz = M_PI * (1 + pow(Jz, 2));

		return boltz_exp * jacR * jacTheta1 * jacTheta2 * jacPhi * jacpR * jacpTheta1 * jacpTheta2 * jacpPhi * jacJx * jacJy * jacJz;
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

int main( int argc, char* argv[] )
{
	// initialize MPI
	MPI_Init( &argc, &argv ); 

	int rank;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	// initializing auxiliary potential values
	potinit();

	// for main process
	if ( rank == 0 )
	{
		cout << "--- Computing FULL EQCONST for N2N2 (Avoird potential) --- " << endl;
	}

	const double LTEMP = 100.0;
	const double HTEMP = 350.0;
	const double STEP = 10.0;

	clock_t full_clock = clock();
	clock_t cycle_clock;

	ofstream file;
	file.open( "full_1_0.txt" );

	hep::mpi_vegas_callback<double>( hep::mpi_vegas_verbose_callback<double> );

	for ( double TEMP = LTEMP; TEMP <= HTEMP; TEMP += STEP )
	{
		cycle_clock = clock();

		auto integrand = bind( integrand_, _1, TEMP );
		
		auto results = hep::mpi_vegas(
			MPI_COMM_WORLD,
			hep::make_integrand<double>( integrand, 11 ),
			std::vector<std::size_t>(20, 3e5)
		);

		auto result = hep::cumulative_result0(results.begin() + 1, results.end());
		double chi_square_dof = hep::chi_square_dof0(results.begin() + 1, results.end());

		if ( rank == 0 )
		{
			cout << ">> Cumulative result: " << endl;
			cout << ">> N = " << result.calls() << "; I = " << result.value() << " +- " << result.error() << endl;
			cout << ">> chi^2/dof = " << chi_square_dof << endl;

			cout << "Time needed: " << ( clock() - cycle_clock ) / (double)(CLOCKS_PER_SEC) << "s" << endl;
		}

		double Qtr_complex = calculate_qtr(COMPLEX_MOLARMASS, TEMP );
		double Qtr_N2 = calculate_qtr(N2_MOLARMASS, TEMP );

		// seems to be right
		double Qrot_N2 = 4 * pow(M_PI, 2) * BOLTZCONST * TEMP / PLANKCONST2 * 0.5 * N_MOLARMASS * DA * pow(N2_LENGTH, 2);
		double Q_N2 = Qtr_N2 * Qrot_N2;
   		

		double eqconst = AVOGADRO / ( UGASCONST * TEMP ) * Qtr_complex / pow(Q_N2, 2) * result.value() * PATOATM / (16 * pow(M_PI, 5)) / 2; 

		if ( rank == 0 )
		{
			cout << "Temperature: " << TEMP << "; EQCONST: " << eqconst << endl << endl;

			file << setprecision(5) << TEMP << " " << setprecision(10) << eqconst << endl;
		}
	}

	file.close();

	if ( rank == 0 )
	{
		cout << "Total time elapsed: " << ( clock() - full_clock ) / (double) CLOCKS_PER_SEC << "s" << endl;
	}

	MPI_Finalize();

	return 0;
}

