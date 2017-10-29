#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

#include "ar_he_pes.h"

// Warn about use of deprecated functions
#define GNUPLOT_DEPRECATE_WARN
#include "gnuplot-iostream.h"

const double HE_MASS = 4.00260325413;
const double AR_MASS = 39.9623831237; 
const double PROTON_TO_ELECTRON_RATIO = 1836.15267389; 

const double MU = HE_MASS * AR_MASS / ( HE_MASS + AR_MASS ) * PROTON_TO_ELECTRON_RATIO; 

// similar to python's linspace
std::vector<double> linspace( const double min, const double max, const int npoints )
{
	const double step = ( max - min ) / ( npoints - 1 );
	
	std::vector<double> res;
	for ( double temp = min; temp <= max; temp += step )
	{
		res.push_back( temp );
	}

	return res;
}

void plot_signal( Gnuplot &gp, std::vector<double> &freqs, std::vector<double> &v, const double lbound = 0.0, const double ubound = 10.0 )
{
	std::vector< std::pair<double, double>> signal;

	for ( int i = 0; i < freqs.size(); i++ )
	{
		signal.push_back( std::make_pair( freqs[i], v[i] ));
	}
	
	std::ostringstream strs;
	strs << lbound;
	std::string lbound_str = strs.str();

	// emptying sstream
	strs.str("");

	strs << ubound;
	std::string ubound_str = strs.str();

	gp << "set xrange [" << lbound_str << ":" << ubound_str << "]\n";
		
	gp << "plot '-' with lines title 'signal'\n";
	gp.send1d( signal );

	gp.flush();
}

void plot_potential( void )
{
	Gnuplot gp;

	const int npoints = 1000;
	std::vector<double> x = linspace( 0.5, 12.0, npoints );
	std::vector<double> pot;

	for ( int i = 0; i < npoints; i++ )
	{
		pot.push_back( ar_he_pot(x[i]) );
	}	

	plot_signal( gp, x, pot, 5.5, 12.0 );
}

class MakeTrajectory
{
	private:
		double energy;
		double j;
	
		std::vector<double> t;
		double rmax = 40.0;
		double rmin;
		
		double trapezoid( const double &lowBound, const double &upBound, const int &nIntervals );
		double Rintegrand( const double &r );

		double bisection_function( const double &r );
		double find_rmin( double a, double b, double eps, int maxSteps );

	public:
		MakeTrajectory( double energy, double j );

		void integrate( std::vector<double> &t,
						std::vector<double> &r );
		
};

MakeTrajectory::MakeTrajectory( double energy, double j )
{
	this->energy = energy;
	this->j = j;	

	const double eps = 1e-6;
	const int maxSteps = 100;
	this->rmin = find_rmin( 0.1, 15.0, eps, maxSteps );

	std::cout << "Rmin: " << this->rmin << std::endl;
}

double MakeTrajectory::bisection_function( const double &r )
{
	double expr = 2.0 / MU * ( this->energy - ar_he_pot(r) ) - pow( this->j / MU / r, 2 );
	return expr;
}

double MakeTrajectory::find_rmin( double a, double b, const double eps, const int maxSteps )
{
	double c;
	if ( bisection_function(a) * bisection_function(b) <= 0 )
	{
		int iter = 1;

		do {
			c = 0.5 * (a + b);
			
			double bfa = bisection_function(a);
			double bfc = bisection_function(c);
			
			//std::cout << "iter: " << iter << 
					//"; a = " << a << 
					//"; b = " << b << 
					//"; f(c) = " << bfc << 
					//"; |a - b| = " << fabs(a - b) << std::endl;

			if ( bfa * bfc > 0 )
			{
				a = c;
			}
			else if ( bfa * bfc < 0 )
			{
				b = c;
			}
			
			iter++;
		} while ( fabs(a - b) >= eps && iter <= maxSteps );

		return c;
	}
	else
	{
		std::cout << "Invalid interval!" << std::endl;

		double bfa = bisection_function( a );
		double bfb = bisection_function( b );
		std::cout << "left-side value = " << bfa << std::endl;
		std::cout << "right-side value = " << bfb << std::endl;
	}
}

double MakeTrajectory::Rintegrand( const double &r )
{
	double expr = 2.0 / MU * (this->energy - ar_he_pot(r) ) - pow(this->j, 2)/ (pow(MU, 2) * pow(r, 2));
	return 1.0 / sqrt( expr );
}

double MakeTrajectory::trapezoid( const double &lowBound, const double &upBound, const int &nIntervals )
{
	double res = 0;
	double step = ( upBound - lowBound ) / ( nIntervals - 1 );

	for ( double x = lowBound + step; x <= upBound - step; x += step )
	{
		res += 2 * Rintegrand( x );
	}

	res += Rintegrand( lowBound );
	res += Rintegrand( upBound );

	double coeff = ( upBound - lowBound ) / ( 2 * nIntervals );
	res *= coeff;

	return res;	
}

void MakeTrajectory::integrate( std::vector<double> &t,
								std::vector<double> &r )
{
	int steps = 50;
	double rstep = (this->rmax - this->rmin) / steps;
	std::cout << "rstep: " << rstep << std::endl;

	double eps = 1e-6;
	double ubound;
	double res;
	double curr_t;

	bool marker = false;

	for ( double lbound = this->rmax - rstep; lbound >= this->rmin - rstep; lbound -= rstep )
	{
		ubound = lbound + rstep;

		const int nIntervals = 500000;

		if ( fabs(lbound - this->rmin) < 0.1 )
		{
			lbound += eps;
			marker = true;
		}
		
		res = trapezoid( lbound, ubound, nIntervals );
		
		if ( marker == true )
		{
			lbound -= eps;
		}

		if ( t.size() == 0 )
		{
			t.push_back( res );
		}	
		else
		{
			curr_t = t.end()[-1] + res;
			std::cout << "t: " << curr_t << "; r: " << lbound << std::endl;
			t.push_back( curr_t );
		}

		r.push_back( lbound );
	}
}

int main()
{
	MakeTrajectory traj( 0.0003869542836545503, 22.2353 ); 
	
	std::vector<double> t;
	std::vector<double> r;
	
	traj.integrate( t, r );

	return 0;
}
