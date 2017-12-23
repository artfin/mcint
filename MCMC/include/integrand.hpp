#pragma once

#include "hep/mc.hpp"
#include "limits.hpp"

#include <iostream>
#include <vector>
#include <string>
#include <functional>

using std::cout;
using std::endl;
using std::function;
using std::tuple;

class Integrand 
{
private:
	function<double(hep::mc_point<double> const& x)> f;
	Limits limits;

	int argn;
	double default_value;
public:

	Integrand( ) 
	{
		cout << "(integrand) empty constructor" << endl;
	}

	Integrand( function<double(hep::mc_point<double> const& x)> f ) : f(f) 
	{
		//cout << "(integrand) function constructor" << endl;
	}

	Integrand( function<double(hep::mc_point<double> const& x)> f, Limits limits ) : f(f), limits(limits)
	{
		//cout << "(integrand) function, Limits constructor" << endl;
	}
	
	double dim( void ) { return limits.limits.size(); }	
	Limits* set_limits( void ) { return &limits; }
	void set_argn( int n ) { argn = n; }
	void set_default_value( double v ) { default_value = v; };
	void show_default_value( void ) { cout << "default_value: " << default_value << endl; } 
	void show_limits( void ) { limits.show_limits(); }

	// returns value of function with some predefined value
	// function if transformed to the n-dimensional cube
	// returned value is multiplied by the value of corresponding jacobian
	double operator()( hep::mc_point<double> const& x )
	{
		std::vector<double> v = x.point(); 
		std::vector<double>::iterator it = v.begin();
		v.insert( it + argn, 1, default_value );

		double jac = 1;

		// transforming to [0, 1]^n - cube	
		for ( size_t i = 0; i < v.size(); i++ )
		{
			if ( i == argn ) 
				continue;

			if ( limits.limits[i].chvarType == Limit::chvarTypes::INFINF )
			{
				v[i] = tan( M_PI * (v[i] - 0.5) );
				jac *= M_PI * (1 + v[i] * v[i]);
			}
			if ( limits.limits[i].chvarType == Limit::chvarTypes::FINITE )
			{
				v[i] = ( limits.limits[i].ub - limits.limits[i].lb ) * v[i];
				jac *= ( limits.limits[i].ub - limits.limits[i].lb );
			}
			if ( limits.limits[i].chvarType == Limit::chvarTypes::ZEROINF )
			{
				v[i] = tan( M_PI / 2 * v[i] );
				jac *= M_PI / 2 * (1 + v[i] * v[i]);
			}
		}

		hep::mc_point<double> p{ v };
		return f( p ) * jac;		
	}
};
