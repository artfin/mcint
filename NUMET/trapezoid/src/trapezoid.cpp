#include <iostream>
#include "trapezoid.h"

numericalIntegration::numericalIntegration( const int &nIntervals )
{
	this->nIntervals = nIntervals;
}

double numericalIntegration::integrate( double (*pFun)(const double &x), const double &lowBound, const double &upBound )
{
	double res = 0;
	double step = ( upBound - lowBound ) / ( this->nIntervals - 1 );

	for ( double x = lowBound + step; x <= upBound - step; x += step )
	{
		res += 2 * (*pFun)( x );
	}

	res += pFun( lowBound );
	res += pFun( upBound );

	double coeff = ( upBound - lowBound ) / ( 2 * this->nIntervals );
	res *= coeff;

	return res;	
}
