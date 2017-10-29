#pragma once

class numericalIntegration
{
	private:
		int nIntervals;

	public:
		numericalIntegration( const int &nIntervals);
		double integrate( double (*pFun)(const double &x), const double &lowBound, const double &upBound );
	
};
