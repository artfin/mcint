#include <boost/math/special_functions/gamma.hpp>

#include <iostream> 
#include <ctime>

using namespace std;

int main()
{
	clock_t start = clock();
	double x = boost::math::gamma_p( 0.1, 0.1 );
	cout << "x: " << x << endl;

	cout << "Time elapsed: " << ( clock() - start ) / (double) CLOCKS_PER_SEC* 1000 << "ms" << endl;

	return 0;
}
