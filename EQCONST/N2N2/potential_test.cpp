#include <iostream>
#include <chrono>
#include <ctime>

extern "C" void potinit( void );
extern "C" void potn2n2( double* rr, double* theta1, double* theta2, double* phi, double* res );

using namespace std;

int main()
{
	potinit();

	int n = 1000;

	double R = 5.0;
	double theta1 = 0.1;
	double theta2 = 0.1;
	double phi = 0.1;

	auto startTime = chrono::high_resolution_clock::now();
	
	for ( int i = 0; i < n; i++ )
	{
		double potential_value;
		potn2n2( &R, &theta1, &theta2, &phi, &potential_value );
	}

	auto endTime = chrono::high_resolution_clock::now();

	cout << "Potential calls: " << n << "; Time elapsed: " << chrono::duration_cast<chrono::milliseconds>( endTime - startTime ).count() << "ms" << endl;

	return 0;
}
