#include "limits.hpp"
#include <iostream>

using namespace std;

int main( )
{
	Limits limits;
	limits.add_limit( 0, 1.0, 2.0 );
	limits.add_limit( 1, "-inf", "+inf" );
	limits.add_limit( 2, "-inf", 0.0 );
	limits.add_limit( 3, 0.0, "+inf" );

	limits.show_limits();

	return 0;
}
