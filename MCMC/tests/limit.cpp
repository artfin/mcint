#include "limit.hpp"

#include <iostream>
#include <string>

using namespace std;

int main()
{
	Limit limit1{ 1, 0.0, "inf" };

	if ( limit1.lbType == limit1.limitTypes::STRING )
		cout << "(limit) string lb_str = " << limit1.lb_str << endl;

	if ( limit1.lbType == limit1.limitTypes::DOUBLE )
		cout << "(limit) double lb = " << limit1.lb << endl;

	if ( limit1.ubType == limit1.limitTypes::STRING )
		cout << "(limit) string ub_str = " << limit1.ub_str << endl;

	if ( limit1.ubType == limit1.limitTypes::DOUBLE )
		cout << "(limit) double ub = " << limit1.ub << endl;
	
	Limit limit2{ 1, "-inf", "inf" }; 

	if ( limit2.lbType == limit2.limitTypes::STRING )
		cout << "(limit) string lb_str = " << limit2.lb_str << endl;

	if ( limit2.lbType == limit2.limitTypes::DOUBLE )
		cout << "(limit) double lb = " << limit2.lb << endl;

	if ( limit2.ubType == limit2.limitTypes::STRING )
		cout << "(limit) string ub_str = " << limit2.ub_str << endl;

	if ( limit2.ubType == limit2.limitTypes::DOUBLE )
		cout << "(limit) double ub = " << limit2.ub << endl;

	Limit limit3{ 1, 0.0, 50.0 };

	if ( limit3.lbType == limit3.limitTypes::STRING )
		cout << "(limit) string lb_str = " << limit3.lb_str << endl;

	if ( limit3.lbType == limit3.limitTypes::DOUBLE )
		cout << "(limit) double lb = " << limit3.lb << endl;

	if ( limit3.ubType == limit3.limitTypes::STRING )
		cout << "(limit) string ub_str = " << limit3.ub_str << endl;

	if ( limit3.ubType == limit3.limitTypes::DOUBLE )
		cout << "(limit) double ub = " << limit3.ub << endl;

	Limit limit4{ 1, "-inf", 0.0 };

	if ( limit4.lbType == limit4.limitTypes::STRING )
		cout << "(limit) string lb_str = " << limit4.lb_str << endl;

	if ( limit4.lbType == limit4.limitTypes::DOUBLE )
		cout << "(limit) double lb = " << limit4.lb << endl;

	if ( limit4.ubType == limit4.limitTypes::STRING )
		cout << "(limit) string ub_str = " << limit4.ub_str << endl;

	if ( limit4.ubType == limit4.limitTypes::DOUBLE )
		cout << "(limit) double ub = " << limit4.ub << endl;

	return 0;
}

