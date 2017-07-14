#include <iostream>

extern "C" void print_hi(void);

using namespace std;

int main()
{
	print_hi();
	cout << "Hello from C++" << endl;

	return 0;
}
