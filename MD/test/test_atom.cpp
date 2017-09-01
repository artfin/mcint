#include "atom.hpp"

using namespace std;

int main()
{
    Atom a1( 1.0 );

    Atom a2( 1.0, 0.0, 0.0, 0.0 );

    cout << "Setting Maxwellian velocity for a2: "<< endl;
    a2.setMaxwellVelocity( 300.0 );
    a2.printVelocity();

    return 0;
}
