#include <iostream>

#include "vector.hpp"
#include "random.hpp"

class Atom
{
    public:
        Atom( double m );
        Atom( double m, double x, double y, double z );

        ~Atom();

    private:
        double mass;

        __vector__ pos;
        __vector__ vel;
        __vector__ force;
};

Atom::Atom( double m )
{
    mass = m;
    std::cout << "Creating atom; m = " << mass << std::endl;
} 

Atom::Atom( double m, double x, double y, double z )
{
    mass = m;
    pos.setCoords( x, y, z );

    std::cout << "Creating atom; m = " << mass << std::endl;
    std::cout << "Coords: " << std::endl;
    pos.print(); 
}

~Atom() 
{
}

// setting the velocity according to Maxwell-Boltzmann distribution
void Atom::setMaxwellVel( double temperature )
{
    double sigma = sqrt( BOLTZCONST * temperature / mass );
}

}

int main()
{
    Atom a1( 1.0 );
    Atom a2( 1.0, 0.0, 0.0, 0.0 );

    return 0;
}


