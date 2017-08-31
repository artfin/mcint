#include "atom.hpp"

void Atom::setMass( double m )
{
    mass = m;
}

double Atom::getMass( void )
{
    return mass;
}

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

void Atom::printCoords( void )
{
    std::cout << "Atom coordinates: " << std::endl;
    pos.print();
}

void Atom::printVelocity( void )
{
    std::cout << "Atom velocity: " << std::endl;
    vel.print();
}

void Atom::printForce( void )
{
    std::cout << "Atom force: " << std::endl;
    force.print();
}

// setting the velocity according to Maxwell-Boltzmann distribution
void Atom::setMaxwellVelocity ( double temperature )
{
    // random Gaussian value with mean = 0, sigma = sqrt( k * T / m )
    double sigma = sqrt( BOLTZCONST * temperature / mass );
    vel.randomGaussian( 0, sigma );
}


