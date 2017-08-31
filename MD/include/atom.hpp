#pragma once

#include <iostream>
#include "vector.hpp"
#include "random.hpp"

// probably need separate file with constants
const double BOLTZCONST = 1e-10;

class Atom
{
    public:
        Atom( double m );
        Atom( double m, double x, double y, double z );

        void setMass( double m );
        double getMass( void );

        void printCoords( void );
        void printVelocity( void );
        void printForce( void );

        void setMaxwellVelocity( double temperature );

    private:
        double mass;

        __vector__ pos;
        __vector__ vel;
        __vector__ force;
};
  

