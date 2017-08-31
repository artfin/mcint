#include <iostream>

class Atom
{
    public:
        double getMass( void );
        void setMass( double m );

        Atom( double m );

    private:
        double mass;
};
