#pragma once

#include <iostream>
#include <math.h>
#include <vector>
#include <assert.h>
#include <random>
#include "random.hpp"

class __vector__
{
    public:
        void print( void );
        void setDimension( int d );
        void setCoords( std::vector<double> c );
        void setCoords( double x, double y, double z );
        void normalize ( void );

        void randomUniform( double scale = 1.0 );
        void randomGaussian( double mean, double sigma );
            
        int getDimension( void );
        std::vector<double> getCoords( void );

        double getLength( void );
        __vector__ power ( const double n );
        double dist( const __vector__ &other );

        // constructors
        __vector__( void );
        __vector__( int d );
        __vector__( double x, double y, double z );
        __vector__( std::vector<double> c );

        // operator overloading
        // assignment
        __vector__ & operator=( const __vector__ &rhs );

        // destructive compound assignment
        __vector__ & operator+=( const __vector__ &rhs );
        __vector__ & operator-=( const __vector__ &rhs );

        __vector__ & operator+=( const double &rhs );
        __vector__ & operator-=( const double &rhs );
        __vector__ & operator*=( const double &rhs );
        __vector__ & operator/=( const double &rhs );


        // binary operators
        const __vector__ operator+( const __vector__ &other ) const;
        const __vector__ operator-( const __vector__ &other ) const;

        // scalar product
        const double operator*( const __vector__ &other ) const;

        const __vector__ operator+( const double &rhs ) const;
        const __vector__ operator-( const double &rhs ) const;
        const __vector__ operator*( const double &rhs ) const;
        const __vector__ operator/( const double &rhs ) const;

        // destructor
        ~__vector__();

    private:
        int dim;
        std::vector<double> coords;
};

