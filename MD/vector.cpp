#include <iostream>
#include <vector>
#include <math.h>
#include <assert.h>

using namespace std;

class __vector__
{
    public:
        void print( void );
        void setDimension( int d );
        void setCoords( vector<double> c );

        int getDimension( void );
        vector<double> getCoords( void );

        double getLength( void );

        // constructors
        __vector__( int d );
        __vector__( double x, double y, double z );
        __vector__( vector<double> c );

        // overload + operator: adding two vectors
        __vector__ operator+( const __vector__ &v )
        {
            assert( v.dim == this->dim && "Vectors' lengths should be equal!" );

            __vector__ res( this->dim );
            for ( size_t i = 0; i < this->dim; i++ )
            {
                res.coords[i] = this->coords[i] + v.coords[i];
            }

            return res;
        }

        // overload + operator: adding double
        __vector__ operator+( const double a )
        {
            __vector__ res( this->dim );
            for ( size_t i = 0; i < this->dim; i++ )
            {
                res.coords[i] = this->coords[i] + a;
            }

            return res;
        }

        // overload - operator: subtracting two vectors
        __vector__ operator-( const __vector__ &v )
        {
            assert( v.dim == this->dim && "Vectors' lengths should be equal!" );

            __vector__ res( this-> dim );
            for ( size_t i = 0; i < this->dim; i++ )
            {
                res.coords[i] = this->coords[i] - v.coords[i];
            }

            return res;
        }

        // overload - operator: subracting double
        __vector__ operator-( const double a )
        {
            __vector__ res ( this-> dim );
            for ( size_t i = 0; i < this->dim; i++ )
            {
                res.coords[i] = this->coords[i] - a;
            }

            return res;
        }

        // overload * operator: scalar product
        double operator*( const __vector__ &v )
        {
            assert( v.dim == this->dim && "Vectors' lengths should be equal!" );

            double res = 0;
            for ( size_t i = 0; i < v.dim; i++ )
            {
                res += this->coords[i] * v.coords[i];
            }

            return res;
        }

        // overload * operator: component-wise multiplication by double
        __vector__ operator* ( const double a )
        {
            __vector__ res ( this-> dim );
            for ( size_t i = 0; i < this->dim; i++ )
            {
                res.coords[i] = this->coords[i] * a;
            }

            return res;
        }

        bool operator=( const __vector__ &v )
        {
            assert( v.dim == this->dim && "Vectors' lengths should be equal!" );

            for ( size_t i = 0; i < this->dim; i++ )
            {
                if ( v.coords[i] != this->coords[i] )
                {
                    return false;
                }
            }

            return true;
        }

        // component-wise power
        __vector__ power ( const double n )
        {
            assert( n > 0 && "Negative powers are unsafe." );

            __vector__ res ( this->dim );
            for ( size_t i = 0; i < this->dim; i++ )
            {
                res.coords[i] = pow( this->coords[i], n );
                cout << "pushing " << pow(this->coords[i], n) << endl;
            } 

            return res;
        }

        // destructor
        ~__vector__();

    private:
        int dim;
        vector<double> coords;
};

// constructor 1:
// creating empty 3-vector
__vector__::__vector__( int d ) : dim( d )
{
    // cout << "Invoking constructor type 1: __vector__ is being created..." << endl;

    for ( int i = 0; i < d; i++ )
    {
        coords.push_back( 0 );
    }
}

// constructor 2:
// creating vector with given coordinates
__vector__::__vector__( double x, double y, double z )
{
    // cout << "Invoking constructor type 2: __vector__ is being created..." << endl;
    
    coords.push_back(x);
    coords.push_back(y);
    coords.push_back(z);
    
    dim = 3;
}

// constructor 3:
// creating vector with given vector of coordinates
__vector__::__vector__( vector<double> c )
{
    // cout << "Invoking constructor 3: __vector__ is being created..." << endl;

    coords = c;
    dim = c.size();
}

// destructor
__vector__::~__vector__( void )
{
}

// printer
void __vector__::print( void )
{
    cout << "Vector dimension: " << dim << "; components: " << coords[0] << ", " << coords[1] << ", " << coords[2] << endl;
}

void __vector__::setCoords( vector<double> c )
{
    for ( size_t i = 0; i < c.size(); i++ )
    {
        coords[i] = c[i];
    }
}

vector<double> __vector__::getCoords( void )
{
    return coords;
}

void __vector__::setDimension( int d )
{
    dim = d;
}

int __vector__::getDimension( void )
{
    return dim;
}

double __vector__::getLength( void )
{
    double len = 0;
    for ( size_t i = 0; i < dim; i++ )
    {
        len += pow(coords[i], 2);
    }

    len = pow(len, 0.5);
    return len;
}

int main()
{
    cout << "Initializing vector:" << endl; 
    __vector__ v2( 2.0, 1.0, 1.5 );
    v2.print();

    cout << endl << endl;

    cout << "Resetting vector coordinates: " << endl;
    vector<double> c {0.0, 3.0, 4.0 };
    v2.setCoords( c );
    v2.print();
    
    cout << endl << endl;

    double l = v2.getLength();
    cout << "Calculating vector length: " << l << endl;
    
    cout << endl << endl;

    vector<double> v {1.0, 1.0, 1.0};
    __vector__ v1 ( v );
    __vector__ v3 ( 3 );

    v3 = v1 + v2;
    v1.print();
    v2.print();
    cout << "Result of adding: " << endl;
    v3.print();

    cout << endl << endl;
    
    v3 = v1 - v2;
    v1.print();
    v2.print();
    cout << "Result of substracting: " << endl;
    v3.print();

    cout << endl << endl;

    double res = v1 * v2;
    v1.print();
    v2.print();
    cout << "Scalar product: " << res << endl;

    cout << endl << endl;

    __vector__ v4 = v1.power( 2 ); 
    v1.print();
    cout << "Square: " << endl;
    v4.print();

    return 0;
}
