#include "vector.hpp"

// constructor 0:
// creating empty 3-vector
__vector__::__vector__ ( void )
{
    dim = 3;
    for ( size_t i = 0; i < dim; i++ )
    {
        coords.push_back( 0 );
    }
}


// constructor 1:
// creating empty d-vector
__vector__::__vector__( int d ) : dim( d )
{
    for ( size_t i = 0; i < d; i++ )
    {
        coords.push_back( 0 );
    }
}

// constructor 2:
// creating vector with given coordinates
__vector__::__vector__( double x, double y, double z )
{
    coords.push_back(x);
    coords.push_back(y);
    coords.push_back(z);
    
    dim = 3;
}

// constructor 3:
// creating vector with given vector of coordinates
__vector__::__vector__( std::vector<double> c )
{
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
    std::cout << "Vector dimension: " << dim << "; components: " << coords[0] << ", " << coords[1] << ", " << coords[2] << std::endl;
}

void __vector__::setCoords( std::vector<double> c )
{
    for ( size_t i = 0; i < c.size(); i++ )
    {
        coords[i] = c[i];
    }
}

void __vector__::setCoords( double x, double y, double z )
{
    assert( dim == 3 && "The vector's size is not 3!" );

    coords[0] = x;
    coords[1] = y;
    coords[2] = z;
}

std::vector<double> __vector__::getCoords( void )
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

void __vector__::normalize ( void )
{
    double len = getLength();
    assert( len != 0 && "Can't normalize zero-length vector!" );

    for ( size_t i = 0; i < this->dim; i++ )
    {
        this->coords[i] /= len;
    }
}

void __vector__::fill_random( double scale ) 
{
    // mt19937 -- Mersenne Twister pseudo-random generator of 32-bit numbers with astate size of 19937 bits
    std::random_device rd; // obtain a random number
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_real_distribution<> distr( 0, scale ); // defines the range

    for ( size_t i = 0; i < this->dim; i++ )
    {
        this->coords[i] = distr( eng );
    } 
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

__vector__ & __vector__::operator=( const __vector__ &rhs )
{
    // checking for self-assignment
    if ( this != &rhs )
    {
        this->dim = rhs.dim;
        this->coords = rhs.coords;
        return *this;
    }

    return *this;
}

__vector__ & __vector__::operator+=( const __vector__ &rhs )
{
    assert( this->dim == rhs.dim && "Vectors' lengths should be equal!");
    
    for ( size_t i = 0; i < this->dim; i++ )
    {
        this->coords[i] += rhs.coords[i];
    }
    
    return *this;
}

__vector__ & __vector__::operator-=( const __vector__ &rhs )
{
    assert( this->dim == rhs.dim && "Vectors' lengths should be equal!");

    for ( size_t i = 0; i < this->dim; i++ )
    {
        this->coords[i] -= rhs.coords[i];
    }

    return *this;
}

__vector__ & __vector__::operator+=( const double &rhs )
{
    for ( size_t i = 0; i < this->dim; i++ )
    {
        this->coords[i] += rhs;
    }

    return *this;
}

__vector__ & __vector__::operator-=( const double &rhs )
{
    for ( size_t i = 0; i < this->dim; i++ )
    {
        this->coords[i] -= rhs;
    }

    return *this;
}

__vector__ & __vector__::operator*=( const double &rhs )
{
    for ( size_t i = 0; i < this->dim; i++ )
    {
        this->coords[i] *= rhs;
    }
    
    return *this;
}

__vector__ & __vector__::operator/=( const double &rhs )
{
    for ( size_t i = 0; i < this->dim; i++ )
    {
        this->coords[i] /= rhs;
    }

    return *this;
}

__vector__ __vector__::power ( const double n )
{
    assert( n > 0 && "Negative powers are unsafe!" );
    
    __vector__ r ( this-> dim ); 
    for ( size_t i = 0; i < this->dim; i++ )
    {
        r.coords[i] = pow( this->coords[i], n );
    }

    return r;
}

double __vector__::dist( const __vector__ &other )
{
    assert( this->dim == other.dim && "Vectors' lengths should be equal!" );

    double r = 0;
    for ( size_t i = 0; i < this->dim; i++ )
    {
        r += pow( this->coords[i] - other.coords[i], 2 ); 
    }

    return pow( r, 0.5 );
}

const double __vector__::operator*( const __vector__ &other ) const
{
    assert( this->dim == other.dim && "Vectors' lengths should be equal!" );

    double r = 0;
    for ( size_t i = 0; i < this->dim; i++ )
    {
        r += this->coords[i] * other.coords[i];
    }

    return r;
}

const __vector__ __vector__::operator+( const __vector__ &other ) const
{
    return __vector__( *this ) += other;
}

const __vector__ __vector__::operator-( const __vector__ &other ) const
{
    return __vector__( *this ) -= other;
}

const __vector__ __vector__::operator+( const double &other ) const
{
    return __vector__( *this ) += other;
}

const __vector__ __vector__::operator-( const double &other ) const
{
    return __vector__( *this ) -= other;
}

const __vector__ __vector__::operator*( const double &other ) const
{
    return __vector__( *this ) *= other;
}

const __vector__ __vector__::operator/( const double &other ) const
{
    return __vector__( *this ) /= other;
}

