#include <iostream>

using namespace std;

class MyClass
{
    public:
        double getA( void );
        void setA ( double var );
       
        MyClass( void ); 
        MyClass( double var );

        // operator overloading:
        // assignment
        MyClass & operator=( const MyClass &rhs );

        // destructive compound assignment
        MyClass & operator+=( const MyClass &rhs );
        MyClass & operator-=( const MyClass &rhs );
        MyClass & operator*=( const MyClass &rhs );
        MyClass & operator/=( const MyClass &rhs );

        // binary operators
        const MyClass operator+( const MyClass &other ) const;
        const MyClass operator-( const MyClass &other ) const;
        const MyClass operator*( const MyClass &other ) const;
        const MyClass operator/( const MyClass &other ) const;

        bool operator==( const MyClass &other ) const;
        bool operator!=( const MyClass &other ) const;

    private:
        double a;
};

// assignment of class would mean setting "a" variable
MyClass & MyClass::operator=( const MyClass &rhs )
{
    // checking for self-assignment
    if ( this != &rhs )
    {
        this->a = rhs.a;
        return *this;
    }

    return *this;
}

// compound assignment operators
MyClass & MyClass::operator+=( const MyClass &rhs )
{
    this->a += rhs.a;
    return *this;
}

MyClass & MyClass::operator-=( const MyClass &rhs )
{
    this->a -= rhs.a;
    return *this;
}

MyClass & MyClass::operator*=( const MyClass &rhs )
{
    this->a *= rhs.a;
    return *this;
}

MyClass & MyClass::operator/=( const MyClass &rhs )
{
    this->a /= rhs.a;
    return *this;
}

const MyClass MyClass::operator+( const MyClass &other ) const
{
    // MyClass( *this ) -- makes a copy of myself
    // and then destructive compound adding
    return MyClass( *this ) += other;
}

const MyClass MyClass::operator-( const MyClass &other ) const
{
    return MyClass( *this ) -= other;
}

const MyClass MyClass::operator*( const MyClass &other ) const
{
    return MyClass( *this ) *= other;
}

const MyClass MyClass::operator/( const MyClass &other ) const 
{
    return MyClass( *this ) /= other;
}

bool MyClass::operator==( const MyClass &other ) const 
{
    return this->a == other.a; 
}

bool MyClass::operator!=( const MyClass &other ) const
{
    return !(*this == other );
}

MyClass::MyClass( void )
{
}

MyClass::MyClass( double var )
{
    a = var;
}

double MyClass::getA( void )
{
    return a;
}

void MyClass::setA( double var )
{
    a = var;
}

int main( void )
{
    MyClass mc1( 5.0 );
    MyClass mc2( 7.0 );
    cout << "Creating mc1: mc1.a = " << mc1.getA() << endl;
    cout << "Creating mc2: mc2.a = " << mc2.getA() << endl;

    cout << "Creating empty mc3 class." << endl << endl;
    MyClass mc3;
    
    cout << "mc3 = mc2: " << endl;
    mc3 = mc2;
    cout << "mc3.a: " << mc3.getA() << endl << endl; 
    
    cout << "mc3 += mc1: " << endl;
    mc3 += mc1;
    cout << "mc3.a: " << mc3.getA() << endl << endl;

    cout << "mc3 -= mc1: " << endl;
    mc3 -= mc1;
    cout << "mc3.a: " << mc3.getA() << endl << endl;

    cout << "mc3 *= mc2: " << endl;
    mc3 *= mc2;
    cout <<" mc3.a: " << mc3.getA() << endl << endl;
    
    mc3 = mc1 + mc2;
    cout << "mc3 = mc1 + mc2:" << endl;
    cout << "mc3.a: " << mc3.getA() << endl << endl;

    mc3 = mc1 - mc2;
    cout << "mc3 = mc1 - mc2:" << endl;
    cout << "mc3.a: " << mc3.getA() << endl << endl;
    
    mc3 = mc1 * mc2;
    cout << "mc3 = mc1 * mc2: " << endl;
    cout << "mc3.a: " << mc3.getA() << endl << endl;

    mc3 = mc1 / mc2;
    cout << "mc3 = mc1 / mc2: " << endl;
    cout << "mc3.a: " << mc3.getA() << endl << endl;

    bool r = mc1 == mc3;
    bool nr = mc1 != mc3;
    cout << "mc1 == mc3: " << r << "; mc1 != mc3: " << nr << endl << endl;
    
    mc1 = mc3;
    r = mc1 == mc3;
    nr = mc1 != mc3;
    cout << "Assigning mc1 to mc3. mc1 == mc3: " << r << "; mc1 != mc3: " << nr << endl << endl;
        
    return 0;
}
