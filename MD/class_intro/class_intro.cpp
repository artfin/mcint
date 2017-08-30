#include <iostream>

using namespace std;

class Box {

    double height;
    double length;
    double breadth;

    public:
        // member function declaration
        double getVolume( void );
       
        double getLength( void );
        double getHeight( void );
        double getBreadth( void );

        void setLength( double len );
        void setBreadth( double bre );
        void setHeight( double hei );
};

// using scope resolution operator ::
double Box::getVolume( void )
{
   return length * breadth * height; 
}

double Box::getLength( void )
{
    return length;
}

double Box::getHeight( void )
{
    return height;
}

double Box::getBreadth( void )
{
    return breadth;
}

void Box::setLength( double len )
{
    length = len;
}

void Box::setBreadth( double bre )
{
    breadth = bre;
}

void Box::setHeight( double hei )
{
    height = hei;
}
int main()
{
    Box box1;

    box1.setLength(6.0);
    box1.setBreadth(5.0);
    box1.setHeight(7.0);    

    double l = box1.getLength();
    cout << "Getting private class parameter using getter-function: " << l << endl;

    double volume1 = box1.getVolume();
    cout << "Calling volume_inside: " << volume1 << endl;

    return 0;
}
