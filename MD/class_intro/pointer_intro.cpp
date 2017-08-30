#include <iostream>

using namespace std;

class Box
{
    public:
        // constructor definition
        Box( double l = 2.0, double b = 2.0, double h = 2.0 )
        {
            cout << "Constructor called. " << endl;
            length = l;
            breadth = b;
            height = h;
        }

        double Volume()
        {
            return length * breadth * height;
        }

    private:
        double length;
        double breadth;
        double height;
};

int main()
{
    Box box1( 3.3, 1.2, 1.5 );
    Box box2( 8.5, 6.0, 2.0 );

    Box *ptrBox; // pointer to a class
    
    // saving the address of the first object
    ptrBox = &box1;

    // accessing a member using member access operator ->
    cout << "Volume of box1: " << ptrBox->Volume() << endl;

    ptrBox = &box2;

    cout << "Volume of box2: " << ptrBox->Volume() << endl;

    return 0;
}
