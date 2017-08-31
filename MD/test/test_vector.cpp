#include "vector.hpp"

using namespace std;

int main()
{
    // creating vector v1 with empty coords
    __vector__ v1 ( 3 );
    __vector__ v3 ( 3 );
    __vector__ v2 ( 1.0, 1.0, 2.0 );

    cout << "Assigning v1 = v2: " << endl;
    v1 = v2;
    v1.print();
    cout << endl;

    cout << "v1 += v2: " << endl;
    v1 += v2;
    v1.print();
    cout << endl;

    cout << "v1 -= v2: " << endl;
    v1 -= v2;
    v1.print();
    cout << endl;

    cout << "v1 += 5: " << endl;
    v1 += 5;
    v1.print();
    cout << endl;

    cout << "v2 -= 4: " << endl;
    v2 -= 4;
    v2.print();
    cout << endl;

    cout << "v1 *= 2: " << endl;
    v1 *= 2;
    v1.print();
    cout << endl;

    cout << "v1 /= 2: " << endl;
    v1 /= 2;
    v1.print();
    cout << endl;

    cout << "v3 = v1 + v2: " << endl;
    v3 = v1 + v2;
    v3.print();
    cout << endl;

    cout << "v3 = v1 - v2: " << endl;
    v3 = v1 - v2;
    v3.print();
    cout << endl;

    cout << "v3 = v1 + 5: " << endl;
    v3 = v1 + 5;
    v3.print();
    cout << endl;
    
    cout << "v3 = v1 * 3: " << endl;
    v3 = v1 * 3;
    v3.print();
    cout << endl;

    cout << "v3 = v1 / 2: " << endl;
    v3 = v1 / 2;
    v3.print();
    cout << endl;

    double r = v1 * v2;
    cout << "r = v1 * v2: " << r << endl << endl;

    cout << "v3 = power( v1, 2 ): " << endl;
    v3 = v1.power( 2 ); 
    v3.print();
    cout << endl;

    v1.setCoords(0, 0, 0);
    v2.setCoords(2, 2, 2);
    r = v1.dist( v2 );
    cout << "v1: (0, 0, 0)" << endl;
    cout << "v2: (2, 2, 2)" << endl;
    cout << "r = v1.dist( v2 ): " << r << endl << endl;

    cout << "Normalizing v2: " << endl;
    v2.normalize();
    v2.print();
    cout << endl;

    cout << "v2.fill_random(): " << endl;
    v2.fill_random();
    v2.print();
    cout << endl;

    return 0;
}
