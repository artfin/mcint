#include <iostream>

using namespace std;

class Line {
    public:
        void setLength( double len );
        void setWidth( double wid );

        double getLength( void );
        double getWidth( void );
        
        Line( double len, double wid ); // This is the constructor!
        ~Line(); // This is the destructor! 

    private:
        double length;
        double width;
};

Line::~Line( void )
{
    cout << "Object is being deleted..." << endl;
}

// member functions definitions including constructor
Line::Line( double len, double wid ): length(len), width(wid)
{
    cout << "Object is being constructed: " << endl;
    cout << ">> length = " << len << endl;
    cout << ">> width = " << wid << endl;
}

void Line::setLength( double len )
{
    length = len;
}

void Line::setWidth( double wid )
{
    width = wid;
}

double Line::getLength( void )
{
    return length;
}

double Line::getWidth( void )
{
    return width;
}

int main()
{
    Line line( 1.0, 1.0 );

    cout << "Initial value of length: " << line.getLength() << endl;
    cout << "Initial value of width: " << line.getWidth() << endl;
    
    // reset line length and width
    line.setLength( 6.0 );
    line.setWidth( 5.0 );
    cout << "Resetting member class parameters using setter functions..." << endl;
    cout << "New length of the line: " << line.getLength() << endl;
    cout << "New width of the line: " << line.getWidth() << endl;

    return 0;
}
