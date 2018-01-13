#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>

using namespace std;

const double ABS_ERROR = 10e-9;
const double REL_ERROR = 10e-9;

template<typename Functor>
double midpnt( Functor func, double a, double b, int n )
{
    double x, sum, h, dh;
    static double s;
    int it = 1, j;

    if ( n == 1 )
        return (s = (b - a) * func(0.5 * (a + b)));
    else
    {
        // number of iterations == 3**(n - 1)
        for ( j = 1; j < n - 1; j++ ) it *= 3;
        
        h = (b - a)/(3.0 * it);
        dh = h + h; 
        x = a + 0.5 * h;

        sum = 0.0;
        for ( j = 1; j <= it; j++ )
        {
            sum += func(x);
            x += dh;
            sum += func(x);
            x += h;
        }
        s = (s + (b-a)*sum/it)/3.0;

        return s;
    }
}

int main( int argc, char* argv[] )
{
    if ( argc != 6 )
    {
        cout << "USAGE: ./... (int) a0, a1, a2, a3, a4" << endl;
        exit( 1 );
    } 

    double a0 = atof( argv[1] );
    double a1 = atof( argv[2] );
    double a2 = atof( argv[3] );
    double a3 = atof( argv[4] );
    double a4 = atof( argv[5] );

    auto integrand = [&](double y){ 
        double ypi = M_PI / 2 * y;
        double x = tan( ypi );
        double cos_value = cos( ypi );
        double p = a0 + a1 * x + a2 * x*x + a3 * x*x*x + a4 * x*x*x*x; 
        return M_PI/2/p/cos_value/cos_value;
    };

    clock_t start = clock();

    double res, last_res;
    int i = 1;
    double abs_error, rel_error;
    while (true)
    {
        res = midpnt( integrand, -1.0, 1.0, i );
        if ( res == 1 ) last_res = res;
        
        abs_error = abs( res - last_res );
        rel_error = abs_error / res;

        cout << setprecision(10) << "res: " << res
             << "; abs_error: " << abs_error
             << "; rel_error: " << rel_error << endl;
        
        last_res = res;
        i++;

        if ( abs_error < ABS_ERROR || rel_error < REL_ERROR ) 
          break;  
    }

    cout << "Time elapsed: " << (clock() - start) / (double) CLOCKS_PER_SEC << " s" << endl;

    return 0;
}
