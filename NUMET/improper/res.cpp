#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include <iomanip>
#include <vector>

using namespace std;

const double ABS_ERROR = 10e-9;
const double REL_ERROR = 10e-9;

const double LBOUND = -1.0;
const double UBOUND = 1.0;
const double LEN = 2.0;

template<typename Functor>
double midpnt( Functor func, const int& n )
{
    double x, sum, h, dh;
    static double s;
    int it = 1, j;

    if ( n == 1 )
        return (s = LEN * func(0));
    else
    {
        // number of iterations == 3**(n - 1)
        for ( j = 1; j < n - 1; j++ ) it *= 3;

        h = LEN/(3.0 * it);
        dh = h + h; 
        x = LBOUND + 0.5 * h;

        sum = 0.0;
        for ( j = 1; j <= it; j++ )
        {
            sum += func(x);
            x += dh;
            sum += func(x);
            x += h;
        }
        s = (s + LEN*sum/it)/3.0;

        return s;
    }
}

int main( int argc, char* argv[] )
{
    int coeffs[5];
    int index = 0, tmp;
    while( (index < 5) && scanf("%d ", &tmp) == 1)
    {
        coeffs[index] = tmp;
        index++;
    }

    auto integrand = [&](const double y){ 
        double ypi = M_PI / 2 * y;
        double sin_value = sin( ypi );
        double cos_value = cos( ypi );
        double x = sin_value / cos_value;
        double x2 = x * x;
        double x3 = x2 * x;
        double x4 = x2 * x2;
        double p = coeffs[0] + coeffs[1] * x + coeffs[2] * x2 + coeffs[3] * x3 + coeffs[4] * x4; 
        return M_PI / 2 / p / cos_value / cos_value;
    };

    double res, last_res;
    int i = 1;
    double abs_error, rel_error;
    cout << setprecision(10);
    while ( true )
    {
        res = midpnt( integrand, i );
        if ( i == 1 ) last_res = res;
        
        abs_error = abs( res - last_res );
        rel_error = abs_error / res;
        
        last_res = res;
        i++;

        //cout << "i: " << i
        //     << "; res: " << res << endl;

        if ( abs_error < ABS_ERROR && rel_error < REL_ERROR && i > 2 ) 
            break;  
    }

    cout << setprecision(10) << res << endl;

    return 0;
}
