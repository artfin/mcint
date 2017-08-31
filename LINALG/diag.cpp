#include <math.h>

#include <iostream>

using namespace std;

// why static ? 
static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

float pythag( float a, float b )
{
    float absa = fabs(a);
    float absb = fabs(b);

    if ( absa > absb )
        return absa * sqrt( 1.0 + SQR( absb / absa ));
    else
        return ( absb == 0.0 ? 0.0 : absb * sqrt( 1.0 + SQR( absa / absb )));
}

void tqli( float d[], float e[], int n )
/*
 * Input:
 *      d[ 1 .. n ] -- containts the diagonal elements of the tridiagonal matrix
 *          on output, it returns the eigenvalues
 *      e[ 1 .. n ] -- inputs the subdiagonal elements of the tridiagonal matrix, with
 *          on output e is destroyed
 */
{
    float pythag( float a, float b );
    int i, m, l, iter, k;
    float s, r, p, g, f, dd, c, b;

    // it is convenient to renumber the elements of e
    for ( i = 2; i <= n; i++ )
    {
        e[i - 1] = e[i];
    }
    e[n] = 0.0;
    
    for ( int l = 1; l <= n; l++ )
    {
        iter = 0;
        
        do
        {
            for ( int m = 1; m <= n - 1; m++ )
            {
                dd = fabs( d[m] ) + fabs( d[m + 1] );
                if ( (float) ( fabs( e[m] ) + dd ) == dd ) 
                {
                    break;
                }
            }

            if ( m != l )
            {
                if ( iter++ == 30 )
                {
                    cout << "Too many iterations in TQLI!" << endl;
                }
                
                // form shift
                g = ( d[l + 1] - d[l] ) / ( 2.0 * e[l] );
                r = pythag( g, 1.0 );
                
                // this is d_m - k_s
                g = d[m] - d[l] + e[l] / (g + SIGN( r, g ));
                s = c = 1.0;
                p = 0.0;
                
                // a plane rotation as in the original QL, followed by Givens rotation to restore 
                // tridiagonal form
                for ( i = m-1; i >= l; i-- )
                {
                    f = s * e[i];
                    b = c * e[i];
                    e[i + 1] = (r = pythag( f, g ));

                    if ( r == 0.0 )
                    {
                        d[i + 1] -= p;
                        e[m] = 0.0;
                        break;
                    }

                    s = f / r;
                    c = g / r;
                    g = d[i + 1] - p;
                    r = ( d[i] - g ) * s + 2.0 * c * b;
                    d[i + 1] = g + ( p = s * r );
                    g = c * r - b; 
                }  
           
                if ( r == 0.0 && i >= l ) 
                {
                    continue;
                }

                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        } while ( m != l );
    }
}

int main()
{
    cout << "hey" << endl;

    // diagonal elements
    float d[3] = 

    return 0;
}

