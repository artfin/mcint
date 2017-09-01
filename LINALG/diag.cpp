#include <math.h>

#include <iostream>
#include <random>

#include <ctime>

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

void tqli( double *d, double *e, int n )
{
   register int   m,l,iter,i,k;
   double         s,r,p,g,f,dd,c,b;

   for(i = 1; i < n; i++) e[i-1] = e[i];
     e[n] = 0.0;
   for(l = 0; l < n; l++) {
      iter = 0;
      do {
         for(m = l; m < n-1; m++) {
            dd = fabs(d[m]) + fabs(d[m+1]);
            if((double)(fabs(e[m])+dd) == dd) break;
         }
         if(m != l) {
            if(iter++ == 50) {
               printf("\n\nToo many iterations in tqli.\n");
               exit(1);
            }
            g = (d[l+1] - d[l])/(2.0 * e[l]);
            r = pythag(g,1.0);
            g = d[m]-d[l]+e[l]/(g+SIGN(r,g));
            s = c = 1.0;
            p = 0.0;
            for(i = m-1; i >= l; i--) {
               f      = s * e[i];
               b      = c*e[i];
               e[i+1] = (r=pythag(f,g));
               if(r == 0.0) {
                  d[i+1] -= p;
                  e[m]    = 0.0;
                  break;
               }
               s      = f/r;
               c      = g/r;
               g      = d[i+1] - p;
               r      = (d[i] - g) * s + 2.0 * c * b;
               d[i+1] = g + (p = s * r);
               g      = c * r - b;
            
            	/*  
               for(k = 0; k < n; k++) {
                  f         = z[k][i+1];
                  z[k][i+1] = s * z[k][i] + c * f;
                  z[k][i]   = c * z[k][i] - s * f;
               } /* end k-loop 
               */
            
            } /* end i-loop */
            if(r == 0.0 && i >= l) continue;
            d[l] -= p;
            e[l]  = g;
            e[m]  = 0.0;
         } /* end if-loop for m != 1 */
      } while(m != l);
   } /* end l-loop */
} /* End: function tqli(), (C) Copr. 1986-92 Numerical Recipes Software )%. */
   
int main()
{
	const int n = 5;
    // diagonal elements
    double d[ n ] = { 4.5, 4.125, 4.0, 4.125, 4.5 };
	// parallel diagonal elements
	double e[ n ] = { 0.0, -2.0, -2.0, -2.0, -2.0 };	

	tqli( d, e, n );

	for ( int i = 0; i < n; i++ )
	{
		cout << "d[" << i << "]: " << d[i] << endl;
	}

	random_device rd1;
	random_device rd2;
	mt19937 eng1( rd1() );
	mt19937 eng2( rd2() );
	uniform_real_distribution<> distr1( 100, 101 );
	uniform_real_distribution<> distr2( 5, 6 );

	const int SIZE = 50000;
	double d2[ SIZE ];
	double e2[ SIZE ];
	for ( int i = 0; i < SIZE; i++ )
	{
		d2[i] = distr1( eng1 );
		e2[i] = distr2( eng2 );
	}

	cout << "Size of input tridiagonal matrix: " << SIZE << endl;

	clock_t start = clock();
	tqli( d2, e2, SIZE );
	cout << "Time elapsed: " << ( clock() - start ) / (double) CLOCKS_PER_SEC << "s" << endl;

	
/*
Eigenvalues:
0: 0.674115
1: 2.30373
2: 4.33237
3: 6.32127
4: 7.61851
*/
    return 0;
}

