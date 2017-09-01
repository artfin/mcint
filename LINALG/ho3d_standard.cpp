/*
  Solves the one-particle Schrodinger equation
  for a potential specified in function
  potential(). 
*/
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "lib.h"
using namespace  std;
// output file as global variable
ofstream ofile;  

// function declarations 
void  **matrix(int, int, int);
void free_matrix(void **);
void initialise(double&, double&, double&, int&, int&) ;
double potential(double, double);
int comp(const double *, const double *);
void output(double, double, int, double *);
void tqli(double *, double *, int, double **);
double pythag(double, double);

// Begin main function
int main(int argc, char* argv[])
{
  int       i, j, max_step, orb_l;
  double    r_min, r_max, step, const_1, const_2, orb_factor, 
    *e, *d, *w, *r, **z, gamma;
  char *outfilename;
  // Read in output file, abort if there are too few command-line arguments
  if( argc <= 1 ){
    cout << "Bad Usage: " << argv[0] << 
      " read also output file on same line" << endl;
    exit(1);
  }
  else{
    outfilename=argv[1];
  }
  ofile.open(outfilename); 
  //   Read in data 
  initialise(r_min, r_max, gamma, orb_l, max_step);
  // initialise constants
  step    = (r_max - r_min) / max_step; 
  const_1 =  1.0 / (step * step);
  const_2 = 0.5 / (step * step);
  orb_factor = orb_l * (orb_l + 1);
  
  // local memory for r and the potential w[r] 
  r = new double[max_step + 1];
  w = new double[max_step + 1];
  for(i = 0; i <= max_step; i++) {
    r[i] = r_min + i * step;
    w[i] = potential(r[i], gamma) + 0.5*orb_factor / (r[i] * r[i]);
  }
  // local memory for the diagonalization process 
  d = new double[max_step];    // diagonal elements 
  e = new double[max_step];    // tri-diagonal off-diagonal elements 
  z = (double **) matrix(max_step, max_step, sizeof(double));
  for(i = 0; i < max_step; i++) {
    d[i]    = const_1 + w[i + 1];
    e[i]    = const_2;
    z[i][i] = 1.0;
    for(j = i + 1; j < max_step; j++)  {
      z[i][j] = 1.0;
    }
  }
  // diagonalize and obtain eigenvalues
  tqli(d, e, max_step - 1, z);      
  // Sort eigenvalues as an ascending series 
  qsort(d,(UL) max_step - 1,sizeof(double),
         (int(*)(const void *,const void *))comp);
  // send results to ouput file
  output(r_min , r_max, max_step, d);
  delete [] r; delete [] w; delete [] e; delete [] d; 
  free_matrix((void **) z); // free memory
  ofile.close();  // close 
  return 0;
} // End: function main() 

/*
  The function potential()
  calculates and return the value of the 
  potential for a given argument x.
  The potential here is for the hydrogen atom
*/        

double potential(double x, double gamma)
{
  return 0.5*x*x;

} // End: function potential()  

/*
  The function   int comp()                  
  is a utility function for the library function qsort()
  to sort double numbers after increasing values.
*/       

int comp(const double *val_1, const double *val_2)
{
  if((*val_1) <= (*val_2))       return -1;
  else  if((*val_1) > (*val_2))  return +1;
  else                     return  0; 
} // End: function comp() 


void initialise(double& r_min, double& r_max, double& gamma, int& orb_l, int& max_step) 
{
  cout << "Min vakues of R = ";
  cin >> r_min;
  cout << "Max value of R = ";
  cin >> r_max;
  cout << "Value of gamma = ";
  cin >> gamma;
  cout << "Orbital momentum = ";
  cin >> orb_l;
  cout << "Number of steps = ";
  cin >> max_step;
}  // end of function initialise   




void output(double r_min , double r_max, int max_step, double *d)
{
  int i;
  ofile << "RESULTS:" << endl;
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile <<"R_min = " << setw(15) << setprecision(8) << r_min << endl;  
  ofile <<"R_max = " << setw(15) << setprecision(8) << r_max << endl;  
  ofile <<"Number of steps = " << setw(15) << max_step << endl;  
  ofile << "Five lowest eigenvalues:" << endl;
  for(i = 0; i < 5; i++) {
    ofile << setw(15) << setprecision(8) << d[i] << endl;
  }
}  // end of function output         



    /*
    ** The function
    **                 tqli()
    ** determine eigenvalues and eigenvectors of a real symmetric
    ** tri-diagonal matrix, or a real, symmetric matrix previously
    ** reduced by function tred2[] to tri-diagonal form. On input,
    ** d[] contains the diagonal element and e[] the sub-diagonal
    ** of the tri-diagonal matrix. On output d[] contains the
    ** eigenvalues and  e[] is destroyed. If eigenvectors are
    ** desired z[][] on input contains the identity matrix. If
    ** eigenvectors of a matrix reduced by tred2() are required,
    ** then z[][] on input is the matrix output from tred2().
    ** On output, the k'th column returns the normalized eigenvector
    ** corresponding to d[k]. 
    ** The function is modified from the version in Numerical recipe.
    */

void tqli(double *d, double *e, int n, double **z)
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
            if(iter++ == 30) {
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
               for(k = 0; k < n; k++) {
                  f         = z[k][i+1];
                  z[k][i+1] = s * z[k][i] + c * f;
                  z[k][i]   = c * z[k][i] - s * f;
               } /* end k-loop */
            } /* end i-loop */
            if(r == 0.0 && i >= l) continue;
            d[l] -= p;
            e[l]  = g;
            e[m]  = 0.0;
         } /* end if-loop for m != 1 */
      } while(m != l);
   } /* end l-loop */
} /* End: function tqli(), (C) Copr. 1986-92 Numerical Recipes Software )%. */
   
double pythag(double a, double b)
{
  double absa,absb;
  absa=fabs(a);
  absb=fabs(b);
  if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}
// End: function pythag(), (C) Copr. 1986-92 Numerical Recipes Software )%.

  /*
   * The function                             
   *      void  **matrix()                    
   * reserves dynamic memory for a two-dimensional matrix 
   * using the C++ command new . No initialization of the elements. 
   * Input data:                      
   *  int row      - number of  rows          
   *  int col      - number of columns        
   *  int num_bytes- number of bytes for each 
   *                 element                  
   * Returns a void  **pointer to the reserved memory location.                                
   */

void **matrix(int row, int col, int num_bytes)
  {
  int      i, num;
  char     **pointer, *ptr;

  pointer = new(nothrow) char* [row];
  if(!pointer) {
    cout << "Exception handling: Memory allocation failed";
    cout << " for "<< row << "row addresses !" << endl;
    return NULL;
  }
  i = (row * col * num_bytes)/sizeof(char);
  pointer[0] = new(nothrow) char [i];
  if(!pointer[0]) {
    cout << "Exception handling: Memory allocation failed";
    cout << " for address to " << i << " characters !" << endl;
    return NULL;
  }
  ptr = pointer[0];
  num = col * num_bytes;
  for(i = 0; i < row; i++, ptr += num )   {
    pointer[i] = ptr; 
  }

  return  (void **)pointer;

  } // end: function void **matrix()

    /*
     * The function                         
     *      void free_matrix()              
     * releases the memory reserved by the function matrix() 
     *for the two-dimensional matrix[][] 
     * Input data:                          
     *  void far **matr - pointer to the matrix
     */

void free_matrix(void **matr)
{

  delete [] (char *) matr[0];
  delete [] matr;

}  // End:  function free_matrix() 

