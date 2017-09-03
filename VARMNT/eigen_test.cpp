#include <iostream>

#include <Eigen/Core>
#include <SymEigsSolver.h>

#include <ctime>

using namespace std;
using namespace Eigen;
using namespace Spectra;

int main()
{
    int SIZE = 10000;
    cout << "Creating matrix of size " << SIZE << "x" << SIZE << endl;

    MatrixXd A = MatrixXd::Random( SIZE, SIZE );
    MatrixXd M = A + A.transpose();

    // construct matrix operation object using the wrapper class
    DenseSymMatProd<double> op(M);

    // construct eigen solver object, requesting the smallest three eigenvalues
    SymEigsSolver<double, SMALLEST_MAGN, DenseSymMatProd<double> > eigs(&op, 3, 6);
    
    cout << "Computing 3 smallest eigenvalues..." << endl;
    
    clock_t start = clock();

    // initialize and compute
    eigs.init();
    int ncovs = eigs.compute();

    // retrieve results
    VectorXd evalues;
    if ( eigs.info() == SUCCESSFUL )
    {
        evalues = eigs.eigenvalues();
    }

    cout << "Eigenvalues found: " << endl << evalues << endl;
    cout << "Time elapsed: " << (double) (clock() - start) / CLOCKS_PER_SEC << "s" << endl;

    return 0;
}
