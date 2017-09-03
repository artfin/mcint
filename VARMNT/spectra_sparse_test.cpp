#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <GenEigsSolver.h>
#include <MatOp/SparseGenMatProd.h>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace Spectra;

int main()
{
    // a band matrix with 1 on the main diagonal, 2 on the below-main subdiagonal, and 3 on the above-main subdiagonal
    
    const int n = 10;
    SparseMatrix<double> M(n, n);
    M.reserve( VectorXi::Constant(n, 3) );

    for ( int i = 0; i < n; i++ )
    {
        M.insert(i, i) = 1.0;
        if ( i > 0 )
        {
            M.insert( i - 1, i ) = 3.0;
        }
        if ( i < n - 1 )
        {
            M.insert( i + 1, i ) = 2.0;
        }
    }

    // construct matrix operation object
    SparseGenMatProd<double> op(M);

    // construct eigen solver object, request 3 smallest eigenvalues
    GenEigsSolver< double, SMALLEST_MAGN, SparseGenMatProd<double> > eigs(&op, 3, 6);

    eigs.init();
    int nconv = eigs.compute();

    VectorXcd evalues;
    if ( eigs.info() == SUCCESSFUL )
    {
        evalues = eigs.eigenvalues();
    }

    cout << "Eigenvalues found: " << endl << evalues << endl;
    
    return 0;
}
