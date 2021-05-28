#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(BH, bigmemory)]]
#include <bigmemory/Matrixaccessor.hpp>
#include <numeric>
// [[Rcpp::export]]
IntegerVector mutualRank(XPtr<BigMatrix> rankMat) {
    // iterators
    int k = 0;
    int i, j;
    // variables
    double rij, rji;
    // Create the matrix accessor
    MatrixAccessor<int> ma(*rankMat);
    // matrix dimensions
    int rows = rankMat->nrow();
    int cols = rankMat->ncol();
    // output
    IntegerVector mr( (rows * (rows - 1)) / 2 );

    for(i=0; i < rows - 1; i++){
        for(j = i + 1; j < cols; j++){
            // bigmemory matrix is accessed like this: [col][row]
            rij = ma[j][i]; // rank of j in i's list
            rji = ma[i][j]; // rank of i in j's list
            // geometric mean integer version to four decimals precision
            mr[k++] = (int) (std::sqrt(rij * rji) * 10000);
        }
    }
    return mr;
}
