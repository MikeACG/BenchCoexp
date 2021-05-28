// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

inline int findBaseIndex(std::size_t x, int n) {
    int havePassed = x * n;
    for(std::size_t i = 1; i <= x; i++){
        havePassed -= i;
    }
    return havePassed;
}

struct AndSim : public Worker {

    const RMatrix<int> M;
    int n, m;

    RVector<int> outv;

    AndSim(IntegerMatrix M, int n, int m, IntegerVector outv)
        : M(M), n(n), m(m), outv(outv) {}
    
    void operator()(std::size_t begin, std::size_t end) {
        int k = findBaseIndex(begin, n);
        for(std::size_t x = begin; x < end; x++){
            for(std::size_t y = x + 1; y < n; y++){
                for(int i = 0; i < m; i++){
                    if (M(x, i) == 1 && M(y, i) == 1) outv[k]++;
                }
                k++;
            }
        }
    }
};

// [[Rcpp::export]]
IntegerVector termSimil(IntegerMatrix M) {
    int n = M.nrow();
    int m = M.ncol();
    IntegerVector outv( (n * (n - 1)) / 2 );
    AndSim andSim(M, n, m, outv);
    parallelFor(0, n, andSim);
    return outv;
}
