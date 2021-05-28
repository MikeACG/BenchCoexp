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

struct StatGW : public Worker {

    const RMatrix<int> M;
    int n, m;

    RVector<int> outv;

    StatGW(IntegerMatrix M, int n, int m, IntegerVector outv)
        : M(M), n(n), m(m), outv(outv) {}
    
    void operator()(std::size_t begin, std::size_t end) {
        int k = findBaseIndex(begin, n);
        int i, j, obs;
        double corr = 0.5;
        double expec, gtmp, g;
        std::vector<int> rowsums(4);
        std::vector<int> colsums(4);
        std::vector<std::vector<int>> table(4, std::vector<int>(4));
        for(std::size_t x = begin; x < end; x++){
            for(std::size_t y = x + 1; y < n; y++){

                for (i = 0; i < m; i++) {
                    table[M(x, i)][M(y, i)]++;
                    rowsums[M(x, i)]++;
                    colsums[M(y, i)]++;
                }

                g = 0;
                for (i = 1; i < 4; i++) {
                    for (j = 1; j < 4; j++) {
                        obs = table[i][j];
                        expec = (double) (rowsums[i] * colsums[j]) / (double) m;
                        if (expec < corr) expec = corr;
                        gtmp = (obs > 0) ? obs * std::log(obs / expec) : 0;
                        if (gtmp > 4 && expec < 1) {
                            expec = 1;
                            gtmp = obs * std::log(obs / expec);
                        }
                        g += gtmp;
                        table[i][j] = 0;
                    }
                }

                outv[k++] = (int) (2 * g * 100);
                for (i = 1; i < 4; i++) rowsums[i] = colsums[i] = 0; 
            }
        }
    }
};

// [[Rcpp::export]]
IntegerVector statG(IntegerMatrix M) {
    int n = M.nrow();
    int m = M.ncol();
    IntegerVector outv( (n * (n - 1)) / 2 );
    StatGW statGW(M, n, m, outv);
    parallelFor(0, n, statGW);
    return outv;
}
