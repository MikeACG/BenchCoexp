// [[Rcpp::depends(RcppParallel)]]
#include <Rcpp.h>
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

struct BenchmarkShort : public Worker {

    RVector<int> ca, ct, opps;
    int k, n;

    RMatrix<int> outm;

    BenchmarkShort(IntegerVector ca, IntegerVector ct, IntegerVector opps, int k, int n, IntegerMatrix outm)
        : ca(ca), ct(ct), opps(opps), k(k), n(n), outm(outm) {}
    
    void operator()(std::size_t begin, std::size_t end) {
        for (std::size_t o = begin; o < end; o++) {
            for (int i = 0; i < n; i++) {
                if (ca[i] <= opps[o]) {
                    if (ct[i] >= k) {
                        outm(o + o, 0)++;
                    } else {
                        outm(o + o, 1)++;
                    }
                } else {
                    if (ct[i] >= k) {
                        outm(o + o + 1, 0)++;
                    } else {
                        outm(o + o + 1, 1)++;
                    }
                }
            }
        }
    }
};

// [[Rcpp::export]]
IntegerMatrix benchShort(IntegerVector ca, IntegerVector ct, IntegerVector opps, int k) {
    int n = ca.length();
    int nsteps = opps.length();
    IntegerMatrix outm(2 * nsteps, 2);
    BenchmarkShort benchmarkShort(ca, ct, opps, k, n, outm);
    parallelFor(0, nsteps, benchmarkShort);
    return outm;
}
