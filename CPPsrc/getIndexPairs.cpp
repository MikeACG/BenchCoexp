#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector getIndexPairs(int idx, IntegerVector triVec, int n) {
    int consecStart, i, j;
    IntegerVector out(n);

    // find position where consecutive pairs for idx start
    consecStart = (idx * n); // if all indices had passed for all before idx
    //for (i = 1; i <= idx; i++) {
    //    consecStart -= i; // for the first index that has passed in reality only n - 1 pairs passed, for the 2nd n -2, etc.
    //}
    
    // find pairs left behind (non-consecutive)
    
    // of every index left behind, find where they're consecutives start
    // consecutive pairs of first index left behind always start at 0, so start from 1
    //for (i = 1; i < idx; i++) {
    //    out[i] = out[i - 1] + (n - i); // next index after pairs that pass for i
    //}

    // refactor last two loops in 1, doesnt matter if you modify out when i = idx
    for (i = 1; i <= idx; i++) {
        consecStart -= i;
        out[i] = out[i - 1] + (n - i);
    }

    // pinpoint index of idx as a pair of each of the left behind indices
    for (i = 0; i < idx; i++) {
        out[i] += idx - i - 1;
    }
    
    // fill the rest of the vector with the consecutive pairs indices, skip own pair
    for (j = i + 1; j < n; j++){
        out[j] = consecStart++;
    }

    // get the pairs
    for (i = 0; i < n; i++) {
        out[i] = triVec[out[i]];
    }

    // own pair is forced to 0
    out[idx] = 0;

    return out;
}
