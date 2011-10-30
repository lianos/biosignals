#ifndef __RCPP_BIOSIGNALS_WRAP_H__
#define __RCPP_BIOSIGNALS_WRAP_H__

namespace Rcpp {
    
    // Creates an IRanges::Rle object out of biosignals::Rle
    template <typename T>
    SEXP wrap(const biosignals::Rle<T>& rle) {
        S4 ans = S4("Rle");
        ans.slot("values") = wrap(rle.values);
        ans.slot("lengths") = wrap(rle.lengths);
        return ans;
    }
    
}

#endif