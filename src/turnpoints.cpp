#if 0

#include "turnpoints.h"

#include <cstdlib>
#include <deque>
#include <iostream>
#include <vector>


/**
 * Returns a pair representing (minima, maxima) vectors
 * with the position of each from x
 */
template <typename T>
std::pair< std::vector<int>, std::vector<int> >
local_turnpoints(std::vector<T> &x, double threshold, int K) {
    std::vector<int> mins;
    std::vector<int> maxs;
    std::pair< std::vector<int>, std::vector<int> > result;

    if (x.size() == 0) {
        return std::make_pair(mins, maxs);
    }

    T cmin = x[0];
    T cmax = x[0];
    T cval_max, cval_min;
    int fidx_min = 0, bidx_min = 0, cidx_min = 0;
    int fidx_max = 0, bidx_max = 0, cidx_max = 0;

    for (int i = 0; i < x.size(); i++) {
        // look for minima
        update_window_stats(x, bidx_max, fidx_max, cidx_max, cmax, K, i, true);
        cval_min = x[i];
        if (SIGN(cval_min) != SIGN(cmin)) {
            cval_min = -1 * cval_min;
        }
        if (cval_min / cmax >= threshold)
        // look for maxima
        update_window_stats(x, bidx_min, fidx_min, cidx_min, cmin, K, i, false);


        cval_max = x[i];
        if (SIGN(cval_max) != SIGN(cmax)) {
            cval_max = -1 * cval_max;
        }

    }

    result = std::make_pair(mins, maxs);
    return result;
}


// ----------------------------------------------------------------- R interface


SEXP
Rlocal_turnpoints(SEXP x_, SEXP threshold_, SEXP wlength_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    double threshold = Rcpp::as<double>(threshold_);
    int wlength = Rcpp::as<int>(wlength_);

    return R_NilValue;
END_RCPP
}

#endif
