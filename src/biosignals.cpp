#include "biosignals.h"

SEXP biosignals_version(SEXP single_) {
    bool single = Rcpp::as<bool>(single);
    
    if (single) {
        return Rcpp::wrap(10000 * BIOSIGNALS_MAJOR_VERSION +
                          100 * BIOSIGNALS_MINOR_VERSION +
                          BIOSIGNALS_PATCH);
    }
    
    return Rcpp::IntegerVector::create(
                _["major"] = BIOSIGNALS_MAJOR_VERSION,
                _["minor"] = BIOSIGNALS_MINOR_VERSION,
                _["patch"] = BIOSIGNALS_PATCH);
}


// ----------------------------------------------------------------------------
// Convolutions


// ----------------------------------------------------------------------------
// vectors
SEXP
Rzero_crossings(SEXP x_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as< std::vector<double> >(x_);
    std::vector<int> z = zero_crossings(x);
    // convert to R coords
    for (int i = 0; i < z.size(); i++) z[i]++;
    return Rcpp::wrap(z);
END_RCPP
}

SEXP
Rsliding_max(SEXP x_, SEXP k_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    int k = Rcpp::as<int>(k_);
    return Rcpp::wrap(sliding_window_maximum(x, k));
END_RCPP
}

SEXP
Rsliding_min(SEXP x_, SEXP k_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    int k = Rcpp::as<int>(k_);
    return Rcpp::wrap(sliding_window_minimum(x, k));
END_RCPP
}

SEXP
Rcoverage_quantiles(SEXP x_, SEXP starts_, SEXP ends_, SEXP breaks_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as< std::vector<double> >(x_);
    std::vector<int> starts = Rcpp::as< std::vector<int> >(starts_);
    std::vector<int> ends = Rcpp::as< std::vector<int> >(ends_);
    std::vector<double> breaks = Rcpp::as< std::vector<double> >(breaks_);
    std::vector< std::vector<int> > result;
    Rcpp::NumericMatrix out(starts.size(), breaks.size());
    int i,j;

    // R --> C indexing
    for (i = 0; i < starts.size(); i++) starts[i] -= 1;
    for (i = 0; i < ends.size(); i++) ends[i] -= 1;

    result = coverage_quantiles(x, starts, ends, breaks);

    for (i = 0; i < starts.size(); i++) {
        for (j = 0; j < breaks.size(); j++) {
            // C --> R indexing
            out(i, j) = result[i][j] + 1;
        }
    }

    return Rcpp::wrap(out);
END_RCPP
}
