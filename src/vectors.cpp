#include "vectors.h"
#include "biosignals/vectors/utils.h"

std::pair< std::vector<int>, std::vector<double> >
as_rle(std::vector<double> &vals, double eps) {
    // std::vector<T> values;
    std::vector<double> values;
    std::vector<int> lengths;
    std::pair< std::vector<int>, std::vector<double> > out;
    
    // T current_val, tmp_val;
    double current_val, tmp_val;
    int current_len;
    
    if (vals.size() == 0) {
        out.first = lengths;
        out.second = values;
        return out;
    }
    
    current_val = vals[0];
    current_len = 1;
    
    for (int i = 1; i < vals.size(); i++) {
        tmp_val = vals[i];
        if (ALMOST_ZERO_EPS(tmp_val - current_val, eps)) {
            current_len++;
        } else {
            lengths.push_back(current_len);
            values.push_back(current_val);
            current_len = 1;
            current_val = tmp_val;
        }
    }
    
    if (values.back() != current_val) {
        values.push_back(current_val);
        lengths.push_back(current_len);
    }
    
    out.first = lengths;
    out.second = values;
    return out;
}

std::vector<double>
expand_rle(std::vector<int> lengths, std::vector<double> &vals) {
    int length = 0;
    int i,j,len,sofar;
    // T val;
    double val;
    
    for (i = 0; i < lengths.size(); i++) {
        length += lengths[i];
    }
    
    // std::vector<T> out = std::vector<T>(length);
    std::vector<double> out = std::vector<double>(length);
    sofar = 0;
    for (i = 0; i < lengths.size(); i++) {
        val = vals[i];
        len = lengths[i];
        for (j = 0; j < len; j++) {
            out[sofar++] = val;
        }
    }
    
    return out;
}

// ----------------------------------------------------------------- Rinterface
SEXP
Rzero_crossings(SEXP x_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as< std::vector<double> >(x_);
    std::vector<int> z = biosignals::zero_crossings(x);

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
    return Rcpp::wrap(biosignals::sliding_window_maximum(x, k));
END_RCPP
}

SEXP
Rsliding_min(SEXP x_, SEXP k_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    int k = Rcpp::as<int>(k_);
    return Rcpp::wrap(biosignals::sliding_window_minimum(x, k));
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

    result = biosignals::coverage_quantiles(x, starts, ends, breaks);

    for (i = 0; i < starts.size(); i++) {
        for (j = 0; j < breaks.size(); j++) {
            // C --> R indexing
            out(i, j) = result[i][j] + 1;
        }
    }

    return Rcpp::wrap(out);
END_RCPP
}
