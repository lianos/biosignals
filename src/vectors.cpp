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

// ----------------------------------------------------------------- Rinterface

SEXP
Ras_rle(SEXP vals_, SEXP eps_) {
BEGIN_RCPP    
    // switch(TYPEOF(x)) {
    // case LGLSXP:
    //     PROTECT(ans = Rle_logical_constructor(x, counts));
    //     break;
    // case INTSXP:
    //     PROTECT(ans = Rle_integer_constructor(x, counts));
    //     break;
    // case REALSXP:
    //     PROTECT(ans = Rle_real_constructor(x, counts));
    //     break;
    // case CPLXSXP:
    //     PROTECT(ans = Rle_complex_constructor(x, counts));
    //     break;
    // case STRSXP:
    //     PROTECT(ans = Rle_string_constructor(x, counts));
    //     break;
    // case RAWSXP:
    //     PROTECT(ans = Rle_raw_constructor(x, counts));
    //     break;
    // default:
    //     error("Rle computation of these types is not implemented");
    // }
    SEXP ret;
    double eps = Rcpp::as<double>(eps_);
    std::vector<double> vals = Rcpp::as< std::vector<double> >(vals_);
    std::pair< std::vector<int>, std::vector<double> > out;
    out = as_rle(vals, eps);

    ret = Rcpp::List::create(Rcpp::Named("lengths", Rcpp::wrap(out.first)),
                             Rcpp::Named("values", Rcpp::wrap(out.second)));
    return ret;
END_RCPP
}

SEXP
Rexpand_rle(SEXP lengths_, SEXP vals_) {
BEGIN_RCPP
    std::vector<int> lengths = Rcpp::as< std::vector<int> >(lengths_);
    std::vector<double> vals = Rcpp::as< std::vector<double> >(vals_);
    biosignals::Rle<double> rle = biosignals::Rle<double>(vals, lengths);
    std::vector<double> out = rle.expand();
    return Rcpp::wrap(out);
END_RCPP
}


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

