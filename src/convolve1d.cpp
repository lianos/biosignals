#include "convolve1d.h"
#include "biosignals/convolve1d.h"

SEXP Rconvolve_1d(SEXP x_, SEXP kernel_, SEXP rescale_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    std::vector<double> kernel = Rcpp::as<std::vector<double> >(kernel_);
    bool rescale = Rcpp::as<bool>(rescale_);
    std::vector<double> *result = biosignals::convolve_1d(x, kernel, rescale);
    return Rcpp::wrap(*result);
END_RCPP
}

SEXP Rfencepost_convolve_1d(SEXP x_, SEXP kernel_, SEXP starts_, SEXP ends_,
                            SEXP rescale_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as< std::vector<double> >(x_);
    std::vector<double> kernel = Rcpp::as< std::vector<double> >(kernel_);
    Rcpp::IntegerVector starts(starts_);
    Rcpp::IntegerVector ends(ends_);
    bool rescale = Rcpp::as<bool>(rescale_);
    Rcpp::NumericVector out(x.size());
    int start,end;

    for (int i = 0; i < starts.size(); i++) {
        start = starts[i] - 1; // R -> C indexing
        end = ends[i] - 1;
        std::vector<double> tmp = std::vector<double>(end - start + 1, 0.0);
        std::copy(x.begin() + start, x.begin() + end, tmp.begin());
        std::vector<double> *conv = biosignals::convolve_1d(tmp, kernel); //, rescale);

        std::copy(conv->begin(), conv->end(), out.begin() + start);
        delete conv;
    }

    return Rcpp::wrap(out);
END_RCPP
}

#if 0
SEXP
Rconvolve_rle(SEXP x_, SEXP kernel_, SEXP starts_, SEXP widths_,
              SEXP rescale_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    std::vector<double> kernel = Rcpp::as<std::vector<double> >(kernel_);
    std::vector<int> starts = Rcpp::as<std::vector<int> >(starts_);
    std::vector<int> widths = Rcpp::as<std::vector<int> >(widths_);
    bool rescale = Rcpp::as<bool>(rescale_);
        
    SEXP bounds, rstart, rwidth;
    SEXP out = NEW_LIST(starts.size());
    
    PROTECT(rstart = NEW_INTEGER(1));
    PROTECT(rwidth = NEW_INTEGER(1));
    
    for (int i = 0; i < starts.size(); i++) {
        INTEGER(rstart)[0] = starts[i];
        INTEGER(rwidth)[0] = widths[i];
        
        // Returns a list with $values and $lengths
        PROTECT(bounds = Rle_seqselect(x_, rstart, rwidth));
        std::vector<double> values(VECTOR_ELT(bounds, 0));
        std::vector<int> lengths(VECTOR_ELT(bounds, 1));
        std::vector<double> this_x = expand_rle(lengths, values);
        std::vector<double> *result = convolve_1d(&x, &kernel, rescale)
        SET_VECTOR_ELT(out, i, Rcpp::wrap(*result));
        UNPROTECT(1);
    }
    
    UNPROTECT(2);
    return out;
END_RCPP
}
#endif
