#include "convolve1d.h"
#include "biosignals/convolve1d.h"
#include "IRanges_interface.h"

SEXP Rconvolve_1d(SEXP x_, SEXP kernel_, SEXP rescale_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    std::vector<double> kernel = Rcpp::as<std::vector<double> >(kernel_);
    bool rescale = Rcpp::as<bool>(rescale_);
    // std::vector<double> *result = biosignals::convolve_1d(x, kernel, rescale);
    // return Rcpp::wrap(*result);
    std::vector<double> result = biosignals::convolve_1d(x, kernel, rescale);
    return Rcpp::wrap(result);
END_RCPP
}

/**
 * Convolves the specified reanges (starts, widths) of an Rle.
 * 
 * This function utilized the IRanges,seqselect_Rle function, which
 * assumes indices are R-based (ie, start at 1) -- all position arithmetic
 * here is done in 1's index. Try not to get confused.
 */
SEXP
Rconvolve_rle(SEXP x_, SEXP kernel_, SEXP starts_, SEXP widths_,
              SEXP rescale_, SEXP eps_) {
    using namespace biosignals;
BEGIN_RCPP
    biosignals::Rle<double> x = Rcpp::as< Rle<double> >(x_);
    std::vector<double> kernel = Rcpp::as< std::vector<double> >(kernel_);
    std::vector<int> starts = Rcpp::as< std::vector<int> >(starts_);
    std::vector<int> widths = Rcpp::as< std::vector<int> >(widths_);
    double eps = Rcpp::as<double>(eps_);
    bool rescale = Rcpp::as<bool>(rescale_);
    int start, width, at = 1;
    SEXP tmp_rle;
    
    std::vector<double> vals = std::vector<double>();
    vals.reserve(x.values.size());
    std::vector<int> lens = std::vector<int>();
    lens.reserve(x.values.size());
    Rle<double> rle;
    
    for (int i = 0; i < starts.size(); i++) {
        start = starts[i];
        if (at < start) {
            // Fill in the regions that are not convolved with values from x
            width = start - at;
            PROTECT(tmp_rle = seqselect_Rle(x_, &at, &width, 1));
            rle = Rcpp::as< Rle<double> >(tmp_rle);
            UNPROTECT(1);
            
            vals.reserve(vals.size() +
                std::distance(rle.values.begin(), rle.values.end()));
            vals.insert(vals.end(), rle.values.begin(), rle.values.end());
            lens.reserve(lens.size() +
                std::distance(rle.lengths.begin(), rle.lengths.end()));
            lens.insert(lens.end(), rle.lengths.begin(), rle.lengths.end());
            at = std::accumulate(rle.lengths.begin(), rle.lengths.end(), at);
        }
        width = widths[i];
        
        PROTECT(tmp_rle = seqselect_Rle(x_, &start, &width, 1));
        rle = Rcpp::as< Rle<double> >(tmp_rle);
        UNPROTECT(1);
        
        rle = convolve_1d(rle, kernel, rescale);
        
        vals.insert(vals.end(), rle.values.begin(), rle.values.end());
        lens.insert(lens.end(), rle.lengths.begin(), rle.lengths.end());
        at += width;
    }
    
    if (at < x.size()) {
        width = x.size() - at + 1;
        PROTECT(tmp_rle = seqselect_Rle(x_, &at, &width, 1));
        rle = Rcpp::as< Rle<double> >(tmp_rle);
        UNPROTECT(1);
        
        vals.reserve(vals.size() +
            std::distance(rle.values.begin(), rle.values.end()));
        vals.insert(vals.end(), rle.values.begin(), rle.values.end());
        lens.reserve(lens.size() +
            std::distance(rle.lengths.begin(), rle.lengths.end()));
        lens.insert(lens.end(), rle.lengths.begin(), rle.lengths.end());
    }
    
    Rcpp::S4 ans("Rle");
    ans.slot("values") = Rcpp::wrap(vals);
    ans.slot("lengths") = Rcpp::wrap(lens);
    return ans;
END_RCPP
}
