#include "convolve1d.h"
#include "vectors.h"

// TODO: You can use two separate vectors to smooth over padded heads/tails
//       and incorporate those into main convolution so you don't have to
//       erase from the head of `result` on the "way out" of this function.
// TODO: Convert convolve code to use FFTs
std::vector<double> *
convolve_1d(std::vector<double> *x, std::vector<double> *kernel, bool rescale) {
    double max_x, max_result, scale_factor;
    int klen = kernel->size();
    int n = x->size() + (2 * klen);
    int N = n + kernel->size(); // used to pad head/tail
    int i,j;
    double val;

    std::vector<double> *result = new std::vector<double>(N, 0.0);
    max_x = abs((*x)[0]);
    max_result = -1 * std::numeric_limits<double>::max();

    // Dodge head edge effects
    for (i = 0; i < kernel->size(); i++) {
        for (j = 0; j < kernel->size(); j++) {
            val = (*result)[i+j] + (*x)[0] * (*kernel)[j];
            (*result)[i+j] = val;
        }
    }

    // Main convolution
    for (i = klen; i < x->size() + klen; i++) {
        if (abs((*x)[i - klen]) > max_x) {
            max_x = abs((*x)[i - klen]);
        }
        for (j = 0; j < klen; j++) {
            val = (*result)[i+j] + (*x)[i - klen] * (*kernel)[j];
            (*result)[i+j] = val;
            if (abs(val) > max_result) {
                max_result = abs(val);
            }
        }
    }

    // Dodge tail edge effects
    for (i = x->size() + klen; i < n; i++) {
        for (j = 0; j < klen; j++) {
            val = (*result)[i+j] + (*x)[x->size() - 1] * (*kernel)[j];
            (*result)[i+j] = val;
        }
    }

    // remove head/tail padding
    int naxe = 1.5 * klen;
    result->erase(result->begin(), result->begin() + naxe);
    result->erase(result->begin() + x->size(), result->end());

    if (ALMOST_ZERO(max_x) || ALMOST_ZERO(max_result)) {
        scale_factor = 0;
    } else {
        scale_factor = max_x / max_result;
    }

    if (rescale) {
        for (i = 0; i < result->size(); i++) {
            (*result)[i] *= scale_factor;
        }
    }

    return result;
}

SEXP Rconvolve_1d(SEXP x_, SEXP kernel_, SEXP rescale_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    std::vector<double> kernel = Rcpp::as<std::vector<double> >(kernel_);
    bool rescale = Rcpp::as<bool>(rescale_);
    std::vector<double> *result = convolve_1d(&x, &kernel, rescale);
    return Rcpp::wrap(*result);
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
        std::vector<double> *conv = convolve_1d(&tmp, &kernel, rescale);

        std::copy(conv->begin(), conv->end(), out.begin() + start);
        delete conv;
    }

    return Rcpp::wrap(out);
END_RCPP
}
