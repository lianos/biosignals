#include "biosignals.h"

SEXP biosignals_version(SEXP single_) {
    bool single = Rcpp::as<bool>(single_);
    
    if (single) {
        return Rcpp::wrap(10000 * BIOSIGNALS_MAJOR_VERSION +
                          100 * BIOSIGNALS_MINOR_VERSION +
                          BIOSIGNALS_PATCH);
    }
    
    return Rcpp::IntegerVector::create(
                Rcpp::_["major"] = BIOSIGNALS_MAJOR_VERSION,
                Rcpp::_["minor"] = BIOSIGNALS_MINOR_VERSION,
                Rcpp::_["patch"] = BIOSIGNALS_PATCH);
}


