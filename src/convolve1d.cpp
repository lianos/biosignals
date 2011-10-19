#include "convolve1d.h"

// TODO: You can use two separate vectors to smooth over padded heads/tails
//       and incorporate those into main convolution so you don't have to
//       erase from the head of `result` on the "way out" of this function.
std::vector<double> *
convolve_1d(std::vector<double> *x, std::vector<double> *kernel, bool rescale) {
  double max_x, max_result, scale_factor;
  int klen = kernel->size();
  int n = x->size() + (2 * klen);
  int N = n + kernel->size(); // used to pad head/tail
  int i,j;
  double val;
  
  std::vector<double> *result = new std::vector<double>(N, 0.0);
  Rprintf("result vector of size %d is allocated\n", result->size());
  max_x = (*x)[0];
  max_result = 0;

  // Dodge head edge effects
  for (i = 0; i < kernel->size(); i++) {
    for (j = 0; j < kernel->size(); j++) {
      val = (*result)[i+j] + (*x)[0] * (*kernel)[j];
      (*result)[i+j] = val;
    }
  }
  
  // Main convolution
  for (i = klen; i < x->size() + klen; i++) {
    for (j = 0; j < klen; j++) {
      val = (*result)[i+j] + (*x)[i - klen] * (*kernel)[j];
      (*result)[i+j] = val;
      if (val > max_result) {
        max_result = val;
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
  
  scale_factor = max_x / max_result;
  
  /* remove head/tail padding */
  int naxe = 1.5 * klen;
  result->erase(result->begin(), result->begin() + naxe);
  result->erase(result->begin() + x->size(), result->end());
  return result;
}

SEXP Rconvolve_1d(SEXP x_, SEXP kernel_, SEXP rescale_) {
  std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
  std::vector<double> kernel = Rcpp::as<std::vector<double> >(kernel_);
  bool rescale = Rcpp::as<bool>(rescale_);
  std::vector<double> *result = convolve_1d(&x, &kernel, rescale);
  return Rcpp::wrap(*result);
}

SEXP Rfencepost_convolve_1d(SEXP x_, SEXP kernel_, SEXP starts_, SEXP ends_) {
  std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
  std::vector<double> kernel = Rcpp::as<std::vector<double> >(kernel_);
  Rcpp::IntegerVector starts(starts_);
  Rcpp::IntegerVector ends(ends_);
  
  return R_NilValue;
}
