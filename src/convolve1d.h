#ifndef __EDGE1D_H__
#define __EDGE1D_H__

#include "biosignals.h"

// ------------------------------------------------------------------- Internal
std::vector<double> *
convolve_1d(std::vector<double> x, std::vector<double> kernel, bool rescale);

// TODO: To convolve using FFTW, see:
// http://hellmund.blogspot.com/2011/02/example-using-rcpp-inline-and-fftw.html

// ------------------------------------------------------------------ Interface
RcppExport SEXP
Rconvolve_1d(SEXP x_, SEXP kernel_, SEXP rescale_);

RcppExport SEXP
Rfencepost_convolve_1d(SEXP x_, SEXP kernel_, SEXP starts_, SEXP ends_,
                       SEXP rescale_);

RcppExport SEXP
Rconvolve_rle(SEXP x_, SEXP kernel_, SEXP starts_, SEXP widths_,
              SEXP rescale_);
#endif