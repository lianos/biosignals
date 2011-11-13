#ifndef __CONVOLVE1D_H__
#define __CONVOLVE1D_H__

#include "RcppBiosignals.h"

RcppExport SEXP
Rconvolve_1d(SEXP x_, SEXP kernel_, SEXP rescale_);

RcppExport SEXP
Rfencepost_convolve_1d(SEXP x_, SEXP kernel_, SEXP starts_, SEXP ends_,
                       SEXP rescale_);

RcppExport SEXP
Rconvolve_rle(SEXP x_, SEXP kernel_, SEXP starts_, SEXP widths_,
              SEXP rescale_, SEXP eps_);
#endif
