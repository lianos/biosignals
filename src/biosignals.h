#ifndef __BIOSIGNALS_H__
#define __BIOSIGNALS_H__

// #include <RcppCommon.h>
// // Override wrap and as?
// #include <Rcpp.h>
// 
// #define SIGN(x) ((x >= 0) ? 1 : -1)
// #define ALMOST_ZERO(x) ((-1e-8 < x) && (x < 1e-8))
// #define ALMOST_ZERO_EPS(x,eps) ((-eps < x) && (x < eps))
// #define ALMOST_EQ(x,y) ((SIGN(x-y) * (x-y)) < 1e-8)
// #define ALMOST_EQ_EPS(x,y,eps) ((SIGN(x-y) * (x-y)) < eps)

#include "RcppBiosignals.h"
#include "convolve1d.h"
#include "vectors.h"

RcppExport SEXP biosignals_version(SEXP single_);


#endif
