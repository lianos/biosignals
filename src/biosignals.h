#ifndef __BIOSIGNALS_H__
#define __BIOSIGNALS_H__

#include <RcppCommon.h>
// Override wrap and as?
#include <Rcpp.h>

#define SIGN(x) ((x >= 0) ? 1 : -1)
#define ALMOST_ZERO(x) ((-1e-8 < x) && (x < 1e-8))
#define ALMOST_EQ(x,y) ((SIGN(x-y) * (x-y)) < 1e-8)

#endif
