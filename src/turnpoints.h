#if 0
#ifndef __TURNPOINTS_H__
#define __TURNPOINTS_H__

// #include "biosignals.h"
#include "biosignals/common.h"
#include "vectors.h"

template <typename T>
std::pair< std::vector<int>, std::vector<int> >
local_turnpoints(std::vector<T> &x, double threshold, int K);

RcppExport SEXP
Rturnpoints(SEXP x_, SEXP threshold_, SEXP wlength_);

#endif
#endif
