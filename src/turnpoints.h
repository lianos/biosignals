#ifndef __TURNPOINTS_H__
#define __TURNPOINTS_H__

#include "biosignals.h"

template <typename T>
std::vector<T> sliding_window_minimum(std::vector<T> & ARR, int K);

template <typename T>
std::vector<T> sliding_window_maximum(std::vector<T> & ARR, int K);

RcppExport SEXP Rsliding_max(SEXP x_, SEXP k_);
RcppExport SEXP Rsliding_min(SEXP x_, SEXP k_);

RcppExport SEXP
turnpoints(SEXP x_, SEXP threshold_, SEXP wlength_);

#endif
