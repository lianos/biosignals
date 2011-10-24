#ifndef __VECTORS_H__
#define __VECTORS_H__

#include "biosignals.h"

template <typename T>
std::vector<int>
zero_crossings(std::vector<T> &x, int start=0, int end=0);

template <typename T>
std::vector<T>
sliding_window_minimum(std::vector<T> & ARR, int K);

template <typename T>
std::vector<T>
sliding_window_maximum(std::vector<T> & ARR, int K);

template <typename T>
std::vector< std::vector<int> >
coverage_quantiles(std::vector<T> &x, std::vector<int> starts,
                   std::vector<int> ends, std::vector<double> percentiles);

// ----------------------------------------------------------------- Rinterface

RcppExport SEXP
Rsliding_max(SEXP x_, SEXP k_);

RcppExport SEXP
Rsliding_min(SEXP x_, SEXP k_);

RcppExport SEXP
Rzero_crossings(SEXP x_);

RcppExport SEXP
Rcoverage_quantiles(SEXP x_, SEXP starts_, SEXP ends_, SEXP breaks_);

#endif

// --------------------------------------------------------------------- scratch
// sliding window implementations are courtesy of:
// https://github.com/keegancsmith/Sliding-Window-Minimum
// template <typename T>
// std::vector<T> sliding_window_minimum(std::vector<T> & ARR, int K) {
//     // pair<int, int> represents the pair (ARR[i], i)
//     std::deque< std::pair<T, int> > window;
//     std::vector<T> idxs(ARR.size());
//     for (int i = 0; i < ARR.size(); i++) {
//         while (!window.empty() && window.back().first >= ARR[i]) {
//             window.pop_back();
//         }
//         window.push_back(std::make_pair(ARR[i], i));
//         while (window.front().second <= i - K) {
//             window.pop_front();
//         }
//         idxs[i] = window.front().first;
//     }
//     return idxs;
// }
//
// template <typename T>
// std::vector<T> sliding_window_maximum(std::vector<T> & ARR, int K) {
//     // pair<int, int> represents the pair (ARR[i], i)
//     std::deque< std::pair<T, int> > window;
//     std::vector<T> idxs(ARR.size());
//     for (int i = 0; i < ARR.size(); i++) {
//         while (!window.empty() && window.back().first <= ARR[i]) {
//             window.pop_back();
//         }
//         window.push_back(std::make_pair(ARR[i], i));
//         while (window.front().second >= i - K) {
//             window.pop_front();
//         }
//         idxs[i] = window.front().first;
//     }
//     return idxs;
// }

