#ifndef __VECTORS_H__
#define __VECTORS_H__

#include "RcppBiosignals.h"

// ----------------------------------------------------------------- Rinterface

RcppExport SEXP
Rsliding_max(SEXP x_, SEXP k_);

RcppExport SEXP
Rsliding_min(SEXP x_, SEXP k_);

RcppExport SEXP
Rzero_crossings(SEXP x_);

RcppExport SEXP
Rcoverage_quantiles(SEXP x_, SEXP starts_, SEXP ends_, SEXP breaks_);

RcppExport SEXP
Ras_rle(SEXP vals_, SEXP eps_);

RcppExport SEXP
Rexpand_rle(SEXP lengths_, SEXP vals_);

RcppExport SEXP
Rexpand_rle_S4(SEXP rle_);

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

