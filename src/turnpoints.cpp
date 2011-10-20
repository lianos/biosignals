#include "turnpoints.h"

#include <cstdlib>
#include <deque>
#include <iostream>
#include <vector>

void update_min_window();
void update_max_window();


template <typename T>
std::vector<T> sliding_window_minimum(std::vector<T> &ARR, int K) {
    std::vector<T> idx(ARR.size());
    if (idx.size() == 0) {
        return idx;
    }

    T cval = ARR[0];
    int fidx = 0, bidx = 0, cidx = 0;

    for (int i = 0; i < ARR.size(); i++) {
        // Add elements at the frontier
        while ((fidx - i <= K) && fidx < ARR.size()) {
            if (ARR[fidx] < cval) {
                cval = ARR[fidx];
                cidx = fidx;
            }
            fidx++;
        }
        fidx--;

        // Remove elements from the tail
        while ((i - bidx > K) && (bidx <= ARR.size() - K)) {
            bidx++;
            // if the max is at bidx, we need to find the next max
            // this can go past position i
            if (bidx > cidx) {
                cval = ARR[bidx];
                cidx = bidx;
                for (int j = bidx + 1; j <= fidx; j++) {
                    if (ARR[j] < cval) {
                        cval = ARR[j];
                        cidx = j;
                    }
                }
            }
        }

        idx[i] = cval;
    }

    return idx;
}

template <typename T>
std::vector<T> sliding_window_maximum(std::vector<T> &ARR, int K) {
    std::vector<T> idx(ARR.size());
    if (idx.size() == 0) {
        return idx;
    }

    T cval = ARR[0];
    int fidx = 0, bidx = 0, cidx = 0;

    for (int i = 0; i < ARR.size(); i++) {
        // Add elements at the frontier
        while ((fidx - i <= K) && fidx < ARR.size()) {
            if (ARR[fidx] > cval) {
                cval = ARR[fidx];
                cidx = fidx;
            }
            fidx++;
        }
        fidx--;

        // Remove elements from the tail
        while ((i - bidx > K) && (bidx <= ARR.size() - K)) {
            bidx++;
            // if the max is at bidx, we need to find the next max
            // this can go past position i
            if (bidx > cidx) {
                cval = ARR[bidx];
                cidx = bidx;
                for (int j = bidx + 1; j <= fidx; j++) {
                    if (ARR[j] > cval) {
                        cval = ARR[j];
                        cidx = j;
                    }
                }
            }
        }

        idx[i] = cval;
    }

    return idx;
}

/**
 * Returns a pair representing (minima, maxima) vectors
 * with the position of each from x
 */
template <typename T>
std::pair< std::vector<int>, std::vector<int> >
local_turnpoints(std::vector<T> &x, int K) {
    std::vector<T> idx(ARR.size());
    if (idx.size() == 0) {
        return idx;
    }

    T cmin = ARR[0];
    T cmax = ARR[0];
    int fidx_min = 0, bidx_min = 0, cidx_min = 0;
    int fidx_max = 0, bidx_max = 0, cidx_max = 0;

    for (int i = 0; i < ARR.size(); i++) {
        // Add elements at the frontier
        while ((fidx - i <= K) && fidx < ARR.size()) {
            if (ARR[fidx] > cval) {
                cval = ARR[fidx];
                cidx = fidx;
            }
            fidx++;
        }
        fidx--;

        // Remove elements from the tail
        while ((i - bidx > K) && (bidx <= ARR.size() - K)) {
            bidx++;
            // if the max is at bidx, we need to find the next max
            // this can go past position i
            if (bidx > cidx) {
                cval = ARR[bidx];
                cidx = bidx;
                for (int j = bidx + 1; j <= fidx; j++) {
                    if (ARR[j] > cval) {
                        cval = ARR[j];
                        cidx = j;
                    }
                }
            }
        }

        idx[i] = cval;
    }
}


// ----------------------------------------------------------------- R interface
SEXP Rsliding_max(SEXP x_, SEXP k_) {
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    int k = Rcpp::as<int>(k_);
    return Rcpp::wrap(sliding_window_maximum(x, k));
}

SEXP Rsliding_min(SEXP x_, SEXP k_) {
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    int k = Rcpp::as<int>(k_);
    return Rcpp::wrap(sliding_window_minimum(x, k));
}


SEXP
Rlocal_turnpoints(SEXP x_, SEXP threshold_, SEXP wlength_) {
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    double threshold = Rcpp::as<double>(threshold_);
    int wlength = Rcpp::as<int>(wlength_);
    
    return R_NilValue;
}


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

