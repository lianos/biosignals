#include "turnpoints.h"

#include <cstdlib>
#include <deque>
#include <iostream>
#include <vector>

template <typename T>
void update_window_stats(std::vector<T> &x, int &back, int &front,
                         int &current, T &val, int win_length, int i,
                         bool min);

template <typename T>
void update_window_stats(std::vector<T> &x, int &back, int &front,
                         int &current, T &val, int win_length, int i,
                         bool is_min) {
    int j;
    while ((front - i <= win_length) && (front < x.size())) {
        // Add elements to the fronteir
        if ((is_min && x[front] < val) || (!is_min && x[front] > val)) {
            current = front;
            val = x[current];
        }
        front++;
    }
    front--;
        
    // remove elements from the tail
    while ((i - back > win_length) && (back <= x.size() - win_length)) {
        back++;
        if (back > current) {
            // The previous min/max just fell outside of the window. Run
            // through the remainder of this window to find the max up
            // through its end
            current = back;
            val = x[current];
            for (j = back + 1; j <= front; j++) {
                if ((is_min && x[j] < val) || (!is_min && x[j] > val)) {
                    current = front;
                    val = x[current];
                }
            }
            
        } // previous min/max fell outside of window
    } // removing elements from tail
}

template <typename T>
std::vector<T> sliding_window_minimum(std::vector<T> &x, int K) {
    std::vector<T> idx(x.size());
    if (idx.size() == 0) {
        return idx;
    }

    T cval = x[0];
    int fidx = 0, bidx = 0, cidx = 0;

    for (int i = 0; i < x.size(); i++) {
        update_window_stats(x, bidx, fidx, cidx, cval, K, i, true);
        idx[i] = cval;
    }

    return idx;
}

template <typename T>
std::vector<T> sliding_window_maximum(std::vector<T> &x, int K) {
    std::vector<T> idx(x.size());
    if (idx.size() == 0) {
        return idx;
    }

    T cval = x[0];
    int fidx = 0, bidx = 0, cidx = 0;

    for (int i = 0; i < x.size(); i++) {
        update_window_stats(x, bidx, fidx, cidx, cval, K, i, false);
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
    std::vector<int> mins;
    std::vector<int> maxs;
    std::pair< std::vector<int>, std::vector<int> > result;
    
    if (x.size() == 0) {
        return std::make_pair(mins, maxs);
    }

    T cmin = x[0];
    T cmax = x[0];
    int fidx_min = 0, bidx_min = 0, cidx_min = 0;
    int fidx_max = 0, bidx_max = 0, cidx_max = 0;

    for (int i = 0; i < x.size(); i++) {
    }
    
    result = std::make_pair(mins, maxs);
    return result;
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

