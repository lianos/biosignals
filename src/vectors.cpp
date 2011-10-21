#include "vectors.h"

template <typename T>
std::vector<int>
zero_crossings(std::vector<T> &x, int start, int end) {
    std::vector<int> zeroes;
    if (end <= 0 || end >= x.size()) {
        end = x.size() - 1;
    }
    if (start >= end) {
        start = end - 1;
    }
    if (start < 0) {
        start = 0;
    }
    for (int i = start; i < end; i++) {
        if ((x[i] >= 0 && x[i+1] < 0) || (x[i] <= 0) && (x[i+1] > 0)) {
            zeroes.push_back(i + 1);
        }
    }
    return zeroes;
}

template <typename T>
void
update_window_stats(std::vector<T> &x, int &back, int &front, int &current, 
                    T &val, int win_length, int i, bool min);

/**
 * Helper function to keep the value of the local min/max in current window
 * 
 * Also slides window along vector so that the max/min is always correct for
 * a window that is (i - win_length:i + win_length:i)
 */
template <typename T>
void 
update_window_stats(std::vector<T> &x, int &back, int &front, int &current,
                    T &val, int win_length, int i, bool is_min) {
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
std::vector<T>
sliding_window_minimum(std::vector<T> &x, int K) {
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
std::vector<T>
sliding_window_maximum(std::vector<T> &x, int K) {
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

// ----------------------------------------------------------------- Rinterface

SEXP
Rzero_crossings(SEXP x_) {
    std::vector<double> x = Rcpp::as< std::vector<double> >(x_);
    std::vector<int> z = zero_crossings(x);
    // convert to R coords
    for (int i = 0; i < z.size(); i++) z[i]++;
    return Rcpp::wrap(z);
}

SEXP
Rsliding_max(SEXP x_, SEXP k_) {
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    int k = Rcpp::as<int>(k_);
    return Rcpp::wrap(sliding_window_maximum(x, k));
}

SEXP
Rsliding_min(SEXP x_, SEXP k_) {
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    int k = Rcpp::as<int>(k_);
    return Rcpp::wrap(sliding_window_minimum(x, k));
}

