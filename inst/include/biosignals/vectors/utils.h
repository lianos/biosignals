#ifndef __BS_VECTORS_UTILS_H__
#define __BS_VECTORS_UTILS_H__

#include <vector>
#include <cmath>
#include "biosignals/macros.h"

namespace biosignals {

/**
 * Find positions in x where the "function" crosses 0
 */
template <typename T>
std::vector<int>
zero_crossings(std::vector<T> &x, int start=0, int end=-1);

template <typename T>
std::vector<int>
zero_crossings(std::vector<T> &x, int start, int end) {
    std::vector<int> zeroes;
    int push;

    if (end < 0 || end >= x.size()) {
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
            push = (abs(x[i]) < abs(x[i + 1])) ? i : i + 1;
            zeroes.push_back(push);
        }
    }
    return zeroes;
}

// Hide this from view
namespace sliding {


/**
 * Helper function to keep the value of the local min/max in current window
 *
 * Also slides window along vector so that the max/min is always correct for
 * a window that is (i - win_length):(i + win_length:i)
 */
template <typename T>
void
update_window_stats(std::vector<T> &x, int &back, int &front, int &current,
                    T &val, int win_length, int i, bool min);

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
                    current = j;
                    val = x[current];
                }
            }
        } // previous min/max fell outside of window
    } // removing elements from tail
}

}

/**
 * Locate the minimum value at point x including the +/- K neighborhood
 * around it
 */
template <typename T>
std::vector<T>
sliding_window_minimum(std::vector<T> & ARR, int K);

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
        sliding::update_window_stats(x, bidx, fidx, cidx, cval, K, i, true);
        idx[i] = cval;
    }

    return idx;
}


/**
 * Locate the maximum value at point x including the +/- K neighborhood
 * around it
 */
template <typename T>
std::vector<T>
sliding_window_maximum(std::vector<T> & ARR, int K);

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
        sliding::update_window_stats(x, bidx, fidx, cidx, cval, K, i, false);
        idx[i] = cval;
    }

    return idx;
}


/**
 * Assume x to be a coverage vector, identify the regions of x that indicate
 * where the percentiles (indicated in percentiles) start (from left to right)
 */
template <typename T>
std::vector< std::vector<int> >
coverage_quantiles(std::vector<T> &x,
                   std::vector<int> starts,
                   std::vector<int> ends,
                   std::vector<double> percentiles);

template <typename T>
std::vector< std::vector<int> >
coverage_quantiles(std::vector<T> &x,
                   std::vector<int> starts,
                   std::vector<int> ends,
                   std::vector<double> percentiles) {
    int np = percentiles.size();
    int nranges = starts.size();
    int istart;
    int iend;
    int iwidth;
    double value;
    std::vector<double> tmpvals;
    double total;
    double current_break;
    double cumsum;
    int i, j, nfound, tmp;

    std::vector< std::vector<int> > output(nranges, std::vector<int>(np));

    for (i = 0; i < nranges; i++) {
        tmpvals.clear();
        istart = starts[i];
        iend = ends[i];
        iwidth = iend - istart + 1;
        total = 0;

        for (j = 0; j < iwidth; j++) {
            total += x[istart + j];
        }

        current_break = percentiles[0];
        j = 0;
        nfound = 0;
        cumsum = 0;

        while ((j < iwidth) && (nfound < np)) {
            cumsum += x[istart + j];
            if (current_break < (cumsum / total)) {
                output[i][nfound] = j;
                if (++nfound < np) {
                    current_break = percentiles[nfound];
                }
            }
            j++;
        }

        if (nfound < np) {
            tmp = (nfound == 0) ? 0 : output[i][nfound - 1];
            for (j = nfound; j < np; j++) {
                output[i][j] = tmp;
            }
        }
    }

    return output;
}

}
#endif