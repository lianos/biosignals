#include "vectors.h"

template <typename T>
std::vector<int>
zero_crossings(std::vector<T> &x, int start, int end) {
    std::vector<int> zeroes;
    int push;

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
            push = (abs(x[i]) < abs(x[i + 1])) ? i : i + 1;
            zeroes.push_back(push);
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
                    current = j;
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

template <typename T>
std::vector< std::vector<int> >
coverage_quantiles(std::vector<T> &x, std::vector<int> starts,
                   std::vector<int> ends, std::vector<double> percentiles) {
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

// ----------------------------------------------------------------- Rinterface

SEXP
Rzero_crossings(SEXP x_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as< std::vector<double> >(x_);
    std::vector<int> z = zero_crossings(x);
    // convert to R coords
    for (int i = 0; i < z.size(); i++) z[i]++;
    return Rcpp::wrap(z);
END_RCPP
}

SEXP
Rsliding_max(SEXP x_, SEXP k_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    int k = Rcpp::as<int>(k_);
    return Rcpp::wrap(sliding_window_maximum(x, k));
END_RCPP
}

SEXP
Rsliding_min(SEXP x_, SEXP k_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
    int k = Rcpp::as<int>(k_);
    return Rcpp::wrap(sliding_window_minimum(x, k));
END_RCPP
}

SEXP
Rcoverage_quantiles(SEXP x_, SEXP starts_, SEXP ends_, SEXP breaks_) {
BEGIN_RCPP
    std::vector<double> x = Rcpp::as< std::vector<double> >(x_);
    std::vector<int> starts = Rcpp::as< std::vector<int> >(starts_);
    std::vector<int> ends = Rcpp::as< std::vector<int> >(ends_);
    std::vector<double> breaks = Rcpp::as< std::vector<double> >(breaks_);
    std::vector< std::vector<int> > result;
    Rcpp::NumericMatrix out(starts.size(), breaks.size());
    int i,j;

    // R --> C indexing
    for (i = 0; i < starts.size(); i++) starts[i] -= 1;
    for (i = 0; i < ends.size(); i++) ends[i] -= 1;

    result = coverage_quantiles(x, starts, ends, breaks);

    for (i = 0; i < starts.size(); i++) {
        for (j = 0; j < breaks.size(); j++) {
            // C --> R indexing
            out(i, j) = result[i][j] + 1;
        }
    }

    return Rcpp::wrap(out);
END_RCPP
}
