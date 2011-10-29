#include "vectors.h"

// template<typename T>
// std::pair< std::vector<int>, std::vector<T> >
std::pair< std::vector<int>, std::vector<double> >
as_rle(std::vector<double> &vals, double eps) {
    // std::vector<T> values;
    std::vector<double> values;
    std::vector<int> lengths;
    std::pair< std::vector<int>, std::vector<double> > out;
    
    // T current_val, tmp_val;
    double current_val, tmp_val;
    int current_len;
    
    if (vals.size() == 0) {
        out.first = lengths;
        out.second = values;
        return out;
    }
    
    current_val = vals[0];
    current_len = 1;
    
    for (int i = 1; i < vals.size(); i++) {
        tmp_val = vals[i];
        if (ALMOST_ZERO_EPS(tmp_val - current_val, eps)) {
            current_len++;
        } else {
            lengths.push_back(current_len);
            values.push_back(current_val);
            current_len = 1;
            current_val = tmp_val;
        }
    }
    
    if (values.back() != current_val) {
        values.push_back(current_val);
        lengths.push_back(current_len);
    }
    
    out.first = lengths;
    out.second = values;
    return out;
}

// template<typename T>
// std::vector<T>
std::vector<double>
expand_rle(std::vector<int> lengths, std::vector<double> &vals) {
    int length = 0;
    int i,j,len,sofar;
    // T val;
    double val;
    
    for (i = 0; i < lengths.size(); i++) {
        length += lengths[i];
    }
    
    // std::vector<T> out = std::vector<T>(length);
    std::vector<double> out = std::vector<double>(length);
    sofar = 0;
    for (i = 0; i < lengths.size(); i++) {
        val = vals[i];
        len = lengths[i];
        for (j = 0; j < len; j++) {
            out[sofar++] = val;
        }
    }
    
    return out;
}


// template <typename T>
// std::vector<int>
// zero_crossings(std::vector<T> &x, int start, int end) {
//     std::vector<int> zeroes;
//     int push;
// 
//     if (end <= 0 || end >= x.size()) {
//         end = x.size() - 1;
//     }
//     if (start >= end) {
//         start = end - 1;
//     }
//     if (start < 0) {
//         start = 0;
//     }
//     for (int i = start; i < end; i++) {
//         if ((x[i] >= 0 && x[i+1] < 0) || (x[i] <= 0) && (x[i+1] > 0)) {
//             push = (abs(x[i]) < abs(x[i + 1])) ? i : i + 1;
//             zeroes.push_back(push);
//         }
//     }
//     return zeroes;
// }
// 
// template <typename T>
// void
// update_window_stats(std::vector<T> &x, int &back, int &front, int &current,
//                     T &val, int win_length, int i, bool min);

/**
 * Helper function to keep the value of the local min/max in current window
 *
 * Also slides window along vector so that the max/min is always correct for
 * a window that is (i - win_length:i + win_length:i)
 */
// template <typename T>
// void
// update_window_stats(std::vector<T> &x, int &back, int &front, int &current,
//                     T &val, int win_length, int i, bool is_min) {
//     int j;
//     while ((front - i <= win_length) && (front < x.size())) {
//         // Add elements to the fronteir
//         if ((is_min && x[front] < val) || (!is_min && x[front] > val)) {
//             current = front;
//             val = x[current];
//         }
//         front++;
//     }
//     front--;
// 
//     // remove elements from the tail
//     while ((i - back > win_length) && (back <= x.size() - win_length)) {
//         back++;
//         if (back > current) {
//             // The previous min/max just fell outside of the window. Run
//             // through the remainder of this window to find the max up
//             // through its end
//             current = back;
//             val = x[current];
//             for (j = back + 1; j <= front; j++) {
//                 if ((is_min && x[j] < val) || (!is_min && x[j] > val)) {
//                     current = j;
//                     val = x[current];
//                 }
//             }
// 
//         } // previous min/max fell outside of window
//     } // removing elements from tail
// }
// 
// template <typename T>
// std::vector<T>
// sliding_window_minimum(std::vector<T> &x, int K) {
//     std::vector<T> idx(x.size());
//     if (idx.size() == 0) {
//         return idx;
//     }
// 
//     T cval = x[0];
//     int fidx = 0, bidx = 0, cidx = 0;
// 
//     for (int i = 0; i < x.size(); i++) {
//         update_window_stats(x, bidx, fidx, cidx, cval, K, i, true);
//         idx[i] = cval;
//     }
// 
//     return idx;
// }
// 
// template <typename T>
// std::vector<T>
// sliding_window_maximum(std::vector<T> &x, int K) {
//     std::vector<T> idx(x.size());
//     if (idx.size() == 0) {
//         return idx;
//     }
// 
//     T cval = x[0];
//     int fidx = 0, bidx = 0, cidx = 0;
// 
//     for (int i = 0; i < x.size(); i++) {
//         update_window_stats(x, bidx, fidx, cidx, cval, K, i, false);
//         idx[i] = cval;
//     }
// 
//     return idx;
// }
// 
// template <typename T>
// std::vector< std::vector<int> >
// coverage_quantiles(std::vector<T> &x, std::vector<int> starts,
//                    std::vector<int> ends, std::vector<double> percentiles) {
//     int np = percentiles.size();
//     int nranges = starts.size();
//     int istart;
//     int iend;
//     int iwidth;
//     double value;
//     std::vector<double> tmpvals;
//     double total;
//     double current_break;
//     double cumsum;
//     int i, j, nfound, tmp;
// 
//     std::vector< std::vector<int> > output(nranges, std::vector<int>(np));
// 
//     for (i = 0; i < nranges; i++) {
//         tmpvals.clear();
//         istart = starts[i];
//         iend = ends[i];
//         iwidth = iend - istart + 1;
//         total = 0;
// 
//         for (j = 0; j < iwidth; j++) {
//             total += x[istart + j];
//         }
// 
//         current_break = percentiles[0];
//         j = 0;
//         nfound = 0;
//         cumsum = 0;
// 
//         while ((j < iwidth) && (nfound < np)) {
//             cumsum += x[istart + j];
//             if (current_break < (cumsum / total)) {
//                 output[i][nfound] = j;
//                 if (++nfound < np) {
//                     current_break = percentiles[nfound];
//                 }
//             }
//             j++;
//         }
// 
//         if (nfound < np) {
//             tmp = (nfound == 0) ? 0 : output[i][nfound - 1];
//             for (j = nfound; j < np; j++) {
//                 output[i][j] = tmp;
//             }
//         }
//     }
// 
//     return output;
// }

// ----------------------------------------------------------------- Rinterface
// SEXP
// Ras_rle(SEXP vals_, SEXP eps_) {
// BEGIN_RCPP    
//     // switch(TYPEOF(x)) {
//     // case LGLSXP:
//     //     PROTECT(ans = Rle_logical_constructor(x, counts));
//     //     break;
//     // case INTSXP:
//     //     PROTECT(ans = Rle_integer_constructor(x, counts));
//     //     break;
//     // case REALSXP:
//     //     PROTECT(ans = Rle_real_constructor(x, counts));
//     //     break;
//     // case CPLXSXP:
//     //     PROTECT(ans = Rle_complex_constructor(x, counts));
//     //     break;
//     // case STRSXP:
//     //     PROTECT(ans = Rle_string_constructor(x, counts));
//     //     break;
//     // case RAWSXP:
//     //     PROTECT(ans = Rle_raw_constructor(x, counts));
//     //     break;
//     // default:
//     //     error("Rle computation of these types is not implemented");
//     // }
//     SEXP ret;
//     double eps = Rcpp::as<double>(eps_);
//     std::vector<double> vals = Rcpp::as< std::vector<double> >(vals_);
//     std::pair< std::vector<int>, std::vector<double> > out;
//     out = as_rle(vals, eps);
// 
//     ret = Rcpp::List::create(Rcpp::Named("lengths", Rcpp::wrap(out.first)),
//                              Rcpp::Named("values", Rcpp::wrap(out.second)));
//     return ret;
// END_RCPP
// }
// 
// SEXP
// Rexpand_rle(SEXP lengths_, SEXP vals_) {
// BEGIN_RCPP
//     std::vector<int> lengths = Rcpp::as< std::vector<int> >(lengths_);
//     std::vector<double> vals = Rcpp::as< std::vector<double> >(vals_);
//     std::vector<double> out = expand_rle(lengths, vals);
//     
//     return Rcpp::wrap(out);
// END_RCPP
// }

// SEXP
// Ras_rle(SEXP vals_, SEXP eps_) {
//     double eps = Rcpp::as<double>(eps_);
//     std::vector<double> vals = Rcpp::as< std::vector<double> >(vals_);
//     biosignals::Rle<double> rle(vals, eps);
//     return rle.asSEXP();
// }
// 
// 
// SEXP
// Rexpand_rle(SEXP lengths_, SEXP vals_) {
//     
// }

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
