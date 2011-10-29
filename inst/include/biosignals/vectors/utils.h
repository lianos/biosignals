#ifndef __BS_VECTORS_UTILS_H__
#define __BS_VECTORS_UTILS_H__

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

#endif