#ifndef __BS_CONVOLVED1D_H__
#define __BS_CONVOLVED1D_H__

std::vector<double> *
convolve_1d(const std::vector<double> &x,
            const std::vector<double> &kernel, bool rescale);

std::vector<double> *
convolve_1d_inbounds(const std::vector<double> &x,
                     const std::vector<double> &kernel,
                     const std::vector<int> &starts,
                     const std::vector<int> &ends, bool rescale,
                     double init = 0.0);

#endif