#ifndef __BS_CONVOLVED1D_H__
#define __BS_CONVOLVED1D_H__

#include <vector>
#include "biosignals/vectors/Rle.h"

namespace biosignals {

std::vector<double>
convolve_1d(std::vector<double> const& x,
            std::vector<double> const& kernel,
            bool rescale=true);

Rle<double>
convolve_1d(Rle<double> &rle,
            std::vector<double> const &kernel,
            bool rescale=true,
            double eps=1e-6);

std::vector<double>
convolve_1d_inbounds(std::vector<double> const& x,
                     std::vector<double> const& kernel,
                     std::vector<int> const& starts,
                     std::vector<int> const& ends,
                     bool rescale, double init=0.0);

}
#endif