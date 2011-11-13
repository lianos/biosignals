#ifndef __BS_CONVOLVED1D_H__
#define __BS_CONVOLVED1D_H__

#include <vector>

namespace biosignals {

std::vector<double>
convolve_1d(std::vector<double> const& x,
            std::vector<double> const& kernel,
            bool rescale=true);

std::vector<double>
convolve_1d_inbounds(std::vector<double> const& x,
                     std::vector<double> const& kernel,
                     std::vector<int> const& starts,
                     std::vector<int> const& ends,
                     bool rescale, double init=0.0);

}
#endif