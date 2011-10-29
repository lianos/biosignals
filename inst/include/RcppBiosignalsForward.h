#ifndef __RCPP_BIOSIGNALS_FORWARD_H__
#define __RCPP_BIOSIGNALS_FORWARD_H__

#include <RcppCommon.h>
#include <Rconfig.h>
#include "RcppBiosignalsConfig.h"

// ----------------------------------------------------------------------------
// include all? biosignals headers
#include "biosignals/macros.h"

#include "biosignals/vectors/utils.h"
#include "biosignals/vectors/Rle.h"

#include "biosignals/convolve1d.h"

// ----------------------------------------------------------------------------

// forward declarations (whatever that means)
namespace Rcpp {
    // support wrap
    template <typename T> SEXP wrap(const biosignals::Rle<T>&);
    
    namespace traits {
        // support for as
        template <typename T> class Exporter< biosignals::Rle<T> >;
    } // trais namespace
} // Rcpp namespace


#endif