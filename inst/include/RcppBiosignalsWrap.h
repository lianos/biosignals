#ifndef __RCPP_BIOSIGNALS_WRAP_H__
#define __RCPP_BIOSIGNALS_WRAP_H__

namespace Rcpp {
    
    // Creates an IRanges::Rle object out of biosignals::Rle
    template <typename T>
    SEXP wrap(const biosignals::Rle<T>& rle) {
        S4 ans = S4("Rle");
        ans.slot("values") = wrap(rle.values);
        ans.slot("lengths") = wrap(rle.lengths);
        return ans;
    }
    
    
    namespace traits {
        // Support Rcpp::as< Rle<double> >
        template<typename T>
        class Exporter< biosignals::Rle<T> > {
        public:
            Exporter(SEXP x) {
                // const int RTYPE = ::Rcpp::traits::r_sexptype_traits<T>::rtype;
                // if (TYPEOF(x) != RTYPE)
                //     throw std::invalid_argument("Wrong R type for mapped vector");
                ::Rcpp::S4 sx = ::Rcpp::S4(x);
                values = ::Rcpp::as< std::vector<T> >(sx.slot("values"));
                lengths = ::Rcpp::as< std::vector<int> >(sx.slot("lengths"));
            }
            
            biosignals::Rle<T> get() {
                return biosignals::Rle<T>(values, lengths);
            }
            
        protected:
            std::vector<T> values;
            std::vector<int> lengths;
        };
    }
}

#endif