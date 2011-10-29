#ifndef __BS__RLE_H__
#define __BS__RLE_H__

#include <vector>
#include "biosignals/macros.h"

namespace biosignals {

template <class T> class Rle {
public:
    
    Rle(std::vector<T> const& vals, double eps=1e-6) {
        this->values = std::vector<T>();
        this->lengths = std::vector<int>();
        
        T current_val, tmp_val;
        int current_len;

        if (vals.size() > 0) {
            current_val = vals[0];
            current_len = 1;

            for (int i = 1; i < vals.size(); i++) {
                tmp_val = vals[i];
                if (ALMOST_ZERO_EPS(tmp_val - current_val, eps)) {
                    current_len++;
                } else {
                    this->lengths.push_back(current_len);
                    this->values.push_back(current_val);
                    current_len = 1;
                    current_val = tmp_val;
                }
            }

            if (this->values.back() != current_val) {
                this->values.push_back(current_val);
                this->lengths.push_back(current_len);
            }
        }
    }
    
    Rle(std::vector<T> const& values, std::vector<int> const& lengths) {
        this->values = values;
        this->lengths = lengths;
    }
    
    std::vector<T> expand() {
        int length = 0;
        int i,j,len,sofar;
        
        T val;

        for (i = 0; i < this->lengths.size(); i++) {
            length += this->lengths[i];
        }

        std::vector<T> out = std::vector<T>(length);
        
        sofar = 0;
        for (i = 0; i < lengths.size(); i++) {
            val = this->values[i];
            len = this->lengths[i];
            for (j = 0; j < len; j++) {
                out[sofar++] = val;
            }
        }

        return out;
    }
    
    // ~Rle() {
    //     delete this->values;
    //     delete this->lengths;
    // }
    
    // T opretator[](int32_t index) {
    //     
    // }
    
    // TODO: Rip of Rle_seqselect here
    // Rle<T> select(int32_t start, int32_t width = -1) {
    //     
    // }
    
    // TODO: Move this to Rcpp::as<Rle>(wut) functionality
    // SEXP asSEXP() {
    //     SEXP ans;
    //     PROTECT(ans = NEW_OBJECT(MAKE_CLASS("Rle")));
    //     SET_SLOT(ans, install("values"), Rcpp::wrap(this->values));
    //     SET_SLOT(ans, install("lengths"), Rcpp::wrap(this->lengths));
    //     return ans;
    // }
    
//protected:

    std::vector<T> values;
    std::vector<int> lengths;
};

}

#endif