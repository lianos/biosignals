#ifndef __BS__RLE_H__
#define __BS__RLE_H__

#include <vector>
#include "biosignals/macros.h"

namespace biosignals {

template <class T> class Rle {
public:
    
    Rle(double eps=1e-6) {
        this->eps = eps;
        this->values = std::vector<T>(0);
        this->lengths = std::vector<int>(0);
        this->length = 0;
    }
    
    Rle(std::vector<T> const& vals, double eps=1e-6) {
        this->values = std::vector<T>();
        this->lengths = std::vector<int>();
        
        this->length = vals.size();
        this->eps = eps;
        
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
            if (this->values.size() == 0 /* vals ?: rep(x, y) or length 1*/
                || (this->values.back() != current_val)) {
                this->values.push_back(current_val);
                this->lengths.push_back(current_len);
            }
        }
    }
    
    Rle(std::vector<T> const& values, std::vector<int> const& lengths,
        double eps=1e-6) {
        this->values = values;
        this->lengths = lengths;
        this->length = std::accumulate(lengths.begin(), lengths.end(), 0);
        this->eps = eps;
    }
    
    std::vector<T> expand() {
        int length = 0;
        int i,j,len,sofar;
        
        T val;

        for (i = 0; i < this->lengths.size(); i++) {
            length += this->lengths[i];
        }

        // std::vector<T> out = std::vector<T>(length);
        std::vector<T> out(length);
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
    
    size_t size() {
        return this->length;
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
    
    
//protected:

    std::vector<T> values;
    std::vector<int> lengths;
    size_t length;
    double eps;
};

}

#endif