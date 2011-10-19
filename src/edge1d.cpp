#include "edge1d.h"

std::vector<double> *
convolve_1d(std::vector<double> *x, std::vector<double> *kernel, bool rescale) {
  double max_x, max_result, scale_factor;
  int klen = kernel->size();
  int n = x->size() + (2 * klen);
  int N = n + kernel->size(); // used to pad head/tail
  int i,j;
  double val;
  
  std::vector<double> *result = new std::vector<double>(N, 0.0);
  Rprintf("result vector of size %d is allocated\n", result->size());
  max_x = (*x)[0];
  max_result = 0;

  // Dodge head edge effects
  for (i = 0; i < kernel->size(); i++) {
    for (j = 0; j < kernel->size(); j++) {
      if ((i + j) >= N) Rprintf("[head] out of bounds i,j : %d,%d\n", i, j);
      val = (*result)[i+j] + (*x)[0] * (*kernel)[j];
      // Rprintf("Head %d = %.2f\n", i + j, val);
      (*result)[i+j] = val;
    }
  }
  
  // Main convolution
  for (i = klen; i < x->size() + klen; i++) {
    for (j = 0; j < klen; j++) {
      if ((i + j) >= N) Rprintf("[main] out of bounds i,j : %d,%d\n", i, j);
      val = (*result)[i+j] + (*x)[i - klen] * (*kernel)[j];
      // Rprintf("Body %d = %.2f\n", i + j, val);
      (*result)[i+j] = val;
      if (val > max_result) {
        max_result = val;
      }
    }
  }
  
  // Dodge tail edge effects
  for (i = x->size() + klen; i < n; i++) {
    for (j = 0; j < klen; j++) {
      if ((i + j) >= N) {
        Rprintf("[tail] out of bounds i,j : %d,%d\n", i, j);
        break;
      }
      val = (*result)[i+j] + (*x)[x->size() - 1] * (*kernel)[j];
      // Rprintf("Tail %d = %.2f\n", i + j, val);
      (*result)[i+j] = val;
    }
  }
  
  scale_factor = max_x / max_result;
  
  /* remove tails */
  int naxe = 1.5 * klen;
  Rprintf("sizse before %d\n", result->size());
  result->erase(result->begin(), result->begin() + naxe);
  Rprintf("sizse after %d\n", result->size());
  result->erase(result->begin() + x->size(), result->end());
  return result;
}

// std::vector<double>
// convolve_1d(std::vector<double> x, std::vector<double> kernel, bool rescale) {
//   double max_x, max_result, scale_factor;
//   int n = x.size() + kernel.size();
//   int N = x.size() + kernel.size() + 2 * kernel.size(); // used to pad head/tail
//   int i,j;
//   double val;
//   
//   std::vector<double> result;
//   result.resize(N);
//   
//   for (i = 0; i < N; i++) {
//     result[i] = 0.0;
//     if ((i < x.size()) && (i == 0 || max_x < x[i])) {
//       max_x = x[i];
//     }
//   }
//   
//   max_x = x[0];
//   max_result = 0;
//   
//   // Dodge head edge effects
//   for (i = 0; i < kernel.size(); i++) {
//     for (j = 0; j < kernel.size(); j++) {
//       val = result[i+j] + x[0] * kernel[j];
//       // Rprintf("Head %d = %.2f\n", i + j, val);
//       result[i+j] = val;
//     }
//   }
//   
//   // Main convolution
//   for (i = kernel.size(); i < x.size(); i++) {
//     for (j = 0; j < kernel.size(); j++) {
//       val = result[i+j] + x[i] * kernel[j];
//       // Rprintf("Body %d = %.2f\n", i + j, val);
//       result[i+j] = val;
//       if (val > max_result) {
//         // Rprintf("[main] New max: (%d,%d : %d) = %.2f\n", i, j, i+j, val);
//         max_result = val;
//       }
//     }
//   }
//   
//   // Dodge tail edge effects
//   for (i = x.size() + kernel.size(); i < result.size(); i++) {
//     for (j = 0; j < kernel.size(); j++) {
//       val = result[i+j] + x[x.size() - 1] * kernel[j];
//       // Rprintf("Tail %d = %.2f\n", i + j, val);
//       result[i+j] = val;
//     }
//   }
// 
//   /* remove tails */
//   int naxe = kernel.size();
//   Rprintf("Size of result before trim: %d\n", result.size());
//   result.erase(result.begin(), result.begin() + naxe);
//   result.erase(result.begin() + x.size(), result.end());
//   Rprintf("Size of result after trim: %d\n", result.size());
//   
//   Rprintf("Max value input: %.2f\n", max_x);
//   Rprintf("Max value output: %.2f\n", max_result);
//   
//   if (rescale) {
//     scale_factor = max_x / max_result;
//     for (i = 0; i < result.size(); i++) {
//       Rprintf("scaling %d: %.2f -> %.2f\n", i, result[i], result[i] * scale_factor);
//       result[i] *= scale_factor;
//     }
//     Rprintf("done scaing\n");
//   }
//   
//   return result;
// }

SEXP Rconvolve_1d(SEXP x_, SEXP kernel_, SEXP rescale_) {
  std::vector<double> x = Rcpp::as<std::vector<double> >(x_);
  std::vector<double> kernel = Rcpp::as<std::vector<double> >(kernel_);
  bool rescale = Rcpp::as<bool>(rescale_);
  std::vector<double> *result = convolve_1d(&x, &kernel, rescale);
  
  // Rprintf("copying output\n");
  // std::vector<double> out(result->size());
  // for (int i = 0; i < result->size(); i++) {
  //   Rprintf("copying %d:%.2f\n", i, (*result)[i]);
  //   out[i] = (*result)[i];
  // }
  // // std::copy(result->begin(), result->end(), out.begin());
  // Rprintf("sending back to R\n");
  // return Rcpp::wrap(out);
  return Rcpp::wrap(*result);
}

SEXP fencepost_convolve_1d(SEXP x_, SEXP kernel_, SEXP starts_, SEXP ends_) {
  return R_NilValue;
}
