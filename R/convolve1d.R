convolve1d <- function(x, kernel, bandwidth=20, rescale=TRUE) {
  # x <- c(rep(x[1], length(kernel)), x, rep(tail(x, 1), length(kernel)))
  ret <- .Call("Rconvolve_1d", as.numeric(x), as.numeric(kernel),
               rescale)
  ret
}

# convolve1d <- cxxfunction(signature(x_='numeric', kernel_='numeric'), body="
# Rcpp::NumericVector x(x_);
# Rcpp::NumericVector kernel(kernel_);
# int n = x.size() + kernel.size() - 1;
# Rcpp::NumericVector result(n);
# int i,j;
# 
# for (i = 0; i < n; i++) result[i] = 0.0;
# for (i = 0; i < x.size(); i++) {
#   for (j = 0; j < kernel.size(); j++) {
#     result[i+j] += x[i] * kernel[j];
#   }
# }
# 
# return result;
# ", plugin="Rcpp")
