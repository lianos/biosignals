## This function was taken from:
## http://hellmund.blogspot.com/2011/02/example-using-rcpp-inline-and-fftw.html
## 
## This isn't correct

require(inline)

plug <- Rcpp:::Rcpp.plugin.maker(
  include.before="#include <fftw3.h>",
  # libs=paste(sprintf("-L%s -lRcpp", dirname(Rcpp:::LdFlags())),
  #            "-Wl,-rpath,/usr/local/lib/R/site-library/Rcpp/lib",
  #            "-lfftw3 -lm -L/usr/lib")
  libs="-lfftw3"
)
registerPlugin("FFTWconv", plug )

convFFTW <- cxxfunction(
signature(xIn = "numeric", yIn = "numeric"),
body = '
Rcpp::NumericVector x(xIn);
Rcpp::NumericVector y(yIn);
int nx=x.size();
int ny=y.size();
int n=nx+ny-1;
Rcpp:NumericVector retPar(n,0.0);

double in_1[n];
double out_2[n];

fftw_complex *in_2,*out_1,*irf_fft,*y_fft;
in_2=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
out_1=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);

fftw_plan forward;
fftw_plan backward;

forward=fftw_plan_dft_r2c_1d(n, in_1, out_1, FFTW_ESTIMATE);
backward=fftw_plan_dft_c2r_1d(n, in_2, out_2, FFTW_ESTIMATE);

for (int i=0; i<(n-nx); i++) {
  in_1[i] = 0.0;
}

for (int i=(n-nx); i fftw_execute(forward);
irf_fft=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
for(int i=0;i for(int j=0;j<2;j++) irf_fft[i][j]=out_1[i][j];

for(int i=0;i for(int i=ny;i fftw_execute(forward);
y_fft=(fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
for(int i=0;i for(int j=0;j<2;j++) y_fft[i][j]=out_1[i][j];

for(int i=0;i in_2[i][0]=irf_fft[i][0]*y_fft[i][0]+irf_fft[i][1]*y_fft[i][1];
in_2[i][1]=irf_fft[i][1]*y_fft[i][0]-irf_fft[i][0]*y_fft[i][1];
}

fftw_execute(backward);

for(int i=0;i 
fftw_destroy_plan(forward);
fftw_destroy_plan(backward);

return retPar;
', plugin="FFTWconv")

convolve_fftw <- function(x,y) {
 convFFTW(x,rev(y))
}
