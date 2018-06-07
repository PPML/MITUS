
#include <Rcpp.h>

#include "dist_orig_func.h"

using namespace Rcpp;

Rcpp::NumericVector dist_orig_func(
  Rcpp::NumericVector V,
  int ag,
  int tb,
  int lt,
  int im,
  int nm,
  int rg,
  int na
){
  double temp_mat[im][nm];
  double mat_sum;
  double dist_orig[im][nm];
  Rcpp::NumericVector dist_orig_v(nm*im);
  Rcpp::NumericVector dist_i_v(nm*im);

mat_sum=0;
for (int i=0; i<im; i++){
  for (int j=0; j<nm; j++){
temp_mat[i][j]=0;
dist_orig[i][j]=0;
  }
}
for (int i=0; i<(im*nm); i++){
dist_orig_v(i)=0;
dist_i_v(i)=0;
}
for(int l=0; l<im; l++) {
  for(int m=0; m<nm; m++) {
    for(int i=0; i<ag; i++) {
      for(int j=0; j<tb; j++) {
        for(int k=0; k<lt; k++) {
          for(int n=0; n<rg; n++) {
            for(int p=0; p<na; p++) {
              temp_mat[l][m]  += V(i+j*11+k*66+l*132+m*528+k*2112+p*4224);
    } } } } }
    mat_sum +=  temp_mat[l][m];
} }
for(int i=0; i<im; i++) {
  for(int j=0; j<nm; j++) {
    dist_orig[i][j]  = temp_mat[i][j]/mat_sum; // determine the proportions
    dist_orig_v((i)+(j*im)) = dist_orig[i][j];
  } }


for (int i=0; i<(nm*im); i++){
  dist_i_v(i)= dist_orig_v(i);
}

return dist_i_v;

}
