#include "f1.h"
SEXP f1(){
  using namespace Rcpp ;

  NumericVector x(1);
  x = 100;
  return x;
}
