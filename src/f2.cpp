#include "f1.h"
#include "f2.h"

SEXP f2(){
  using namespace Rcpp;
  NumericVector x = f1();
  x[0] = x[0] + 1;
  return x;
}
