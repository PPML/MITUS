#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]


//////////////////// CODE FOR THE OUTPUTS OF THE MODEL ////////////////////
////////////////////       FORMAT FOLLOWS BELOW        ////////////////////
/// LOCATION | YR | AGE | NATIV | SCEN # | OUTPUT | VALUE | CI LOW | CI HIGH ///

//////////////////// ONLY PULL RESULTS FROM JUNE OUTPUT ///////////////////
if(m==6) {
  Outputs[loc][yr][ag][na][scen][var][ci_low][ci_low]
