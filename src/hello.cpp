#include <Rcpp.h>
#include "hello.h"

typedef void (*func_type)(int);



void run_forloop(func_type myfunction)
{
  for(int i = 0; i < 2; i++)
  {
    myfunction(i);

  }


}

void dummyfunc(int x)
{

  Rcpp::Rcout << x << "\n";
}

void hello_world(const double *input, double *output)
{
  Rcpp::Rcout << "input[0] " << input[0] << " input[1] " << input[1] << "\n";
  output[0] = input[1];
  output[1] = input[0];

  run_forloop(&dummyfunc);
}
