

# STEP 1: make sure you have the work chain all set-up. This is best described on the web.
# If you have a mac it should be easy, if  PC it will be painful. Sorry.

# Then....

library(Rcpp)

# Code model

cppFunction('
  List modC(
    int                 a_scalar_integer,
    NumericMatrix       a_matrix,
    std::vector<double> a_vector,
    double              a_scalar_real
      ) { 
    /// objects created internally -- declare type, name, and dimension
    /// note that doubles can have multple dimensions
    double        a_vector2[a_scalar_integer];
    double        a_matrix2[a_matrix.nrow()][a_matrix.ncol()];
    double        an_array[5][5][5]; 

    /// There are big speed gains to operative in primitive types, such as doubles and integers
    /// therefore we will convert our NumericMatrix a_matrix to a array of doubles (a_matrix2)
    /// not in C++ that indices start at zero!

    for(int i=0; i<a_matrix.nrow(); i++) { 
      for(int j=0; j<a_matrix.ncol(); j++) { 
        a_matrix2[i][j] = a_matrix(i,j); 
      }
    }

    /// Now lets do some stuff.. set all cells in an_array equal to 2.0
    for(int i=0; i<5; i++) { 
      for(int j=0; j<5; j++) {
        for(int k=0; k<5; k++) {
          an_array[i][j][k] = 2.0;
        }
      }
    }
    /// Then add three if i=j=k, else multiply by 2
    for(int i=0; i<5; i++) { 
      for(int j=0; j<5; j++) {
        for(int k=0; k<5; k++) {
          if( i==j & j==k ) {
            an_array[i][j][k] += 3.0;
          } else
            an_array[i][j][k] *= 2.0;
          }
        }
      }

    /// And finally...
    for(int i=0; i<a_scalar_integer; i++) { 
      a_vector2[i] = i;
    }

    /// NOW export results
    /// Can export array of doubles, so convert back to NumericMatrix and NumericVector...
    NumericVector        a_vector3(a_scalar_integer); 
    NumericMatrix        a_matrix3(a_matrix.nrow(),a_matrix.ncol()); 

    for(int i=0; i<a_scalar_integer; i++) { 
      a_vector3(i) = a_vector2[i];
    }

    for(int i=0; i<a_matrix.nrow(); i++) { 
      for(int j=0; j<a_matrix.ncol(); j++) { 
        a_matrix3(i,j) = a_matrix2[i][j]; 
      }
    }

    return Rcpp::List::create(
      Rcpp::Named("a_matrix3") = a_matrix3,
      Rcpp::Named("a_vector3") = a_vector3);
    }')


# Test model
modC(a_scalar_integer = 7,
     a_matrix = matrix(1,3,4),
     a_vector = 1:7,
     a_scalar_real = 5)

