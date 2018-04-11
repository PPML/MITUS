//////input variables
#include <RcppArmadillo.h>
#include <algorithm>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat mat_mult(
    arma::mat& a,
    arma::mat& b) {
  return(a % b);
}


//[[Rcpp::export]]
Rcpp::List rebal(
  double              V1[11][6][2][4][4][2][3],
  Rcpp::NumericMatrix dist_gen,
  arma::mat           dist_new,
  arma::mat           can_go,
  arma::mat           did_go,
  arma::mat           dist_goal,
  arma::vec           dist_goal_v,
  arma::vec           dist_orig_v,
  arma::vec           diff_i_v
){

  ////////////////////////////////////////////////////////////////////////////////
  ////////    BELOW IS A LIST OF THE VARIABLES CREATED INTERNALLY IN MODEL   /////
  ////////////////////////////////////////////////////////////////////////////////
  int N; int m; int p; int r; int c; int m2; int p2;
  double        temp;
  arma::rowvec  temp_vec(16);
  arma::mat     temp_mat(16,16);
  arma::mat     temp_mat2(16,16);
  arma::mat     trans_mat(16,16);
  arma::mat     trans_mat_tot(16,16);
  arma::vec     dist_i_v(16);
  arma::mat     dist_orig(4,4);
  double        dist_new_fin[4][4];
  double        frc;
  double        pop_t;
  double        sse;

  ///////////////////////////////////////////////////////////////////////////////
  ///////                            INITIALIZE                             /////
  ///////////////////////////////////////////////////////////////////////////////

  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      dist_new_fin[i][j] = 0;
    } }

//////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// REBALANCE THE POPULATION //////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////

  ////// need to define the current distribution of persons across the RG at this timestep
////// first calculate the total population at this time step
pop_t =0;
////// second calculate the distribution of population across the two risk factors
for(int ag=0; ag<11; ag++) {
  for(int tb=0; tb<6; tb++) {
    for(int lt=0; lt<2; lt++) {
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++) {
              pop_t += V1[ag][tb][lt][im][nm][rg][na];
              dist_orig(nm,im) += V1[ag][tb][lt][im][nm][rg][na]; //sum across risk factors
              dist_orig(nm,im)  = dist_orig(nm,im)/pop_t; // determine the proportions

              ///insert error that distribution does not sum to one and then go from there;
              //Rcpp::Rcout <<"dist_orig at"<< im << "and"<< nm << "is" << dist_orig[nm][im]<< "\n";
            } } } } } } }


for(int n=0; n<N; n++){

  /////// CALCULATE DISTANCE FROM CURRENT DISTRIBUTION TO GOAL DISTRIBUTION /////
    for (int i=0; i<sizeof(dist_i_v); i++){
      diff_i_v[i] = dist_i_v[i] - dist_goal_v[i];
    }
  //////////                  CREATE TRANSITION MATRIX                    ////////
    for (int r=0; r<16; r++){
      for (int c=0; c<16; c++){
        temp=diff_i_v[r]-diff_i_v[c];
        trans_mat[r,c] = can_go[r,c]*(std::max(0.0,temp));
        Rcpp::Rcout << "initial trans_mat is" << trans_mat;
      }
    }
  //////////                ADJUST THE TRANSITION MATRIX                  ////////
    //////////   1ST SCALE UP RATES, 2ND MAKE SURE DOES NOT SUM OVER 1    //////////
    frc = 0.1;  // approach seems quite sensitive to this value, = fraction of change to

    for(int i=0; i<16; i++){
      for(int j=0; j<16; j++){
        temp_vec      = sum(trans_mat, 0); ///needs to be the column sums
        temp_mat[i,j] =  trans_mat[i,j] / dist_i_v[i]*frc;
        temp_mat[i,j] =  temp_mat [i,j] / (std::max(1.0, temp_vec[i]));
      }}
    Rcpp::Rcout << "temp trans_mat is" << trans_mat;

    temp_mat2 = diagmat(1-(sum(temp_mat,1)));

    //////////                      FINALIZE TRANS_MAT                    //////////
      for(int i=0; i<16; i++){
        for(int j=0; j<16; j++){
          if (i != j) {
            trans_mat[i,j] = temp_mat[i,j];
          } else {trans_mat_tot[i,j]=temp_mat2[i,j];
          } } }
    //////////                RECORD ABSOLUTE TRANSITIONS                 //////////
      for(int i=0; i<16; i++){
        for(int j=0; j<16; j++){
          if (i==j){
            did_go[i,j]=0;
          } else {
            did_go[i,j] += dist_i_v[i]*trans_mat[i,j];
          }
          //////////               UPDATE THE DISTRIBUTION VECTOR             ////////////
            dist_i_v[i] = dist_i_v[i]*trans_mat[i,j];
        }}

    //////////                    NOW UPDATE IN ONE STEP                 ///////////
      temp_mat=did_go;

    for(int i=0; i<16; i++){
      for(int j=0; j<16; j++){
        temp_mat[i,j] = did_go[i,j] / dist_orig_v[i,j];
      } }
    Rcpp::Rcout<< "temp_mat is" << temp;

    ///Use the diagmat functin of RcppArmadillo to create a diagonal matrix
    ///with the vector of values defined below.

    trans_mat_tot = arma::diagmat(1-sum(trans_mat_tot, 1));

    ///Replace non-diagonal values with the values from temp above

    for(int i=0; i<16; i++){
      for(int j=0; j<16; j++){
        if (i != j) {
          trans_mat_tot[i,j] = temp_mat[i,j];
        } else {trans_mat_tot[i,j]=trans_mat_tot[i,j];
        } } }

    Rcpp::Rcout<< "trans_mat_tot is" << trans_mat_tot;

    //////////           NOW FINALLY UPDATE THE DISTRIBUTION           ///////////
      for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
          dist_new[i,j] = dist_orig[i,j];
        } }
    for (int m=0; m<4; m++){
      for (int p=0; p<4; p++){
        for (int m2=0; m2<4; m2++){
          for (int p2=0; p2<4; p2++){
            ////scalar multiplication (not matrix multiplication)
            dist_new[m+1,p+1] += mat_mult(dist_orig[m2+1,p2+1], trans_mat_tot[1+m2+p2*4,1+m+p*4]);
          } } } }
    for(int i=0; i<16; i++){
      for(int j=0; j<16; j++){
        temp_mat[i,j] = dist_goal[i,j];
        sse = arma::accu(pow(dist_goal[i,j] - dist_new[i,j], 2)) /
        arma::accu(pow(dist_goal[i,j] - dist_orig[i,j], 2));
      }}
    Rcpp::Rcout << "sse is" << sse;
}

/////////////APPLY THE NEW DISTRIBUTION TO THE POPULATION /////////////////////
  /////////////revert the arma::vec to a double;
for (int im=0; im<4; im++){
  for (int nm=0; nm<4; nm++){
    dist_new_fin[nm][im] = dist_new[nm,im];
  } }
for(int ag=0; ag<11; ag++) {
  for(int tb=0; tb<6; tb++) {
    for(int lt=0; lt<2; lt++){
      for (int im=0; im<4; im++){
        for (int nm=0; nm<4; nm++){
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++){
              V1[ag][tb][lt][im][nm][rg][na] = V1[ag][tb][lt][im][nm][rg][na]*dist_new_fin[nm][im];
            } } } } } } }

////////// RETURN THE NEW DISTRIBUTION AND THE NEWLY DISTRIBUTED POPULATION
return
  Rcpp::List::create(
    Rcpp::Named("V1") = V1,
    Rcpp::Named("new_dist") = dist_new_fin
  );
}
