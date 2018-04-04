//////input variables
#include <Rcpp.h>
//#include <algorithm>

using namespace Rcpp;


//[[Rcpp::export]]
Rcpp::List cRebal(
  Rcpp::NumericMatrix                    V1,
  Rcpp::NumericMatrix       dist_new,
  Rcpp::NumericMatrix       can_go,
  Rcpp::NumericMatrix       did_go,
  Rcpp::NumericMatrix       dist_goal,
  std::vector<double>       dist_goal_v,
  std::vector<double>       dist_orig_v,
  std::vector<double>       diff_i_v
) {

  ////////////////////////////////////////////////////////////////////////////////
  ////////    BELOW IS A LIST OF THE VARIABLES CREATED INTERNALLY IN MODEL   /////
  ////////////////////////////////////////////////////////////////////////////////
  int N; int m; int p; int r; int c; int m2; int p2;
  double        temp;
//  double   dist_genN[dist_gen.nrow()][dist_gen.ncol()];
  double   dist_newN[dist_new.nrow()][dist_new.ncol()];
  double   can_goN[can_go.nrow()][can_go.ncol()];
  double   did_goN[did_go.nrow()][did_go.ncol()];
  double   dist_goalN[dist_goal.nrow()][dist_goal.ncol()];
  double   temp_vec[16];
  double     temp_mat[16][16];
  double    temp_mat2[16][16];
  double     trans_mat[16][16];
  double    trans_mat_tot[16][16];
  double    dist_i_v[16];
  double    dist_orig[4][4];
 Rcpp::NumericMatrix    dist_new_fin(4,4);
 Rcpp::NumericMatrix V1_fin(4,4);
  double    V1N[4][4];
  double        frc;
  double        pop_t;
  double        sse;
  double colsum[16];
  double rowsum[16];

  ///////////////////////////////////////////////////////////////////////////////
  ///////                            INITIALIZE                             /////
  ///////////////////////////////////////////////////////////////////////////////
  N=30;
  pop_t =0;

  // for(int i=0; i<dist_gen.nrow(); i++) {
  //   for(int j=0; j<dist_gen.ncol(); j++) {
  //     dist_genN[i][j] = dist_gen(i,j);
  //   } }

  for(int i=0; i<dist_new.nrow(); i++) {
    for(int j=0; j<dist_new.ncol(); j++) {
      dist_newN[i][j] = dist_new(i,j);
    } }

  for(int i=0; i<dist_goal.nrow(); i++) {
    for(int j=0; j<dist_goal.ncol(); j++) {
      dist_goalN[i][j] = dist_new(i,j);
    } }

  for(int i=0; i<can_go.nrow(); i++) {
    for(int j=0; j<can_go.ncol(); j++) {
      can_goN[i][j] = can_go(i,j);
    } }

  for(int i=0; i<did_go.nrow(); i++) {
    for(int j=0; j<did_go.ncol(); j++) {
      did_goN[i][j] = did_go(i,j);
    } }

  for(int i=0; i<V1.nrow(); i++) {
    for(int j=0; j<V1.ncol(); j++) {
      V1N[i][j] = V1(i,j);
    } }

  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      dist_orig[i][j] = 0;
    } }

  for(int i=0; i<16; i++) {
    for(int j=0; j<16; j++) {
      temp_mat[i][j] = 0;
      temp_mat2[i][j] = 0;
      trans_mat[i][j] = 0;
      trans_mat_tot[i][j] = 0;
    } }

  for(int i=0; i<16; i++) {
    dist_i_v[i]=0;
    temp_vec[i]=0;
    colsum[i]=0;
    rowsum[i]=0;
    }

//////////////////////////////////////////////////////////////////////////////////
  ////////////////////////// REBALANCE THE POPULATION //////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////

  ////// need to define the current distribution of persons across the RG at this timestep
////// first calculate the total population at this time step

////// second calculate the distribution of population across the two risk factors
// for(int ag=0; ag<11; ag++) {
//   for(int tb=0; tb<6; tb++) {
//     for(int lt=0; lt<2; lt++) {
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++) {
          // for(int rg=0; rg<2; rg++) {
          //   for(int na=0; na<3; na++) {
              pop_t += V1N[im][nm];
              dist_orig[nm][im] = V1N[im][nm];
              // pop_t             += V1[ag][tb][lt][im][nm][rg][na];
              // dist_orig[nm][im] += V1[ag][tb][lt][im][nm][rg][na]; //sum across risk factors
              dist_orig[nm][im]  = dist_orig[nm][im]/pop_t; // determine the proportions

              ///insert error that distribution does not sum to one and then go from there;
              //Rcpp::Rcout <<"dist_orig at"<< im << "and"<< nm << "is" << dist_orig[nm][im]<< "\n";
            } }
      // } } } } }


for(int n=0; n<N; n++){

  /////// CALCULATE DISTANCE FROM CURRENT DISTRIBUTION TO GOAL DISTRIBUTION /////
    for (int i=0; i<sizeof(dist_i_v); i++){
      diff_i_v[i] = dist_i_v[i] - dist_goal_v[i];
    }
  //////////                  CREATE TRANSITION MATRIX                    ////////
    for (int r=0; r<16; r++){
      for (int c=0; c<16; c++){
        temp=diff_i_v[r]-diff_i_v[c];
        trans_mat[r][c] = can_goN[r][c]*(std::max(0.0,temp));
        Rcpp::Rcout << "initial trans_mat is" << trans_mat;
      }
    }
  //////////                ADJUST THE TRANSITION MATRIX                  ////////
    //////////   1ST SCALE UP RATES, 2ND MAKE SURE DOES NOT SUM OVER 1    //////////
    frc = 0.1;  // approach seems quite sensitive to this value, = fraction of change to

    for(int j=0; j<16; j++){
        colsum[j]=0;
      for(int i=0; i<16; i++){
        colsum[j] += trans_mat[i][j];
     ///needs to be the column sums
        temp_mat[i][j] =  trans_mat[i][j] / dist_i_v[i]*frc;
        temp_mat[i][j] =  temp_mat [i][j] / (std::max(1.0, colsum[j]));
      }}
    Rcpp::Rcout << "temp trans_mat is" << trans_mat;

    //////////                      FINALIZE TRANS_MAT                    //////////
      for(int i=0; i<16; i++){
        rowsum[i]=0;
        for(int j=0; j<16; j++){
          rowsum[i]+=temp_mat[i][j];

          if (i==j){
            temp_mat2[i][j]=rowsum[i];
          } else {
            temp_mat2[i][j]=0;
          }

          if (i != j) {
            trans_mat[i][j] = temp_mat[i][j];
          } else {trans_mat_tot[i][j]=temp_mat2[i][j];
          } } }
    //////////                RECORD ABSOLUTE TRANSITIONS                 //////////
      for(int i=0; i<16; i++){
        for(int j=0; j<16; j++){
          if (i==j){
            did_goN[i][j]=0;
          } else {
            did_goN[i][j] += dist_i_v[i]*trans_mat[i][j];
          }
          //////////               UPDATE THE DISTRIBUTION VECTOR             ////////////
            dist_i_v[i] = dist_i_v[i]*trans_mat[i][j];
        }}
} //end of N loop
    //////////                    NOW UPDATE IN ONE STEP                 ///////////
    for(int i=0; i<16; i++){
      for(int j=0; j<16; j++){
        temp_mat[i][j] = did_goN[i][j];
        temp_mat[i][j] = did_goN[i][j] / dist_orig_v[i];
      } }
    Rcpp::Rcout<< "temp_mat is" << temp;

    ///Use the diagmat functin of RcppArmadillo to create a diagonal matrix
    ///with the vector of values defined below.
    for(int i=0; i<16; i++){
      rowsum[i]=0;
      for(int j=0; j<16; j++){
    rowsum[i] +=trans_mat_tot[i][j]; //rowsum
        if (i==j){
          trans_mat_tot[i][j]=1-rowsum[i];
        } else {
          trans_mat_tot[i][j]=0;
        }
      } }
    ///Replace non-diagonal values with the values from temp above

    for(int i=0; i<16; i++){
      for(int j=0; j<16; j++){
        if (i != j) {
          trans_mat_tot[i][j] = temp_mat[i][j];
        } else {trans_mat_tot[i][j]=trans_mat_tot[i][j];
        } } }

    Rcpp::Rcout<< "trans_mat_tot is" << trans_mat_tot;

    //////////           NOW FINALLY UPDATE THE DISTRIBUTION           ///////////
      for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
          dist_newN[i][j] = dist_orig[i][j];
        } }
    for (int m=0; m<4; m++){
      for (int p=0; p<4; p++){
        for (int m2=0; m2<4; m2++){
          for (int p2=0; p2<4; p2++){
            ////scalar multiplication (not matrix multiplication)
            dist_newN[m+1][p+1] += dist_orig[m2+1][p2+1] * trans_mat_tot[1+m2+p2*4][1+m+p*4];
          } } } }
    // for(int i=0; i<16; i++){
    //   for(int j=0; j<16; j++){
    //     temp_mat[i][j] = dist_goal[i][j];
    //     sum1=
    //     sse = arma::accu(pow(dist_goal[i][j] - dist_new[i][j], 2)) /
    //       arma::accu(pow(dist_goal[i][j] - dist_orig[i][j], 2));
    //   }}
    // Rcpp::Rcout << "sse is" << sse;


/////////////APPLY THE NEW DISTRIBUTION TO THE POPULATION /////////////////////
  /////////////revert the arma::vec to a double;
  for (int im=0; im<4; im++){
    for (int nm=0; nm<4; nm++){
      dist_new_fin(nm,im) = dist_newN[nm][im];
    } }

// for(int ag=0; ag<11; ag++) {
//   for(int tb=0; tb<6; tb++) {
//     for(int lt=0; lt<2; lt++){
      for (int im=0; im<4; im++){
        for (int nm=0; nm<4; nm++){
          // for(int rg=0; rg<2; rg++) {
          //   for(int na=0; na<3; na++){
              V1N[im][nm] = V1N[im][nm]*dist_newN[nm][im];
            } }

      for (int im=0; im<4; im++){
        for (int nm=0; nm<4; nm++){
          V1_fin(nm,im) = V1N[nm][im];
        } }
///return stuff
return
  Rcpp::List::create(
    Rcpp::Named("V1") = V1_fin,
    Rcpp::Named("new_dist") = dist_new_fin
  );

}
