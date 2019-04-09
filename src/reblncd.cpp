#include <Rcpp.h>

using namespace Rcpp;

//'@title reblncd
//'@description this function calculates a matrix of transition probabilities for rebalancing in the model.
//'@param IP dataframe of formatted parameters
//'@return trans_mat_tot_ages
//[[Rcpp::export]]

Rcpp::NumericMatrix reblncd(
    Rcpp::NumericMatrix mubt,
    Rcpp::NumericMatrix can_go,
    double RRmuHR,
    std::vector<double> RRmuRF,
    std::vector<double> HRdist,
    std::vector<double> dist_gen_v,
    std::vector<double> adj_fact
){
  double mubtN[11];
  double dist_i_v[4];
  double diff_i_v[4];
  double dist_t1_v[4];
  double temp_vec[4];
  double row_sum[4];
  double can_goN[4][4];
  double did_go[4][4];
  double trans_mat[4][4];
  double trans_mat_tot[4][4];
  double trans_mat_tot_ag[4][44];
  Rcpp::NumericMatrix trans_mat_tot_ages(4,44);
  double frc;
  int mat_sum; int N;


  /////Initialize the above objects

  for(int i=0;i<4;i++){
    dist_i_v[i]=0;
    row_sum[i]=0;
    temp_vec[i]=0;
    for (int j=0;j<4;j++){
      can_goN[i][j]=can_go(i,j);
      did_go[i][j]=0;
      trans_mat[i][j]=0;
      trans_mat_tot[i][j]=0;
    }
    for (int k=0; k<44; k++){
      trans_mat_tot_ag[i][k]=0;
    }
  }

  for (int i=0; i<11; i++){
    mubtN[i]=mubt(0,i);
  }

  frc=0.01;
  mat_sum=0;
  N=0;


  ////////////////////////////////////////////////////////////////////////////////////////

  //'Open the Age Loop to Calculate Age Specific Transition Matrices
  for(int ag=0; ag<11; ag++){
    //' reset the appropriate age specific variables
    for(int i=0;i<4;i++){
      dist_t1_v[i]=0;
      dist_i_v[i]=0;
      for (int j=0;j<4;j++){
        did_go[i][j]=0;
        trans_mat_tot[i][j]=0;
      } }
    mat_sum=0;
    //' Calculate the dist_t1_v
    for (int nm=0; nm<4; nm++){

        // if ((ag<9) & ((RRmuRF[nm]*RRmuHR)<5)){
          dist_t1_v[nm]=dist_gen_v[nm]*(1-((mubtN[ag]*RRmuRF[nm]*RRmuHR*HRdist[ag])));

        // } else {
        //   dist_t1_v[nm+im*4]=dist_gen_v[nm+im*4]*(1-((mubtN[ag]*5*HRdist[ag])));
        //
        // }

      }


    for (int i=0; i<4; i++){
      // Rcpp::Rcout <<"dist_t1_v at ag = "<< ag << "and index = "<< i << " is "<<  dist_t1_v[i]<< "\n";
      dist_i_v[i]=dist_t1_v[i];
    }
    //' Open the iteration loop

    N=100;

    for (int n=0; n<N; n++){
      //'calculate the difference between dist_gen and current dist
      for (int i=0; i<4; i++){
        diff_i_v[i] = dist_i_v[i] - dist_gen_v[i];
        // Rcpp::Rcout <<"diff_i_v at ag = "<< ag << "and index = "<< i << " is "<<  diff_i_v[i]<< "\n";

      }
      //'create the transition matrix
      //'initialize to zero for each loop
      for (int r=0; r<4; r++){
        for (int c=0; c<4; c++){
          trans_mat[r][c] = 0;
        } }
      //'Set the Value of trans_mat to the max of 0 and the difference between row & column
      for (int r=0; r<4; r++){
        for (int c=0; c<4; c++){
          if ((diff_i_v[r]-diff_i_v[c]) > 0.0) {
            trans_mat[r][c] = can_goN[r][c]*(diff_i_v[r]-diff_i_v[c]);
          }
        } }

      //'Adjust the Transition Matrix
      //'1st scale up rates, 2nd make sure that it does not sum over one
      frc = 0.01;  // approach seems quite sensitive to this value, = fraction of change to
      for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
          trans_mat[i][j] =  (trans_mat[i][j] /(dist_i_v[i]+1e-200))*frc; //makes this number bigger
        } }

      //'Calculate Row Sums
      for(int i=0; i<4; i++){
        row_sum[i]=0; //reset row sum
        for(int j=0; j<4; j++){
          row_sum[i]+=trans_mat[i][j];
        } }
      ////////This is the step that is failing!
      for(int i=0; i<4; i++){
        if ((1-row_sum[i])<0.0){ //max of 1 and sum(trans_mat)
          for(int j=0; j<4; j++){
            trans_mat[i][j] =  trans_mat[i][j] / (row_sum[i]+1e-200); //dividing by the row sum guarantees that the new row sum does not exceed 1
          }
        } }

      //'Finalize trans_mat
      for(int i=0; i<4; i++){
        row_sum[i]=0; //reset row sum
        for(int j=0; j<4; j++){
          row_sum[i]+=trans_mat[i][j]; //calculate the row sum of temp mat above
        } }

      // for(int i=0; i<4; i++){
      //
      //   if(row_sum[i]>1){
      //     row_sum[i]=1;
      //   } }

      //this step is what causes the negative
      for(int i=0; i<4; i++){
        if (row_sum[i]<0) {
          Rcpp::Rcout<<"row sum is negative at ag = " << ag<< " and i = "<< i <<" and n = "<< n <<"\n";
        }
        if (row_sum[i]>1) {
          Rcpp::Rcout<<"row sum is too large at ag = " << ag<< " and i = "<< i <<" and n = "<< n <<"rowsum is = " << 1-row_sum[i]<<"\n";
        }
        for(int j=0; j<4; j++){
          if (i==j){
            trans_mat[i][j]=(1-row_sum[i]);
          }
        } }

      for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
          if (trans_mat[i][j]<0) {
            Rcpp::Rcout<<"trans_mat is negative at ag = " << ag<< " and i = "<< i <<" and j = "<< j <<" and n = "<< n <<"\n";
          } } }
      //'Record Absolute Transitions
      for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
          // if (trans_mat[i][j]<0) {
          //   Rcpp::Rcout<<"trans_mat is negative at ag = " << ag<< " and i = "<< i <<" and j = "<< j <<" and n = "<< n <<"\n";
          // }

          did_go[i][j] += (dist_i_v[i]*trans_mat[i][j]);
        } }



      for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
          if (i==j){
            did_go[i][j]=0;
          }
        } }
      //'Update the distribution vector
      //'it is supposed to be matrix multiplication
      for(int c=0; c<4; c++){
        temp_vec[c] = 0;
      }
      for(int r=0; r<4; r++){
        for(int c=0; c<4; c++){
          temp_vec[c] +=dist_i_v[r]*trans_mat[r][c];
        } }
      for(int c=0; c<4; c++){
        dist_i_v[c] = temp_vec[c];
      }
    }//end of iteration loop
    //'create trans_mat_tot
    for(int i=0; i<4; i++){
      for(int j=0; j<4; j++){
        trans_mat_tot[i][j] = did_go[i][j] / (dist_gen_v[i]);
      } }
    for(int i=0; i<4; i++){
      row_sum[i]=0;
      for(int j=0; j<4; j++){
        row_sum[i] +=trans_mat_tot[i][j]*adj_fact[ag]; ///calculate the number of total transitions from a state
      } }
    for(int i=0; i<4; i++){
      if (row_sum[i]>1){
        Rcpp::Rcout<<"row sum 2 is too large at ag = " << ag<< " and i = "<< i  << "& Rs= "<< row_sum[i]<<"\n";
      }
      for(int j=0; j<4; j++){
        if (i==j){
          trans_mat_tot[i][j]=(1-row_sum[i]); //the remainder of those that did not transition to another state
        } } }
    for(int i=0; i<4; i++){
      for(int j=0; j<4; j++){
        trans_mat_tot[i][j]=trans_mat_tot[i][j];//*adj_fact[ag];
      } }


    for(int i=0; i<4; i++){
      for(int j=0; j<4; j++){
        // Rcpp::Rcout<<(4*ag)+j<<"\n";
        // Rcpp::Rcout<<"ag="<<ag<< "i="<< i <<"j="<< j << trans_mat_tot[i][j] <<"\n";
        trans_mat_tot_ag[i][(4*ag)+j]=trans_mat_tot[i][j];
      } }
  } //end of age loop

  for(int i=0; i<4; i++){
    for(int j=0; j<44; j++){
      trans_mat_tot_ages(i,j)= trans_mat_tot_ag[i][j];
    } }

  return trans_mat_tot_ages;
}//end of function
