
////////THIS CODE ALLOWS FOR THE DISTRIBUTION OF THE POPULATION AMONG
/////// THE TWO NEW GENERIC RISK DIMENSIONS

NumericMatrix trans_mat
double dist_i_v
double dist_orig_v
double diff_i_v
double can_go
int N  /// number of iterations to run the model
int r
int c
int i
int n
double frc

/////////       FUNCTIONS USED BELOW ///////////

std::vector


//////////        INITIALIZE         //////////

for(int i=0; i<trans_mat.nrow(); i++) {
    for(int j=0; j<trans_mat.ncol(); j++) {
        trans_matN[i][j] = trans_mat(i,j);
    } }

for(int n=0; n<N; n++){
/////// CALCULATE DISTANCE FROM CURRENT DISTRIBUTION TO GOAL DISTRIBUTION //////

    diff_i_v = dist_i_v - dist_goal_v

////////// CREATE TRANSITION MATRIX //////////
    
for (int r=0; r<16; r++){
    for (int c=0; c<16; c++){
        trans_mat[r,c] = can_go[r,c]*max(0,(diff_i_v[r]-diff_i_v[c]))
    }
}

////////// ADJUST THE TRANSITION MATRIX
////////// 1ST SCALE UP RATES, 2ND MAKE SURE DOES NOT SUM OVER 1 //////////
frc = 0.1  // approach seems quite sensitive to this value, = fraction of change to

for(int i=0; i<16; i++){
    trans_matN[i,] =  trans_matN[i,] / dist_i_v[i]*frc
    trans_matN[i,] =  trans_matN[i,] / max(1.0,sum(trans_matN[i,])) // should be dist0_v?
}
   
//////////  FINALIZE TRANS_MAT
    for(int i=0; i<sizeof(trans_matN); i++){
        for(int j=0; j<sizeof(trans_matN); j++){
        if (i == j)
            diagonal[i] = trans_matN[i,j]
        } }
}
