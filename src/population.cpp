////////CODE FOR THE DISTRIBUTION OF THE POPULATION AMONG
/////// THE TWO NEW GENERIC RISK DIMENSIONS

using namespace Rcpp;

//[[Rcpp::export]]


NumericMatrix trans_mat
NumericMatrix trans_matN
NumericMatrix trans_mat_tot

std::vector<double> dist_i_v
std::vector<double> dist_orig_v
std::vector<double> dist_goal_v
std::vector<double> dist_new
std::vector<double> dist_orig
std::vector<double> diff_i_v
std::vector<double> can_go
std::vector<double> did_go

int N  /// number of iterations to run the model
int r;
int c;
int i;
int n;
int m;
int p;

double frc;
////////////////////////////////////////////////////////////////////////////////
/////////                    FUNCTIONS USED BELOW                    ///////////
////////////////////////////////////////////////////////////////////////////////

/////////           FUNCTION TO RETURN DIAGONAL OF A MATRIX            /////////
std::vector<double> diagonal (matrix, diagonal){

    std::vector<double> diagonal;

    for(int i=0; i<sizeof(matrix); i++){
        for(int j=0; j<sizeof(matrix); j++){
            if (i == j)
                diagonal[i] = matrix[i,j];
                
    return diagonal;
            
        }
    }
}
//////////           FUNCTION TO CALCULATE ROW SUMS OF MATRIX         //////////
std::vector<double> row_sum(matrix){
    int row_sum;
    
    for(int i=0; i<sizeof(matrix); i++){
        for(int j=0; j<sizeof(matrix); j++){
            
            row_sum += matrix[i,j];
            
            return row_sum;
        }
    }
}
//////////      FUNCTION TO CALCULATE SUM OF ALL VALUES IN MATRIX     //////////
int sum (matrix) {
    
    int sum;

    for(int i=0; i<sizeof(matrix); i++){
        for(int j=0; j<sizeof(matrix); j++){
            
            sum += matrix[i,j];
            
            return sum;
        }
    }
}
//////////                         INITIALIZE                         //////////

    for(int i=0; i<trans_mat.nrow(); i++) {
        for(int j=0; j<trans_mat.ncol(); j++) {
            trans_matN[i][j] = trans_mat(i,j);
    } }

for(int n=0; n<N; n++){
/////// CALCULATE DISTANCE FROM CURRENT DISTRIBUTION TO GOAL DISTRIBUTION  /////
        diff_i_v = dist_i_v - dist_goal_v;

//////////                  CREATE TRANSITION MATRIX                    ////////
    for (int r=0; r<16; r++){
        for (int c=0; c<16; c++){
            trans_mat[r,c] = can_go[r,c]*max(0,(diff_i_v[r]-diff_i_v[c]));
        }
    }

//////////                ADJUST THE TRANSITION MATRIX                  ////////
//////////   1ST SCALE UP RATES, 2ND MAKE SURE DOES NOT SUM OVER 1    //////////
    frc = 0.1;  // approach seems quite sensitive to this value, = fraction of change to

    for(int i=0; i<16; i++){
        trans_matN[i,] =  trans_matN[i,] / dist_i_v[i]*frc;
        trans_matN[i,] =  trans_matN[i,] / std::max(1.0,sum(trans_matN[i,])); // should be dist0_v?
    }
    
//////////                      FINALIZE TRANS_MAT                    //////////
            diagonal(trans_matN) = 1-row_sums(trans_matN);
    
//////////                RECORD ABSOLUTE TRANSITIONS                 //////////
    for(i=0; i<16; i++){
        did_go[i,] += dist_i_v[i]*trans_mat[i,];
    }
    diag(did_go) = 0
//////////               UPDATE THE DISTRIBUTION VECTOR             ////////////
    dist_i_v = dist_i_v*trans_matN;

//////////                    NOW UPDATE IN ONE STEP                 ///////////
    trans_mat_tot = did_go;
    
    for(i=0; i<16; i++){
        trans_mat_tot[i,] = did_go[i,] / dist_orig_v[i,];
    }

    diagonal(trans_mat_tot) = 1 - row_sums(trans_mat_tot);
    
//////////           NOW FINALLY UPDATE THE DISTRIBUTION           ///////////
    dist_new = dist_orig

    dist_new [,] = 0;
    for (m=0; m<4; m++){
        for (p=0; p<4; p++){
            for (m2=0; m2<4; m2++){
                for (p2=0; p2<4; p2++){
                    dist_new[m+1,p+1] += dist_orig[m2+1,p2+1]*trans_mat_tot[1+m2+p2*4,1+m+p*4];
                }
            }
        }
    }

    matrix_sum(pow(dist_goal - dist_new, 2)) /
            matrix_sum(pow(dist_goal - dist_orig, 2));
}
