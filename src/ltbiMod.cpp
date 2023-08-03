#include "ltbiMod.h"
SEXP ltbiMod(
  int           s,
  double  	    V0,
  double  	    V1,
  double        rLtScrtN,
  double        LtTxParN,
  double        LtDxPar_ltN,
  double        LtDxPar_noltN,
  double        temp,
  double        temp2,
  double        temp3,
  double        temp4,
  double        rr_ltbi,
  double        rrTestLrNoTb,
  double        rrTestHr,
  double        VLdx,
  double        VLtest,
  double        VLtemp,
  double        VLDxtemp,
  double        Vinf
){
  using namespace Rcpp ;
  ///////////////////////////////////////////////////////////////////////////////
  /////                  LTBI SCREENING AND TLTBI INITIATION                /////
  ///////////////////////////////////////////////////////////////////////////////
  ///// BECAUSE WE CALCULATE THE NUMBER OF SUSCEPTIBLES BASED ON EXISTING   /////
  ///// POPULATION SIZES THE ORDER OF BASELINE VS. TTT SCREENING MATTERS    /////
  ///// BELOW WE EXPLORE HAVING THE TTT FIRST THEN THE BASELINE
  ////////////////////////////////////////////////////////////////////////////////
  /////                     BASE LINE VALUE CALCULATIONS                     /////
  ////////////////////////////////////////////////////////////////////////////////
  ///// calculate the rates of tb positive for the basecase and the rate of  /////
  ///// TB negatives these are equal to the rate of screening combined with  /////
  ///// the sensitivity or specificity in the base case                      /////
  ///// tempHR is adjustment to add extra screening for HR in param          /////
  ////////////////////////////////////////////////////////////////////////////////

  ///// Variables created in the model
  // double        rTbP;
  // double        rTbN;
  // int tempHR = 1;
  // int tempNTB = 1;
  // for(int rg=0; rg<2; rg++) {
  //   for(int na=0; na<3; na++) {
  //     for(int ag=0; ag<11; ag++) {
  //       ///// US BORN, LOW RISK  ////////////////////////////////////////////////
  //       if(rg==0 && na==0) {
  //         rTbP = rLtScrtN[s][0]*LtDxPar_ltN[0][s];
  //         rTbN = rLtScrtN[s][0]*LtDxPar_noltN[0][s];
  //         tempHR=1;
  //         tempNTB = rrTestLrNoTb;
  //       }
  //       ///// US BORN, HIGH RISK  ///////////////////////////////////////////////
  //       if(rg==1 && na==0) {
  //         rTbP = rLtScrtN[s][1]*LtDxPar_ltN[1][s];
  //         rTbN = rLtScrtN[s][1]*LtDxPar_noltN[1][s];
  //         tempHR = rrTestHr;
  //         tempNTB = 1;
  //       }
  //       ///// Young NUS (under 5)  ////////////////////////////////////////////
  //       if(rg==0 && na > 0 && ag==0) {
  //         rTbP = rLtScrtN[s][0]*LtDxPar_ltN[2][s];
  //         rTbN = rLtScrtN[s][0]*LtDxPar_noltN[2][s];
  //         tempHR = 1;
  //         tempNTB = 1;
  //       }
  //       ///// NON US BORN  ////////////////////////////////////////////////////
  //       if(rg==0 && na > 0 && ag > 0) {
  //         rTbP = rLtScrtN[s][0]*LtDxPar_ltN[3][s];
  //         rTbN = rLtScrtN[s][0]*LtDxPar_noltN[3][s];
  //         tempHR = 1;
  //         tempNTB = 1;
  //       }
  //       ///// NON US BORN, HIGH RISK  //////////////////////////////////////////
  //       if(rg==1 && na >0) {
  //         rTbP = rLtScrtN[s][1]*LtDxPar_ltN[4][s];
  //         rTbN = rLtScrtN[s][1]*LtDxPar_noltN[4][s];
  //         tempHR = rrTestHr;
  //         tempNTB = 1;
        // }
  // }}}
  // return V1;
}
