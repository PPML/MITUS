#include <Rcpp.h>
using namespace Rcpp;
//'@name cSim_demo_ag
//'@description a simple demographic model with inputs of cSim
//'@param nYrs number of years to run the model.
//'@param nRes number of results of the model
//'@param InitPop Initial Population matrix
//'@param Birthst Births over time
//'@param ImmTot Immigration over time
//'@param mubt background mortality over time
//'@param ag_den denominator used in the aging process
//'@return a list of outputs
//[[Rcpp::export]]

Rcpp::List cSim_demo_ag(
    int                       nYrs,
    int                       nRes,
    Rcpp::NumericMatrix       InitPop,
    std::vector<double>       Birthst,
    Rcpp::NumericMatrix       ImmNon,
    Rcpp::NumericMatrix       ImmLat,
    Rcpp::NumericMatrix       ImmAct,
    Rcpp::NumericMatrix       ImmFst,
    Rcpp::NumericMatrix       mubt,
    Rcpp::NumericMatrix       ag_den


){
  int           s;
  double        temp;
  double        temp2;
  double        InitPopN[InitPop.nrow()][InitPop.ncol()];
  double        ImmTotN[ImmNon.nrow()][ImmNon.ncol()];
  double        mubtN[mubt.nrow()][mubt.ncol()];
  double        Outputs[nYrs][nRes];
  Rcpp::NumericMatrix Outputs2(nYrs,nRes);
  double        ag_denN[ag_den.nrow()][ag_den.ncol()];
  double   V0[11];
  double   V1[11];
  double   VMort[11];

///////////////////////////////////////////////////////////////////////////////
///////                            INITIALIZE                             /////
///////////////////////////////////////////////////////////////////////////////
for(int i=0; i<InitPop.nrow(); i++) {
  for(int j=0; j<InitPop.ncol(); j++) {
    InitPopN[i][j] = InitPop(i,j);
  } }

for(int i=0; i<ImmNon.nrow(); i++) {
  for(int j=0; j<ImmNon.ncol(); j++) {
    ImmTotN[i][j] = ImmNon(i,j)+ImmLat(i,j)+ImmFst(i,j)+ImmAct(i,j);
  } }

for(int i=0; i<mubt.nrow(); i++) {
  for(int j=0; j<mubt.ncol(); j++) {
    mubtN[i][j] = mubt(i,j);

  } }
for(int i=0; i<ag_den.nrow(); i++) {
  for(int j=0; j<ag_den.ncol(); j++) {
    ag_denN[i][j] = ag_den(i,j);
  } }

for(int ag=0; ag<11; ag++) {
              V0[ag] = 0;
              V1[ag] = 0;
              VMort[ag] = 0;
  }
for(int i=0; i<nYrs; i++) {
  for(int j=0; j<nRes; j++) {
    Outputs[i][j] = 0;
  }  }
temp2=0; temp=0;

////////////////////////////////////////////////////////////////////////////////
//////                             StatList                                /////
//////                             BURN IN                                 /////
//////                           Populate model                            /////
////////////////////////////////////////////////////////////////////////////////
// set the vector equal to initial population
for(int ag=0; ag<11; ag++) {
  V0[ag] = InitPopN[ag][0]+InitPopN[ag][1];
}
for(int ag=0; ag<11; ag++) {
  V1[ag] = V0[ag];
}
//Main Model Loop
for(int y=0; y<nYrs; y++) {
  /////////////////////////////////MONTH LOOP////////////////////////////////////
  for(int m=0; m<12; m++) {
    /////////////////////CREATE A COUNTER OF MONTHS SINCE START////////////////////
    s = y*12+m;
// Births
  V1[0] += Birthst[s];
// Immigration
for(int ag=0; ag<11; ag++) {
  V1[ag] += ImmTotN[s][ag];
}
// Emmigration
  // for(int ag=0; ag<11; ag++) {
  //   V1[ag] = V1[ag]*.995;
  // }
// Mortality
for(int ag=0; ag<11; ag++) {
  VMort[ag] = V0[ag]*mubtN[s][ag];

  // VMort[ag] = ((1-p_HR)*V0[ag]*mubtN[s][ag])+(p_HR*V0[ag]*mubtN[s][ag]*RRmuHR[1]);
}
// for(int ag=9; ag<11; ag++) {
//   VMort[ag] = V0[ag]*mubtN[s][ag]*1.2;
//
//   // VMort[ag] = ((1-p_HR)*V0[ag]*mubtN[s][ag])+(p_HR*V0[ag]*mubtN[s][ag]*RRmuHR[1]);
// }
for(int ag=0; ag<11; ag++) {
V1[ag]   -= VMort[ag];
  // Rcpp::Rcout << "for ag =" << ag << "and t = "<< s<< "mort is" << VMort[ag] << "\n";
}




// Aging
for(int ag=0; ag<10; ag++) {
/////          IF AGE > 4, IT TAKES 120 MONTHS TO LEAVE AGE GROUP          /////
//   if(ag>0) {
//     temp2 = 120;
// /////          IF AGE < 4, IT TAKES 60 MONTHS TO LEAVE AGE GROUP           /////
//   } else {
//     temp2 = 60;
//   }
  temp=V0[ag]/ag_denN[s][ag];
  V1[ag]-=temp;
  V1[ag+1]+=temp;
}

// FILL IN THE OUTPUTS
if(m==6) {
  for(int ag=0; ag<11; ag++) {
    Outputs[y][0]      = y+1950;  // Year
    Outputs[y][1]     += V1[ag];
    Outputs[y][2+ag]   = V1[ag];
    Outputs[y][13+ag]  = VMort[ag]*12;
}
}
// Update V0 to V1
for(int ag=0; ag<11; ag++) {
  V0[ag] = V1[ag];
  // VMort[ag] =0;
}
  } } //end of year and month loops

// RETURN STUFF
for(int i=0; i<nYrs; i++) {
  for(int j=0; j<nRes; j++) {
    Outputs2(i,j)  = Outputs[i][j];
  }  }

return
  Rcpp::List::create(
    Rcpp::Named("Outputs") = Outputs2
  );
}
