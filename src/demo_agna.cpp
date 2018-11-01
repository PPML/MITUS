#include <Rcpp.h>
using namespace Rcpp;
//'@name cSim_demo_agna
//'@description a simple demographic model with inputs of cSim
//'@param nYrs number of years to run the model.
//'@param nRes number of results of the model
//'@param InitPop Initial Population matrix
//'@param Birthst Births over time
//'@param ImmTot Immigration over time
//'@param mubt background mortality over time
//'@return a list of outputs
//[[Rcpp::export]]

Rcpp::List cSim_demo_agna(
int                       nYrs,
int                       nRes,
Rcpp::NumericMatrix       InitPop,
std::vector<double>       Birthst,
Rcpp::NumericMatrix       ImmNon,
Rcpp::NumericMatrix       ImmLat,
Rcpp::NumericMatrix       ImmAct,
Rcpp::NumericMatrix       ImmFst,
std::vector<double>       rEmmigFB,
Rcpp::NumericMatrix       mubt
){
int           s;
double        temp;
double        temp2;
double        InitPopN[InitPop.nrow()][InitPop.ncol()];
double        ImmTotN[ImmNon.nrow()][ImmNon.ncol()];
double        mubtN[mubt.nrow()][mubt.ncol()];
double        Outputs[nYrs][nRes];
Rcpp::NumericMatrix Outputs2(nYrs,nRes);
double   V0[11][3];
double   V1[11][3];
double   VMort[11][3];

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

for(int ag=0; ag<11; ag++) {
  for(int na=0; na<3; na++) {

V0[ag][na] = 0;
V1[ag][na] = 0;
VMort[ag][na] = 0;
  } }
temp2=0; temp=0;

////////////////////////////////////////////////////////////////////////////////
//////                             StatList                                /////
//////                             BURN IN                                 /////
//////                           Populate model                            /////
////////////////////////////////////////////////////////////////////////////////
// set the vector equal to initial population
for(int ag=0; ag<11; ag++) {
V0[ag][0] = InitPopN[ag][0];
V0[ag][2] = InitPopN[ag][1];
}
for(int ag=0; ag<11; ag++) {
  for(int na=0; na<3; na++) {
V1[ag][na] = V0[ag][na];
} }
//Main Model Loop
for(int y=0; y<nYrs; y++) {
/////////////////////////////////MONTH LOOP////////////////////////////////////
for(int m=0; m<12; m++) {
/////////////////////CREATE A COUNTER OF MONTHS SINCE START////////////////////
s = y*12+m;

// Births
V1[0][0] += Birthst[s];
// Immigration
for(int ag=0; ag<11; ag++) {
V1[ag][1] += ImmTotN[s][ag];
  }
//Emmigration
// for(int ag=0; ag<11; ag++) {
//   V1[ag][1] -= V0[ag][1]*rEmmigFB[0];
//   V1[ag][2] -= V0[ag][2]*rEmmigFB[1];
// }
//long term immigrant
for(int ag=0; ag<11; ag++) {
  temp=V0[ag][1]/24;
  V1[ag][1] -= temp;
  V1[ag][2] += temp;
}
// Mortality
for(int ag=0; ag<11; ag++) {
  for(int na=0; na<3; na++) {
VMort[ag][na] += V0[ag][na]*mubtN[s][ag];
V1[ag][na]   -= VMort[ag][na];
if(s==50){

Rcpp::Rcout << "for ag =" << ag << "and t = "<< s<< "mort is" << VMort[ag][na] << "\n";

  } }}
// for(int ag=0; ag<11; ag++) {
//   for(int na=0; na<3; na++) {
//
// V1[ag][na]   -= VMort[ag][na];
// // Rcpp::Rcout << "for ag =" << ag << "and t = "<< s<< "mort is" << VMort[ag] << "\n";
//   } }

// Aging
for(int ag=0; ag<10; ag++) {
  /////          IF AGE > 4, IT TAKES 120 MONTHS TO LEAVE AGE GROUP          /////
  if(ag>0) {
    temp2 = 120;
    /////          IF AGE < 4, IT TAKES 60 MONTHS TO LEAVE AGE GROUP           /////
  } else {
    temp2 = 60;
  }
  for(int na=0; na<3; na++) {
temp=V0[ag][na]/temp2;
V1[ag][na]-=temp;
V1[ag+1][na]+=temp;
} }

// FILL IN THE OUTPUTS
if(m==6) {
for(int ag=0; ag<11; ag++) {
  for(int na=0; na<3; na++) {
Outputs[y][0]      = y+1950;  // Year
Outputs[y][1]     += V1[ag][na];
Outputs[y][2+ag]  += V1[ag][na];
Outputs[y][13+ag]  = VMort[ag][na]*12;
Outputs[y][24+na] += V1[ag][na];
  }} }
// Update V0 to V1
for(int ag=0; ag<11; ag++) {
  for(int na=0; na<3; na++) {
V0[ag][na] = V1[ag][na];
// VMort[ag] =0;
  } }
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
