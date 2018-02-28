#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
List cSim(
          int                 nYrs,      //number of years to run the model
          int                 nRes,      // results
          NumericMatrix       rDxt,      // rate of diagnosis over time
          std::vector<double> TxQualt,   // treatment quality over time
          NumericMatrix       InitPop,   // initial population matrix
          NumericMatrix       Mpfast,    // matrix of probability of fast TB progression
          std::vector<double> ExogInf,   // exogenous infection risk
          NumericMatrix       MpfastPI,  // matrix of probability of fast TB progression in those with partial immunity from prior infection
          NumericMatrix       Mrslow,    // matrix of the rate of slow TB progression
          std::vector<double> rrSlowFB,  // rate ratios that are applied to the rate of slow progression for foreign born population
          double              rfast,     // rate of fast TB progression
          double              RRcurDef,  // rate ratio of cure given treatment default
          double              rSlfCur,   // rate of self cure
          double              p_HR,
          NumericMatrix       dist,
          NumericMatrix       vTMort,    // matrix of TB mortality
          NumericMatrix       vRFMort,    // matrix of RF mortality
          std::vector<double> RRmuHR,
          double              muTbRF,    //factor for comorbidity btw TB and non-TB,
          std::vector<double> Birthst,   // vector of absolute births over time
          NumericMatrix       HrEntEx,
          NumericMatrix       ImmNon,
          NumericMatrix       ImmLat,
          NumericMatrix       ImmAct,
          NumericMatrix       ImmFst,
          NumericMatrix       mubt,       //background mortality over time
          NumericMatrix       RelInf,     //relative infectiousness
          std::vector<double> RelInfRg,   //relative infectiousness by risk group
          std::vector<double> Vmix,       // vector of mixing parameters (sigmas)
          NumericMatrix       rRFt,        //rate of risk factor (population) of interest over time
          std::vector<double> rEmmigFB,    // rate of emmigration among the foreign born
          std::vector<double> TxVec,       //vector of TB treatment parameters
          double              TunTxMort,   //tuning for treatment mortality
          std::vector<double> rDeft,       // rate of treatment default over time
          std::vector<double> rLtScrt,
          std::vector<double> LtTxPar,     // latent treatment parameters
          NumericMatrix       LtDxPar,     // latent diagnosis parameters
          std::vector<double> RRdxAge,     // rate ratios for diagnosis with respect to age
          double              rRecov,      //rate from latent slow to partially immune TB
          double              pImmScen,    // lack of reactivitiy to IGRA for Sp
          std::vector<double>  EarlyTrend,  // TB natural history parameter
          std::vector<double> dLtt,        //latent diagnosis over time
          std::vector<double> pReTx,
          std::vector<double> NixTrans,
          NumericMatrix       EffLtX,
          double              EffLt
) {
    ////////////////////////////////////////////////////////////////////////////////
    ////////    BELOW IS A LIST OF THE VARIABLES CREATED INTERNALLY IN MODEL   /////
    ////////////////////////////////////////////////////////////////////////////////
    int           ti;
    int           ti2;
    int           s;
    int           extrV[10];
    double        ExogInfN[ExogInf.size()];
    double        InitPopN[InitPop.nrow()][InitPop.ncol()];
    double        InitPopZ[InitPop.nrow()][InitPop.ncol()];
    double        MpfastN[Mpfast.nrow()][Mpfast.ncol()];
    double        MpslowN[Mpfast.nrow()][Mpfast.ncol()];
    double        MpfastPIN[MpfastPI.nrow()][MpfastPI.ncol()];
    double        MpslowPIN[MpfastPI.nrow()][MpfastPI.ncol()];
    double        MrslowN[Mrslow.nrow()][Mrslow.ncol()];
    double        vTMortN[vTMort.nrow()][vTMort.ncol()];
    double        vRFMortN[vRFMort.nrow()][vRFMort.ncol()];
   // double        vIsxtoIsyN[vIsxtoIsy.nrow()][vIsxtoIsy.ncol()];
   // double        vNmxtoNmyN[vNmxtoNmy.nrow()][vNmxtoNmy.ncol()];
   // double        vrgxtorgyN[vrgxtorgy.nrow()][vrgxtorgy.ncol()];
    double        distN[dist.nrow()][dist.ncol()];
    double        HrEntExN[HrEntEx.nrow()][HrEntEx.ncol()];
    double        ImmNonN[ImmNon.nrow()][ImmNon.ncol()];
    double        ImmLatN[ImmLat.nrow()][ImmLat.ncol()];
    double        ImmFstN[ImmFst.nrow()][ImmFst.ncol()];
    double        ImmActN[ImmAct.nrow()][ImmAct.ncol()];
    double        TBImm[11][ImmAct.nrow()][2];
    double        mubtN[mubt.nrow()][mubt.ncol()];
    double        RelInfN[6];
  //  double        rIntvInitN[rIntvInit.nrow()][rIntvInit.ncol()];
    double        rDxtN[rDxt.nrow()][rDxt.ncol()];
    double        TxVecN[TxVec.size()];
    double        LtDxParN[LtDxPar.nrow()][LtDxPar.ncol()];
    double        EffLtXN[EffLtX.size()];
    double        TxVecZ[6];
    double        temp;
    double        temp2;
    double        temp3;
    double        temp4;
    double        temp4V[11][5];
    double        rTbP;
    double        rTbN;
    double        Outputs[nYrs][nRes];
    double        V0[11][6][2][4][4][2][3];
    double        V1[11][6][2][4][4][2][3];
    double        VMort[11][6][2][4][4][2][3];
    double        Vdx[11][6][2][4][4][2][3];
    double        VNkl[2][2];  ///HIGH AND LOW RISK, NATIVITY
    double        VGjkl[2][2]; ///HIGH AND LOW RISK, NATIVITY
    double        Vjaf[4];     ///BY NUMBER OF MIXING GROUPS
    double        VLjkl[2][2];  ///HIGH AND LOW RISK, NATIVITY
    NumericMatrix Outputs2(nYrs,nRes);

///////////////////////////////////////////////////////////////////////////////
///////                            INITIALIZE                             /////
///////////////////////////////////////////////////////////////////////////////
for(int i=0; i<InitPop.nrow(); i++) {
        for(int j=0; j<InitPop.ncol(); j++) {
            InitPopN[i][j] = InitPop(i,j);
        } }
for(int i=0; i<ExogInf.size(); i++) {
            ExogInfN[i] = ExogInf[i];
        }
for(int i=0; i<Mpfast.nrow(); i++) {
        for(int j=0; j<Mpfast.ncol(); j++) {
            MpfastN[i][j]     = Mpfast(i,j);
            MpfastPIN[i][j]   = MpfastPI(i,j);
            MpslowN[i][j]     = 1-Mpfast(i,j);
            MpslowPIN[i][j]   = 1-MpfastPI(i,j);
            MrslowN[i][j]     = Mrslow(i,j);
        } }
for(int i=0; i<LtDxPar.nrow(); i++) {
        for(int j=0; j<LtDxPar.ncol(); j++) {
            LtDxParN[i][j] = LtDxPar(i,j);
        } }
for(int i=0; i<EffLtX.size(); i++) {
            EffLtXN[i] = EffLtX[i];
        }
for(int i=0; i<vTMort.nrow(); i++) {
        for(int j=0; j<vTMort.ncol(); j++) {
            vTMortN[i][j] = vTMort(i,j);
        } }
///this was HIV mortality temporarily replaced with RF of interest
///CHECK PARAM FILE; THIS WILL BECOME A VECTOR
for(int i=0; i<vRFMort.nrow(); i++) {
        for(int j=0; j<vRFMort.ncol(); j++) {
            vRFMortN[i][j] = vRFMort(i,j);
        } }
for(int i=0; i<ImmNon.nrow(); i++) {
        for(int j=0; j<ImmNon.ncol(); j++) {
            ImmNonN[i][j] = ImmNon(i,j);
        } }
for(int i=0; i<ImmLat.nrow(); i++) {
        for(int j=0; j<ImmLat.ncol(); j++) {
            ImmLatN[i][j] = ImmLat(i,j);
        } }
for(int i=0; i<ImmFst.nrow(); i++) {
        for(int j=0; j<ImmFst.ncol(); j++) {
            ImmFstN[i][j] = ImmFst(i,j);
        } }
for(int i=0; i<ImmAct.nrow(); i++) {
        for(int j=0; j<ImmAct.ncol(); j++) {
            ImmActN[i][j] = ImmAct(i,j);
        } }
////////do we want this to be dependent on treatment history (latent treatment)
for(int ag=0; ag<11; ag++) {
      for(int s=0; s< ImmAct.nrow(); s++) {
        TBImm[ag][s][0] = ImmActN[s][ag];
        TBImm[ag][s][1] = ImmFstN[s][ag];
      } }
for(int i=0; i<mubt.nrow(); i++) {
        for(int j=0; j<mubt.ncol(); j++) {
            mubtN[i][j] = mubt(i,j);
        } }
for(int i=0; i<RelInf.size(); i++) {
            RelInfN[i] = RelInf[i];
        }
//////was the rate of HIV and has been replaced as rate of risk factor of interest////
//    for(int i=0; i<rRFt.nrow(); i++) {
//        for(int j=0; j<rRFt.ncol(); j++) {
//            rRFtN[i][j][0] = rRFt(i,j);
//            rRFtN[i][j][1] = rRFt(i,j)*RFHrPar;
//        } }
for(int i=0; i< TxVec.size(); i++) {
            TxVecN[i] = TxVec[i];
        }
for(int i=0; i<rDxt.nrow(); i++) {
        for(int j=0; j<rDxt.ncol(); j++) {
            rDxtN[i][j] = rDxt(i,j);
        } }
for(int i=0; i<dist.nrow(); i++) {
      for(int j=0; j<dist.ncol(); j++) {
        distN[i][j] = dist(i,j);
      } }
for(int i=0; i<HrEntEx.nrow(); i++) {
        for(int j=0; j<HrEntEx.ncol(); j++) {
            HrEntExN[i][j] = HrEntEx(i,j);
        } }
for(int i=0; i<5; i++) {
            TxVecZ[i] = 0.0;
    }
for(int i=0; i<nYrs; i++) {
    for(int j=0; j<nRes; j++) {
        Outputs[i][j] = 0;
    }  }
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++){
            for(int im=0; im<4; im++){
                for(int nm=0; nm<4; nm++){
                    for(int rg=0; rg<2; rg++) {
                        for(int na=0; na<2; na++){
                            V0[ag][tb][lt][im][nm][rg][na]    = 0;
                            V1[ag][tb][lt][im][nm][rg][na]    = 0;
                            VMort[ag][tb][lt][im][nm][rg][na] = 0;
                            Vdx[ag][tb][lt][im][nm][rg][na]   = 0;
} } } } } } }
for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
        VNkl[i][j] = 0;
} }
for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
        VGjkl[i][j] = 0;
        VLjkl[i][j] = 0;
} }
///effective contact rates
for(int i=0; i<4; i++) {
    Vjaf[i]= 0;
}
for(int i=0; i<5; i++) {
    extrV[i*2] = extrV[i*2+1] = i;
}
for(int ag=0; ag<11; ag++) {
    for(int im=0; im<5; im++) {
        temp4V[ag][im] = (1-pow(1-MrslowN[ag][im]-rRecov,24.0))*MrslowN[ag][im];
} }

////////////////////////////////////////////////////////////////////////////////
///////                  UPDATING TREATMENT METERS                        //////
////////     THIS DIFFERENT TO MAIN MODEL DUE TO SIMPLIFIED OUTCOMES      //////
////////////////////////////////////////////////////////////////////////////////
//////// TREATMENT EFFICACY UPDATED FOR TREATMENT QUALITY //////////////////////
    TxVecZ[1] = TxVecN[1]*TxQualt[0];
///////// RATE OF TREATMENT EXIT TO CURE (LS) //////////////////////////////////
    TxVecZ[2] = TxVecN[0]*TxVecZ[1] + rDeft[0]*TxVecZ[1]*RRcurDef;
//////// RATE OF TREATMENT EXIT TO FAILURE (IN/IP) /////////////////////////////
    TxVecZ[3] = TxVecN[0]*(1-TxVecZ[1]) + rDeft[0]*(1-TxVecZ[1]*RRcurDef);

////////////////////////////////////////////////////////////////////////////////
//////                             StatList                                /////
//////                             BURN IN                                 /////
//////                           Populate model                            /////
////////////////////////////////////////////////////////////////////////////////
///////         SOURCE IN THE CODE TO DISTRIBUTE THE POPULATION            /////
for(int ag=0; ag<11; ag++) {
    for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++) {
              for (int i=0; i<4; i++) {
                for (int j=0; j<4; j++) {
////////////////////        UNINFECTED/SUSCEPTIBLE POP /////////////////////////
      V0[ag][0][0][i][j][0][0] = InitPopN[ag][0]*0.40*distN[i][j]*(1-p_HR); //low risk US born
      V0[ag][0][0][i][j][1][0] = InitPopN[ag][0]*0.40*distN[i][j]*(p_HR); //high risk US born
      V0[ag][0][0][i][j][0][1] = InitPopN[ag][1]*0.40*distN[i][j]*(1-p_HR); //low risk non-US born
      V0[ag][0][0][i][j][1][1] = InitPopN[ag][1]*0.40*distN[i][j]*(p_HR);  //high risk non-US born
/////////////////////////   LATENT SLOW INFECTED POP  //////////////////////////
      V0[ag][2][0][i][j][0][0] = InitPopN[ag][0]*0.60*distN[i][j]*(1-p_HR);
      V0[ag][2][0][i][j][1][0] = InitPopN[ag][0]*0.60*distN[i][j]*(p_HR);
      V0[ag][2][0][i][j][0][1] = InitPopN[ag][1]*0.60*distN[i][j]*(1-p_HR);
      V0[ag][2][0][i][j][1][1] = InitPopN[ag][1]*0.60*distN[i][j]*(p_HR);
  } } } } }
//////create a 2nd array with same dimensions as V0 & populate w/ same values//
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++){
                    for (int na=0; na<3; na++){
                        V1[ag][tb][0][im][nm][rg][na]  = V1[ag][tb][0][im][nm][rg][na];
} } } } } }
////////////////////////RUN THE MODEL FOR 3000 MONTHS /////////////////////////
for(int m=0; m<3001; m++) {
/////////////////////////////////START BURN IN//////////////////////////////////
////////////////////////////////////BIRTHS//////////////////////////////////////
///////////USE DISTRIBUTION TO POPULATE THE MODEL ACROSS RISK GROUPS////////////
for (int i=0; i<4; i++){
    for (int j=0;j<4; j++){
          V0[0][0][0][i][j][0][0]  += Birthst[0]*distN[i][j]*(1-p_HR);
          V0[0][0][0][i][j][1][0]  += Birthst[0]*distN[i][j]*(p_HR);
} }
//////////////////////////////////IMMIGRATION///////////////////////////////////

for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
      for (int i=0; i<4; i++){
        for (int j=0; j<4; j++){
          V1[ag][0][0][i][j][0][1]   += ImmNonN[0][ag]*distN[i][j]*(1-p_HR);  // NO TB, low risk
          V1[ag][0][0][i][j][1][1]   += ImmNonN[0][ag]*distN[i][j]*(p_HR);    // NO TB, high risk

          V1[ag][2][0][i][j][0][1]   += ImmLatN[0][ag]*distN[i][j]*(1-p_HR); // LATENT SLOW TB, low risk
          V1[ag][2][0][i][j][1][1]   += ImmLatN[0][ag]*distN[i][j]*(p_HR);   // LATENT SLOW TB, high risk

          V1[ag][3][0][i][j][0][1]   += (TBImm[ag][0][1])*distN[i][j]*(1-p_HR);   // LATENT FAST, low risk
          V1[ag][3][0][i][j][1][1]   += (TBImm[ag][0][1])*distN[i][j]*(p_HR);   // LATENT FAST, high risk

          V1[ag][4][0][i][j][0][1]   += (TBImm[ag][0][0])*distN[i][j]*(1-p_HR);   //ACTIVE TB, low risk
          V1[ag][4][0][i][j][1][1]   += (TBImm[ag][0][0])*distN[i][j]*(p_HR);   //ACTIVE TB, high risk

} } } }
///////////////////////////////EMMIGRATION//////////////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++){
          for(int rg=0; rg<2; rg++){
            V1[ag][tb][0][im][nm][rg][1]  -= V0[ag][tb][0][im][nm][rg][1]*rEmmigFB[0];   // FB1
            V1[ag][tb][0][im][nm][rg][2]  -= V0[ag][tb][0][im][nm][rg][2]*rEmmigFB[1];   // FB2
          } } } } }
////////////////////////////////MORTALITY//////////////////////////////////////
for(int ag=0; ag<11; ag++) {
  for(int im=0; im<4; im++) {
    for(int nm=0; nm<4; nm++) {
      for(int rg=0; rg<2; rg++) {
        for(int im=0; im<4; im++) {
          for(int na=0; na<3; na++){
            for(int tb=0; tb<5; tb++) {
                 V1[ag][tb][0][im][nm][rg][na]  -= V0[ag][tb][0][im][nm][rg][na]*(mubtN[0][ag]*RRmuHR[rg]+vTMortN[ag][tb]);
            }
////////////////          MORTALITY WITH TB TREATMENT         ////////////////////
                 V1[ag][5 ][0][0][0][rg][na]  -= V0[ag][5 ][0][0][0][rg][na]*(mubtN[0][ag]*RRmuHR[rg]+vTMortN[ag][5 ]*pow(1.0-TxVecZ[1],TunTxMort)); //check the mortality in param
  } } } } } }
/////////////////////////////////////AGING///////////////////////////////////////
for(int ag=0; ag<10; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int rg=0; rg<2; rg++) {
/////          IF AGE > 4, IT TAKES 120 MONTHS TO LEAVE AGE GROUP          /////
                if(ag>0) {
                    temp2 = 120;
/////          IF AGE < 4, IT TAKES 60 MONTHS TO LEAVE AGE GROUP           /////
                } else {
                    temp2 = 60;
                }
            temp = V0[ag][tb][0][0][0][0][rg]/temp2;
                   V1[ag  ][tb][0][0][0][0][rg]  -= temp;
                   V1[ag+1][tb][0][0][0][0][rg]  += temp;
} } }
///////////////////////// NEW FB -> ESTABLISHED FB ///////////////////////////////
///////////////////////// TWO YEARS FOR TRANSITION ///////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
      for(int rg=0; rg<6; rg++) {
        temp = V0[ag][tb][0][0][0][rg][1] / 24;
        V1[ag][tb][0][0][0][rg][1]  -= temp;
        V1[ag][tb][0][0][0][rg][2]  += temp;
      } } }
//////////////////////////// HIGH-RISK ENTRY/EXIT ////////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        temp  = V0[ag][tb][0][0][0][0][0]*HrEntExN[ag][0];
        temp2 = V0[ag][tb][0][0][0][0][1]*HrEntExN[ag][1];
//THESE CODES WERE UPDATED, BUT REMAIN ALMOST THE SAME
        V1[ag][tb][0][0][0][0][0]  += temp2-temp;
        V1[ag][tb][0][0][0][0][1]  += temp-temp2;
} }
////////////////////////////  TRANSMISSION RISK  ////////////////////////////////
for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
        VNkl [i][j] = 0;
        VGjkl[i][j] = 0;   // set to zero
    } }
// Step 1
// DO WE WANT TO LOOP OVER THE DISTRIBUTION OF RISK FACTORS?
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++){
///////// RISK FACTOR FREE & LOW RISK US BORN
            VNkl[0][0]  += V0[ag][tb][lt][0][0][0][0];
///////// RISK FACTOR FREE & HIGH RISK US BORN
            VNkl[0][1]  += V0[ag][tb][lt][0][0][1][0];
///////// RISK FACTOR FREE & LOW RISK NON US BORN
            VNkl[1][0]  += V0[ag][tb][lt][0][0][0][1] + V0[ag][tb][lt][0][0][0][2];
///////// RISK FACTOR FREE & HIGH RISK NON US BORN
            VNkl[1][1]  += V0[ag][tb][lt][0][0][1][1] + V0[ag][tb][lt][0][0][1][2];
} } }
// Step 2  (active TB)
// DO WE WANT TO LOOP OVER THE DISTRIBUTION OF RISK FACTORS?
for(int ag=0; ag<11; ag++) {
/////////  LOW RISK US BORN
    VGjkl[0][0]  +=  V0[ag][4][0][0][0][0][0]                            *RelInfN[4];
///////// HIGH RISK US BORN
    VGjkl[1][0]  +=  V0[ag][4][0][0][0][1][1]                            *RelInfN[4];
///////// LOW RISK NON US BORN
    VGjkl[0][1]  += (V0[ag][4][0][0][0][0][1] + V0[ag][4][0][0][0][0][2])*RelInfN[4];
/////////  HIGH RISK NON US BORN
    VGjkl[1][1]  += (V0[ag][4][0][0][0][1][1] + V0[ag][4][0][0][0][1][2])*RelInfN[4];
}
// Step 2 (treated TB)
// No contribution to force of infection

// Step 3
Vjaf[0]  = (RelInfRg[0]*VGjkl[0][0]         +       //LOW RISK US BORN
            RelInfRg[1]*VGjkl[1][0]*Vmix[0] +       //HIGH RISK US BORN
            RelInfRg[0]*VGjkl[0][1]*Vmix[1])+       //LOW RISK FOREIGN BORN
            RelInfRg[1]*VGjkl[1][1]*Vmix[0]*Vmix[1] //HIGH RISK FOREIGN BORN
        /  (RelInfRg[0]*VNkl[0][0]                    +
            RelInfRg[1]*VNkl[1][0]*Vmix[0]            +
            RelInfRg[0]*VNkl[0][1]*Vmix[1]            +
            RelInfRg[1]*VNkl[0][1]*Vmix[0]*Vmix[1]    +
            1e-12);

Vjaf[1]  = (Vmix[0]*VGjkl[1][1] + Vmix[1]*VGjkl[1][0]) / (Vmix[0]*VNkl[1][1] +VNkl[1][0]+1e-12);

Vjaf[2]  = (Vmix[1]*VGjkl[1][1] + Vmix[0]*VGjkl[0][1]) / (Vmix[1]*VNkl[1][1] + VNkl[0][1]+1e-12);

Vjaf[3]  = VGjkl[1][1] / (VNkl[1][1]+1e-12);

// Step 4
/// LOW RISK US BORN
VLjkl[0 ][0 ]  = RelInfRg[0]*Vjaf[0];
///////// HIGH RISK US BORN
VLjkl[1 ][0 ]  = RelInfRg[1]*((1-Vmix[0])*Vjaf[1]+Vmix[0]*Vjaf[0]);
///////// LOW RISK NON US BORN
VLjkl[0 ][1 ]  = RelInfRg[0]*((1-Vmix[1])*Vjaf[2]+Vmix[1]*Vjaf[0]) + ExogInfN[0];
///////// HIGH RISK NON US BORN
VLjkl[1 ][1 ]  = RelInfRg[1]*((1-Vmix[0])*(1-Vmix[1])*Vjaf[2]+Vmix[0]*Vmix[1]*Vjaf[0])+ExogInfN[0];
///////////////////////////////INFECTION///////////////////////////////////////
///////////////////////for all age groups, risk groups/////////////////////////
///////INFECTION IS CALCULATED WITH THE FORCE OF INFECTION BY RISK GROUP///////
/////// THE TOTAL NUMBER OF INFECTED THEN ENTER BOTH THE LATENT SLOW &  ///////
///////& LATENT FAST TB STATES DEPENDENT ON PROBABILITY OF FAST LATENCY ///////
///////    PEOPLE ARE REMOVED FROM SUSCEPTIBLE, PARTIALLY IMMUNE &      ///////
///////                   LATENT SLOW STATES                            ///////
///////////////////////////////////////////////////////////////////////////////
for(int rg=0; rg<2; rg++) {
  for (int na=0; na<3; na++){
    for(int ag=0; ag<11; ag++) {
///////////////////////////////   SUCEPTIBLE  /////////////////////////////////
    temp = V0[ag][0][0][0][0][rg][na]*VLjkl[rg][na]*EarlyTrend[m];
//////////////////////////// REMOVE FROM SUSCEPTIBLE //////////////////////////
    V1[ag][0][0][0][0][rg][na]  -= temp;
//////////////////////////////// LATENT TB SLOW ///////////////////////////////
    V1[ag][2][0][0][0][rg][na]  += temp*MpslowN[ag][0];
//////////////////////////////// LATENT TB FAST ///////////////////////////////
    V1[ag][3][0][0][0][rg][na]  += temp*MpfastN[ag][0];
///////////////////////////////////////////////////////////////////////////////

/////////////////////////////// SUPER-INFECTION SP ////////////////////////////
    temp = V0[ag][1][0][0][0][rg][na]*VLjkl[rg][na];
    V1[ag][1][0][0][0][rg][na] -= temp;
    V1[ag][2][0][0][0][rg][na] += temp*MpslowPIN[ag][0];
    V1[ag][3][0][0][0][rg][na] += temp*MpfastPIN[ag][0];
///////////////////////////////////////////////////////////////////////////////

/////////////////////////////// SUPER-INFECTION LS ////////////////////////////
    temp = V0[ag][2][0][0][0][rg][na]*VLjkl[rg][na];
    V1[ag][2][0][0][0][rg][na]  -= temp;
    V1[ag][2][0][0][0][rg][na]  += temp*MpslowPIN[ag][0];
    V1[ag][3][0][0][0][rg][na]  += temp*MpfastPIN[ag][0];
    } } }
////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////BREAKDOWN///////////////////////////////////
///////////////////////for all age groups, risk groups/////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int rg=0; rg<2; rg++) {
        for (int na=0; na<3; na++){
        temp  = V0[ag][2][0][0][0][rg][na]*MrslowN[ag][0]*rrSlowFB[rg];
        temp2 = V0[ag][3][0][0][0][rg][na]*rfast;
        V1[ag][2][0][0][0][rg][na]  -= temp;
        V1[ag][3][0][0][0][rg][na]  -= temp2;
        V1[ag][4][0][0][0][rg][na]  += (temp+temp2);
} } }
///////////////////////////     LATENT SLOW TO SAFE /////////////////////////////
///////////////////////for all age groups, risk groups///////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int rg=0; rg<2; rg++) {
      for (int na=0; na<3; na++){
        temp  = V0[ag][2][0][0][0][rg][na]*rRecov;
        V1[ag][2][0][0][0][rg][na]  -= temp;
        V1[ag][1][0][0][0][rg][na]  += temp;
 } } }

////////////////////////////////// SELF CURE/////////////////////////////////////
/////////////////for all age groups, risk groups, only TB 4 /////////////////////
for(int ag=0; ag<11; ag++) {
    for(int rg=0; rg<2; rg++) {
        for (int na=0; na<3; na++){
        temp = V0[ag][4][0][0][0][rg][na]*rSlfCur;
               V1[ag][4][0][0][0][rg][na]  -= temp;
               V1[ag][2][0][0][0][rg][na]  += temp;
} } }
////////////////////TB DIAGNOSIS AND TX INITIATION PUBLIC ///////////////////////
///////////////////////for all age groups, living cond///////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int rg=0; rg<2; rg++) {
      for(int na=0; na<3; na++) {
//       if(rg!=1) {
//            ti = 0;
//     } else {
//            ti = 1;
//     }
////////replace ti with rg?
    temp  = V0[ag][4 ][0][0][0][rg][na]*rDxtN[0][rg]/RRdxAge[ag]/EarlyTrend[m];
            V1[ag][4 ][0][0][0][rg][na]  -= temp;
            V1[ag][5 ][0][0][0][rg][na]  += temp;
      } } }
///////////////////////////   TREATMENT OUTCOMES    /////////////////////////////
///////////////////////for all age groups, risk groups///////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int rg=0; rg<2; rg++) {
        for(int na=0; na<3; na++){
      if(rg!=1) {
         ti2 = 0;
    } else {
         ti2 = 1;
    }
//////////CURES//////////////////////////////////////////////////////////////////
            temp  = V0[ag][5][0][0][0][rg][na]*TxVecZ[2];
            V1[ag][5][0][0][0][rg][na]  -= temp;
            V1[ag][2][0][0][0][rg][na]  += temp;
//////////FAILURES(INCLUDING TREATMENT DEFAULT)//////////////////////////////////
            temp  = V0[ag][5][0][0][0][rg][na]*TxVecZ[3]; ///check this line
                    V1[ag][5][0][0][0][rg][na]  -= temp;
                    V1[ag][4][0][0][0][rg][na]  += temp;
} } }
///////////////////////////////RESET POPULATION SIZE/////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int i=0; i<2; i++) {
        InitPopZ[ag][i] = 0;
} }
for(int ag=0; ag<11; ag++) {
  for(int tb=0; tb<6; tb++) {
////////////////NEED TO UPDATE THESE FOR NEW RISK GROUPS ///////////////////////
    InitPopZ[ag][0]  += V1[ag][tb][0][0][0][0][0];
    InitPopZ[ag][1]  += V1[ag][tb][0][0][0][0][1]+V1[ag][tb][0][0][0][0][2];
} }
for(int ag=0; ag<11; ag++) {
    for(int i=0; i<2; i++) {                        // factor for pop size reset
        InitPopZ[ag][i] = InitPopN[ag][i]/(InitPopZ[ag][i]+1e-12);
} }
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {                    // reset pop to InitPop
////reset under the assumption that InitPopZis [age][nativity]--CHECK ASAP
        V1[ag][tb][0][0][0][0][0]  = V1[ag][tb][0][0][0][0][0]*InitPopZ[ag][0]+
                                     V1[ag][tb][0][0][0][0][0]*InitPopZ[ag][1];
        V1[ag][tb][0][0][0][0][1]  = V1[ag][tb][0][0][0][0][1]*InitPopZ[ag][0]+
                                     V1[ag][tb][0][0][0][0][1]*InitPopZ[ag][1];
for(int rg=0; rg<2; rg++) {
  for (int na=0; na<3; na++){
    V0[ag][tb][0][0][0][rg][na]  = V1[ag][tb][0][0][0][rg][na];
} } } }

///////////////////////////////////////////////////////////////////////////////
} ///////////////////////////END BURN IN///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
NumericVector  CheckV0(14784);
    for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++){
                for(int im=0; im<4; im++){
                    for(int nm=0; nm<4; nm++){
                        for(int rg=0; rg<2; rg++) {
                            for(int na=0; na<2; na++) {
CheckV0(ag+tb*11+lt*66+lt*132+im*528+nm*2112+rg*4928+na*12672) = V1[ag][tb][lt][im][nm][rg][na];
 } } } } } } }
///////////////////////////////////////////////////////////////////////////////
/////////////////////////////BEGIN MODEL///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////////////////////CREATE TIME LOOPS FOR MODEL RUN///////////////////////////
//////////////////////////////////YEAR LOOP////////////////////////////////////
for(int y=0; y<nYrs; y++) {
/////////////////////////////////MONTH LOOP////////////////////////////////////
for(int m=0; m<12; m++) {
/////////////////////CREATE A COUNTER OF MONTHS SINCE START////////////////////
s = y*12+m;
/////////////////////////UPDATING TREATMENT PARAMETERS/////////////////////////
///// TxMatZ: 0=completion rate, 1 = tx success, 2 = RATE OF EXIT TO CURE /////
///// 3 = RATE OF EXIT TO ACTIVE TB, 4 = RATE OF EXIT TO RETREATMENT      /////
///// 5 = PROBABILITY OF TREATMENT COMPLETION                           ///////
///////////////////////////////////////////////////////////////////////////////
//////// TREATMENT EFFICACY UPDATED FOR TREATMENT QUALITY //////////////////////
TxVecZ[1] = TxVecN[1]*TxQualt[s];
///////// RATE OF TREATMENT EXIT TO CURE (LS) //////////////////////////////////
TxVecZ[2] = TxVecN[0]*TxVecZ[1] + rDeft[s]*TxVecZ[1]*RRcurDef;
///////// RATE OF TREATMENT EXIT TO ACTIVE TB //////////////////////////////////
////need to check the pReTx

TxVecZ[3] = TxVecN[0]*(1-TxVecZ[1])*(1-pReTx[s]) + rDeft[s]*(1-TxVecZ[1])*RRcurDef*(1-pReTx[s]);
///////// RATE OF TREATMENT EXIT TO RE TREATMENT////////////////////////////////
TxVecZ[4] = TxVecN[0]*(1-TxVecZ[1])*(pReTx[s]) + rDeft[s]*(1-TxVecZ[1])*RRcurDef*pReTx[s];
////////////////////// P(TREATMENT COMPLETION) /////////////////////////////////
TxVecZ[5] = TxVecN[0]*(1-(1.0-TxVecZ[1])*pReTx[s]);
/////////////////////////////////////BIRTHS//////////////////////////////////////
for (int i=0; i<4; i++){
  for (int j=0; j<4; j++){
/////LOW RISK GROUP BIRTHS//////////////////////////////////////////////////////
    V0[0][0][0][i][j][0][0]  += Birthst[s]*distN[i][j]*(1-p_HR);
/////HIGH RISK GROUP BIRTHS//////////////////////////////////////////////////////
    V0[0][0][0][i][j][1][0]  += Birthst[s]*distN[i][j]*(p_HR); //check this
} }
///////////////////////////////// IMMIGRATION ///////////////////////////////////

for(int ag=0; ag<11; ag++) {
  for (int i=0; i<4; i++){
    for (int j=0;j<4; j++){
////not concerned with previous latent treatment so we can assume all tx naive///
        V1[ag][0][0][i][j][0][1]   += ImmNonN[s][ag]*distN[i][j]*(1-p_HR);  // NO TB, low risk
        V1[ag][0][0][i][j][1][1]   += ImmNonN[s][ag]*distN[i][j]*(p_HR);    // NO TB, high risk

//      V1[ag][1][0][i][j][0][1]   += ImmPIN[s][ag]*dist[i,j]*(1-p_HR);  // previously treated, low risk
//      V1[ag][1][0][i][j][1][1]   += ImmPIN[s][ag]*dist[i,j]*(p_HR);    // previously treated, high risk

        V1[ag][2][0][i][j][0][1]   += ImmLatN[s][ag]*distN[i][j]*(1-p_HR); // LATENT SLOW TB, low risk
        V1[ag][2][0][i][j][1][1]   += ImmLatN[s][ag]*distN[i][j]*(p_HR);   // LATENT SLOW TB, high risk

        V1[ag][3][0][i][j][0][1]   += ImmFstN[s][ag]*distN[i][j]*(1-p_HR);   // LATENT FAST, low risk
        V1[ag][3][0][i][j][1][1]   += ImmFstN[s][ag]*distN[i][j]*(p_HR);   // LATENT FAST, high risk

        V1[ag][4][0][i][j][0][1]   += ImmActN[s][ag]*distN[i][j]*(1-p_HR);   //ACTIVE TB, low risk
        V1[ag][4][0][i][j][1][1]   += ImmActN[s][ag]*distN[i][j]*(p_HR);   //ACTIVE TB, high risk

      } } }

/////////////////////////////////  EMMIGRATION ///////////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++){
            for(int im=0; im<4; im++){
                for(int nm=0; nm<4; nm++){
                    for(int rg=0; rg<2; rg++) {
                        V1[ag][tb][lt][im][nm][rg][1]  -= V0[ag][tb][lt][im][nm][rg][1]*rEmmigFB[0];      // FB1
                        V1[ag][tb][lt][im][nm][rg][2]  -= V0[ag][tb][lt][im][nm][rg][2]*rEmmigFB[1];      // FB2
} } } } } }
/////////////////////////////////  MORTALITY  ///////////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2; lt++){
        for(int nm=0; nm<4 ; nm++) {
            for(int im=0; im<4 ; im++) {
                if(im==3) {
                    temp = muTbRF;
                } else { temp = 0;
                }
              for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
////////////////////////UNINFECTED, SUSCEPTIBLE//////////////////////////////////
    VMort[ag][0 ][lt][im][nm][rg][na]  = V0[ag][0 ][lt][im][nm][rg][na]*
                                 (mubtN[s][ag]*RRmuHR[rg]+vRFMortN[ag][nm] );
////////////////////////UNINFECTED, PART. IMMUNE/////////////////////////////////
    VMort[ag][1 ][lt][im][nm][rg][na]  = V0[ag][1 ][lt][im][nm][rg][na]*
                                  (mubtN[s][ag]*RRmuHR[rg]+vRFMortN[ag][nm] );
////////////////////////    LATENT TB SLOW      /////////////////////////////////
    VMort[ag][2 ][lt][im][nm][rg][na]  = V0[ag][2 ][lt][im][nm][rg][na]*
                                (mubtN[s][ag]*RRmuHR[rg]+vRFMortN[ag][nm] );
////////////////////////    LATENT TB FAST      /////////////////////////////////
    VMort[ag][3 ][lt][im][nm][rg][na]  = V0[ag][3 ][lt][im][nm][rg][na]*
                      (mubtN[s][ag]*RRmuHR[rg]+vRFMortN[ag][nm] );
////////////////////////      ACTIVE TB         /////////////////////////////////
    VMort[ag][4 ][lt][im][nm][rg][na]  = V0[ag][4 ][lt][im][nm][rg][na]*
              (mubtN[s][ag]*RRmuHR[rg]+vRFMortN[ag][nm]+vTMortN[ag][4 ]+temp );
////////////////////////    TB TREATMENT        /////////////////////////////////
    VMort[ag][5 ][lt][im][nm][rg][na]  = V0[ag][5 ][lt][im][nm][rg][na]*
    (mubtN[s][ag]*RRmuHR[rg]+vRFMortN[ag][nm]+(vTMortN[ag][5 ]+temp)*pow(1.0-TxVecZ[1],TunTxMort));
/////////////// UPDATE THE PRIMARY VECTOR BY REMOVING MORTALITY /////////////////
          for(int tb=0; tb<6; tb++) {
             V1[ag][tb][lt][im][nm][rg][na]  -= VMort[ag][tb][lt][im][nm][rg][na];
          }
  } } } } } }
/////////////////////////////////////AGING///////////////////////////////////////
for(int ag=0; ag<1; ag++) {
/////          IF AGE > 4, IT TAKES 120 MONTHS TO LEAVE AGE GROUP          /////
    if(ag>0) {
            temp2 = 120;
/////          IF AGE < 4, IT TAKES 60 MONTHS TO LEAVE AGE GROUP           /////
    } else {
            temp2 = 60;
    }
for(int tb=0; tb<6; tb++) {
    for(int lt=0; lt<2; lt++){
        for (int im=0; im<4; im++){
            for (int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++){
                            temp = V0[ag  ][tb][lt][im][nm][rg][na]/temp2;
                                   V1[ag  ][tb][lt][im][nm][rg][na]  -= temp;
                                   V1[ag+1][tb][lt][im][nm][rg][na]  += temp;
} } } } } } }
///////////////////////// NEW FB -> ESTABLISHED FB ///////////////////////////////
///////////////////////// TWO YEARS FOR TRANSITION ///////////////////////////////
for(int ag=0; ag<11; ag++){
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++){
            for (int im=0; im<4; im++){
                for (int nm=0; nm<4; nm++){
                    for(int rg=0; rg<2; rg++) {
                        temp = V0[ag][tb][lt][im][nm][rg][1]/24;
                        V1[ag][tb][lt][im][nm][rg][1]  -= temp;
                        V1[ag][tb][lt][im][nm][rg][2]  += temp;
                    } } } } } }
//////////////////////////// HIGH-RISK ENTRY/EXIT ////////////////////////////////
///////////////////REVIEW AND UPDATE PROPERLY; JUST EDITED INDEXES////////////////
//////////////////////////// HIGH-RISK ENTRY/EXIT ////////////////////////////////
for(int ag=0; ag<11; ag++) {
  for(int tb=0; tb<6; tb++) {
    temp  = V0[ag][tb][0][0][0][0][0]*HrEntExN[ag][0];
    temp2 = V0[ag][tb][0][0][0][0][1]*HrEntExN[ag][1];
    //THESE CODES WERE UPDATED, BUT REMAIN ALMOST THE SAME
    V1[ag][tb][0][0][0][0][0]  += temp2-temp;
    V1[ag][tb][0][0][0][0][1]  += temp-temp2;
  } }
////////////////////////////  TRANSMISSION RISK  ////////////////////////////////
for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
        VNkl[i][j] = 0; // set to zero
        VGjkl[i][j] = 0; // set to zero
} }
////////////////////////////          Step 1         ////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++){
                for(int nm=0; nm<4; nm++){
////////////////////////////    LOW RISK, US BORN    ////////////////////////////
                    VNkl[0][0]  += V0[ag][tb][lt][im][nm][0][0];
////////////////////////////   HIGH RISK, US BORN    ////////////////////////////
                    VNkl[1][0]  += V0[ag][tb][lt][im][nm][1][0];
////////////////////////////  LOW RISK, NON US BORN  ////////////////////////////
                    VNkl[2][0]  += V0[ag][tb][lt][im][nm][0][1] + V0[ag][tb][lt][im][nm][0][2];
//////////////////////////// HIGH RISK, NON US BORN  ////////////////////////////
                    VNkl[3][0]  += V0[ag][tb][lt][im][nm][1][1] + V0[ag][tb][lt][im][nm][1][2];
} } } } }
/////////////////////////// /   Step 2  (ACTIVE TB)   ////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2 ; lt++) {
        for(int im=0; im<4 ; im++) {
            for(int nm=0; nm<4; nm++){
////////////////////////////    LOW RISK, US BORN    ////////////////////////////
                VGjkl[0][0]  +=  V0[ag][4][lt][im][nm][0][0]                               *RelInfN[4];
////////////////////////////   HIGH RISK, US BORN    ////////////////////////////
                VGjkl[1][0]  +=  V0[ag][4][lt][im][nm][1][0]                               *RelInfN[4];
////////////////////////////  LOW RISK, NON US BORN  ////////////////////////////
                VGjkl[0][1]  += (V0[ag][4][lt][im][nm][0][1] + V0[ag][4][lt][im][nm][0][2])*RelInfN[4];
//////////////////////////// HIGH RISK, NON US BORN  ////////////////////////////
                VGjkl[1][1]  += (V0[ag][4][lt][im][nm][1][1] + V0[ag][4][lt][im][nm][1][2])*RelInfN[4];
} } } }
////////////////////////////   Step 2 (TREATED TB)   ////////////////////////////
////////////////////  No contribution to force of infection  ////////////////////

////////////////////////////          Step 3         ////////////////////////////
Vjaf[0]  = (RelInfRg[0]*VGjkl[0][0]         +       //LOW RISK US BORN
            RelInfRg[1]*VGjkl[1][0]*Vmix[0] +       //HIGH RISK US BORN
            RelInfRg[0]*VGjkl[0][1]*Vmix[1])+       //LOW RISK FOREIGN BORN
            RelInfRg[1]*VGjkl[1][1]*Vmix[0]*Vmix[1] //HIGH RISK FOREIGN BORN
        /  (RelInfRg[0]*VNkl[0][0]                    +
            RelInfRg[1]*VNkl[1][0]*Vmix[0]            +
            RelInfRg[0]*VNkl[0][1]*Vmix[1]            +
            RelInfRg[1]*VNkl[0][1]*Vmix[0]*Vmix[1]    +
            1e-12);

Vjaf[1]  = (Vmix[0]*VGjkl[1][1] + Vmix[1]*VGjkl[1][0]) / (Vmix[0]*VNkl[1][1] +VNkl[1][0]+1e-12);

Vjaf[2]  = (Vmix[1]*VGjkl[1][1] + Vmix[0]*VGjkl[0][1]) / (Vmix[1]*VNkl[1][1] + VNkl[0][1]+1e-12);

Vjaf[3]  = VGjkl[1][1] / (VNkl[1][1]+1e-12);

// Step 4
//CREATE LAMBDA FORCE OF INFECTION
/// LOW RISK US BORN
VLjkl[0 ][0 ]  = RelInfRg[0]*Vjaf[0];
///////// HIGH RISK US BORN
VLjkl[1 ][0 ]  = RelInfRg[1]*((1-Vmix[0])*Vjaf[1]+Vmix[0]*Vjaf[0]);
///////// LOW RISK NON US BORN
VLjkl[0 ][1 ]  = RelInfRg[0]*((1-Vmix[1])*Vjaf[2]+Vmix[1]*Vjaf[0])+ExogInfN[0];
///////// HIGH RISK NON US BORN
VLjkl[1 ][1 ]  = RelInfRg[1]*((1-Vmix[0])*(1-Vmix[1])*Vjaf[2]+Vmix[0]*Vmix[1]*Vjaf[0])+ExogInfN[0];

///////////////////////////////    INFECTION   /////////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2 ; lt++) {
        for(int im=0; im<4 ; im++) {
            for(int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++){
                    for(int na=0; na<3; na++){
///////////////////////////////   SUCEPTIBLE  /////////////////////////////////
  temp = V0[ag][0][lt][im][nm][rg][na]*VLjkl[rg][na]*NixTrans[s];  // Su
         V1[ag][0][lt][im][nm][rg][na]  -= temp;
         V1[ag][2][lt][im][nm][rg][na]  += temp*MpslowN[ag][im];
         V1[ag][3][lt][im][nm][rg][na]  += temp*MpfastN[ag][im];
///////////////////////////////   SUCEPTIBLE, PI  /////////////////////////////////
  temp = V0[ag][1][lt][im][nm][rg][na]*VLjkl[rg][na]*NixTrans[s];  // Sp
         V1[ag][1][lt][im][nm][rg][na]  -= temp;
         V1[ag][2][lt][im][nm][rg][na]  += temp*MpslowPIN[ag][im];
         V1[ag][3][lt][im][nm][rg][na]  += temp*MpfastPIN[ag][im];
/////////////////SUPER INFECTION LATENT TB SLOW ///////////////////////////////
  temp = V0[ag][2][lt][im][nm][rg][na]*VLjkl[rg][na]*NixTrans[s];  // Ls
         V1[ag][2][lt][im][nm][rg][na] -= temp;
         V1[ag][2][lt][im][nm][rg][na] += temp*MpslowPIN[ag][im]/2;
         V1[ag][2][lt][im][nm][rg][na] += temp*MpslowPIN[ag][im]/2;
         V1[ag][3][lt][im][nm][rg][na] += temp*MpfastPIN[ag][im];
} } } } } }
///////////////////////////////////BREAKDOWN///////////////////////////////////
///////////////////////for all age groups, risk groups/////////////////////////
for(int ag=0; ag<11; ag++) {
  for(int lt=0; lt<2 ; lt++) {
    for(int im=0; im<4 ; im++) {
      for(int nm=0; nm<4; nm++) {
        for(int rg=0; rg<2; rg++) {
          for(int na=0; na<3; na++) {
            temp  = V0[ag][2][lt][im][nm][rg][na]*MrslowN[ag][im]*rrSlowFB[na];  // Latent Slow
            temp2 = V0[ag][3][lt][im][nm][rg][na]*rfast;  // Latent Fast

            V1[ag][2][lt][im][nm][rg][na]  -= temp; //REMOVE FROM LATENT SLOW
            V1[ag][3][lt][im][nm][rg][na]  -= temp2; //REMOVE FROM LATENT FAST
            V1[ag][4][lt][im][nm][rg][na]  += (temp+temp2); //PLACE IN ACTIVE DISEASE
///ASK ABOUT THIS CHUNK///     // Tl progression if INH resistant
///if they were in the latent treatment chunk, but the treatment failed bc INH resistant
////progress to active disease
////NEED TO BRAINSTORM HOW THIS WILL WORK;
 //           temp =  V0[ag][2][dr][tx][im][rg][na]*MrslowN[ag][im]*rrSlowFB[rg]*(1-EffLtXN[s][dr]);
   //         temp2 = V0[ag][3][dr][tx][im][rg][na]*MrslowN[ag][im]*rrSlowFB[rg]*(1-EffLtXN[s][dr]);
     //               V1[ag][6][dr][tx][im][rg][na]  -= temp;
       //             V1[ag][4][dr][tx][im][rg][na]  += (temp+temp2);

} } } } } }
///////////////////////////   LATENT SLOW TO SAFE   /////////////////////////////
for(int ag=0; ag<11; ag++) {
  for(int lt=0; lt<2 ; lt++) {
    for(int im=0; im<4 ; im++) {
      for(int nm=0; nm<4; nm++) {
        for(int rg=0; rg<2; rg++) {
          for(int na=0; na<3; na++) {
            temp  = V0[ag][2][lt][im][nm][rg][na]*rRecov;  // Ls
                    V1[ag][2][lt][im][nm][rg][na]  -= temp;
                    V1[ag][1][lt][im][nm][rg][na]  += temp;
} } } } } }
////////////////////////////////// SELF CURE/////////////////////////////////////
for(int ag=0; ag<11; ag++) {
  for(int lt=0; lt<2 ; lt++) {
    for(int im=0; im<4 ; im++) {
       for(int nm=0; nm<4; nm++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++) {
              temp  = V0[ag][4 ][lt][im][nm][rg][na]*rSlfCur;
              V1[ag][4 ][lt][im][nm][rg][na]  -= temp;
              V1[ag][2 ][lt][im][nm][rg][na]  += temp;
} } } } } }
/// LTBI SCREENING AND TLTBI INITIATION /// only for no previous TB or LTBI tx
for(int nm=0; nm<4; nm++) {
for(int im=0; im<4; im++) {
  for(int rg=0; rg<2; rg++) {
    for(int na=0; na<3; na++) {
////////////// US BORN, LOW RISK  //////////////////
      if(im==0 & rg==0 & na==0) {
        rTbP = rLtScrt[s]*LtDxParN[0][0];
        rTbN = rLtScrt[s]*LtDxParN[0][1];
      }
////////////// US BORN, HIGH RISK  /////////////////
      if(im==0 & rg==1 & na==0 ) {
        rTbP = rLtScrt[s]*LtDxParN[1][0];
        rTbN = rLtScrt[s]*LtDxParN[1][1];
      }
//////////// NON US BORN  ////////////////
      if(im==0 & na>0) {
        rTbP = rLtScrt[s]*LtDxParN[2][0];
        rTbN = rLtScrt[s]*LtDxParN[2][1];
      }
for(int ag=0; ag<11; ag++) {
// Have LTBI
//    temp  = V0[ag][2][0][im][nm][rg][na]*rTbP;
//    temp2 = V0[ag][3][0][im][nm][rg][na]*rTbP;
//    V1[ag][2][0][im][nm][rg][na]  -= temp;
//    V1[ag][3][0][im][nm][rg][na]  -= temp2;
//    V1[ag][6][0][im][nm][rg][na]  += temp+temp2;
// Dont have LTBI
    temp  = V0[ag][0][0][im][nm][rg][na]*rTbN;
    temp2 = V0[ag][1][0][im][nm][rg][na]*rTbN;
    V1[ag][0][0][im][nm][rg][na]  -= temp;
    V1[ag][1][0][im][nm][rg][na]  -= temp2;
///////moving to latent tx experienced as in last model -- is this correct?
    V1[ag][0][1][im][nm][rg][na]  += temp;
    V1[ag][1][1][im][nm][rg][na]  += temp2;
} } } } }
/// TLTBI: TX COMPLETION + DEFAULT /// only need to consider tx naive compartment
for(int ag=0; ag<11; ag++) {
  for(int im=0; im<4 ; im++) {
    for(int nm=0; nm<4; nm++) {
      for(int rg=0; rg<2; rg++) {
        for(int na=0; na<3; na++) {
          temp  = V0[ag][2][0][im][nm][rg][na]*rTbP*dLtt[s]; // tx completion
          temp2 = V0[ag][3][0][im][nm][rg][na]*rTbP*dLtt[s]; // tx completion
          temp3 = V0[ag][2][0][im][nm][rg][na]*rTbP*LtTxPar[1]; // default
          temp4 = V0[ag][3][0][im][nm][rg][na]*rTbP*LtTxPar[1]; // default
                  V1[ag][2][0][im][nm][rg][na]  -= temp+temp3; //remove from latent slow
                  V1[ag][3][0][im][nm][rg][na]  -= temp2+temp4;  //remove from latent fast
                  V1[ag][1][1][im][nm][rg][na]  += (temp+temp2)*EffLtXN[s]*EffLt; //exit to cure
                  V1[ag][2][1][im][nm][rg][na]  += (temp+temp2)*(1-EffLtXN[s]*EffLt);  //tx comp fail to latent slow
                  V1[ag][2][0][im][nm][rg][na]  += (temp3+temp4); //latent tx default to latent slow
} } } } }
///////////////////// TB DIAGNOSIS AND TX INITIATION  /////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2; lt++) {
        for(int im=0; im<4 ; im++) {
            for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                        if(rg!=1) {
                            ti = 0;
                        } else { ti = 1;
                        }
                        temp  = V0[ag][4][lt][im][nm][rg][na]*rDxtN[s][ti  ]/RRdxAge[ag];
                        V1[ag][4][lt][im][nm][rg][na]      -= temp;
                        V1[ag][5][lt][im][nm][rg][na]      += temp;
                        Vdx[ag][4][lt][im][nm][rg][na]      = temp;
} } } } } }
//////////////////// RISK FACTOR OF INTEREST INCIDENCE //////////////////////
//for(int ag=0; ag<11; ag++) {
//    for(int tb=0; tb<6; tb++) {
//        for(int lt=0; lt<2; lt++) {
//            for(int rg=0; rg<2; rg++) {
//                for(int na=0; na<3; na++) {
//                    if(rg!=1) {
//                        temp = V0[ag][tb][lt][0][0][rg][na]*rRFtN[s][ag][0];
//                    } else {
//                        temp = V0[ag][tb][lt][0][0][rg][na]*rRFtN[s][ag][1];
//                    }
//                    V1[ag][tb][lt][0][0][rg][na]  -= temp;
//                    V1[ag][tb][lt][1][1][rg][na]  += temp;
//} } } } }
//////////////////////// TB TREATMENT OUTCOMES /////////////////////////////
for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++) {
            for(int nm=0; nm<4; nm++) {
              for(int im=0; im<4; im++) {
               for(int rg=0; rg<2; rg++) {
                   for(int na=0; na<3; na++) {
                    // Cures back to Ls state
                    temp=V0[ag][5][lt][im][nm][rg][na]*TxVecZ[2];
                    V0[ag][5][lt][im][nm][rg][na]  -= temp;
                    V0[ag][2][lt][im][nm][rg][na]  += temp;
                    ///// EXIT TO ACTIVE DISEASE //////
                    temp=V0[ag][5][lt][im][nm][rg][na]*TxVecZ[3];
                    V0[ag][5][lt][im][nm][rg][na]  -= temp;
                    V0[ag][4][lt][im][nm][rg][na]  += temp;
                    ///// EXIT TO TB RETREATMENT //////
                    temp=V0[ag][5][lt][im][nm][rg][na]*TxVecZ[4];
                    V0[ag][5][lt][im][nm][rg][na]  -= temp;
                    V0[ag][5][lt][im][nm][rg][na]  += temp;
                   } } } } } }
///////////////////////////////////////////////////////////////////////////////
/////////////////////////    FILL RESULTS TABLE    ////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////////////////////////     MID-YEAR RESULTS     ////////////////////////////
if(m==6) {
/////////////////////////////////// YEAR //////////////////////////////////////
Outputs[y][0]      = y+1950;  // Year
////////////////    COUNTS BY TOTAL, AGE, TB, RF, AND RG     //////////////////
for(int ag=0; ag<11; ag++) {
  for(int tb=0; tb<6; tb++) {
    for(int lt=0; lt<2; lt++) {
        for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                        Outputs[y][1    ] += V1[ag][tb][lt][im][nm][rg][na];   // N_ALL
                        Outputs[y][2 +ag] += V1[ag][tb][lt][im][nm][rg][na];   // N_ by age (11)
                        Outputs[y][13+tb] += V1[ag][tb][lt][im][nm][rg][na];   // N_ by tb (6)
                        Outputs[y][19+im] += V1[ag][tb][lt][im][nm][rg][na];   // N_ by im (4)
                        Outputs[y][23+nm] += V1[ag][tb][lt][im][nm][rg][na];   // N_ by nm (4)
                        Outputs[y][27+nm] += V1[ag][tb][lt][im][nm][rg][na];   // N_ by rg (2)
                        Outputs[y][29+nm] += V1[ag][tb][lt][im][nm][rg][na];   // N_ by na (3)
} } } } } } }
////////////////////    COUNTS BY NATIVITY AND AGE    ////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                    for(int rg=0; rg<2; rg++) {
                        Outputs[y][33+ag] += V1[ag][tb][lt][im][nm][rg][0]; // N_ by age and US (11)
                        Outputs[y][44+ag] += V1[ag][tb][lt][im][nm][rg][1]+V1[ag][tb][lt][im][nm][rg][2];   // N_ by age and FB (11)
} } } } } }
/////////////////    COUNTS BY NATIVITY, AGE, & LTBI    //////////////////////
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2; lt++) {
        for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                    Outputs[y][55+ag] += (V1[ag][1][lt][im][nm][rg][0])*(1-pImmScen)+
                                          V1[ag][2][lt][im][nm][rg][0]+V1[ag][3][lt][im][nm][rg][0];   // N_ by age and US (11) LATENT INFECTION

                    Outputs[y][66+ag] += (V1[ag][1][lt][im][nm][rg][1]+V1[ag][1][lt][im][nm][rg][2])*(1-pImmScen)+
                        V1[ag][2][lt][im][nm][rg][1]+V1[ag][2][lt][im][nm][rg][2]+
                        V1[ag][3][lt][im][nm][rg][1]+V1[ag][3][lt][im][nm][rg][2]; // N_ by age and FB (11) LATENT INFECTION
} } } } }
/////////////////////RISK FACTOR OF INTEREST COUNT BY AGE/////////////////////
///////will need to be updated; if im>1 | nm>1 then i=1, else i=0 ////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                    for(int rg=0; rg<2; rg++) {
                        for(int na=0; na<3; na++) {
                            Outputs[y][77+ag] += V1[ag][tb][lt][im][nm][rg][na];   // N_RF by age (11)
} } } } } } }
///////////// TB MORTALITY COUNT BY AGE, RISK FACTOR OF INTEREST///////////////
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2; lt++) {
        for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                        if(im>0) { ti = 11; } else { ti = 0; }
                        if(im==3) { temp = muTbRF; } else { temp = 0; }
                        Outputs[y][88+ag+ti]  += V0[ag][4 ][lt][im][nm][rg][na]*(vTMortN[ag][4 ]+temp);
                        Outputs[y][88+ag+ti]  += V0[ag][5 ][lt][im][nm][rg][na]*(vTMortN[ag][5 ]+temp)*pow(1.0-TxVecZ[1],TunTxMort); ///check txmatz
} } } } } }
////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
for(int i=88; i<110; i++) { Outputs[y][i] = Outputs[y][i]*12; }

// rf MORTALITY BY AGE
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                    for(int rg=0; rg<2; rg++) {
                        for(int na=0; na<3; na++) {
                            Outputs[y][110+ag]  += V0[ag][tb][lt][im][nm][rg][na]*vRFMortN[ag][im];
                        } } } } } } }
////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
for(int i=110; i<121; i++) { Outputs[y][i] = Outputs[y][i]*12;  }
///////////////////////    TOTAL MORTALITY BY AGE    /////////////////////////
for(int ag=0; ag<11; ag++) {
for(int tb=0; tb<6; tb++) {
for(int lt=0; lt<2; lt++) {
for(int im=0; im<4; im++) {
for(int nm=0; nm<4; nm++) {
for(int rg=0; rg<2; rg++) {
for(int na=0; na<3; na++) {
Outputs[y][121+ag]  += VMort[ag][tb][lt][im][nm][rg][na];
} } } } } } }
////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
for(int i=121; i<132; i++) { Outputs[y][i] = Outputs[y][i]*12; }
///////////////////////     TB TREATMENT OUTCOMES    /////////////////////////
for(int ag=0; ag<11; ag++) {
  for(int lt=0; lt<2; lt++) {
    for(int im=0; im<4; im++) {
      for(int nm=0; nm<4; nm++) {
        for(int rg=0; rg<2; rg++) {
          for(int na=0; na<3; na++) {
//if(rg!=1) { temp = rDeft[s]; } else { temp = rDeftH[s]; }
//////updated tx mat z 22 to tx mat z 5
Outputs[y][132]  += V0[ag][5 ][lt][im][nm][rg][na]*TxVecZ[5]; // tx completion  check this!!!

Outputs[y][133]  += V0[ag][5 ][lt][im][nm][rg][na]*rDeft[s]; // tx discontinuation
Outputs[y][134]  += VMort[ag][5 ][lt][im][nm][rg][na];// tx mort
} } } } } }
for(int i=132; i<135; i++) { Outputs[y][i] = Outputs[y][i]*12; }

// NOTIFICATIONS
        for(int ag=0; ag<11; ag++) {
            for(int lt=0; lt<2; lt++) {
                for(int im=0; im<4; im++) {
                    for(int nm=0; nm<4; nm++) {
                        for(int rg=0; rg<2; rg++) {
                            for(int na=0; na<3; na++) {
                                Outputs[y][135   ] += Vdx[ag][4 ][lt][im][nm][rg][na];   // All dx (1)
                                Outputs[y][136+ag] += Vdx[ag][4 ][lt][im][nm][rg][na];   // dx by age (11)
                                        if(lt<0) {
                                Outputs[y][147   ] += Vdx[ag][4 ][lt][im][nm][rg][na];   // dx by tx history - naive (2)
                                        } else {
                                Outputs[y][148   ] += Vdx[ag][4 ][lt][im][nm][rg][na]; } // dx by tx history - experienced (2)
                                        if(im>0) {
                                Outputs[y][149   ] += Vdx[ag][4 ][lt][im][nm][rg][na];   // dx HIV pos (1)
                                      } else {
                                Outputs[y][150   ] += Vdx[ag][4 ][lt][im][nm][rg][na]; }  // dx HIV neg (1)
                                Outputs[y][151+rg] += Vdx[ag][4 ][lt][im][nm][rg][na];   // N_ by rg (2)
                            } } } } } }
                    for(int i=135; i<153; i++) { Outputs[y][i] = Outputs[y][i]*12; }
/// TLTBI INITS ///
for(int im=0; im<4; im++) {
  for(int nm=0; nm<4; nm++) {
    for(int rg=0; rg<2; rg++) {
      for(int na=0; na<3; na++) {
    if(im==0 & rg==0 ) { rTbP = rLtScrt[s]*LtDxParN[0][0]; rTbN = rLtScrt[s]*LtDxParN[0][1]; }
    if(im==0 & rg==1 ) { rTbP = rLtScrt[s]*LtDxParN[1][0]; rTbN = rLtScrt[s]*LtDxParN[1][1]; }
    if(im>0  & rg==0 ) { rTbP = rLtScrt[s]*LtDxParN[3][0]; rTbN = rLtScrt[s]*LtDxParN[3][1]; }
    if(im>0  & rg==1 ) { rTbP = rLtScrt[s]*LtDxParN[4][0]; rTbN = rLtScrt[s]*LtDxParN[4][1]; }
    for(int ag=0; ag<11; ag++) {
        Outputs[y][175] += (V0[ag][3 ][0 ][im][nm][rg][na]+V0[ag][2 ][0 ][im][nm][rg][na])*rTbP +
                           (V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*rTbN; //all init
                                                                if(na>0) {
                                                                   Outputs[y][176] += (V0[ag][3 ][0 ][im][nm][rg][na]+V0[ag][2 ][0 ][im][nm][rg][na])*rTbP +
                                                                    (V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*rTbN; } // FB inits
                                                                if(rg==1) {
                                                                    Outputs[y][177] +=  (V0[ag][3 ][0 ][im][nm][rg][na]+V0[ag][2 ][0 ][im][nm][rg][na])*rTbP +
                                                                    (V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*rTbN; } // high risk inits
                                                                if(im>0) {
                                                                    Outputs[y][178] += (V0[ag][3 ][0 ][im][nm][rg][na]+V0[ag][2 ][0 ][im][nm][rg][na])*rTbP +
                                                                    (V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*rTbN; } // RF inits

                                                                    Outputs[y][179] += (V0[ag][3 ][0 ][im][nm][rg][na]+V0[ag][2 ][0 ][im][nm][rg][na])*rTbP; // inits with LTBI
    } } } } }
for(int i=175; i<180; i++) { Outputs[y][i] = Outputs[y][i]*12; } // annualize

/// TB INCIDENCE, BY ALL VS RECENT  ///
// By recency (<2 years) == all immediate, 1-(1-rfast)^24 x all Lf


// BREAKDOWN from Ls, Lf
temp3 = (1-pow(1-rfast,24.0))*rfast;

for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2; lt++) {
        for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                        temp2 = V0[ag][2 ][lt][im][nm][rg][na]*temp4V[ag][im] + V0[ag][2 ][lt][im][nm][rg][na]*temp3;// Progression from recent infection
                        temp = V0[ag][2 ][lt][im][nm][rg][na]*MrslowN[ag][im] + V0[ag][3 ][lt][im][nm][rg][na]*rfast; // All progression
                        Outputs[y][180      ] += temp ;   // all incidence
                        Outputs[y][180+17   ] += temp2;   // all incidence, recent infection
                        Outputs[y][181+ag   ] += temp ;   // incidence by age
                        Outputs[y][181+ag+17] += temp2;   // incidence by age, recent infection
                        if(na<1) {
                            Outputs[y][192      ] += temp ;   //  incidence, US born
                            Outputs[y][192+17   ] += temp2;   //  incidence, US born, recent infection
                    } else {
                        Outputs[y][193      ] += temp ;   //  incidence, FB
                        Outputs[y][193+17   ] += temp2;  //  incidence, FB, recent infection
                      if(na==2) {
                          Outputs[y][194      ] += temp ;   //  incidence, FB2 born
                          Outputs[y][194+17   ] += temp2; } } //  incidence, FB2 born, recent infection
                      if(rg==1) {
                         Outputs[y][195      ] += temp ;   //  incidence, HR
                         Outputs[y][195+17   ] += temp2; } //  incidence, HR, recent infection
                      if(im>0) {
                            Outputs[y][196      ] += temp ;   //  incidence, HIV pos
                            Outputs[y][196+17   ] += temp2; } //  incidence, HIV pos, recent infection

                    } } } } } }

for(int i=180; i<214; i++) { Outputs[y][i] = Outputs[y][i]*12; } // annualize

// NOTIFICATIONS, dead at diagnosis
for(int nm=0; nm<4 ; nm++) {
    if(nm > 2) {
        temp = muTbRF;
    } else {
        temp = 0;
    }
    for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++){
            for (int im=0; im<4; im++){
                for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++){
                        temp2 = V0[ag][4 ][lt][im][nm][rg][na]*(vTMortN[ag][4]+temp);
                        Outputs[y][214   ] += temp2;   // All dx (1)
                        Outputs[y][215+ag] += temp2;   // dx by age (11)
                        if(lt<1) {
                            Outputs[y][226   ] += temp2;   // dx by tx history - naive (2)
                        } else {
                            Outputs[y][227   ] += temp2; } // dx by tx history - experienced (2)
                        if(im>0) {
                            Outputs[y][228   ] += temp2;   // dx HIV pos (1)
                        } else {
                            Outputs[y][229   ] += temp2; }  // dx HIV neg (1)
                            Outputs[y][230+rg] += temp2;   // N_ by rg (2)
} } } } } }
for(int i=214; i<232; i++) { Outputs[y][i] = Outputs[y][i]*12; }

// NOTIFICATIONS US
    for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++){
            for (int im=0; im<4; im++){
                 for (int nm=0; nm<4; nm++){
                     for(int rg=0; rg<2; rg++) {
                         for(int na=0; na<3; na++){
                            Outputs[y][232+ag] += Vdx[ag][4 ][lt][im][nm][rg][na];   // dx by age (11)
} } } } } }
   for(int i=232; i<243; i++) { Outputs[y][i] = Outputs[y][i]*12; }
// NOTIFICATIONS US, dead at diagnosis
for(int nm=0; nm<4 ; nm++) {
  for (int im=0; im<4; im++){
if(im > 2) { temp = muTbRF;  } else { temp = 0;  }
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2; lt++){
        for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++){
                temp2 = V0[ag][4 ][lt][im][nm][rg][na]*(vTMortN[ag][4]+temp); //vTMort[ag][tb]
                Outputs[y][243+ag] += temp2;   // dx by age (11)
} } } } } }
for(int i=243; i<254; i++) { Outputs[y][i] = Outputs[y][i]*12; }

// NOTIFICATIONS, dead at diagnosis  HIV_NEGATIVE
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2; lt++){
        for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++){
                temp2 = V0[ag][4 ][lt][0][0][rg][na]*(vTMortN[ag][4]+temp);
                Outputs[y][254+ag] += temp2;   // dx by age (11)
} } } }
////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
for(int i=254; i<265; i++) { Outputs[y][i] = Outputs[y][i]*12; }
// TOTAL MORTALITY BY AGE, HAVE HIV
    for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++){
                for (int im=0; im<4; im++){
                    for (int nm=0; nm<4; nm++){
                        for(int rg=0; rg<2; rg++) {
                            for(int na=0; na<3; na++){
                                    Outputs[y][267+ag]  += VMort[ag][tb][lt][im][nm][rg][na];
                            } } } } } } }
////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
for(int i=267; i<278; i++) { Outputs[y][i] = Outputs[y][i]*12; }

/////////////////////  TOTAL MORTALITY BY AGE, HAVE TB   /////////////////////
    for(int ag=0; ag<11; ag++) {
            for(int lt=0; lt<2; lt++){
                for (int im=0; im<4; im++){
                    for (int nm=0; nm<4; nm++){
                        for(int rg=0; rg<2; rg++) {
                            for(int na=0; na<3; na++){
                                Outputs[y][278+ag]  += VMort[ag][4][lt][im][nm][rg][na];
                                Outputs[y][278+ag]  += VMort[ag][5][lt][im][nm][rg][na];

                        } } } } } }
////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
for(int i=278; i<289; i++) { Outputs[y][i] = Outputs[y][i]*12; }

// COUNTS TB BY US/FB
        for(int ag=0; ag<11; ag++){
            for(int tb=0; tb<6; tb++) {
                for(int lt=0; lt<2; lt++){
                    for (int im=0; im<4; im++){
                        for (int nm=0; nm<4; nm++){
                            for(int rg=0; rg<2; rg++) {
                                for(int na=0; na<3; na++){
                                    if(na<1) {
                                        Outputs[y][289] += V1[ag][2][lt][im][nm][rg][na];
                                        Outputs[y][290] += V1[ag][3][lt][im][nm][rg][na];
                                        Outputs[y][291] += V1[ag][4][lt][im][nm][rg][na];
                                        Outputs[y][292] += V1[ag][5][lt][im][nm][rg][na];
                                    } else {
                                        Outputs[y][293] += V1[ag][2][lt][im][nm][rg][na];
                                        Outputs[y][294] += V1[ag][3][lt][im][nm][rg][na];
                                        Outputs[y][295] += V1[ag][4][lt][im][nm][rg][na];
                                        Outputs[y][296] += V1[ag][5][lt][im][nm][rg][na];
                                    } } } } } } } }
// Force of infection
  Outputs[y][297] += VLjkl[0][0];
  Outputs[y][298] += VLjkl[0][1];
  Outputs[y][299] += VLjkl[1][0];
  Outputs[y][300] += VLjkl[1][1];
////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
for(int i=297; i<301; i++) { Outputs[y][i] = Outputs[y][i]*12; }

///  NEW INFECTIONS + SUPER INFECTIOPN ///
for(int ag=0; ag<11; ag++) {
  for(int lt=0; lt<2 ; lt++) {
    for(int im=0; im<4 ; im++) {
      for(int nm=0; nm<4 ; nm++) {
      for(int rg=0; rg<2; rg++) {
        for(int na=0; na<3; na++){
  Outputs[y][303+rg] += (V0[ag][0][lt][im][nm][rg][na]+V0[ag][1][lt][im][nm][rg][na]+V0[ag][2][lt][im][nm][rg][na])*VLjkl[rg][na]*NixTrans[s];
        } } } } } }
////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
for(int i=303; i<307; i++) { Outputs[y][i] = Outputs[y][i]*12; }

////////////////////// TB MORTALITY BY NATIVITY //////////////////////////////
for(int ag=0; ag<11; ag++){
    for(int lt=0; lt<2; lt++){
        for (int im=0; im<4; im++){
            for (int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++){
                    if(rg>1) { ti = 1; } else { ti = 0; }
                    if(im==3) { temp = muTbRF; } else { temp = 0; }
                        Outputs[y][307+ti]  += V0[ag][4][lt][im][nm][rg][na]*(vTMortN[ag][4 ]+temp);
                        Outputs[y][307+ti]  += V0[ag][5][lt][im][nm][rg][na]*(vTMortN[ag][5 ]+temp)*pow(1.0-TxVecZ[1],TunTxMort); //check the TxMatZ portion
                            } } } } } }
////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
for(int i=307; i<309; i++) { Outputs[y][i] = Outputs[y][i]*12; }

} ////end of mid-year results bracket
///////////////////////////////////////////////////////////////////////////////////
//////////////////////////////END MIDYEAR RESULTS//////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////                       UPDATE V0 as V1                       ///////////
///////////////////////////////////////////////////////////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++){
            for (int im=0; im<4; im++){
                for (int nm=0; nm<4; nm++){
                    for(int rg=0; rg<2; rg++) {
                        for(int na=0; na<3; na++){
                            V0[ag][tb][lt][im][nm][rg][na] = V1[ag][tb][lt][im][nm][rg][na];
} } } } } } }
} //// end of month loop!//////////////////////////////////////////////////////////
} //// end of year loop!///////////////////////////////////////////////////////////
////////////////   RE-CHECK THAT NO STATES HAVE BECOME NEGATIVE    ////////////////
NumericVector  CheckV(14784);
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++){
            for(int im=0; im<4; im++){
                for(int nm=0; nm<4; nm++){
                    for(int rg=0; rg<2; rg++) {
                        for(int na=0; na<2; na++){
                            CheckV(ag+tb*11+lt*66+im*132+nm*528+rg*2112+na*4224) = V1[ag][tb][lt][im][nm][rg][na];
} } } } } } }
///////////////////////////////////////////////////////////////////////////////////
//////                              RETURN STUFF                              /////
///////////////////////////////////////////////////////////////////////////////////
for(int i=0; i<nYrs; i++) {
    for(int j=0; j<nRes; j++) {
        Outputs2(i,j)  = Outputs[i][j];
}  }

return Rcpp::List::create(
Rcpp::Named("Outputs") = Outputs2,
Rcpp::Named("V1") = CheckV,
Rcpp::Named("V0") = CheckV0
);

}
