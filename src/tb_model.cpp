#include <Rcpp.h>
#include <algorithm>

using namespace Rcpp;

//[[Rcpp::export]]
Rcpp::List cSim(
    int                 nYrs,      //number of years to run the model
    int                 nRes,      // results
    Rcpp::NumericMatrix rDxt,      // rate of diagnosis over time
    std::vector<double> TxQualt,   // treatment quality over time
    Rcpp::NumericMatrix       InitPop,   // initial population matrix
    Rcpp::NumericMatrix       Mpfast,    // matrix of probability of fast TB progression
    std::vector<double> ExogInf,   // exogenous infection risk
    Rcpp::NumericMatrix       MpfastPI,  // matrix of probability of fast TB progression in those with partial immunity from prior infection
    Rcpp::NumericMatrix       Mrslow,    // matrix of the rate of slow TB progression
    std::vector<double> rrSlowFB,  // rate ratios that are applied to the rate of slow progression for foreign born population
    double              rfast,     // rate of fast TB progression
    double              RRcurDef,  // rate ratio of cure given treatment default
    double              rSlfCur,   // rate of self cure
    double              p_HR,
    Rcpp::NumericMatrix       vTMort,    // matrix of TB mortality
    Rcpp::NumericMatrix vRFMort,    // matrix of RF mortality
    std::vector<double> RRmuHR,
    double              muTbRF,    //factor for comorbidity btw TB and non-TB,
    std::vector<double> Birthst,   // vector of absolute births over time
    Rcpp::NumericMatrix       HrEntEx,
    Rcpp::NumericMatrix       ImmNon,
    Rcpp::NumericMatrix       ImmLat,
    Rcpp::NumericMatrix       ImmAct,
    Rcpp::NumericMatrix       ImmFst,
    Rcpp::NumericMatrix       mubt,       //background mortality over time
    std::vector<double> RelInf,     //relative infectiousness
    std::vector<double> RelInfRg,   //relative infectiousness by risk group
    std::vector<double> Vmix,       // vector of mixing parameters (sigmas)
    std::vector<double> rEmmigFB,    // rate of emmigration among the foreign born
    std::vector<double> TxVec,       //vector of TB treatment parameters
    double              TunTxMort,   //tuning for treatment mortality
    std::vector<double> rDeft,       // rate of treatment default over time
    std::vector<double> rLtScrt,
    std::vector<double> LtTxPar,     // latent treatment parameters
    Rcpp::NumericMatrix       LtDxPar,     // latent diagnosis parameters
    std::vector<double> RRdxAge,     // rate ratios for diagnosis with respect to age
    double              rRecov,      //rate from latent slow to partially immune TB
    double              pImmScen,    // lack of reactivitiy to IGRA for Sp
    std::vector<double>  EarlyTrend,  // TB natural history parameter
    std::vector<double> dLtt,        //latent diagnosis over time
    std::vector<double> pReTx,
    std::vector<double> EffLt0,
    double              EffLt,
    std::vector<double> NixTrans,
    Rcpp::NumericMatrix       dist_gen,
    Rcpp::NumericMatrix       can_go,
    std::vector<double>       dist_goal_v,
    std::vector<double>       dist_orig_v,
    std::vector<double>       diff_i_v
) {
  ////////////////////////////////////////////////////////////////////////////////
  ////////    BELOW IS A LIST OF THE VARIABLES CREATED INTERNALLY IN MODEL   /////
  ////////////////////////////////////////////////////////////////////////////////
  int           ti;
  int           s;
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
  double        HrEntExN[HrEntEx.nrow()][HrEntEx.ncol()];
  double        ImmNonN[ImmNon.nrow()][ImmNon.ncol()];
  double        ImmLatN[ImmLat.nrow()][ImmLat.ncol()];
  double        ImmFstN[ImmFst.nrow()][ImmFst.ncol()];
  double        ImmActN[ImmAct.nrow()][ImmAct.ncol()];
  double        TBImm[11][ImmAct.nrow()][2]; //age time
  double        mubtN[mubt.nrow()][mubt.ncol()];
  //  double        rIntvInitN[rIntvInit.nrow()][rIntvInit.ncol()];
  double        rDxtN[rDxt.nrow()][rDxt.ncol()];
  double        LtDxParN[LtDxPar.nrow()][LtDxPar.ncol()];
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
  double        V2[11][6][2][4][4][2][3];
  double        VMort[11][6][2][4][4][2][3];
  double        Vdx[11][6][2][4][4][2][3];
  double        VLdx[11][6][2][4][4][2][3];
  double        VNkl[2][2];  ///HIGH AND LOW RISK, NATIVITY
  double        VGjkl[2][2]; ///HIGH AND LOW RISK, NATIVITY
  double        Vjaf[4];     ///BY NUMBER OF MIXING GROUPS
  double        VLjkl[2][2];  ///HIGH AND LOW RISK, NATIVITY
  int       N;
  double   dist_genN[dist_gen.nrow()][dist_gen.ncol()];
  double   can_goN[can_go.nrow()][can_go.ncol()];
  long double   did_goN[16][16];
  double   temp_vec[16];
  double   temp_mat[4][4];
  // double   temp_mat2[528][528];
  long double     trans_mat[16][16];
  long double    trans_mat_tot[16][16];
  double    dist_i_v[16];
  double    dist_orig[4][4];
  Rcpp::NumericMatrix    trans_mat_fin(16,16);
  Rcpp::NumericVector   dist_i_v_fin(16);
  double        frc;
  double        pop_t;
  // double        sse;
  double rowsum[16];
  double mat_sum;
  int reblnc; int tb_dyn;
  Rcpp::NumericMatrix Outputs2(nYrs,nRes);

  ///////////////////////////////////////////////////////////////////////////////
  ///////                            INITIALIZE                             /////
  ///////////////////////////////////////////////////////////////////////////////
  for(int i=0; i<InitPop.nrow(); i++) {
    for(int j=0; j<InitPop.ncol(); j++) {
      InitPopN[i][j] = InitPop(i,j);
    } }
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
  for(int i=0; i<vTMort.nrow(); i++) {
    for(int j=0; j<vTMort.ncol(); j++) {
      vTMortN[i][j] = vTMort(i,j);
    } }
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
  //////was the rate of HIV and has been replaced as rate of risk factor of interest////
  //    for(int i=0; i<rRFt.nrow(); i++) {
  //        for(int j=0; j<rRFt.ncol(); j++) {
  //            rRFtN[i][j][0] = rRFt(i,j);
  //            rRFtN[i][j][1] = rRFt(i,j)*RFHrPar;
  //        } }
  for(int i=0; i<rDxt.nrow(); i++) {
    for(int j=0; j<rDxt.ncol(); j++) {
      rDxtN[i][j] = rDxt(i,j);
    } }
  for(int i=0; i<dist_gen.nrow(); i++) {
    for(int j=0; j<dist_gen.ncol(); j++) {
      dist_genN[i][j] = dist_gen(i,j);
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
              for(int na=0; na<3; na++){
                V0[ag][tb][lt][im][nm][rg][na]    = 0;
                V1[ag][tb][lt][im][nm][rg][na]    = 0;
                VMort[ag][tb][lt][im][nm][rg][na] = 0;
                Vdx[ag][tb][lt][im][nm][rg][na]   = 0;
                VLdx[ag][tb][lt][im][nm][rg][na]   = 0;
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
  for(int ag=0; ag<11; ag++) {
    for(int im=0; im<5; im++) {
      temp4V[ag][im] = (1-pow(1-MrslowN[ag][im]-rRecov,24.0))*MrslowN[ag][im];
    } }

  for(int i=0; i<can_go.nrow(); i++) {
    for(int j=0; j<can_go.ncol(); j++) {
      can_goN[i][j] = can_go(i,j);
    } }

  for(int i=0; i<16; i++) {
    for(int j=0; j<16; j++) {
      did_goN[i][j] = 0;
    } }

  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      // for(int k=0; k<11; k++) {
      //   for(int l=0; l<3; l++) {
      dist_orig[i][j]= 0;
        } }
  // } }

  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {

      temp_mat[i][j] = 0;
        } }
for(int i=0; i<16; i++) {
  for(int j=0; j<16; j++) {
      // temp_mat2[i][j] = 0;
      trans_mat[i][j] = 0;
      trans_mat_tot[i][j] = 0;
    } }

  for(int i=0; i<16; i++) {
    dist_i_v[i]=0;
    temp_vec[i]=0;
  }
  N=30;
  reblnc=1;
  tb_dyn=0;
  if (tb_dyn != 1){
    for(int ag=0; ag<11; ag++) {
      for(int tb=0; tb<6; tb++) {
        vTMortN[ag][tb] =0;
      } }
    // for(int ag=0; ag<11; ag++) {
    // for(int nm=0; nm<4; nm++) {
    //   vRFMortN[ag][nm] = 0;
    // } }
    // for(int rg=0; rg<2; rg++){RRmuHR[rg]=1; }

    muTbRF =0;
    }
//   for(int i=0; i<16; i++) {
//     for(int j=0; j<16; j++) {
// Rcpp::Rcout << "at time one size" << can_goN[i][j] << "\n"; }}
////////////////////////////////////////////////////////////////////////////////
  ///////                  UPDATING TREATMENT METERS                        //////
  ////////     THIS DIFFERENT TO MAIN MODEL DUE TO SIMPLIFIED OUTCOMES      //////
  ////////////////////////////////////////////////////////////////////////////////
  //////// TREATMENT EFFICACY UPDATED FOR TREATMENT QUALITY //////////////////////
  TxVecZ[1] = TxVec[1]*TxQualt[0];
  ///////// RATE OF TREATMENT EXIT TO CURE (LS) //////////////////////////////////
  TxVecZ[2] = TxVec[0]*TxVecZ[1] + rDeft[0]*TxVecZ[1]*RRcurDef;
  //////// RATE OF TREATMENT EXIT TO FAILURE (IN/IP) /////////////////////////////
  TxVecZ[3] = TxVec[0]*(1-TxVecZ[1]) + rDeft[0]*(1-TxVecZ[1]*RRcurDef);

  ////////////////////////////////////////////////////////////////////////////////
  //////                             StatList                                /////
  //////                             BURN IN                                 /////
  //////                           Populate model                            /////
  ////////////////////////////////////////////////////////////////////////////////
  ///////         SOURCE IN THE CODE TO DISTRIBUTE THE POPULATION            /////
  for(int ag=0; ag<11; ag++) {
    for(int im=0; im<4; im++) {
      for(int nm=0; nm<4; nm++) {
        ////////////////////        UNINFECTED/SUSCEPTIBLE POP /////////////////////////
        V0[ag][0][0][im][nm][0][0] = InitPopN[ag][0]*0.40*(1-p_HR)*dist_genN[nm][im]; //low risk US born
        V0[ag][0][0][im][nm][1][0] = InitPopN[ag][0]*0.40*(p_HR)*dist_genN[nm][im]; //high risk US born

        V0[ag][0][0][im][nm][0][2] = InitPopN[ag][1]*0.40*(1-p_HR)*dist_genN[nm][im]; //low risk non-US born
        V0[ag][0][0][im][nm][1][2] = InitPopN[ag][1]*0.40*(p_HR)*dist_genN[nm][im];  //high risk non-US born
        /////////////////////////   LATENT SLOW INFECTED POP  //////////////////////////
        V0[ag][2][0][im][nm][0][0] = InitPopN[ag][0]*0.60*(1-p_HR)*dist_genN[nm][im];
        V0[ag][2][0][im][nm][1][0] = InitPopN[ag][0]*0.60*(p_HR)*dist_genN[nm][im];

        V0[ag][2][0][im][nm][0][2] = InitPopN[ag][1]*0.60*(1-p_HR)*dist_genN[nm][im];
        V0[ag][2][0][im][nm][1][2] = InitPopN[ag][1]*0.60*(p_HR)*dist_genN[nm][im];
      } } }

  //////create a 2nd array with same dimensions as V0 & populate w/ same values//
  for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++){
          for(int rg=0; rg<2; rg++){
            for (int na=0; na<3; na++){
              V1[ag][tb][0][im][nm][rg][na]  = V0[ag][tb][0][im][nm][rg][na];
            } } } } } }

  ////////////////////////RUN THE MODEL FOR 3000 MONTHS /////////////////////////
  for(int m=0; m<3001; m++) {
    /////////////////////////////////START BURN IN//////////////////////////////////
    ////////////////////////////////////BIRTHS//////////////////////////////////////
    ///////////USE DISTRIBUTION TO POPULATE THE MODEL ACROSS RISK GROUPS////////////
    for(int im=0; im<4; im++) {
      for(int nm=0; nm<4; nm++){
        V1[0][0][0][im][nm][0][0]  += Birthst[0]*dist_genN[nm][im]*(1-p_HR);

        V1[0][0][0][im][nm][1][0]  += Birthst[0]*dist_genN[nm][im]*(p_HR);
      } }
    //////////////////////////////////IMMIGRATION///////////////////////////////////
    /////////////SINCE WE NO LONGER HAVE PREVIOUS TREATMENT AS A STATE SHOULD WE
    /////////////HAVE IMMIGRANTS IMMIGRATE STRAIGHT INTO PARTIALLY IMMUNE STATE
    /////////////QUESTION: CURRENTLY ACTIVE TX CURES TO LATENT DISEASE; IS THIS RIGHT?
    for(int ag=0; ag<11; ag++) {
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++){

            V1[ag][0][0][im][nm][0][1]   += ImmNonN[0][ag]*dist_genN[nm][im]*(1-p_HR);  // NO TB, low risk
            V1[ag][0][0][im][nm][1][1]   += ImmNonN[0][ag]*dist_genN[nm][im]*(p_HR);    // NO TB, high risk

            V1[ag][2][0][im][nm][0][1]   += ImmLatN[0][ag]*dist_genN[nm][im]*(1-p_HR); // LATENT SLOW TB, low risk
            V1[ag][2][0][im][nm][1][1]   += ImmLatN[0][ag]*dist_genN[nm][im]*(p_HR);   // LATENT SLOW TB, high risk

            V1[ag][3][0][im][nm][0][1]   += TBImm[ag][0][1]*dist_genN[nm][im]*(1-p_HR);   // LATENT FAST, low risk
            V1[ag][3][0][im][nm][1][1]   += TBImm[ag][0][1]*dist_genN[nm][im]*(p_HR);   // LATENT FAST, high risk

            V1[ag][4][0][im][nm][0][1]   += TBImm[ag][0][0]*dist_genN[nm][im]*(1-p_HR);   //ACTIVE TB, low risk
            V1[ag][4][0][im][nm][1][1]   += TBImm[ag][0][0]*dist_genN[nm][im]*(p_HR);   //ACTIVE TB, high risk
          } } }
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
            for(int na=0; na<3; na++){
              for(int tb=0; tb<5; tb++) {
                V1[ag][tb][0][im][nm][rg][na]  -= V0[ag][tb][0][im][nm][rg][na]*(mubtN[0][ag]*RRmuHR[rg]+vRFMortN[ag][nm]+vTMortN[ag][tb]);
              }
  ////////////////          MORTALITY WITH TB TREATMENT         ////////////////////
                V1[ag][5 ][0][im][nm][rg][na]  -= V0[ag][5 ][0][im][nm][rg][na]*(mubtN[0][ag]*RRmuHR[rg]+vRFMortN[ag][nm]+vTMortN[ag][5 ]*pow(1.0-TxVecZ[1],TunTxMort)); //check the mortality in param
            } } } } }
    /////////////////////////////////////AGING///////////////////////////////////////
    for(int ag=0; ag<10; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++) {
                /////          IF AGE > 4, IT TAKES 120 MONTHS TO LEAVE AGE GROUP          /////
                if(ag>0) {
                  temp2 = 120;
                  /////          IF AGE < 4, IT TAKES 60 MONTHS TO LEAVE AGE GROUP           /////
                } else {
                  temp2 = 60;
                }
                temp = V0[ag][tb][0][im][nm][rg][na]/temp2;
                V1[ag  ][tb][0][im][nm][rg][na]  -= temp;
                V1[ag+1][tb][0][im][nm][rg][na]  += temp;
              } } } } } }
    ///////////////////////// NEW FB -> ESTABLISHED FB ///////////////////////////////
    ///////////////////////// TWO YEARS FOR TRANSITION ///////////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {
              temp = V0[ag][tb][0][im][nm][rg][1] / 24;
              V1[ag][tb][0][im][nm][rg][1]  -= temp;
              V1[ag][tb][0][im][nm][rg][2]  += temp;
            } } } } }
    //////////////////////////// HIGH-RISK ENTRY/EXIT ////////////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int na=0; na<3; na++) {
              temp  = V0[ag][tb][0][im][nm][0][na]*HrEntExN[ag][0];
              temp2 = V0[ag][tb][0][im][nm][1][na]*HrEntExN[ag][1];
              //THESE CODES WERE UPDATED, BUT REMAIN ALMOST THE SAME
              V1[ag][tb][0][im][nm][0][na]  += temp2-temp;
              V1[ag][tb][0][im][nm][1][na]  += temp-temp2;

            } } } } }
if (tb_dyn==1){
    ///////////////////////////         BREAK DOWN      /////////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int im=0; im<4 ; im++) {
        for(int nm=0; nm<4; nm++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++) {
              temp  =  V0[ag][2][0][im][nm][rg][na]*MrslowN[ag][im]*rrSlowFB[na];  // Ls
              temp2  = V0[ag][3][0][im][nm][rg][na]*rfast;

              V1[ag][2][0][im][nm][rg][na]  -= temp;
              V1[ag][3][0][im][nm][rg][na]  -= temp2;
              V1[ag][4][0][im][nm][rg][na]  += temp+temp2;
            } } } } }
    ///////////////////////////   LATENT SLOW TO SAFE   /////////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int im=0; im<4 ; im++) {
        for(int nm=0; nm<4; nm++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++) {
              temp  = V0[ag][2][0][im][nm][rg][na]*rRecov;  // Ls
              V1[ag][2][0][im][nm][rg][na]  -= temp;
              V1[ag][1][0][im][nm][rg][na]  += temp;
            } } } } }
    ////////////////////////////////// SELF CURE/////////////////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int im=0; im<4 ; im++) {
        for(int nm=0; nm<4; nm++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++) {
              temp  = V0[ag][4 ][0][im][nm][rg][na]*rSlfCur;
              V1[ag][4 ][0][im][nm][rg][na]  -= temp;
              V1[ag][2 ][0][im][nm][rg][na]  += temp;
            } } } } }
    ////////////////////////////  TRANSMISSION RISK  ////////////////////////////////
    for(int i=0; i<2; i++) {
      for(int j=0; j<2; j++) {
        VNkl [i][j] = 0;
        VGjkl[i][j] = 0;   // set to zero
      } }
    // Step 1
    // take total population of mixing groups
    for(int ag=0; ag<11; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            ///////// RISK FACTOR FREE & LOW RISK US BORN
            ///////// ranges from ~8-160 persons per time step (wide variation due to risk groups? )
            VNkl[0][0]  += V0[ag][tb][0][im][nm][0][0];
            ///////// RISK FACTOR FREE & HIGH RISK US BORN
            ///////// ranges under 1 every time step
            VNkl[1][0]  += V0[ag][tb][0][im][nm][1][0];
            ///////// RISK FACTOR FREE & LOW RISK NON US BORN
            ///////// ranges from .5 -10 people between time steps
            VNkl[0][1]  += V0[ag][tb][0][im][nm][0][1] + V0[ag][tb][0][im][nm][0][2];
            ///////// also ranges from .5 -10 people between time steps
            ///////// RISK FACTOR FREE & HIGH RISK NON US BORN
            VNkl[1][1]  += V0[ag][tb][0][im][nm][1][1] + V0[ag][tb][0][im][nm][1][2];
          } } } }
    // Step 2  (active TB)
    // Number of active cases* relative infectiousness
    for(int ag=0; ag<11; ag++) {
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++) {
          /////////  LOW RISK US BORN
          VGjkl[0][0]  +=  V0[ag][4][0][im][nm][0][0]                              *RelInf[4];
          ///////// HIGH RISK US BORN
          VGjkl[1][0]  +=  V0[ag][4][0][im][nm][1][0]                              *RelInf[4];
          ///////// LOW RISK NON US BORN
          VGjkl[0][1]  += (V0[ag][4][0][im][nm][0][1] + V0[ag][4][0][im][nm][0][2])*RelInf[4];
          /////////  HIGH RISK NON US BORN
          VGjkl[1][1]  += (V0[ag][4][0][im][nm][1][1] + V0[ag][4][0][im][nm][1][2])*RelInf[4];
        } } }
    // Step 2 (treated TB)
    // No contribution to force of infection

    // Step 3
    Vjaf[0]  = (RelInfRg[0]*VGjkl[0][0]         +       //LOW RISK US BORN
      RelInfRg[1]*VGjkl[1][0]*Vmix[0] +       //HIGH RISK US BORN
      RelInfRg[0]*VGjkl[0][1]*Vmix[1] +       //LOW RISK FOREIGN BORN
      RelInfRg[1]*VGjkl[1][1]*Vmix[0]*Vmix[1]) //HIGH RISK FOREIGN BORN
      /   (RelInfRg[0]*VNkl[0][0]           +
        RelInfRg[1]*VNkl[1][0]*Vmix[0]   +
        RelInfRg[0]*VNkl[0][1]*Vmix[1]   +
        RelInfRg[1]*VNkl[1][1]*Vmix[0]*Vmix[1] +
        1e-12);

    Vjaf[1]  = (Vmix[0]*VGjkl[1][1] + VGjkl[0][1]) / (Vmix[0]*VNkl[1][1] +VNkl[0][1]+1e-12); //high risk us born population


    Vjaf[2]  = (Vmix[1]*VGjkl[1][1] + VGjkl[1][0]) / (Vmix[1]*VNkl[1][1] + VNkl[1][0]+1e-12); //low risk non us born population



    Vjaf[3]  = VGjkl[1][1] / (VNkl[1][1]+1e-12); //high risk foreign born population



    // Step 4
    /// LOW RISK US BORN
    VLjkl[0 ][0 ]  = RelInfRg[0]*Vjaf[0];
    ///////// HIGH RISK US BORN
    VLjkl[1 ][0 ]  = RelInfRg[1]*(Vmix[0]*Vjaf[1]+(1-Vmix[0])*Vjaf[0]);
    ///////// LOW RISK NON US BORN
    VLjkl[0 ][1 ]  = RelInfRg[0]*(Vmix[1]*Vjaf[2]+(1-Vmix[1])*Vjaf[0]) + ExogInf[0];
    ///////// HIGH RISK NON US BORN
    ///check the use of RelInfRg here as beta, might need to be a combo param but unclear check the old param file
    //beta HR* [sigmaHR* (sigmaFB*Vjaf[3]+(1-sigmaFB)*Vjaf[1])+(1-sigmaHR)*(sigmaFB*Vjaf[2])+(1-sigmaFB)*Vjaf[0]]+ExogInf
    VLjkl[1 ][1 ]  = RelInfRg[1]*(Vmix[0]*(Vmix[1]*Vjaf[3]+(1-Vmix[0])*Vjaf[1])+
      ((1-Vmix[0])*Vmix[1]*Vjaf[2]+(1-Vmix[1])*Vjaf[0])
    )+ExogInf[0];

    ///////////////////////////////INFECTION///////////////////////////////////////
    ///////////////////////for all age groups, risk groups/////////////////////////
    ///////INFECTION IS CALCULATED WITH THE FORCE OF INFECTION BY RISK GROUP///////
    /////// THE TOTAL NUMBER OF INFECTED THEN ENTER BOTH THE LATENT SLOW &  ///////
    ///////& LATENT FAST TB STATES DEPENDENT ON PROBABILITY OF FAST LATENCY ///////
    ///////    PEOPLE ARE REMOVED FROM SUSCEPTIBLE, PARTIALLY IMMUNE &      ///////
    ///////                   LATENT SLOW STATES                            ///////
    ///////////////////////////////////////////////////////////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++) {
          for(int rg=0; rg<2; rg++) {
            for (int na=0; na<3; na++){

              ///////////////////////////////   SUCEPTIBLE  /////////////////////////////////
              temp = V0[ag][0][0][im][nm][rg][na]*EarlyTrend[m]*VLjkl[rg][na];
              //////////////////////////// REMOVE FROM SUSCEPTIBLE //////////////////////////
              V1[ag][0][0][im][nm][rg][na]  -= temp;
              //////////////////////////////// LATENT TB SLOW ///////////////////////////////
              V1[ag][2][0][im][nm][rg][na]  += temp*MpslowN[ag][im];
              //////////////////////////////// LATENT TB FAST ///////////////////////////////
              V1[ag][3][0][im][nm][rg][na]  += temp*MpfastN[ag][im];
              ///////////////////////////////////////////////////////////////////////////////

              /////////////////////////////// SUPER-INFECTION SP ////////////////////////////
              temp = V0[ag][1][0][im][nm][rg][na]*VLjkl[rg][na];
              V1[ag][1][0][im][nm][rg][na] -= temp;
              V1[ag][2][0][im][nm][rg][na] += temp*MpslowPIN[ag][im];
              V1[ag][3][0][im][nm][rg][na] += temp*MpfastPIN[ag][im];
              ///////////////////////////////////////////////////////////////////////////////

              /////////////////////////////// SUPER-INFECTION LS ////////////////////////////
              temp = V0[ag][2][0][im][nm][rg][na]*VLjkl[rg][na];
              V1[ag][2][0][im][nm][rg][na]  -= temp;
              V1[ag][2][0][im][nm][rg][na]  += temp*MpslowPIN[ag][im];
              V1[ag][3][0][im][nm][rg][na]  += temp*MpfastPIN[ag][im];

            } } } } }

    ////////////////////TB DIAGNOSIS AND TX INITIATION PUBLIC ///////////////////////
    ///////////////////////for all age groups, living cond///////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++) {
              temp  = V0[ag][4 ][0][im][nm][rg][na]*rDxtN[0][rg]/RRdxAge[ag]*EarlyTrend[m];
              V1[ag][4 ][0][im][nm][rg][na]  -= temp;
              V1[ag][5 ][0][im][nm][rg][na]  += temp;
            } } } } }
    ///////////////////////////   TREATMENT OUTCOMES    /////////////////////////////
    ///////////////////////for all age groups, risk groups///////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++){
              //////////CURES//////////////////////////////////////////////////////////////////
              temp  = V0[ag][5][0][im][nm][rg][na]*TxVecZ[2];
              V1[ag][5][0][im][nm][rg][na]  -= temp;
              V1[ag][2][0][im][nm][rg][na]  += temp;
              //////////FAILURES(INCLUDING TREATMENT DEFAULT)//////////////////////////////////
              temp  = V0[ag][5][0][im][nm][rg][na]*TxVecZ[3]; ///check this line
              V1[ag][5][0][im][nm][rg][na]  -= temp;
              V1[ag][4][0][im][nm][rg][na]  += temp;
            } } } } }
    } //end of TB loop

    /////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////RESET POPULATION SIZE/////////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int i=0; i<2; i++) {
        InitPopZ[ag][i] = 0; //i is for us vs non us born
      } }
    for(int ag=0; ag<11; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {
              InitPopZ[ag][0]  += V1[ag][tb][0][im][nm][rg][0];

              InitPopZ[ag][1]  += V1[ag][tb][0][im][nm][rg][1]+V1[ag][tb][0][im][nm][rg][2];

            } } } } }
    for(int ag=0; ag<11; ag++) {
      for(int i=0; i<2; i++) {                        // factor for pop size reset
        InitPopZ[ag][i] = InitPopN[ag][i]/(InitPopZ[ag][i]+1e-12);

      } }
    for(int ag=0; ag<11; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {  // reset pop to InitPop
              V1[ag][tb][0][im][nm][rg][0]  = V1[ag][tb][0][im][nm][rg][0]*InitPopZ[ag][0];
              V1[ag][tb][0][im][nm][rg][1]  = V1[ag][tb][0][im][nm][rg][1]*InitPopZ[ag][1];
              V1[ag][tb][0][im][nm][rg][2]  = V1[ag][tb][0][im][nm][rg][2]*InitPopZ[ag][1];

              for (int na=0; na<3; na++){
                V0[ag][tb][0][im][nm][rg][na]  = V1[ag][tb][0][im][nm][rg][na];
              } } } } } }

    ///////////////////////////////////////////////////////////////////////////////
  } ///////////////////////////END BURN IN///////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  Rcpp::NumericVector  CheckV0(12672);
  for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
      for(int lt=0; lt<2; lt++){
        for(int im=0; im<4; im++){
          for(int nm=0; nm<4; nm++){
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++) {
                CheckV0(ag+tb*11+lt*66+im*132+nm*528+rg*2112+na*4224) = V1[ag][tb][lt][im][nm][rg][na];
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
      TxVecZ[1] = TxVec[1]*TxQualt[s];
      ///////// RATE OF TREATMENT EXIT TO CURE (LS) //////////////////////////////////
      TxVecZ[2] = TxVec[0]*TxVecZ[1] + rDeft[s]*TxVecZ[1]*RRcurDef;
      ///////// RATE OF TREATMENT EXIT TO ACTIVE TB //////////////////////////////////
      TxVecZ[3] = TxVec[0]*(1-TxVecZ[1])*(1-pReTx[s]) + rDeft[s]*(1-TxVecZ[1])*RRcurDef*(1-pReTx[s]);
      ///////// RATE OF TREATMENT EXIT TO RE TREATMENT////////////////////////////////
      TxVecZ[4] = TxVec[0]*(1-TxVecZ[1])*(pReTx[s]) + rDeft[s]*(1-TxVecZ[1])*RRcurDef*pReTx[s];
      ////////////////////// P(TREATMENT COMPLETION) /////////////////////////////////
      TxVecZ[5] = TxVec[0]*(1-(1.0-TxVecZ[1])*pReTx[s]);
      /////////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////BIRTHS//////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////////////
      for (int im=0; im<4; im++) {
        for (int nm=0; nm<4; nm++) {
          /////LOW RISK GROUP BIRTHS//////////////////////////////////////////////////////
          V1[0][0][0][im][nm][0][0]  += Birthst[s]*dist_genN[nm][im]*(1-p_HR);
          /////HIGH RISK GROUP BIRTHS//////////////////////////////////////////////////////
          V1[0][0][0][im][nm][1][0]  += Birthst[s]*dist_genN[nm][im]*(p_HR);
        } }
      ///////////////////////////////// IMMIGRATION ///////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        // for (int lt=0; lt<2; lt++) {
          for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {

              V1[ag][0][0][im][nm][0][1]   += ImmNonN[s][ag]*dist_genN[nm][im]*(1-p_HR);  // NO TB, low risk
              V1[ag][0][0][im][nm][1][1]   += ImmNonN[s][ag]*dist_genN[nm][im]*(p_HR);    // NO TB, high risk

              V1[ag][2][0][im][nm][0][1]   += ImmLatN[s][ag]*dist_genN[nm][im]*(1-p_HR); // LATENT SLOW TB, low risk
               V1[ag][2][0][im][nm][1][1]   += ImmLatN[s][ag]*dist_genN[nm][im]*(p_HR);   // LATENT SLOW TB, high risk

              V1[ag][3][0][im][nm][0][1]   += TBImm[ag][s][1]*dist_genN[nm][im]*(1-p_HR);   // LATENT FAST, low risk
               V1[ag][3][0][im][nm][1][1]   += TBImm[ag][s][1]*dist_genN[nm][im]*(p_HR);   // LATENT FAST, high risk

              V1[ag][4][0][im][nm][0][1]   += TBImm[ag][s][0]*dist_genN[nm][im]*(1-p_HR);   //ACTIVE TB, low risk
              V1[ag][4][0][im][nm][1][1]   += TBImm[ag][s][0]*dist_genN[nm][im]*(p_HR);   //ACTIVE TB, high risk
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
          for(int im=0; im<4 ; im++) {
            for(int nm=0; nm<4; nm++) {
              ////we need to decide about TB co-mortality factor
              if(im == 3) {
                temp = muTbRF;
              } else { temp = 0; }
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
                    temp += VMort[ag][tb][lt][im][nm][rg][na];
                    // Rcpp::Rcout << "total mortality at time" << s << "is" << temp << "\n";

                  }

                } } } } } }


      /////////////////////////////////////AGING///////////////////////////////////////
      for(int ag=0; ag<10; ag++) {
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
          for(int lt=0; lt<2; lt++) {
            for (int im=0; im<4; im++) {
              for (int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  temp = V0[ag][tb][lt][im][nm][rg][1]/24;
                  V1[ag][tb][lt][im][nm][rg][1]  -= temp;
                  V1[ag][tb][lt][im][nm][rg][2]  += temp;
                } } } } } }
      //////////////////////////// HIGH-RISK ENTRY/EXIT ////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int na=0; na<3; na++) {
                  temp  = V0[ag][tb][lt][im][nm][0][na]*HrEntExN[ag][0];
                  temp2 = V0[ag][tb][lt][im][nm][1][na]*HrEntExN[ag][1];

                  V1[ag][tb][lt][im][nm][0][na]  += temp2-temp;
                  V1[ag][tb][lt][im][nm][1][na]  += temp-temp2;
                } } } } } }
if (tb_dyn==1){
                  ////////////////////////////  TRANSMISSION RISK  ////////////////////////////////
      for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++) {
          VNkl [i][j] = 0;
          VGjkl[i][j] = 0;   // set to zero
        } }
      // Step 1
      // take total population of mixing groups
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                ///////// RISK FACTOR FREE & LOW RISK US BORN
                VNkl[0][0]  += V0[ag][tb][lt][im][nm][0][0];
                ///////// RISK FACTOR FREE & HIGH RISK US BORN
                VNkl[1][0]  += V0[ag][tb][lt][im][nm][1][0];
                ///////// RISK FACTOR FREE & LOW RISK NON US BORN
                VNkl[0][1]  += V0[ag][tb][lt][im][nm][0][1] + V0[ag][tb][lt][im][nm][0][2];
                ///////// RISK FACTOR FREE & HIGH RISK NON US BORN
                VNkl[1][1]  += V0[ag][tb][lt][im][nm][1][1] + V0[ag][tb][lt][im][nm][1][2];
              } } } } }
      // Step 2  (active TB)
      // Number of active cases* relative infectiousness
      for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++) {
          for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
              /////////  LOW RISK US BORN
              VGjkl[0][0]  +=  V0[ag][4][lt][im][nm][0][0]                              *RelInf[4];
              ///////// HIGH RISK US BORN
              VGjkl[1][0]  +=  V0[ag][4][lt][im][nm][1][0]                              *RelInf[4];
              ///////// LOW RISK NON US BORN
              VGjkl[0][1]  += (V0[ag][4][lt][im][nm][0][1] + V0[ag][4][lt][im][nm][0][2])*RelInf[4];
              /////////  HIGH RISK NON US BORN
              VGjkl[1][1]  += (V0[ag][4][lt][im][nm][1][1] + V0[ag][4][lt][im][nm][1][2])*RelInf[4];
            } } } }
      // Step 2 (treated TB)
      // No contribution to force of infection

      // Step 3
      Vjaf[0]  = (RelInfRg[0]*VGjkl[0][0]          +       //LOW RISK US BORN
        RelInfRg[1]*VGjkl[1][0]*Vmix[0]  +       //HIGH RISK US BORN
        RelInfRg[0]*VGjkl[0][1]*Vmix[1]  +       //LOW RISK FOREIGN BORN
        RelInfRg[1]*VGjkl[1][1]*Vmix[0]*Vmix[1]) //HIGH RISK FOREIGN BORN
        /              (RelInfRg[0]*VNkl[0][0]          +
          RelInfRg[1]*VNkl[1][0]*Vmix[0]  +
          RelInfRg[0]*VNkl[0][1]*Vmix[1]  +
          RelInfRg[1]*VNkl[1][1]*Vmix[0]*Vmix[1] +
          1e-12);

      Vjaf[1]  = (Vmix[0]*VGjkl[1][1] + VGjkl[0][1]) / (Vmix[0]*VNkl[1][1] +VNkl[0][1]+1e-12); //high risk us born population


      Vjaf[2]  = (Vmix[1]*VGjkl[1][1] + VGjkl[1][0]) / (Vmix[1]*VNkl[1][1] + VNkl[1][0]+1e-12); //low risk non us born population

      Vjaf[3]  = VGjkl[1][1] / (VNkl[1][1]+1e-12); //high risk foreign born population

      // Step 4
      /// LOW RISK US BORN
      VLjkl[0 ][0 ]  = RelInfRg[0]*Vjaf[0];
      ///////// HIGH RISK US BORN
      VLjkl[1 ][0 ]  = RelInfRg[1]*(Vmix[0]*Vjaf[1]+(1-Vmix[0])*Vjaf[0]);
      ///////// LOW RISK NON US BORN
      VLjkl[0 ][1 ]  = RelInfRg[0]*(Vmix[1]*Vjaf[2]+(1-Vmix[1])*Vjaf[0]) + ExogInf[0];
      ///////// HIGH RISK NON US BORN
      ///check the use of RelInfRg here as beta, might need to be a combo param but unclear check the old param file
      //beta HR* [sigmaHR* (sigmaFB*Vjaf[3]+(1-sigmaFB)*Vjaf[1])+(1-sigmaHR)*(sigmaFB*Vjaf[2])+(1-sigmaFB)*Vjaf[0]]+ExogInf
      VLjkl[1 ][1 ]  = RelInfRg[1]*(Vmix[0]*(Vmix[1]*Vjaf[3]+(1-Vmix[0])*Vjaf[1])+
        ((1-Vmix[0])*Vmix[1]*Vjaf[2]+(1-Vmix[1])*Vjaf[0])
      )+ExogInf[0];

      ///////////////////////////////INFECTION///////////////////////////////////////
      ///////////////////////for all age groups, risk groups/////////////////////////
      ///////INFECTION IS CALCULATED WITH THE FORCE OF INFECTION BY RISK GROUP///////
      /////// THE TOTAL NUMBER OF INFECTED THEN ENTER BOTH THE LATENT SLOW &  ///////
      ///////& LATENT FAST TB STATES DEPENDENT ON PROBABILITY OF FAST LATENCY ///////
      ///////    PEOPLE ARE REMOVED FROM SUSCEPTIBLE, PARTIALLY IMMUNE &      ///////
      ///////                   LATENT SLOW STATES                            ///////
      ///////////////////////////////////////////////////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++) {
          for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
              for(int rg=0; rg<2; rg++) {
                for (int na=0; na<3; na++){

                  ///////////////////////////////   SUCEPTIBLE  /////////////////////////////////
                  temp = V0[ag][0][lt][im][nm][rg][na]*VLjkl[rg][na]*NixTrans[s];
                  //////////////////////////// REMOVE FROM SUSCEPTIBLE //////////////////////////
                  V1[ag][0][lt][im][nm][rg][na]  -= temp;
                  //////////////////////////////// LATENT TB SLOW ///////////////////////////////
                  V1[ag][2][lt][im][nm][rg][na]  += temp*MpslowN[ag][im];
                  //////////////////////////////// LATENT TB FAST ///////////////////////////////
                  V1[ag][3][lt][im][nm][rg][na]  += temp*MpfastN[ag][im];
                  ///////////////////////////////////////////////////////////////////////////////

                  /////////////////////////////// SUPER-INFECTION SP ////////////////////////////
                  temp = V0[ag][1][lt][im][nm][rg][na]*VLjkl[rg][na]*NixTrans[s];
                  V1[ag][1][lt][im][nm][rg][na] -= temp;
                  V1[ag][2][lt][im][nm][rg][na] += temp*MpslowPIN[ag][im];
                  V1[ag][3][lt][im][nm][rg][na] += temp*MpfastPIN[ag][im];
                  ///////////////////////////////////////////////////////////////////////////////

                  /////////////////////////////// SUPER-INFECTION LS ////////////////////////////
                  temp = V0[ag][2][lt][im][nm][rg][na]*VLjkl[rg][na]*NixTrans[s];
                  V1[ag][2][lt][im][nm][rg][na]  -= temp;
                  V1[ag][2][lt][im][nm][rg][na]  += temp*MpslowPIN[ag][im];
                  V1[ag][3][lt][im][nm][rg][na]  += temp*MpfastPIN[ag][im];

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
                  //////REMOVE FROM LATENT SLOW AND LATENT FAST AND PLACE IN ACTIVE DISEASE
                  V1[ag][2][lt][im][nm][rg][na]  -= temp;  //REMOVE FROM LATENT SLOW
                  V1[ag][3][lt][im][nm][rg][na]  -= temp2; //REMOVE FROM LATENT FAST
                  V1[ag][4][lt][im][nm][rg][na]  += (temp+temp2); //PLACE IN ACTIVE DISEASE
                  //////PROGRESSION OF DISEASE IF LATENT TREATMENT FAILS
                  temp  = V0[ag][2][lt][im][nm][rg][na]*MrslowN[ag][im]*rrSlowFB[na];//*(1-EffLt0[s]);
                  temp2 = V0[ag][3][lt][im][nm][rg][na]*MrslowN[ag][im]*rrSlowFB[na];//*(1-EffLt0[s]);

                  V1[ag][2][lt][im][nm][rg][na]  -= temp;
                  V1[ag][3][lt][im][nm][rg][na]  -= temp2;
                  V1[ag][4][lt][im][nm][rg][na]  += temp+temp2;
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
              /////////DOES OUR GENERIC RISK GROUP AFFECT LT DX PARAMETERS?
              /////////should there be a joint number for foreign born and high risk (more elevated than high risk; particularly for the rTBN BCG)
              ////////////// US BORN, LOW RISK  //////////////////
              if( rg==0 & na==0) {
                rTbP = rLtScrt[s]*LtDxParN[0][0];
                rTbN = rLtScrt[s]*LtDxParN[0][1];
              }
              //////////// NON US BORN  ////////////////
              if(rg==0 & na > 0) {
                rTbP = rLtScrt[s]*LtDxParN[2][0];
                rTbN = rLtScrt[s]*LtDxParN[2][1];
              }
              ////////////// US BORN, HIGH RISK  /////////////////
              if(rg==1) {
                rTbP = rLtScrt[s]*LtDxParN[1][0];
                rTbN = rLtScrt[s]*LtDxParN[1][1];
              }
              for(int ag=0; ag<11; ag++) {
                ////if we want a count of the number of people with LTBI diagnosis we will need to create a new vector
                ////////////// HAVE LTBI
                temp  = V0[ag][2][0][im][nm][rg][na]*rTbP;
                temp2 = V0[ag][3][0][im][nm][rg][na]*rTbP;
                for(int tb=2; tb<3; tb++){
                  VLdx[ag][tb][0][im][nm][rg][na] = temp+temp2;
                }
                ////////////// Dont have LTBI
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
        for(int lt=0; lt<2; lt++) {
          for(int im=0; im<4 ; im++) {
            for(int nm=0; nm<4; nm++) {
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                  temp  = V0[ag][2][lt][im][nm][rg][na]*rTbP*dLtt[s]; // tx completion
                  temp2 = V0[ag][3][lt][im][nm][rg][na]*rTbP*dLtt[s]; // tx completion
                  temp3 = V0[ag][2][lt][im][nm][rg][na]*rTbP*LtTxPar[1]; // default
                  temp4 = V0[ag][3][lt][im][nm][rg][na]*rTbP*LtTxPar[1]; // default

                  V1[ag][2][lt][im][nm][rg][na]  -= temp+temp3; //remove from latent slow
                  V1[ag][3][lt][im][nm][rg][na]  -= temp2+temp4;  //remove from latent fast
                  V1[ag][1][1][im][nm][rg][na]   += (temp+temp2)*EffLt;//*EffLt0[s]; //exit to cure
                  V1[ag][2][1][im][nm][rg][na]   += (temp+temp2)*(1-EffLt);//*(1-EffLt0[s]*EffLt);  //tx comp fail to latent slow
                  ///placed in tx naive for consistency in last model
                  V1[ag][2][0][im][nm][rg][na]  += (temp3+temp4); //latent tx default to latent slow
                } } } } } }
      ///////////////////// TB DIAGNOSIS AND TX INITIATION  /////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++) {
          for(int im=0; im<4 ; im++) {
            for(int nm=0; nm<4; nm++) {
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                  temp  = V0[ag][4][lt][im][nm][rg][na]*rDxtN[s][rg]/RRdxAge[ag];
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
          for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                  // Cures back to Ls state (should this actually be into partially immune state?)
                  temp=V0[ag][5][lt][im][nm][rg][na]*TxVecZ[2];
                  V1[ag][5][lt][im][nm][rg][na]  -= temp;
                  V1[ag][2][lt][im][nm][rg][na]  += temp;
                  ///// EXIT TO ACTIVE DISEASE //////
                  temp=V0[ag][5][lt][im][nm][rg][na]*TxVecZ[3];
                  V1[ag][5][lt][im][nm][rg][na]  -= temp;
                  V1[ag][4][lt][im][nm][rg][na]  += temp;
                  ///// EXIT TO TB RETREATMENT //////
                  temp=V0[ag][5][lt][im][nm][rg][na]*TxVecZ[4];
                  V1[ag][5][lt][im][nm][rg][na]  -= temp;
                  V1[ag][5][lt][im][nm][rg][na]  += temp;
                } } } } } }
      }//end of TB loop
//////////////////////////////////////////////////////////////////////////////////
////////////////////////// REBALANCE THE POPULATION //////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
////only create the distirbution; do not apply until after the outputs are//////
//////created for this timestep.
if (reblnc==1){
////// need to define the current distribution of persons across the RG at this timestep
////// RESET ALL THE VARIABLES
for(int ag=0; ag<11; ag++) {
  for(int na=0; na<3; na++) {
    mat_sum=0;
    for (int i=0; i<4; i++){
      for (int j=0; j<4; j++){
        dist_orig[i][j]=0;
        temp_mat[i][j]=0;
        dist_orig_v[i]=0;
      } }
    for(int im=0; im<4; im++) {
      for(int nm=0; nm<4; nm++) {
        for(int tb=0; tb<6; tb++) {
           for(int lt=0; lt<2; lt++) {
             for(int rg=0; rg<2; rg++) {
             temp_mat[nm][im]  += V1[ag][tb][lt][im][nm][rg][na];
             } } }
        mat_sum +=  temp_mat[nm][im];
      } }
//////insert error that distribution does not sum to one and then go from there;
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++) {
          dist_orig[nm][im]  = temp_mat[nm][im]/mat_sum; // determine the proportions
          dist_orig_v[(nm)+(im*4)] = dist_orig[nm][im];
    // if (std::isnan(dist_orig[nm][im]) > 0  ){
    // Rcpp::Rcout << "@ nm " << nm <<"& im "<< im << "orig is nan @" << s << "\n";}
    //  if (std::isnan(dist_orig_v[(nm)+(im*4)]) > 0){
    // Rcpp::Rcout << "@ nm " << nm <<"& im "<< im << "orig is " << dist_orig_v[(nm)+(im*4)]<< "\n";}
        } }
      // temp=0;
      // // // // /////check that the distributions sum to 1;
      // for (int i=0; i<4; i++){
      //   for (int j=0; j<4; j++){
      //     temp += dist_orig[i][j];
      // // Rcpp::Rcout <<"dist_orig is" << dist_orig[i][j] <<"at i "<< i << "& j " << j << "\n";
      // } }
      //  Rcpp::Rcout <<"sum of dist_orig is" << temp << "at s= "<< s<<"\n";


      //dist_orig[nm][im]+= V1[ag][tb][lt][im][nm][rg][na];

      for (int i=0; i<16; i++){
        dist_i_v[i] = dist_orig_v[i];
      }

      for(int n=0; n<30; n++){
      //   // /////// CALCULATE DISTANCE FROM CURRENT DISTRIBUTION TO GOAL DISTRIBUTION /////
        for (int i=0; i<16; i++){
          diff_i_v[i] = dist_i_v[i] - dist_goal_v[i];
        }
      //
      //   // for (int i=0; i<16; i++){
      //   //   Rcpp::Rcout <<"at s=" << s << "diff_i_v is" <<  diff_i_v[i] << "\n";
      //   // }
      //   // //////////                  CREATE TRANSITION MATRIX                    ////////
        for (int r=0; r<16; r++){
          for (int c=0; c<16; c++){
            trans_mat[r][c] = 0;
          } }
        for (int r=0; r<16; r++){
          for (int c=0; c<16; c++){
            /////max of 0 or (diff_i_v[r]-diff_i_v[c])
            if ((diff_i_v[r]-diff_i_v[c]) > 0) {
              trans_mat[r][c] = can_goN[r][c]*(diff_i_v[r]-diff_i_v[c]);
            } else {
              trans_mat[r][c] = 0; //looks good all 0 & numbers under .001;
            }
            // //             // Rcpp::Rcout << "initial trans_mat is" << trans_mat;
          } }
        // //
        // // //////////                ADJUST THE TRANSITION MATRIX                  ////////
        // // //////////   1ST SCALE UP RATES, 2ND MAKE SURE DOES NOT SUM OVER 1    //////////
        frc = 0.1;  // approach seems quite sensitive to this value, = fraction of change to
        mat_sum=0;
        ////is this correct, idk
        for(int i=0; i<16; i++){
          for(int j=0; j<16; j++){
            trans_mat[i][j] =  trans_mat[i][j] / dist_i_v[i]*frc;
            mat_sum += trans_mat[i][j]; //should this be calculated in separate loop?
            if (mat_sum > 1.0){ //max of 1 and sum(trans_mat)
              trans_mat[i][j] =  trans_mat [i][j] / mat_sum;
            }
          } }
        // // //       Rcpp::Rcout << "temp trans_mat is" << trans_mat;
        // // // //////////                      FINALIZE TRANS_MAT                    //////////
        for(int i=0; i<16; i++){
          rowsum[i]=0; //reset row sum
          for(int j=0; j<16; j++){
            rowsum[i]+=trans_mat[i][j]; //calculate the row sum of temp mat above
            if (i==j){
              trans_mat[i][j]=(1-rowsum[i]); //rowsum is equal to zero?
            }  } }
        // // // //////////                RECORD ABSOLUTE TRANSITIONS                 //////////
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      did_goN[i][j]=0;
    } }
        for(int i=0; i<16; i++){
          for(int j=0; j<16; j++){
            did_goN[i][j] += dist_i_v[i]* trans_mat[i][j]; } }
        for(int i=0; i<16; i++){
          for(int j=0; j<16; j++){
            if (i==j){
              did_goN[i][j]=0;
            }
          } }
  //       // //////////               UPDATE THE DISTRIBUTION VECTOR             ////////////
  //       // //////////           This is supposed to be matrix multiplication   ////////////
        for(int c=0; c<16; c++){
          temp_vec[c] = 0;
        }
        for(int r=0; r<16; r++){
          for(int c=0; c<16; c++){
            temp_vec[c] += dist_i_v[r]*trans_mat[r][c];
          } } //looks good after one iteration; explodes after 30
        for(int c=0; c<16; c++){
          dist_i_v[c] = temp_vec[c];
        }
      } //end of N loop
  //     // //////////                    NOW UPDATE IN ONE STEP                 ///////////
      for(int i=0; i<16; i++){
        for(int j=0; j<16; j++){
          trans_mat_tot[i][j] = did_goN[i][j]; }}
      for(int i=0; i<16; i++){
        for(int j=0; j<16; j++){
          trans_mat_tot[i][j] = did_goN[i][j] / dist_orig_v[i];
        } }
  //     //
      for(int i=0; i<16; i++){
        rowsum[i]=0;
        for(int j=0; j<16; j++){
          rowsum[i] +=trans_mat_tot[i][j]; //rowsum
          if (i==j){
            trans_mat_tot[i][j]=(1-rowsum[i]);
          } else {
            trans_mat_tot[i][j]=trans_mat_tot[i][j];
          }
        } }
  //     //
  //     // //      Rcpp::Rcout<< "trans_mat_tot is" << trans_mat_tot;
  //     // //////////           NOW FINALLY UPDATE THE DISTRIBUTION           ///////////
      for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++){
          for (int im=0; im<4; im++){
            for (int nm=0; nm<4; nm++){
              for(int rg=0; rg<2; rg++) {
                V2[ag][tb][lt][im][nm][rg][na]=0;
              } } } } }

  for(int tb=0; tb<6; tb++) {
    for(int lt=0; lt<2; lt++){
      for (int im=0; im<4; im++){
        for (int nm=0; nm<4; nm++){
          for(int rg=0; rg<2; rg++) {
            for (int m2=0; m2<4; m2++){
              for (int p2=0; p2<4; p2++){
                //         ////scalar multiplication (not matrix multiplication)
                //      dist_newN[m][p] +=  dist_orig[m2][p2] * trans_mat_tot[m2+p2*4][m+p*4];////removed +1 index
                V2[ag][tb][lt][im][nm][rg][na] += V1[ag][tb][lt][im][nm][rg][na] * trans_mat_tot[m2+p2*4][nm+im*4];
              } } } } } } }
  } } //end of age & nativity loops
    } //end of rebalancing loop

if (reblnc==2){
  mat_sum=0;
  for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
      dist_orig[i][j]=0;
      temp_mat[i][j] =0;
    } }
  // dist_orig_v=0;
  // dist_i_v=0;

  ////// first calculate the total population at this time step

  ////// second calculate the distribution of population across the two risk factors
  for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
      for(int lt=0; lt<2; lt++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++) {
                // pop_t             += V1[ag][tb][lt][im][nm][rg][na];
                dist_orig[nm][im]  += V1[ag][tb][lt][im][nm][rg][na];
              }}}}}}}
  for(int im=0; im<4; im++) {
    for(int nm=0; nm<4; nm++) {
      mat_sum+=dist_orig[nm][im] ; } }
  // ///insert error that distribution does not sum to one and then go from there;
  for(int im=0; im<4; im++) {
    for(int nm=0; nm<4; nm++) {
      dist_orig[nm][im]  = dist_orig[nm][im]/mat_sum; // determine the proportions
      dist_orig_v[(nm)+(im*4)] = dist_orig[nm][im]; //removed +1 //still nan's
      //               // if (std::isnan(dist_orig[nm][im]) > 0  ){
      //              // Rcpp::Rcout << "@ nm " << nm <<"& im "<< im << "orig is nan @" << s << "\n";}
      //              //  if (std::isnan(dist_orig_v[(nm)+(im*4)]) > 0){
      //              // Rcpp::Rcout << "@ nm " << nm <<"& im "<< im << "orig is " << dist_orig_v[(nm)+(im*4)]<< "\n";}
    } }
  /////check that the distributions sum to 1;
  // for (int i=0; i<4; i++){
  //   for (int j=0; j<4; j++){
  //
  //     temp += dist_orig[i][j];
  //     Rcpp::Rcout <<"sum of dist_orig is" << temp << "\n";
  //
  //} }

  for (int i=0; i<16; i++){
    dist_i_v[i] = dist_orig_v[i];
  }

  for(int n=0; n<30; n++){
    // /////// CALCULATE DISTANCE FROM CURRENT DISTRIBUTION TO GOAL DISTRIBUTION /////
    for (int i=0; i<16; i++){
      diff_i_v[i] = dist_i_v[i] - dist_goal_v[i];
    }

    // for (int i=0; i<16; i++){
    //   Rcpp::Rcout <<"at s=" << s << "diff_i_v is" <<  diff_i_v[i] << "\n";
    // }
    // //////////                  CREATE TRANSITION MATRIX                    ////////
    for (int r=0; r<16; r++){
      for (int c=0; c<16; c++){
        trans_mat[r][c] = 0;
      } }
    for (int r=0; r<16; r++){
      for (int c=0; c<16; c++){
        /////max of 0 or (diff_i_v[r]-diff_i_v[c])
        if ((diff_i_v[r]-diff_i_v[c]) > 0) {
          trans_mat[r][c] = can_goN[r][c]*(diff_i_v[r]-diff_i_v[c]);
        } else {
          trans_mat[r][c] = can_goN[r][c]*0; //looks good all 0 & numbers under .001;
        }
        // //             // Rcpp::Rcout << "initial trans_mat is" << trans_mat;
      } }
    //
    // //////////                ADJUST THE TRANSITION MATRIX                  ////////
    // //////////   1ST SCALE UP RATES, 2ND MAKE SURE DOES NOT SUM OVER 1    //////////
    frc = 0.1;  // approach seems quite sensitive to this value, = fraction of change to
    mat_sum=0;
    ////is this correct, idk
    for(int i=0; i<16; i++){
      for(int j=0; j<16; j++){
        trans_mat[i][j] =  trans_mat[i][j] / dist_i_v[i]*frc;
        mat_sum += trans_mat[i][j]; //should this be calculated in separate loop?
        if (mat_sum > 1.0){ //max of 1 and sum(trans_mat)
          trans_mat[i][j] =  trans_mat [i][j] / mat_sum;
        } else {
          trans_mat[i][j]=trans_mat[i][j]/1;
        }
      } }
    // // //       Rcpp::Rcout << "temp trans_mat is" << trans_mat;
    // // // //////////                      FINALIZE TRANS_MAT                    //////////
    for(int i=0; i<16; i++){
      rowsum[i]=0; //reset row sum
      for(int j=0; j<16; j++){
        rowsum[i]+=trans_mat[i][j]; //calculate the row sum of temp mat above
        if (i==j){
          trans_mat[i][j]=(1-rowsum[i]); //rowsum is equal to zero?
        } else {
          trans_mat[i][j]=trans_mat[i][j];
        } } }
    // // //////////                RECORD ABSOLUTE TRANSITIONS                 //////////
    for(int i=0; i<16; i++){
      for(int j=0; j<16; j++){
        did_goN[i][j]=0;
      } }
    for(int i=0; i<16; i++){
      for(int j=0; j<16; j++){
        did_goN[i][j] += dist_i_v[i]* trans_mat[i][j]; } }
    for(int i=0; i<16; i++){
      for(int j=0; j<16; j++){
        if (i==j){
          did_goN[i][j]=0;
        } else {
          did_goN[i][j]=did_goN[i][j];
        }
      } }
    // //////////               UPDATE THE DISTRIBUTION VECTOR             ////////////
    // //////////           This is supposed to be matrix multiplication   ////////////
    for(int c=0; c<16; c++){
      temp_vec[c] = 0;
    }
    for(int r=0; r<16; r++){
      for(int c=0; c<16; c++){
        temp_vec[c] += dist_i_v[r]*trans_mat[r][c];
      } } //looks good after one iteration; explodes after 30
    for(int c=0; c<16; c++){
      dist_i_v[c] = temp_vec[c];
    }
  } //end of N loop
  // //////////                    NOW UPDATE IN ONE STEP                 ///////////
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      trans_mat_tot[i][j] = did_goN[i][j]; }}
  for(int i=0; i<16; i++){
    for(int j=0; j<16; j++){
      trans_mat_tot[i][j] = did_goN[i][j] / dist_orig_v[i];
    } }
  //
  for(int i=0; i<16; i++){
    rowsum[i]=0;
    for(int j=0; j<16; j++){
      rowsum[i] +=trans_mat_tot[i][j]; //rowsum
      if (i==j){
        trans_mat_tot[i][j]=(1-rowsum[i]);
      } else {
        trans_mat_tot[i][j]=trans_mat_tot[i][j];
      }
    } }
  //
  // //      Rcpp::Rcout<< "trans_mat_tot is" << trans_mat_tot;
  // //////////           NOW FINALLY UPDATE THE DISTRIBUTION           ///////////
  for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
      for(int lt=0; lt<2; lt++){
        for (int im=0; im<4; im++){
          for (int nm=0; nm<4; nm++){
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++){
                V2[ag][tb][lt][im][nm][rg][na]=0;
              } } } } } } }

  for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
      for(int lt=0; lt<2; lt++){
        for (int im=0; im<4; im++){
          for (int nm=0; nm<4; nm++){
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++){
                for (int m2=0; m2<4; m2++){
                  for (int p2=0; p2<4; p2++){
                    //         ////scalar multiplication (not matrix multiplication)
                    //      dist_newN[m][p] +=  dist_orig[m2][p2] * trans_mat_tot[m2+p2*4][m+p*4];////removed +1 index
                    V2[ag][tb][lt][im][nm][rg][na] += V1[ag][tb][lt][im][nm][rg][na] * trans_mat_tot[m2+p2*4][nm+im*4];
                  } } } } } } } } }
} //end of rebalancing loop w/o ag or nativity
    ///////////////////////////////////////////////////////////////////////////////
    /////////////////////////    FILL RESULTS TABLE    ////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////
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
                    Outputs[y][27+rg] += V1[ag][tb][lt][im][nm][rg][na];   // N_ by rg (2)
                    Outputs[y][29+na] += V1[ag][tb][lt][im][nm][rg][na];   // N_ by na (3)
                  } } } } } } }
       // Rcpp::Rcout << "total population is" << Outputs[y][1]<< "\n";

      ////////////////////    COUNTS BY NATIVITY AND AGE    ////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  Outputs[y][32+ag] += V1[ag][tb][lt][im][nm][rg][0]; // N_ by age and US (11)
                  Outputs[y][43+ag] += V1[ag][tb][lt][im][nm][rg][1]+V1[ag][tb][lt][im][nm][rg][2];   // N_ by age and FB (11)
                } } } } } }
      /////////////////    LATENT COUNTS BY NATIVITY, AGE, & LTBI    //////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++) {
          for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
              for(int rg=0; rg<2; rg++) {
                Outputs[y][54+ag] += (V1[ag][1][lt][im][nm][rg][0])*(1-pImmScen)+
                  V1[ag][2][lt][im][nm][rg][0]+V1[ag][3][lt][im][nm][rg][0];   // N_ by age and US (11) LATENT INFECTION

                Outputs[y][65+ag] += (V1[ag][1][lt][im][nm][rg][1]+V1[ag][1][lt][im][nm][rg][2])*(1-pImmScen)+
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
                    Outputs[y][76+ag] += V1[ag][tb][lt][im][nm][rg][na];   // N_RF by age (11)
                  } } } } } } }
      ///////////// TB MORTALITY COUNT BY AGE, RISK FACTOR OF INTEREST///////////////
      for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++) {
          for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                  if(im>0) {
                    ti = 11;
                  } else { ti = 0; }
                  if(im==3) {
                    temp = muTbRF;
                  } else { temp = 0;
                  }
                  Outputs[y][87+ag+ti]  += V0[ag][4 ][lt][im][nm][rg][na]*(vTMortN[ag][4 ]+temp);
                  Outputs[y][87+ag+ti]  += V0[ag][5 ][lt][im][nm][rg][na]*(vTMortN[ag][5 ]+temp)*pow(1.0-TxVecZ[1],TunTxMort); //last term is approx .6
                } } } } } }
      ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
      for(int i=87; i<109; i++) { Outputs[y][i] = Outputs[y][i]*12; }

      ///////////////////////  RISK FACTOR MORTALITY BY AGE /////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    Outputs[y][109+ag]  += V0[ag][tb][lt][im][nm][rg][na]*vRFMortN[ag][nm];
                  } } } } } } }
      ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
      for(int i=109; i<120; i++) { Outputs[y][i] = Outputs[y][i]*12;  }
      ///////////////////////    TOTAL MORTALITY BY AGE    /////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    Outputs[y][120+ag]  += VMort[ag][tb][lt][im][nm][rg][na];
                  } } } } } } }
      ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
      for(int i=120; i<131; i++) { Outputs[y][i] = Outputs[y][i]*12; }
      ///////////////////////     TB TREATMENT OUTCOMES    /////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++) {
          for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                  //if(rg!=1) { temp = rDeft[s]; } else { temp = rDeftH[s]; }
                  //////updated tx mat z 22 to tx mat z 5
                  Outputs[y][131]  += V0[ag][5 ][lt][im][nm][rg][na]*TxVecZ[5]; // tx completion
                  Outputs[y][132]  += V0[ag][5 ][lt][im][nm][rg][na]*rDeft[s]; // tx discontinuation
                  Outputs[y][133]  += VMort[ag][5 ][lt][im][nm][rg][na];// tx mort
                } } } } } }
      ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
      for(int i=131; i<134; i++) { Outputs[y][i] = Outputs[y][i]*12; }

      // NOTIFICATIONS
      for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++) {
          for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                  Outputs[y][134   ] += Vdx[ag][4 ][lt][im][nm][rg][na];   // All dx (1)
                  Outputs[y][135+ag] += Vdx[ag][4 ][lt][im][nm][rg][na];   // dx by age (11)
                  Outputs[y][146+na] += Vdx[ag][4 ][lt][im][nm][rg][na];   // dx by na (3)
                  // if(im>0) {
                  //   Outputs[y][148   ] += Vdx[ag][4 ][lt][im][nm][rg][na];   // dx HIV pos (1)
                  // } else {
                  //   Outputs[y][149   ] += Vdx[ag][4 ][lt][im][nm][rg][na]; }  // dx HIV neg (1)

                  Outputs[y][149+rg] += Vdx[ag][4 ][lt][im][nm][rg][na];   // N_ by rg (2)
                } } } } } }
      for(int i=134; i<151; i++) { Outputs[y][i] = Outputs[y][i]*12; }
      /// TLTBI INITS ///
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++) {
              if(na==0 & rg==0 ) { rTbP = rLtScrt[s]*LtDxParN[0][0];
                rTbN = rLtScrt[s]*LtDxParN[0][1]; }
              if(na==0 & rg==1 ) { rTbP = rLtScrt[s]*LtDxParN[1][0];
                rTbN = rLtScrt[s]*LtDxParN[1][1]; }
              if(na > 0) { rTbP = rLtScrt[s]*LtDxParN[3][0];
                rTbN = rLtScrt[s]*LtDxParN[3][1]; }
              for(int ag=0; ag<11; ag++) {
                Outputs[y][151] += (V0[ag][3 ][0 ][im][nm][rg][na]+V0[ag][2 ][0 ][im][nm][rg][na])*rTbP +
                  (V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*rTbN; //all init
                if(na>0) {
                  Outputs[y][152] += (V0[ag][3 ][0 ][im][nm][rg][na]+V0[ag][2 ][0 ][im][nm][rg][na])*rTbP +
                    (V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*rTbN; } // FB inits
                if(rg==1) {
                  Outputs[y][153] +=  (V0[ag][3 ][0 ][im][nm][rg][na]+V0[ag][2 ][0 ][im][nm][rg][na])*rTbP +
                    (V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*rTbN; } // high risk inits
                // if(im>0) {
                //   Outputs[y][156] += (V0[ag][3 ][0 ][im][nm][rg][na]+V0[ag][2 ][0 ][im][nm][rg][na])*rTbP +
                //     (V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*rTbN; } // RF inits

                Outputs[y][154] += (V0[ag][3 ][0 ][im][nm][rg][na]+V0[ag][2 ][0 ][im][nm][rg][na])*rTbP; // inits with LTBI
              } } } } }
      for(int i=149; i<155; i++) { Outputs[y][i] = Outputs[y][i]*12; } // annualize

      /// TB INCIDENCE, BY ALL VS RECENT  ///
      // By recency (<2 years) == all immediate, 1-(1-rfast)^24 x all Lf
      // Break down from latent fast and slow
      temp3 = (1-pow(1-rfast,24.0))*rfast;

      for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++) {
          for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                  temp2 = V0[ag][2 ][lt][im][nm][rg][na]*temp4V[ag][im]  + V0[ag][3 ][lt][im][nm][rg][na]*temp3;// Progression from recent infection
                  temp =  V0[ag][2 ][lt][im][nm][rg][na]*MrslowN[ag][im] + V0[ag][3 ][lt][im][nm][rg][na]*rfast; // All progression
                  Outputs[y][155      ] += temp ;   // all incidence
                  Outputs[y][155+16   ] += temp2;   // all incidence, recent infection
                  Outputs[y][156+ag   ] += temp ;   // incidence by age
                  Outputs[y][156+ag+16] += temp2;   // incidence by age, recent infection
                  if(na<1) {
                    Outputs[y][167      ] += temp ;   //  incidence, US born
                    Outputs[y][167+16   ] += temp2;   //  incidence, US born, recent infection
                  } else {
                    Outputs[y][168      ] += temp ;   //  incidence, FB
                    Outputs[y][168+16   ] += temp2;
                  }//  incidence, FB, recent infection
                  if(na==2) {
                    Outputs[y][169      ] += temp ;   //  incidence, FB2 born
                    Outputs[y][169+16   ] += temp2;    //  incidence, FB2 born, recent infection
                  }
                  if(rg==1) {
                    Outputs[y][170      ] += temp ;   //  incidence, HR
                    Outputs[y][170+16   ] += temp2;
                  } //  incidence, HR, recent infection
                  // if(im>0) {
                  //   Outputs[y][175      ] += temp ;   //  incidence, HIV pos
                  //   Outputs[y][175+17   ] += temp2; } //  incidence, HIV pos, recent infection
                } } } } } }
      for(int i=154; i<187; i++) { Outputs[y][i] = Outputs[y][i]*12; } // annualize

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
                  Outputs[y][187   ] += temp2;   // All dx (1)
                  Outputs[y][188+ag] += temp2;   // dx by age (11)
                  Outputs[y][199+na] += temp2;
                  // if(im>0) {
                  //   Outputs[y][228   ] += temp2;   // dx HIV pos (1)
                  // } else {
                  //   Outputs[y][229   ] += temp2; }  // dx HIV neg (1)
                  Outputs[y][202+rg] += temp2;   // N_ by rg (2)
                } } } } } }
      for(int i=187; i<204; i++) { Outputs[y][i] = Outputs[y][i]*12; }

      // NOTIFICATIONS US
      for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++){
          for (int im=0; im<4; im++){
            for (int nm=0; nm<4; nm++){
              for(int rg=0; rg<2; rg++) {
                Outputs[y][204+ag] += Vdx[ag][4 ][lt][im][nm][rg][0];   // dx by age (11)
              } } } } }
      for(int i=204; i<215; i++) { Outputs[y][i] = Outputs[y][i]*12; }
      // NOTIFICATIONS US, dead at diagnosis
      for(int nm=0; nm<4 ; nm++) {
        for (int im=0; im<4; im++){
          if(im > 2) { temp = muTbRF;  } else { temp = 0;  }
          for(int ag=0; ag<11; ag++) {
            for(int lt=0; lt<2; lt++){
              for(int rg=0; rg<2; rg++) {
                temp2 = V0[ag][4 ][lt][im][nm][rg][0]*(vTMortN[ag][4]+temp); //vTMort[ag][tb]
                Outputs[y][215+ag] += temp2;   // dx by age (11)
              } } } } }
      for(int i=215; i<226; i++) { Outputs[y][i] = Outputs[y][i]*12; }

      // // NOTIFICATIONS, dead at diagnosis  HIV_NEGATIVE
      // for(int ag=0; ag<11; ag++) {
      //   for(int lt=0; lt<2; lt++){
      //     for(int rg=0; rg<2; rg++) {
      //       for(int na=0; na<3; na++){
      //         temp2 = V0[ag][4 ][lt][0][0][rg][na]*(vTMortN[ag][4]+temp);
      //         Outputs[y][254+ag] += temp2;   // dx by age (11)
      //       } } } }


      ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
      // for(int i=229; i<240; i++) { Outputs[y][i] = Outputs[y][i]*12; }
      // // TOTAL MORTALITY BY AGE, HAVE HIV
      // for(int ag=0; ag<11; ag++) {
      //   for(int tb=0; tb<6; tb++) {
      //     for(int lt=0; lt<2; lt++){
      //       for (int im=0; im<4; im++){
      //         for (int nm=0; nm<4; nm++){
      //           for(int rg=0; rg<2; rg++) {
      //             for(int na=0; na<3; na++){
      //               Outputs[y][267+ag]  += VMort[ag][tb][lt][im][nm][rg][na];
      //             } } } } } } }

      ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
      // for(int i=241; i<252; i++) { Outputs[y][i] = Outputs[y][i]*12; }

      /////////////////////  TOTAL MORTALITY BY AGE, HAVE TB   /////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++){
          for (int im=0; im<4; im++){
            for (int nm=0; nm<4; nm++){
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++){
                  Outputs[y][226+ag]  += VMort[ag][4][lt][im][nm][rg][na];
                } } } } } }
      ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
      for(int i=226; i<237; i++) { Outputs[y][i] = Outputs[y][i]*12; }

      // COUNTS TB BY US/FB
      for(int ag=0; ag<11; ag++){
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++){
            for (int im=0; im<4; im++){
              for (int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++){
                    if(na<1) {
                      Outputs[y][237] += V1[ag][2][lt][im][nm][rg][na];
                      Outputs[y][238] += V1[ag][3][lt][im][nm][rg][na];
                      Outputs[y][239] += V1[ag][4][lt][im][nm][rg][na];
                    } else {
                      Outputs[y][240] += V1[ag][2][lt][im][nm][rg][na];
                      Outputs[y][241] += V1[ag][3][lt][im][nm][rg][na];
                      Outputs[y][242] += V1[ag][4][lt][im][nm][rg][na];
                    } } } } } } } }
      // Force of infection
      Outputs[y][243] += VLjkl[0][0];
      Outputs[y][244] += VLjkl[1][0];
      Outputs[y][245] += VLjkl[0][1];
      Outputs[y][246] += VLjkl[1][1];
      ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
      for(int i=243; i<247; i++) { Outputs[y][i] = Outputs[y][i]*12; }

      ///  NEW INFECTIONS + SUPER INFECTIOPN ///
      for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2 ; lt++) {
          for(int im=0; im<4 ; im++) {
            for(int nm=0; nm<4 ; nm++) {
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++){
                  Outputs[y][247+rg] += (V0[ag][0][lt][im][nm][rg][na]+V0[ag][1][lt][im][nm][rg][na]+V0[ag][2][lt][im][nm][rg][na])*VLjkl[rg][na]*NixTrans[s];
                  Outputs[y][249+na] += (V0[ag][0][lt][im][nm][rg][na]+V0[ag][1][lt][im][nm][rg][na]+V0[ag][2][lt][im][nm][rg][na])*VLjkl[rg][na]*NixTrans[s];
                } } } } } }
      ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
      for(int i=245; i<252; i++) { Outputs[y][i] = Outputs[y][i]*12; }

      ////////////////////// TB MORTALITY BY NATIVITY //////////////////////////////
      for(int ag=0; ag<11; ag++){
        for(int lt=0; lt<2; lt++){
          for (int im=0; im<4; im++){
            for (int nm=0; nm<4; nm++){
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++){
                  if(na>1) {
                    ti = 1;
                  } else { ti = 0; }
                  if(im > 2) { temp = muTbRF;
                  } else { temp = 0; }
                  Outputs[y][252+ti]  += V0[ag][4][lt][im][nm][rg][na]*(vTMortN[ag][4 ]+temp);
                  Outputs[y][252+ti]  += V0[ag][5][lt][im][nm][rg][na]*(vTMortN[ag][5 ]+temp)*pow(1.0-TxVecZ[1],TunTxMort); //check the TxMatZ portion
                } } } } } }
      ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
      for(int i=252; i<254; i++) { Outputs[y][i] = Outputs[y][i]*12; }
      ///////////////////////    TOTAL MORTALITY BY AGE    /////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    if (na<1){
                      Outputs[y][254+ag]  += VMort[ag][tb][lt][im][nm][rg][na];
                      Outputs[y][276+im]  += VMort[ag][tb][lt][im][nm][rg][na];
                      Outputs[y][284+nm]  += VMort[ag][tb][lt][im][nm][rg][na];
                      Outputs[y][292+rg]  += VMort[ag][tb][lt][im][nm][rg][na];
                    } else {
                      Outputs[y][265+ag]  += VMort[ag][tb][lt][im][nm][rg][na];
                      Outputs[y][280+im]  += VMort[ag][tb][lt][im][nm][rg][na];
                      Outputs[y][288+nm]  += VMort[ag][tb][lt][im][nm][rg][na];
                      Outputs[y][294+rg]  += VMort[ag][tb][lt][im][nm][rg][na];
                    }
                  } } } } } } }
      ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
      for(int i=254; i<296; i++) { Outputs[y][i] = Outputs[y][i]*12; }
///////////////////////  POPULATION  /////////////////////////

      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    if (na<1){
                      Outputs[y][296+im]  += V1[ag][tb][lt][im][nm][rg][na];
                      Outputs[y][304+rg]  += V1[ag][tb][lt][im][nm][rg][na];
                      Outputs[y][308+nm]  += V1[ag][tb][lt][im][nm][rg][na];
                    } else {
                      Outputs[y][300+im]  += V1[ag][tb][lt][im][nm][rg][na];
                      Outputs[y][306+rg]  += V1[ag][tb][lt][im][nm][rg][na];
                      Outputs[y][312+nm]  += V1[ag][tb][lt][im][nm][rg][na];

                    }

                  } } } } } } }

      for(int i=296; i<316; i++) { Outputs[y][i] = Outputs[y][i]*12; }


      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    Outputs[y][316]  += VMort[ag][tb][lt][im][nm][rg][na];
                  } } } } } } }

      for(int i=316; i<317; i++) { Outputs[y][i] = Outputs[y][i]*12; }

    } ////end of mid-year results bracket
    ///////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////END MIDYEAR RESULTS//////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////
    // for(int i=0; i<16; i++){
    //   for(int j=0; j<16; j++){
    //     temp_mat[i][j] = dist_goal[i][j];
    //     sum1=
    //     sse = arma::accu(pow(dist_goal[i][j] - dist_new[i][j], 2)) /
    //       arma::accu(pow(dist_goal[i][j] - dist_orig[i][j], 2));
    //   }}
    // Rcpp::Rcout << "sse is" << sse;

    //////REVERT TO A FORM THAT CAN BE EXPORTED TO R
    // for (int im=0; im<4; im++){
    //   for (int nm=0; nm<4; nm++){
    //     dist_new_fin(nm,im) = dist_newN[nm][im];
    //   } }
    //
    // for (int i=0; i<528; i++){
    //   for (int j=0; j<528; j++){
    //     trans_mat_fin(i,j) = trans_mat[i][j];
    //   } }
    // //
    // for (int i=0; i<528; i++){
    //   dist_i_v_fin(i) = dist_i_v[i];
    // }
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
                  if (reblnc > 0){
                    V0[ag][tb][lt][im][nm][rg][na] = V2[ag][tb][lt][im][nm][rg][na];
                  } else {
                    V0[ag][tb][lt][im][nm][rg][na] = V1[ag][tb][lt][im][nm][rg][na];
                  }
                } } } } } } }

    } //// end of month loop!//////////////////////////////////////////////////////////
  } //// end of year loop!///////////////////////////////////////////////////////////
  ////////////////   RE-CHECK THAT NO STATES HAVE BECOME NEGATIVE    ////////////////
  Rcpp::NumericVector  CheckV(12672);
  for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
      for(int lt=0; lt<2; lt++){
        for(int im=0; im<4; im++){
          for(int nm=0; nm<4; nm++){
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++) {
                CheckV(ag+tb*11+lt*66+im*132+nm*528+rg*2112+na*4224) = V1[ag][tb][lt][im][nm][rg][na];
              } } } } } } }
  ///////////////////////////////////////////////////////////////////////////////////
  //////                              RETURN STUFF                              /////
  ///////////////////////////////////////////////////////////////////////////////////
  for(int i=0; i<nYrs; i++) {
    for(int j=0; j<nRes; j++) {
      Outputs2(i,j)  = Outputs[i][j];
    }  }
  return
    Rcpp::List::create(
      Rcpp::Named("Outputs") = Outputs2,
      Rcpp::Named("dist_i_v") = dist_i_v_fin,
      Rcpp::Named("trans_mat") = trans_mat_fin,
      Rcpp::Named("matsum") = mat_sum,
      Rcpp::Named("pop_t") = pop_t,
      Rcpp::Named("V0") = CheckV0,
      Rcpp::Named("V1") = CheckV
    );

}
