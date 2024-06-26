#include <Rcpp.h>
#include <math.h>
#include <cmath>
using namespace Rcpp;
//'@name national_cSim
//'@description runs a simulation of the tb model
//'@param nYrs number of years to run the model.
//'@param nRes number of results of the model
//'@param rDxt Rate of active TB diagnosis over time
//'@param TxQualt active TB treatment over time
//'@param InitPop Initial Population matrix
//'@param Mpfast Matrix of the probabilities of fast TB progression (age x TB prog risk group)
//'@param ExogInf Exogenous infection risk for non-US born population
//'@param MpfastPI Matrix of the probabilities of fast TB progression w/ partial immunity (age x TB prog risk group)
//'@param Mrslow matrix of the rates of slow progression (age x TB prog risk group)
//'@param rrSlowFB rate of fast TB progression
//'@param RRcurDef Rate Ratio for cure given treatment defaul
//'@param rSlfCur rate of self cure from active TB
//'@param p_HR probability of high risk population @ entry into the model
//'@param vTMort vector of TB mortality rates
//'@param RRmuRF rate ratio of mortality across mortality risk group
//'@param RRmuHR rate ratio of mortality across low/high risk dimension
//'@param RRmuTBPand rate ratio of TB mortality during pandemic
//'@param Birthst Births over time
//'@param HrEntEx Matrix of Entry and Exit rates into the High Risk population
//'@param ImmNon Immigration with no TB
//'@param ImmLat Immigration with Latent TB
//'@param ImmAct Immigration with Active TB
//'@param ImmFst Immigration with Fast Progressing TB
//'@param Int1Test Additional tests for Int1
//'@param Int1Init Additional treatment inits for Int1
//'@param Int1Tx Additional treatment completions for Int1
//'@param SpImmNon Immigration with cured TB -- Scen 4 only
//'@param net_mig_usb net internal migration usb
//'@param net_mig_nusb net internal migration nusb
//'@param mubt background mortality over time
//'@param RelInf beta
//'@param RelInfRg beta based off of risk group
//'@param RRcrAG rate ratio for contact rate by ag
//'@param Vmix 1-sigma
//'@param rEmmigFB rate of emmigration in non-US born population
//'@param TxVec vector of parameters for TB Tx
//'@param TunTxMort Tuning parameter for mortality on TB Tx
//'@param rDeft rate of default from TB treatment over time
//'@param rLtScrt rate of latent screening over time
//'@param ttt_samp_dist probabilities of screening for ttt intervention
//'@param ttt_ltbi how much to increase ltbi vector
//'@param ttt_ltbi_init counterfactual tltbi initiation
//'@param ttt_ltbi_comp counterfactual tltbi completion
//'@param ttt_ltbi_eff counterfactual tltbi effectiveness
//'@param ttt_ltbi_acceot counterfactual ltbi testing
//'@param ttt_ltbi_sens  counterfactual tltbi sensitivity
//'@param ttt_ltbi_spec counterfactual tltbi specificity
//'@param ttt_month when to apply ttt interventions
//'@param ttt_pop_scrn population size to apply the ltbi prev to
//'@param LtDxPar_lt matrix of latent diagnosis parameters
//'@param LtDxPar_nolt matrix of latent diagnosis parameters
//'@param LtTxPar matrix of latent treatment parameters
//'@param RRdxAge vector of rate ratios for TB diagnosis by age
//'@param rRecov rate of recovery from latent slow to safe tb state
//'@param pImmScen lack of reactivitiy to IGRA for Sp
//'@param EarlyTrend ramp down of TB in burn-in
//'@param pReTx probability of re-treatment for TB
//'@param ag_den denominator used in the aging process
//'@param NixTrans reduction of transmission over time
//'@param NixTb reduction of all TB over time
//'@param dist_gen general distribution across tb progression and mort
//'@param trans_mat_tot_ages
//'@return Outputs a list of outputs
//[[Rcpp::export]]

Rcpp::List national_cSim(
    std::vector<int> setup_pars,
    Rcpp::NumericMatrix rDxt,
    std::vector<double> TxQualt,
    Rcpp::NumericMatrix InitPop,
    Rcpp::NumericMatrix Mpfast,
    std::vector<double> ExogInf,
    Rcpp::NumericMatrix MpfastPI,
    Rcpp::NumericMatrix Mrslow,
    std::vector<double> rrSlowFB,
    double              rfast,
    double              RRcurDef,
    double              rSlfCur,
    double              p_HR,
    Rcpp::NumericMatrix vTMort,
    std::vector<double> RRmuRF,
    std::vector<double> RRmuHR,
    std::vector<double> RRmuTBPand,
    std::vector<double> Birthst,
    Rcpp::NumericMatrix HrEntEx,
    Rcpp::NumericMatrix ImmNon,
    Rcpp::NumericMatrix ImmLat,
    Rcpp::NumericMatrix ImmAct,
    Rcpp::NumericMatrix ImmFst,
    Rcpp::NumericMatrix SpImmNon,
    std::vector<double> net_mig_usb,
    std::vector<double> net_mig_nusb,
    Rcpp::NumericMatrix mubt,
    std::vector<double> RelInf,
    std::vector<double> RelInfRg,
    std::vector<double> RRcrAG,
    std::vector<double> Vmix,
    std::vector<double> rEmmigFB,
    std::vector<double> TxVec,
    double              TunTxMort,
    std::vector<double> rDeft,
    Rcpp::NumericMatrix rLtScrt,
    Rcpp::NumericMatrix ttt_samp_dist,
    std::vector<int>    ttt_month,
    double              ttt_pop_scrn,
    std::vector<double> ttt_ltbi,
    int                 ttt_ltbi_accept,
    int                 ttt_ltbi_init,
    int                 ttt_ltbi_comp,
    int                 ttt_ltbi_eff,
    int                 ttt_ltbi_sens,
    int                 ttt_ltbi_spec,
    Rcpp::NumericMatrix LtTxPar,
    Rcpp::NumericMatrix LtDxPar_lt,
    Rcpp::NumericMatrix LtDxPar_nolt,
    double rrTestLrNoTb,
    double rrTestHr,
    Rcpp::NumericMatrix Int1Test,
    Rcpp::NumericMatrix Int1Init,
    Rcpp::NumericMatrix Int1Tx,
    std::vector<double> RRdxAge,
    double              rRecov,
    double              pImmScen,
    std::vector<double> EarlyTrend,
    std::vector<double> pReTx,
    Rcpp::NumericMatrix ag_den,
    std::vector<double> NixTrans,
    std::vector<double> NixTb,
    Rcpp::NumericMatrix dist_gen,
    Rcpp::NumericMatrix trans_mat_tot_ages
) {
  ////////////////////////////////////////////////////////////////////////////////
  ////////    BELOW IS A LIST OF THE VARIABLES CREATED INTERNALLY IN MODEL   /////
  ////////////////////////////////////////////////////////////////////////////////
  int           ti; int n2;
  int           s;
  int           nYrs = setup_pars[0];
  int           nRes = setup_pars[1];
  int           res_m = setup_pars[2];
  double        InitPopN[InitPop.nrow()][InitPop.ncol()];
  double        InitPopZ[InitPop.nrow()][InitPop.ncol()];
  double        MpfastN[Mpfast.nrow()][Mpfast.ncol()];
  double        MpslowN[Mpfast.nrow()][Mpfast.ncol()];
  double        MpfastPIN[MpfastPI.nrow()][MpfastPI.ncol()];
  double        MpslowPIN[MpfastPI.nrow()][MpfastPI.ncol()];
  double        MrslowN[Mrslow.nrow()][Mrslow.ncol()];
  double        vTMortN[vTMort.nrow()][vTMort.ncol()];
  double        ag_denN[ag_den.nrow()][ag_den.ncol()];
  double        HrEntExN[HrEntEx.nrow()][HrEntEx.ncol()];
  double        ImmNonN[ImmNon.nrow()][ImmNon.ncol()];
  double        SpImmNonN[ImmNon.nrow()][ImmNon.ncol()];
  double        ImmLatN[ImmLat.nrow()][ImmLat.ncol()];
  double        ImmFstN[ImmFst.nrow()][ImmFst.ncol()];
  double        ImmActN[ImmAct.nrow()][ImmAct.ncol()];
  double        Int1TestN[Int1Test.nrow()][Int1Test.ncol()];
  double        Int1InitN[Int1Init.nrow()][Int1Init.ncol()];
  double        Int1TxN[Int1Tx.nrow()][Int1Tx.ncol()];
  double        mubtN[mubt.nrow()][mubt.ncol()];
  double        rDxtN[rDxt.nrow()][rDxt.ncol()];
  double        rLtScrtN[rLtScrt.nrow()][rLtScrt.ncol()];
  double        ttt_samp_distN[ttt_samp_dist.nrow()][ttt_samp_dist.ncol()];
  double        ttt_dist[ttt_samp_dist.nrow()][ttt_samp_dist.ncol()];
  double        LtTxParN[LtTxPar.nrow()][LtTxPar.ncol()];
  double        LtDxPar_ltN[LtDxPar_lt.nrow()][LtDxPar_lt.ncol()];
  double        LtDxPar_noltN[LtDxPar_nolt.nrow()][LtDxPar_nolt.ncol()];
  double        TxVecZ[6];
  double        temp;  double temp2; double temp3;
  double        temp4; double temp5; double temp6; double temp11;
  double        temp7; double temp8; double temp9; double temp10;
  double        indextemp;
  double        temp4V[11][5][2];
  double        temp3V[2];
  double        rTbP;
  double        rTbP_norm;
  double        rTbN;
  double        rTbN_norm;
  double        pop_scrn;
  double        rr_ltbi;
  double        Outputs[nYrs][nRes];
  double  	    V0[11][6][2][4][4][2][3];
  double  	    V1[11][6][2][4][4][2][3];
  double  	    VMort[11][6][2][4][4][2][3];
  double       Vdx[11][6][2][4][4][2][3];
  double       VLdx[11][6][2][4][4][2][3];
  double       VLtest[11][6][2][4][4][2][3];
  double        VLtemp[11][6][2][4][4][2][3];
  double        VLDxtemp[11][6][2][4][4][2][3];
  double       Vinf[11][6][2][4][4][2][3];
  double        VAtemp[11][4][4][2][3];
  double            VNkl[2][2];  ///HIGH AND LOW RISK, NATIVITY, AGE
  double            VGkl[2][2]; ///HIGH AND LOW RISK, NATIVITY,AGE
  double            Vjaf[4];     ///BY NUMBER OF MIXING GROUPS, AGE
  double            VLkla[2][2][11];  ///HIGH AND LOW RISK, NATIVITY,AGE
  int               N;
  double        dist_genN[dist_gen.nrow()][dist_gen.ncol()];
  double        temp_vec[4];
  double        temp_mat[4][4];
  double      	temp_mat2[4][4];
  double 	      trans_mat_tot_agesN[trans_mat_tot_ages.nrow()][trans_mat_tot_ages.ncol()];
  double        ttt_test_susc;
  double        ttt_test_PI;
  double        ttt_test_lf;
  double        ttt_test_ls;
  double        ttt_ltbi_sensN;
  double        ttt_ltbi_specN;
  double        ttt_ltbi_initN;
  double        ttt_ltbi_acceptN;
  double        ttt_ltbi_compN;
  double        ttt_ltbi_effN;
  double        ttt_frc_scrn;
  double        mat_sum;
  double        temp_vec2[4];
  int reblnc;
  int tb_dyn;
  double base_diag;
  Rcpp::NumericMatrix Outputs2(nYrs,nRes);
  Rcpp::NumericMatrix dist_mat(4,4);
  double        RRmuRFN[4];
  double        mort_dist[4];
  double temp20 = 0 ; double temp21 = 0; double temp22 = 0; double temp23 = 0;
  double temp19 = 0 ; double temp24 = 0; double temp25 = 0; double temp26 = 0;
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
  for(int i=0; i<ttt_samp_dist.nrow(); i++) {
    for(int j=0; j<ttt_samp_dist.ncol(); j++) {
      ttt_samp_distN[i][j] = ttt_samp_dist(i,j);
      ttt_dist[i][j] = 0;
    } }
  for(int i=0; i<LtDxPar_lt.nrow(); i++) {
    for(int j=0; j<LtDxPar_lt.ncol(); j++) {
      LtDxPar_ltN[i][j] = LtDxPar_lt(i,j);
      LtDxPar_noltN[i][j] = LtDxPar_nolt(i,j);
    } }
  for(int i=0; i<LtTxPar.nrow(); i++) {
    for(int j=0; j<LtTxPar.ncol(); j++) {
      LtTxParN[i][j] = LtTxPar(i,j);
    } }
  for(int i=0; i<vTMort.nrow(); i++) {
    for(int j=0; j<vTMort.ncol(); j++) {
      vTMortN[i][j] = vTMort(i,j);
    } }
  for(int i=0; i<ag_den.nrow(); i++) {
    for(int j=0; j<ag_den.ncol(); j++) {
      ag_denN[i][j] = ag_den(i,j);
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
  for(int i=0; i<SpImmNon.nrow(); i++) {
    for(int j=0; j<SpImmNon.ncol(); j++) {
      SpImmNonN[i][j] = SpImmNon(i,j);
    } }
  for(int i=0; i<Int1Test.nrow(); i++) {
    for(int j=0; j<Int1Test.ncol(); j++) {
      Int1TestN[i][j] = Int1Test(i,j);
      Int1InitN[i][j] = Int1Init(i,j);
      Int1TxN[i][j] = Int1Tx(i,j);
    } }
  for(int i=0; i<mubt.nrow(); i++) {
    for(int j=0; j<mubt.ncol(); j++) {
      mubtN[i][j] = mubt(i,j);
    } }
  for(int i=0; i<rDxt.nrow(); i++) {
    for(int j=0; j<rDxt.ncol(); j++) {
      rDxtN[i][j] = rDxt(i,j);
    } }
  for(int i=0; i<rLtScrt.nrow(); i++) {
    for(int j=0; j<rLtScrt.ncol(); j++) {
      rLtScrtN[i][j] = rLtScrt(i,j);
    } }
  for(int i=0; i<dist_gen.nrow(); i++) {
    for(int j=0; j<dist_gen.ncol(); j++) {
      dist_genN[i][j] =dist_gen(i,j);
    } }
  for(int i=0; i<HrEntEx.nrow(); i++) {
    for(int j=0; j<HrEntEx.ncol(); j++) {
      HrEntExN[i][j] = HrEntEx(i,j);
    } }
  for(int i=0; i<6; i++) {
    TxVecZ[i] = 0.0;
  }
  for(int i=0; i<nYrs; i++) {
    for(int j=0; j<nRes; j++) {
      Outputs[i][j] = 0;
    }  }
  for(int ag=0; ag<11; ag++) {
    for(int im=0; im<4; im++){
      for(int nm=0; nm<4; nm++){
        for(int rg=0; rg<2; rg++) {
          for(int na=0; na<3; na++){
            for(int lt=0; lt<2; lt++){
              for(int tb=0; tb<6; tb++) {
                V0    [ag][tb][lt][im][nm][rg][na] = 0;
                V1    [ag][tb][lt][im][nm][rg][na] = 0;
                VMort [ag][tb][lt][im][nm][rg][na] = 0;
                Vdx   [ag][tb][lt][im][nm][rg][na] = 0;
                VLdx  [ag][tb][lt][im][nm][rg][na] = 0;
                VLtest[ag][tb][lt][im][nm][rg][na] = 0;
                VLtemp[ag][tb][lt][im][nm][rg][na] = 0;
                VLDxtemp[ag][tb][lt][im][nm][rg][na] = 0;
                Vinf[ag][tb][lt][im][nm][rg][na]  = 0;
            }
            VAtemp[ag][im][nm][rg][na] = 0;
          } } } } }}

  for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
      VNkl[i][j] = 0;
      VGkl[i][j] = 0;
      for(int k=0; k<11; k++) {
        VLkla[i][j][k] = 0;
      } } }
  ///effective contact rates
  for(int i=0; i<4; i++) {
    mort_dist[i]=0;
    temp_vec2[i]=0;
    temp_vec[i]=0;
    RRmuRFN[i]=RRmuRF[i];
  }
  for(int rg=0; rg<2; rg++) {
    for(int ag=0; ag<11; ag++) {
      for(int im=0; im<4; im++) {
        temp4V[ag][im][rg] = (1-pow(1-(MrslowN[ag][im])-rRecov,24.0-rDxtN[0][rg]))*(MrslowN[ag][im]);
      } }
    temp3V[rg]=0;
  }

  for(int j=0; j<4; j++) {
    Vjaf[j] = 0;
  }
  for(int i=0; i<trans_mat_tot_ages.nrow(); i++) {
    for(int j=0; j<trans_mat_tot_ages.ncol(); j++) {
      trans_mat_tot_agesN[i][j] = trans_mat_tot_ages(i,j);
    } }
  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      temp_mat2[i][j] = 0;
      temp_mat[i][j] = 0;
    } }
  N=30; indextemp=0;
  reblnc=1;
  tb_dyn=1;
  temp=0; n2=0; int ni=0;
  temp10=0; temp2=0; temp3=0; temp4=0; temp5=0; temp6=0;
  temp7=0; temp8=0; temp9=0; temp10=0; temp11=0;
  ttt_test_susc=0; ttt_test_PI=0; ttt_test_lf=0; ttt_test_ls=0; base_diag=0;
  rr_ltbi=1;
  pop_scrn=0;
  if (tb_dyn != 1){
    for(int ag=0; ag<11; ag++) {
      for(int tb=0; tb<6; tb++) {
        vTMortN[ag][tb] =0;
      } }
  }
  ttt_ltbi_sensN=0; ttt_ltbi_specN=0; ttt_ltbi_acceptN=0;
  ttt_ltbi_initN=0; ttt_ltbi_compN=0; ttt_ltbi_effN=0;
  rTbN_norm=0; rTbP_norm=0;;
  ///// CREATE A COUNTER FOR THE TTT INTERVENTION; THIS IS INTIALIZED TO ZERO /////
  ///// OUTSIDE OF THE TIME LOOP AND WILL BE INCREMENTED AT THE BEGINNING OF //////
  ///// EACH TIME STEP OF THE TTT INTERVENTION; IT IS USED TO ADJUST THE FRC //////
  ///// SCREENED.
  double ttt_m = 0.0;
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
        ////////////////////////  UNINFECTED/SUSCEPTIBLE POP /////////////////////
        ///// Low risk US born ///////////////////////////////////////////////////
        V0[ag][0][0][im][nm][0][0] = InitPopN[ag][0]*0.40*(1-p_HR)*dist_genN[nm][im];
        ///// High risk US born //////////////////////////////////////////////////
        V0[ag][0][0][im][nm][1][0] = InitPopN[ag][0]*0.40*(p_HR)*dist_genN[nm][im];
        ///// Low risk non-US born ///////////////////////////////////////////////
        V0[ag][0][0][im][nm][0][2] = InitPopN[ag][1]*0.40*(1-p_HR)*dist_genN[nm][im];
        ///// High risk non-US born //////////////////////////////////////////////
        V0[ag][0][0][im][nm][1][2] = InitPopN[ag][1]*0.40*(p_HR)*dist_genN[nm][im];
        ///////////////////////// LATENT SLOW INFECTED POP  ///////////////////////
        ///// Low risk US born ////////////////////////////////////////////////////
        V0[ag][2][0][im][nm][0][0] = InitPopN[ag][0]*0.60*(1-p_HR)*dist_genN[nm][im];
        ///// High risk US born  //////////////////////////////////////////////////
        V0[ag][2][0][im][nm][1][0] = InitPopN[ag][0]*0.60*(p_HR)*dist_genN[nm][im];
        ///// Low risk non-US born ////////////////////////////////////////////////
        V0[ag][2][0][im][nm][0][2] = InitPopN[ag][1]*0.60*(1-p_HR)*dist_genN[nm][im];
        ///// High risk non-US born ///////////////////////////////////////////////
        V0[ag][2][0][im][nm][1][2] = InitPopN[ag][1]*0.60*(p_HR)*dist_genN[nm][im];
      } } }

  ////////////////////////////////////////////////////////////////////////////////
  ////// create a 2nd array with same dimensions as V0 & populate w/ same values//
  ////////////////////////////////////////////////////////////////////////////////
  for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++){
          for(int rg=0; rg<2; rg++){
            for (int na=0; na<3; na++){
              V1[ag][tb][0][im][nm][rg][na]  = V0[ag][tb][0][im][nm][rg][na];
            } } } } } }

  ////////////////////////RUN THE MODEL FOR 1201 MONTHS /////////////////////////
  for(int m=0; m<1201; m++) {
    /////////////////////////////////////////////////////////////////////////////
    /////                            START BURN IN                          /////
    /////////////////////////////////////////////////////////////////////////////
    /////                               BIRTHS                              /////
    /////      USE DISTRIBUTION TO POPULATE THE MODEL ACROSS RISK GROUPS    /////
    /////////////////////////////////////////////////////////////////////////////
    for(int im=0; im<4; im++) {
      for(int nm=0; nm<4; nm++){
        V1[0][0][0][im][nm][0][0]  += Birthst[0]*dist_genN[nm][im]*(1-p_HR);
        V1[0][0][0][im][nm][1][0]  += Birthst[0]*dist_genN[nm][im]*(p_HR);
      } }
    /////////////////////////////////////////////////////////////////////////////
    /////                             IMMIGRATION                           /////
    /////////////////////////////////////////////////////////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++){
          ///// NO TB, Low risk  ////////////////////////////////////////////////
          V1[ag][0][0][im][nm][0][1]   += ImmNonN[0][ag]*dist_genN[nm][im]*(1-p_HR);
          ///// NO TB, High risk  ///////////////////////////////////////////////
          V1[ag][0][0][im][nm][1][1]   += ImmNonN[0][ag]*dist_genN[nm][im]*(p_HR);
          ///// LATENT SLOW TB, Low risk  ////////////////////////////////////////
          V1[ag][2][0][im][nm][0][1]   += ImmLatN[0][ag]*dist_genN[nm][im]*(1-p_HR);
          ///// LATENT SLOW TB, High risk  ///////////////////////////////////////
          V1[ag][2][0][im][nm][1][1]   += ImmLatN[0][ag]*dist_genN[nm][im]*(p_HR);
          ///// LATENT FAST TB, Low risk  ////////////////////////////////////////
          V1[ag][3][0][im][nm][0][1]   += ImmFstN[0][ag]*dist_genN[nm][im]*(1-p_HR);
          ///// LATENT FAST TB, High risk  ////////////////////////////////////////
          V1[ag][3][0][im][nm][1][1]   += ImmFstN[0][ag]*dist_genN[nm][im]*(p_HR);
          ///// ACTIVE TB, Low risk  //////////////////////////////////////////////
          V1[ag][4][0][im][nm][0][1]   += ImmActN[0][ag]*dist_genN[nm][im]*(1-p_HR);
          ///// ACTIVE TB, High risk  //////////////////////////////////////////////
          V1[ag][4][0][im][nm][1][1]   += ImmActN[0][ag]*dist_genN[nm][im]*(p_HR);
        } } }
    /////////////////////////////////////////////////////////////////////////////
    /////                          EMMIGRATION                              /////
    /////////////////////////////////////////////////////////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++){
            for(int rg=0; rg<2; rg++){
              ///// In US under 2 years  ////////////////////////////////////////
              V1[ag][tb][0][im][nm][rg][1]  -= V0[ag][tb][0][im][nm][rg][1]*rEmmigFB[0];
              ///// In US over 2 years  /////////////////////////////////////////
              V1[ag][tb][0][im][nm][rg][2]  -= V0[ag][tb][0][im][nm][rg][2]*rEmmigFB[1];
            } } } } }
    /////////////////////////////////////////////////////////////////////////////
    /////                            MORTALITY                              /////
    /////////////////////////////////////////////////////////////////////////////
    temp=0;
    for(int ag=0; ag<11; ag++) {
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++){
              if ((ag<9) | ((mubtN[0][ag]*RRmuRFN[nm]*RRmuHR[rg]) < .5)){
                temp = ((RRmuRFN[nm]*RRmuHR[rg])*mubtN[0][ag]);
              } else {
                temp =  (.5);
              }
              ///// WITHOUT TB TREATMENT ////////////////////////////////////////
              for(int tb=0; tb<5; tb++) {
                V1[ag][tb][0][im][nm][rg][na]  -= (V0[ag][tb][0][im][nm][rg][na]*(temp+vTMortN[ag][tb]));
              }
              ///// WITH TB TREATMENT ///////////////////////////////////////////
              V1[ag][5 ][0][im][nm][rg][na]  -= V0[ag][5 ][0][im][nm][rg][na]*(temp+vTMortN[ag][5 ]*pow(1.0-TxVecZ[1],TunTxMort));
            } } } } }
    /////////////////////////////////////////////////////////////////////////////
    /////                              AGING                                /////
    /////////////////////////////////////////////////////////////////////////////
    for(int ag=0; ag<10; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++) {
                temp = V0[ag][tb][0][im][nm][rg][na]/ag_denN[0][ag];
                V1[ag  ][tb][0][im][nm][rg][na]  -= temp;
                V1[ag+1][tb][0][im][nm][rg][na]  += temp;
              } } } } } }
    /////////////////////////////////////////////////////////////////////////////
    /////                     NEW FB -> ESTABLISHED FB                      /////
    /////                     TWO YEARS FOR TRANSITION                      /////
    /////////////////////////////////////////////////////////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {
              temp = V0[ag][tb][0][im][nm][rg][1] / 24;
              V1[ag][tb][0][im][nm][rg][1]  -= temp;
              V1[ag][tb][0][im][nm][rg][2]  += temp;
            } } } } }
    /////////////////////////////////////////////////////////////////////////////
    /////                        HIGH-RISK ENTRY/EXIT                       /////
    /////////////////////////////////////////////////////////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int na=0; na<3; na++) {
              ///// ENTRY //////////////////////////////////////////////////////
              temp  = V0[ag][tb][0][im][nm][0][na]*HrEntExN[ag][0];
              ///// EXIT ///////////////////////////////////////////////////////
              temp2 = V0[ag][tb][0][im][nm][1][na]*HrEntExN[ag][1];
              V1[ag][tb][0][im][nm][0][na]  += temp2-temp;
              V1[ag][tb][0][im][nm][1][na]  += temp-temp2;
            } } } } }
    /////////////////////////////////////////////////////////////////////////////
    /////                      OPEN TB DYNAMICS LOOP                        /////
    /////////////////////////////////////////////////////////////////////////////
    if (tb_dyn==1){
      ////////////////////////////  TRANSMISSION RISK  ////////////////////////////////
      for(int i=0; i<2; i++){
        for(int j=0; j<2; j++){
          VNkl[i][j]=0;
          VGkl[i][j]=0;
        }
      }
      // Step 1 & 2
      // take total population of mixing groups
      // mixing groups are only risk group and nativity
      // we do not need to stratify by age here because we are only allowing the number of contacts
      // of infectious persons to vary as a constant across all mixing groups; no assortative mixing  RRcrAG[ag]
      // by age.
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
              ///////// LOW RISK US BORN
              VNkl[0][0] += V0[ag][tb][0][im][nm][0][0]
              * RRcrAG[ag] * RelInfRg[0];
              VGkl[0][0] += V0[ag][tb][0][im][nm][0][0]
              * RRcrAG[ag] * RelInfRg[0] * RelInf[tb];
              ///////// HIGH RISK US BORN
              VNkl[1][0] += V0[ag][tb][0][im][nm][1][0]
              * RRcrAG[ag] * RelInfRg[1];
              VGkl[1][0] += V0[ag][tb][0][im][nm][1][0]
              * RRcrAG[ag] * RelInfRg[1] * RelInf[tb];
              ///////// LOW RISK NON US BORN
              VNkl[0][1] += (V0[ag][tb][0][im][nm][0][1] + V0[ag][tb][0][im][nm][0][2])
                * RRcrAG[ag] * RelInfRg[2];
              VGkl[0][1] += (V0[ag][tb][0][im][nm][0][1] + V0[ag][tb][0][im][nm][0][2])
                * RRcrAG[ag] * RelInfRg[2] * RelInf[tb];
              ///////// HIGH RISK NON US BORN
              VNkl[1][1] += (V0[ag][tb][0][im][nm][1][1] + V0[ag][tb][0][im][nm][1][2])
                * RRcrAG[ag] * RelInfRg[3];
              VGkl[1][1] += (V0[ag][tb][0][im][nm][1][1] + V0[ag][tb][0][im][nm][1][2])
                * RRcrAG[ag] * RelInfRg[3] * RelInf[tb];
            } } } }

      // Step 3 (Infected/Total)
      //probability of infection  -- maybe not stratified by age either?
      //for each mixing group
      // calculated as the a weighted average based on the number of effective contacts
      // each nativity and risk group contributes to that mixing group
      // 0 = common pool; 1 = exclusive nusb, 2 = exclusive hr, 3 = exclusive nusb-hr
      // 0 = common pool; 1 = exclusive hr, 2= exclusive nusb , 3 = exclusive nusb-hr

      //Vmix[0] is HR ; Vmix[1] is NUSB
      Vjaf[0] = (VGkl[0][0]         +
        VGkl[1][0]*Vmix[0] +
        VGkl[0][1]*Vmix[1] +
        VGkl[1][1]*Vmix[1]*Vmix[0]) /
          (VNkl[0][0]         +
            VNkl[1][0]*Vmix[0] +
            VNkl[0][1]*Vmix[1] +
            VNkl[1][1]*Vmix[1]*Vmix[0] + 1e-12);

      Vjaf[1] = (VGkl[0][1] + VGkl[1][1]*Vmix[0]) /
        (VNkl[0][1] + VNkl[1][1]*Vmix[0] + 1e-12);  //exclusive HR

      Vjaf[2] = (VGkl[1][0] + VGkl[1][1]*Vmix[1])/ //exclusive NUSB
        (VNkl[1][0] + VNkl[1][1]*Vmix[1] + 1e-12);

      Vjaf[3] = VGkl[1][1] / (VNkl[1][1] + 1e-12); //exclusive HR NUSB

      // Step 4
      // calculate force of infection
      // the number of contacts you have*probability of infection*mixing
      for(int ag=0; ag<11; ag++) {
        /// LOW RISK US BORN
        VLkla[0 ][0 ][ag]  = RelInfRg[0] * RRcrAG[ag] * Vjaf[0] ;
        ///////// HIGH RISK US BORN
        VLkla[1 ][0 ][ag]  = RelInfRg[1] * RRcrAG[ag] * Vjaf[1]*(1-Vmix[0]) +
          RelInfRg[0] * RRcrAG[ag] *  Vjaf[0]*Vmix[0];
        ///////// LOW RISK NON US BORN
        VLkla[0 ][1 ][ag]  = RelInfRg[2] * RRcrAG[ag] * Vjaf[2]*(1-Vmix[1]) +
          RelInfRg[0] * RRcrAG[ag] * Vjaf[0]*Vmix[1] + ExogInf[0];
        ///////// HIGH RISK NON US BORN
        VLkla[1 ][1 ][ag]  = RelInfRg[3] * RRcrAG[ag] *
          ((Vjaf[3] * (1-Vmix[0]) * (1-Vmix[1]) +
          Vjaf[2] *    Vmix[0]  * (1-Vmix[1]) +
          Vjaf[1] * (1-Vmix[0]) *    Vmix[1]  +
          Vjaf[0] *    Vmix[0]  *    Vmix[1]))  +
          ExogInf[0];
      }
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
                if (na==0){
                  n2=na;
                } else {n2=1;}

                ///////////////////////////////   SUCEPTIBLE  /////////////////////////////////
                temp = V0[ag][0][0][im][nm][rg][na]*VLkla[rg][n2][ag]*EarlyTrend[m];
                //////////////////////////// REMOVE FROM SUSCEPTIBLE //////////////////////////
                V1[ag][0][0][im][nm][rg][na]  -= temp;
                // Rcpp::Rcout << "susceptible" << (V1[ag][0][0][im][nm][rg][na]  -= temp) << "age= " << ag << "na " << na << "rg " << rg << "\n";

                //////////////////////////////// LATENT TB SLOW ///////////////////////////////
                V1[ag][2][0][im][nm][rg][na]  += temp*MpslowN[ag][im];
                // Rcpp::Rcout << "latent slow" << (V1[ag][2][0][im][nm][rg][na]  += temp*MpslowN[ag][im]) << "age= " << ag << "na " << na << "rg " << rg << "\n";

                //////////////////////////////// LATENT TB FAST ///////////////////////////////
                V1[ag][3][0][im][nm][rg][na]  += temp*MpfastN[ag][im];
                // Rcpp::Rcout << "latent fast" << (V1[ag][3][0][im][nm][rg][na]  += temp*MpfastN[ag][im]) << "age= " << ag << "na " << na << "rg " << rg << "\n";

                ///////////////////////////////////////////////////////////////////////////////

                /////////////////////////////// SUPER-INFECTION SP ////////////////////////////
                temp = V0[ag][1][0][im][nm][rg][na]*VLkla[rg][n2][ag];
                V1[ag][1][0][im][nm][rg][na] -= temp;
                V1[ag][2][0][im][nm][rg][na] += temp*MpslowPIN[ag][im];
                V1[ag][3][0][im][nm][rg][na] += temp*MpfastPIN[ag][im];
                ///////////////////////////////////////////////////////////////////////////////

                /////////////////////////////// SUPER-INFECTION LS ////////////////////////////
                temp = V0[ag][2][0][im][nm][rg][na]*VLkla[rg][n2][ag];
                V1[ag][2][0][im][nm][rg][na]  -= temp;
                V1[ag][2][0][im][nm][rg][na]  += temp*MpslowPIN[ag][im];
                V1[ag][3][0][im][nm][rg][na]  += temp*MpfastPIN[ag][im];

              } } } } }
      ///////////////////////////////////////////////////////////////////////////////
      /////                             BREAK DOWN                              /////
      ///////////////////////////////////////////////////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int im=0; im<4 ; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++) {
                temp  =  V0[ag][2][0][im][nm][rg][na]*(MrslowN[ag][im])*rrSlowFB[na];  // Ls
                temp2  = V0[ag][3][0][im][nm][rg][na]*(rfast);
                V1[ag][2][0][im][nm][rg][na]  -= temp;
                V1[ag][3][0][im][nm][rg][na]  -= temp2;
                V1[ag][4][0][im][nm][rg][na]  += (temp+temp2);
              } } } } }
      ///////////////////////////////////////////////////////////////////////////////
      /////                         LATENT SLOW TO SAFE                         /////
      ///////////////////////////////////////////////////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int im=0; im<4 ; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++) {
                temp  = V0[ag][2][0][im][nm][rg][na]*rRecov;  // Ls
                V1[ag][2][0][im][nm][rg][na]  -= temp;
                V1[ag][1][0][im][nm][rg][na]  += temp;
              } } } } }
      ///////////////////////////////////////////////////////////////////////////////
      /////                              SELF CURE                              /////
      ///////////////////////////////////////////////////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int im=0; im<4 ; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++) {
                temp  = V0[ag][4 ][0][im][nm][rg][na]*rSlfCur;
                V1[ag][4 ][0][im][nm][rg][na]  -= temp;
                V1[ag][2 ][0][im][nm][rg][na]  += temp;
              } } } } }
      ///////////////////////////////////////////////////////////////////////////////
      /////             TB DIAGNOSIS AND TX INITIATION PUBLIC                   /////
      /////                      for all age groups                             /////
      ///////////////////////////////////////////////////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++) {
                temp  = V0[ag][4 ][0][im][nm][rg][na]*rDxtN[0][rg]/RRdxAge[ag]/EarlyTrend[m];
                V1[ag][4 ][0][im][nm][rg][na]  -= temp;
                V1[ag][5 ][0][im][nm][rg][na]  += temp;
              } } } } }
      /////////////////////////////////////////////////////////////////////////////////
      /////                           TREATMENT OUTCOMES                          /////
      /////                     for all age groups, risk groups                   /////
      /////////////////////////////////////////////////////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++){
                ///// CURES ///////////////////////////////////////////////////////////
                temp  = V0[ag][5][0][im][nm][rg][na]*TxVecZ[2];
                V1[ag][5][0][im][nm][rg][na]  -= temp;
                V1[ag][2][0][im][nm][rg][na]  += temp;
                ///// FAILURES(INCLUDING TREATMENT DEFAULT) ///////////////////////////
                temp  = V0[ag][5][0][im][nm][rg][na]*TxVecZ[3];
                V1[ag][5][0][im][nm][rg][na]  -= temp;
                V1[ag][4][0][im][nm][rg][na]  += temp;
              } } } } }
    } //end of tb dynamics loop
    /////////////////////////////////////////////////////////////////////////////////
    /////                          REBALANCE THE POPULATION                     /////
    /////                          & UPDATE THE DISTRIBUTION                    /////
    /////////////////////////////////////////////////////////////////////////////////
    if (reblnc==1){
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for (int im=0; im<4; im++){
            for (int nm=0; nm<4; nm++){
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++){
                  for (int m2=0; m2<4; m2++){
                    for (int p2=0; p2<4; p2++){
                      temp = V1[ag][tb][0][p2][m2][rg][na]*(trans_mat_tot_agesN[m2+p2*4][(16*(ag+1))-(16-(nm+im*4))]);
                      V1[ag][tb][0][im][nm][rg][na] += temp;
                      V1[ag][tb][0][p2][m2][rg][na] -= temp;

                    } } } } } } } }
    } //end of rebalancing loop
    /////////////////////////////////////////////////////////////////////////////////
    /////                          RESET POPULATION SIZE                        /////
    /////////////////////////////////////////////////////////////////////////////////
    for(int ag=0; ag<11; ag++) {
      for(int i=0; i<2; i++) {
        InitPopZ[ag][i] = 0; //i is for us vs non us born
      } }
    for(int ag=0; ag<11; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int nm=0; nm<4; nm++) {
          for(int im=0; im<4; im++) {
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
        for(int nm=0; nm<4; nm++) {
          for(int im=0; im<4; im++) {
            for(int rg=0; rg<2; rg++) {  // reset pop to InitPop
              V1[ag][tb][0][im][nm][rg][0]  = V1[ag][tb][0][im][nm][rg][0]*InitPopZ[ag][0];
              V1[ag][tb][0][im][nm][rg][1]  = V1[ag][tb][0][im][nm][rg][1]*InitPopZ[ag][1];
              V1[ag][tb][0][im][nm][rg][2]  = V1[ag][tb][0][im][nm][rg][2]*InitPopZ[ag][1];
            } } } } }

    for(int ag=0; ag<11; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int nm=0; nm<4; nm++) {
          for(int im=0; im<4; im++) {
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++){
                V0[ag][tb][0][im][nm][rg][na] = V1[ag][tb][0][im][nm][rg][na];
              } } } } } }
  }
  /////////////////////////////////////////////////////////////////////////////////
  /////                            END BURN IN                                /////
  //////////////////////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////////////////////
  /////              CHECK THE SIZE OF THE POPULATION STATES                  /////
  /////////////////////////////////////////////////////////////////////////////////
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
  /////////////////////////////////////////////////////////////////////////////////
  /////                            BEGIN MODEL                                /////
  /////////////////////////////////////////////////////////////////////////////////
  /////                  CREATE TIME LOOPS FOR MODEL RUN                      /////
  /////                               YEAR LOOP                               /////
  for(int y=0; y<nYrs; y++) {
    /////                            MONTH  LOOP                              /////
    for(int m=0; m<12; m++) {
      /////              CREATE A COUNTER OF MONTHS SINCE START               /////
      s = y*12+m;
      /////////////////////////////////////////////////////////////////////////////////
      /////                      UPDATING TREATMENT PARAMETERS                    /////
      /////  TxMatZ: 0=completion rate, 1 = tx success, 2 = RATE OF EXIT TO CURE  /////
      /////     3 = RATE OF EXIT TO ACTIVE TB, 4 = RATE OF EXIT TO RETREATMENT    /////
      /////                5 = PROBABILITY OF TREATMENT COMPLETION                /////
      /////////////////////////////////////////////////////////////////////////////////
      ///// TREATMENT EFFICACY UPDATED FOR TREATMENT QUALITY //////////////////////////
      TxVecZ[1] = TxVec[1]*TxQualt[s];
      ///// RATE OF TREATMENT EXIT TO CURE (LS) ///////////////////////////////////////
      TxVecZ[2] = TxVec[0]*TxVecZ[1] + rDeft[s]*TxVecZ[1]*RRcurDef;
      ///// RATE OF TREATMENT EXIT TO ACTIVE TB ///////////////////////////////////////
      TxVecZ[3] = TxVec[0]*(1-TxVecZ[1])*(1-pReTx[s]) + rDeft[s]*(1-TxVecZ[1]*RRcurDef);
      ///// RATE OF TREATMENT EXIT TO RE TREATMENT/////////////////////////////////////
      TxVecZ[4] = TxVec[0]*(1-TxVecZ[1])*(pReTx[s]);
      ///// PROBABILITY (TREATMENT COMPLETION) ////////////////////////////////////////
      TxVecZ[5] = TxVec[0]*(1-(1.0-TxVecZ[1])*pReTx[s]);
      /////////////////////////////////////////////////////////////////////////////////
      /////                                 BIRTHS                                /////
      /////        USE DISTRIBUTION TO POPULATE THE MODEL ACROSS RISK GROUPS      /////
      /////////////////////////////////////////////////////////////////////////////////
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++){
          V1[0][0][0][im][nm][0][0]  += Birthst[s]*dist_genN[nm][im]*(1-p_HR);
          V1[0][0][0][im][nm][1][0]  += Birthst[s]*dist_genN[nm][im]*(p_HR);
        } }
      /////////////////////////////////////////////////////////////////////////////////
      /////                               IMMIGRATION                             /////
      /////////////////////////////////////////////////////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++){
            ///// NO TB, Low risk  ////////////////////////////////////////////////////
            V1[ag][0][0][im][nm][0][1]   += ImmNonN[s][ag]*dist_genN[nm][im]*(1-p_HR);
            ///// NO TB, High risk  ///////////////////////////////////////////////////
            V1[ag][0][0][im][nm][1][1]   += ImmNonN[s][ag]*dist_genN[nm][im]*(p_HR);
            ///// NO TB, Low risk  ////////////////////////////////////////////////////
            V1[ag][1][0][im][nm][0][1]   += SpImmNonN[s][ag]*(1-p_HR)*dist_genN[nm][im];
            ///// NO TB, High risk  ///////////////////////////////////////////////////
            V1[ag][1][0][im][nm][1][1]   += SpImmNonN[s][ag]*(p_HR)*dist_genN[nm][im];
            ///// LATENT SLOW TB, Low risk  ///////////////////////////////////////////
            V1[ag][2][0][im][nm][0][1]   += ImmLatN[s][ag]*dist_genN[nm][im]*(1-p_HR);
            ///// LATENT SLOW TB, High risk  //////////////////////////////////////////
            V1[ag][2][0][im][nm][1][1]   += ImmLatN[s][ag]*dist_genN[nm][im]*(p_HR);
            ///// LATENT FAST TB, Low risk  ///////////////////////////////////////////
            V1[ag][3][0][im][nm][0][1]   += ImmFstN[s][ag]*dist_genN[nm][im]*(1-p_HR);
            ///// LATENT FAST TB, High risk  //////////////////////////////////////////
            V1[ag][3][0][im][nm][1][1]   += ImmFstN[s][ag]*dist_genN[nm][im]*(p_HR);
            ///// ACTIVE TB, Low risk  ////////////////////////////////////////////////
            V1[ag][4][0][im][nm][0][1]   += ImmActN[s][ag]*dist_genN[nm][im]*(1-p_HR);
            ///// ACTIVE TB, High risk  ///////////////////////////////////////////////
            V1[ag][4][0][im][nm][1][1]   += ImmActN[s][ag]*dist_genN[nm][im]*(p_HR);
          } } }
      /////////////////////////////////////////////////////////////////////////////////
      /////                            EMMIGRATION                                /////
      /////////////////////////////////////////////////////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++){
            for(int im=0; im<4; im++){
              for(int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++) {
                  ///// In US under 2 years  ////////////////////////////////////////
                  V1[ag][tb][lt][im][nm][rg][1]  -= V0[ag][tb][lt][im][nm][rg][1]*rEmmigFB[0];
                  ///// In US over 2 years  /////////////////////////////////////////
                  V1[ag][tb][lt][im][nm][rg][2]  -= V0[ag][tb][lt][im][nm][rg][2]*rEmmigFB[1];
                } } } } } }
      /////////////////////////////////////////////////////////////////////////////////
      /////                         NET INTERNAL MIGRATION                        /////
      /////////////////////////////////////////////////////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++){
            for(int im=0; im<4; im++){
              for(int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++) {
                  ///// USB ///////////////////////////////////////////////////////////
                  V1[ag][tb][lt][im][nm][rg][0]  += V0[ag][tb][lt][im][nm][rg][0] * net_mig_usb[ag];
                  ///// NUSB in US under 2 years  /////////////////////////////////////
                  V1[ag][tb][lt][im][nm][rg][1]  += V0[ag][tb][lt][im][nm][rg][1] * net_mig_nusb[ag];
                  ///// NUSB in US over 2 years  //////////////////////////////////////
                  V1[ag][tb][lt][im][nm][rg][2]  += V0[ag][tb][lt][im][nm][rg][2] * net_mig_nusb[ag];
                } } } } } }
      /////////////////////////////////////////////////////////////////////////////
      /////                              MORTALITY                            /////
      /////////////////////////////////////////////////////////////////////////////
      temp = 0; temp2 = 0;
      for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++){
          for(int im=0; im<4 ; im++) {
            for(int nm=0; nm<4; nm++) {
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                  if ((ag<9) | ((mubtN[s][ag]*RRmuRFN[nm]*RRmuHR[rg]) < .5)){
                    temp = ((RRmuRFN[nm]*RRmuHR[rg])*mubtN[s][ag]);}
                  else {
                    temp =  .5;
                  }
                  ///// WITHOUT TB  ///////////////////////////////////////////////////
                  for(int tb=0; tb<6; tb++) {
                    if (tb < 4){
                      temp2  = V0[ag][tb][lt][im][nm][rg][na]*temp;
                    } else if (tb == 4) {
                      ///// WITH ACTIVE TB //////////////////////////////////////////////////
                      temp2  = V0[ag][4][lt][im][nm][rg][na]*(temp+vTMortN[ag][4]*RRmuTBPand[s]);
                    } else if (tb == 5) {
                      ///// WITH TB TREATMENT ///////////////////////////////////////////////
                      temp2 = V0[ag][5][lt][im][nm][rg][na]*(temp+vTMortN[ag][5 ]*pow(1.0-TxVecZ[1],TunTxMort)*RRmuTBPand[s]);
                    }
                    ///// UPDATE THE PRIMARY VECTOR BY REMOVING MORTALITY /////////////////
                    V1   [ag][tb][lt][im][nm][rg][na] -= temp2;
                    VMort[ag][tb][lt][im][nm][rg][na] += temp2;
                  } } } } } } }
      /////////////////////////////////////////////////////////////////////////////
      /////                              AGING                                /////
      /////////////////////////////////////////////////////////////////////////////
      for(int ag=0; ag<10; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++){
            for (int im=0; im<4; im++){
              for (int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++){
                    temp = V0[ag  ][tb][lt][im][nm][rg][na]/ag_denN[s][ag];
                    V1[ag  ][tb][lt][im][nm][rg][na]  -= temp;
                    V1[ag+1][tb][lt][im][nm][rg][na]  += temp;
                  } } } } } } }
      /////////////////////////////////////////////////////////////////////////////
      /////                     NEW FB -> ESTABLISHED FB                      /////
      /////                     TWO YEARS FOR TRANSITION                      /////
      /////////////////////////////////////////////////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++){
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  temp = V0[ag][tb][lt][im][nm][rg][1] / 24;
                  V1[ag][tb][lt][im][nm][rg][1]  -= temp;
                  V1[ag][tb][lt][im][nm][rg][2]  += temp;
                } } } } } }
      /////////////////////////////////////////////////////////////////////////////
      /////                        HIGH-RISK ENTRY/EXIT                       /////
      /////////////////////////////////////////////////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int na=0; na<3; na++) {
                  ///// ENTRY /////////////////////////////////////////////////////
                  temp  = V0[ag][tb][lt][im][nm][0][na]*HrEntExN[ag][0];
                  ///// EXIT //////////////////////////////////////////////////////
                  temp2 = V0[ag][tb][lt][im][nm][1][na]*HrEntExN[ag][1];
                  V1[ag][tb][lt][im][nm][0][na]  += (temp2-temp);
                  V1[ag][tb][lt][im][nm][1][na]  += (temp-temp2);
                } } } } } }
      /////////////////////////////////////////////////////////////////////////////
      /////                      OPEN TB DYNAMICS LOOP                        /////
      /////////////////////////////////////////////////////////////////////////////
      if (tb_dyn==1){
        ///////////////////////////////////////////////////////////////////////////
        ///// Step 1 & 2
        ///// take total population of mixing groups
        ///// mixing groups are only risk group and nativity
        ///// we do not need to stratify by age here because we are only allowing
        ///// the number of contacts of infectious persons to vary as a constant
        ///// across all mixing groups; no assortative mixing  RRcrAG[ag] by age.
        ///////////////////////////////////////////////////////////////////////////
        ////////////////////////////  TRANSMISSION RISK  ////////////////////////////////

        // Step 1 & 2
        // take total population of mixing groups
        // mixing groups are only risk group and nativity
        // we do not need to stratify by age here because we are only allowing the number of contacts
        // of infectious persons to vary as a constant across all mixing groups; no assortative mixing  RRcrAG[ag]
        // by age.

        for(int i=0; i<2; i++){
          for(int j=0; j<2; j++){
            VNkl[i][j]=0;
            VGkl[i][j]=0;
          }
        }
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                  ///////// LOW RISK US BORN
                  VNkl[0][0] += V0[ag][tb][lt][im][nm][0][0]
                  * RRcrAG[ag] * RelInfRg[0];
                  VGkl[0][0] += V0[ag][tb][lt][im][nm][0][0]
                  * RRcrAG[ag] * RelInfRg[0] * RelInf[tb];
                  ///////// HIGH RISK US BORN
                  VNkl[1][0] += V0[ag][tb][lt][im][nm][1][0]
                  * RRcrAG[ag] * RelInfRg[1];
                  VGkl[1][0] += V0[ag][tb][lt][im][nm][1][0]
                  * RRcrAG[ag] * RelInfRg[1] * RelInf[tb];
                  ///////// LOW RISK NON US BORN
                  VNkl[0][1] += (V0[ag][tb][lt][im][nm][0][1] + V0[ag][tb][lt][im][nm][0][2])
                    * RRcrAG[ag] * RelInfRg[2];
                  VGkl[0][1] += (V0[ag][tb][lt][im][nm][0][1] + V0[ag][tb][lt][im][nm][0][2])
                    * RRcrAG[ag] * RelInfRg[2] * RelInf[tb];
                  ///////// HIGH RISK NON US BORN
                  VNkl[1][1] += (V0[ag][tb][lt][im][nm][1][1] + V0[ag][tb][lt][im][nm][1][2])
                    * RRcrAG[ag] * RelInfRg[3];
                  VGkl[1][1] += (V0[ag][tb][lt][im][nm][1][1] + V0[ag][tb][lt][im][nm][1][2])
                    * RRcrAG[ag] * RelInfRg[3] * RelInf[tb];
                } } } } }


        // Step 3 (Infected/Total)
        //probability of infection  -- maybe not stratified by age either?
        //for each mixing group
        // calculated as the a weighted average based on the number of effective contacts
        // each nativity and risk group contributes to that mixing group
        // 0 = common pool; 1 = exclusive nusb, 2 = exclusive hr, 3 = exclusive nusb-hr

        Vjaf[0] = (VGkl[0][0]         + VGkl[1][0]*Vmix[0] +
          VGkl[0][1]*Vmix[1] + VGkl[1][1]*Vmix[1]*Vmix[0]) /
            (VNkl[0][0]         + VNkl[1][0]*Vmix[0]+
              VNkl[0][1]*Vmix[1] + VNkl[1][1]*Vmix[1]*Vmix[0] + 1e-12);

        Vjaf[1] = (VGkl[0][1] + VGkl[1][1]*Vmix[0]) /
          (VNkl[0][1] + VNkl[1][1]*Vmix[0] + 1e-12);

        Vjaf[2] = (VGkl[1][0] + VGkl[1][1]*Vmix[1])/
          (VNkl[1][0] + VNkl[1][1]*Vmix[1] + 1e-12);

        Vjaf[3] = VGkl[1][1] / (VNkl[1][1] + 1e-12);

        // Step 4
        // calculate force of infection
        // the number of contacts you have*probability of infection*mixing
        for(int ag=0; ag<11; ag++) {
          /// LOW RISK US BORN
          VLkla[0 ][0 ][ag]  = RelInfRg[0] * RRcrAG[ag] * Vjaf[0] ;
          ///////// HIGH RISK US BORN
          VLkla[1 ][0 ][ag]  = RelInfRg[1] * RRcrAG[ag] * Vjaf[1]*(1-Vmix[0]) +
            RelInfRg[0] * RRcrAG[ag] *  Vjaf[0]*Vmix[0];
          ///////// LOW RISK NON US BORN
          VLkla[0 ][1 ][ag]  = RelInfRg[2] * RRcrAG[ag] * Vjaf[2]*(1-Vmix[1]) +
            RelInfRg[0] * RRcrAG[ag] * Vjaf[0]*Vmix[1] + ExogInf[s];
          ///////// HIGH RISK NON US BORN
          VLkla[1 ][1 ][ag]  =  VLkla[1 ][1 ][ag]  = RelInfRg[3] * RRcrAG[ag] *
            ((Vjaf[3] * (1-Vmix[0]) * (1-Vmix[1]) +
            Vjaf[2] *    Vmix[0]  * (1-Vmix[1]) +
            Vjaf[1] * (1-Vmix[0]) *    Vmix[1]  +
            Vjaf[0] *    Vmix[0]  *    Vmix[1]))  +
            ExogInf[s];
        }

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

                    if (na==0){
                      n2=na;
                    } else {n2=1;}

                    ///////////////////////////////   SUCEPTIBLE  /////////////////////////////////
                    temp = V0[ag][0][lt][im][nm][rg][na]*(VLkla[rg][n2][ag])*NixTrans[s];
                    Vinf[ag][0][lt][im][nm][rg][na] = (V0[ag][0][lt][im][nm][rg][na]*VLkla[rg][n2][ag]*NixTrans[s]);
                    //////////////////////////// REMOVE FROM SUSCEPTIBLE //////////////////////////
                    V1[ag][0][lt][im][nm][rg][na]  -= temp;
                    //////////////////////////////// LATENT TB SLOW ///////////////////////////////
                    V1[ag][2][lt][im][nm][rg][na]  += temp*MpslowN[ag][im];
                    //////////////////////////////// LATENT TB FAST ///////////////////////////////
                    V1[ag][3][lt][im][nm][rg][na]  += temp*MpfastN[ag][im];
                    ///////////////////////////////////////////////////////////////////////////////

                    /////////////////////////////// SUPER-INFECTION SP ////////////////////////////
                    temp = V0[ag][1][lt][im][nm][rg][na]*(VLkla[rg][n2][ag])*NixTrans[s];
                    // if((V0[ag][1][lt][im][nm][rg][na]*VLkla[rg][n2][ag]*NixTrans[s]) > V1[ag][1][lt][im][nm][rg][na] ){
                    //
                    //   Rcpp::Rcout << "sum is "<<  (V0[ag][1][lt][im][nm][rg][na]*VLkla[rg][n2][ag]*NixTrans[s]) << "at ag = " << ag  << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
                    //   Rcpp::Rcout << "sum is "<<  V1[ag][1][lt][im][nm][rg][na] << "at ag = " << ag <<  "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
                    //
                    //   return
                    //     Rcpp::List::create(
                    //       Rcpp::Named("Outputs") = Outputs2
                    //     );
                    // }
                    V1[ag][1][lt][im][nm][rg][na] -= temp;
                    V1[ag][2][lt][im][nm][rg][na] += temp*MpslowPIN[ag][im];
                    V1[ag][3][lt][im][nm][rg][na] += temp*MpfastPIN[ag][im];
                    ///////////////////////////////////////////////////////////////////////////////

                    /////////////////////////////// SUPER-INFECTION LS ////////////////////////////
                    temp = V0[ag][2][lt][im][nm][rg][na]*(VLkla[rg][n2][ag])*NixTrans[s];
                    V1[ag][2][lt][im][nm][rg][na]  -= temp;
                    V1[ag][2][lt][im][nm][rg][na]  += temp*MpslowPIN[ag][im];
                    V1[ag][3][lt][im][nm][rg][na]  += temp*MpfastPIN[ag][im];

                  } } } } } }
          ///////////////////////////////////////////////////////////////////////////////
          /////                             BREAK DOWN                              /////
          ///////////////////////////////////////////////////////////////////////////////
          for(int ag=0; ag<11; ag++) {
            for(int lt=0; lt<2 ; lt++) {
              for(int im=0; im<4 ; im++) {
                for(int nm=0; nm<4; nm++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      temp  = V0[ag][2][lt][im][nm][rg][na]*(MrslowN[ag][im])*rrSlowFB[na];  // Latent Slow
                      temp2 = V0[ag][3][lt][im][nm][rg][na]*(rfast);  // Latent Fast
                      //////REMOVE FROM LATENT SLOW AND LATENT FAST AND PLACE IN ACTIVE DISEASE
                      V1[ag][2][lt][im][nm][rg][na]  -= temp;  //REMOVE FROM LATENT SLOW
                      V1[ag][3][lt][im][nm][rg][na]  -= temp2; //REMOVE FROM LATENT FAST
                      V1[ag][4][lt][im][nm][rg][na]  += (temp+temp2); //PLACE IN ACTIVE DISEASE
                    } } } } } }
          ///////////////////////////////////////////////////////////////////////////////
          /////                         LATENT SLOW TO SAFE                         /////
          ///////////////////////////////////////////////////////////////////////////////
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
          ///////////////////////////////////////////////////////////////////////////////
          /////                              SELF CURE                              /////
          ///////////////////////////////////////////////////////////////////////////////
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
          ///////////////////  //////////////////  //////////////////  ///////////////////
          ///////////////////  LTBI SCREENING AND TLTBI INITIATION     ///////////////////
          ///////////////////  //////////////////  //////////////////  ///////////////////

          ////////////////////////////////////////////////////////////////////////////////
          /////                     BASE LINE VALUE CALCULATIONS                     /////
          ////////////////////////////////////////////////////////////////////////////////
          ///// calculate the rates of tb positive for the basecase and the rate of  /////
          ///// TB negatives these are equal to the rate of screening combined with  /////
          ///// the sensitivity or specificity in the base case                      /////
          ///// tempHR is adjustment to add extra screening for HR in param          /////
          ////////////////////////////////////////////////////////////////////////////////

          int tempHR = 1;
          int tempNTB = 1;
          double tx_eligible_pop = 0;
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++) {
              for(int ag=0; ag<11; ag++) {
                ///// US BORN, LOW RISK  ////////////////////////////////////////////////
                if(rg==0 && na==0) {
                  rTbP = rLtScrtN[s][0]*LtDxPar_ltN[0][s];
                  rTbN = rLtScrtN[s][0]*LtDxPar_noltN[0][s];
                  tempHR=1;
                  tempNTB = rrTestLrNoTb;
                }
                ///// US BORN, HIGH RISK  ///////////////////////////////////////////////
                if(rg==1 && na==0) {
                  rTbP = rLtScrtN[s][1]*LtDxPar_ltN[1][s];
                  rTbN = rLtScrtN[s][1]*LtDxPar_noltN[1][s];
                  tempHR = rrTestHr;
                  tempNTB = 1;
                }
                ///// Young NUS (under 5)  ////////////////////////////////////////////
                if(rg==0 && na > 0 && ag==0) {
                  rTbP = rLtScrtN[s][0]*LtDxPar_ltN[2][s];
                  rTbN = rLtScrtN[s][0]*LtDxPar_noltN[2][s];
                  tempHR = 1;
                  tempNTB = 1;
                }
                ///// NON US BORN  ////////////////////////////////////////////////////
                if(rg==0 && na > 0 && ag > 0) {
                  rTbP = rLtScrtN[s][0]*LtDxPar_ltN[3][s];
                  rTbN = rLtScrtN[s][0]*LtDxPar_noltN[3][s];
                  tempHR = 1;
                  tempNTB = 1;
                }
                ///// NON US BORN, HIGH RISK  //////////////////////////////////////////
                if(rg==1 && na >0) {
                  rTbP = rLtScrtN[s][1]*LtDxPar_ltN[4][s];
                  rTbN = rLtScrtN[s][1]*LtDxPar_noltN[4][s];
                  tempHR = rrTestHr;
                  tempNTB = 1;
                }
                for(int im=0; im<4; im++) {
                  for(int nm=0; nm<4; nm++) {
                    /////////////////////////////////////////////////////////////////////
                    ///// calculate the number of tests                              ////
                    /////////////////////////////////////////////////////////////////////
                    VLtest[ag][0][0][im][nm][rg][na] += V0[ag][0][0][im][nm][rg][na] * rLtScrtN[s][rg] * tempHR * rrTestLrNoTb;
                    VLtest[ag][1][0][im][nm][rg][na] += V0[ag][1][0][im][nm][rg][na] * rLtScrtN[s][rg] * tempHR * rrTestLrNoTb;
                    VLtest[ag][2][0][im][nm][rg][na] += V0[ag][2][0][im][nm][rg][na] * rLtScrtN[s][rg] * tempHR;
                    VLtest[ag][3][0][im][nm][rg][na] += V0[ag][3][0][im][nm][rg][na] * rLtScrtN[s][rg] * tempHR;
                    /////////////////////////////////////////////////////////////////////
                    ///// we also need to store the number of tests in this time step
                    ///// to decrease the susceptibles if there is a TTT intervention.
                    /////////////////////////////////////////////////////////////////////
                    // for (int tb = 0; tb < 4; tb++){
                    //   VLtemp[ag][tb][0][im][nm][rg][na] = VLtest[ag][tb][0][im][nm][rg][na];
                    // }
                    ///////////////////////////////////////////////////////////////////////////////
                    /////                     FALSE LATENT POSITIVES                          /////
                    ///////////////////////////////////////////////////////////////////////////////
                    VLDxtemp[ag][0][0][im][nm][rg][na] = V0[ag][0][0][im][nm][rg][na]*rTbN;
                    VLdx[ag][0][0][im][nm][rg][na] += V0[ag][0][0][im][nm][rg][na]*rTbN;
                    VLDxtemp[ag][1][0][im][nm][rg][na] = V0[ag][1][0][im][nm][rg][na]*rTbN;
                    VLdx[ag][1][0][im][nm][rg][na] += V0[ag][1][0][im][nm][rg][na]*rTbN;
                    ///// Remove from the TB naive and PI states ////////////////////////
                    V1[ag][0][0][im][nm][rg][na]  -= VLDxtemp[ag][0][0][im][nm][rg][na];
                    V1[ag][1][0][im][nm][rg][na]  -= VLDxtemp[ag][1][0][im][nm][rg][na];
                    ///// Move to latent tx experienced /////////////////////////////////
                    V1[ag][0][1][im][nm][rg][na]  += VLDxtemp[ag][0][0][im][nm][rg][na];
                    V1[ag][1][1][im][nm][rg][na]  += VLDxtemp[ag][1][0][im][nm][rg][na];
                    ///////////////////////////////////////////////////////////////////////////////
                    /////                     TRUE LATENT POSITIVES                           /////
                    ///////////////////////////////////////////////////////////////////////////////
                    ///// Update the number of diagnoses
                    VLDxtemp[ag][2][0][im][nm][rg][na] = V0[ag][2][0][im][nm][rg][na] * rTbP;
                    VLdx[ag][2][0][im][nm][rg][na] += V0[ag][2][0][im][nm][rg][na] * rTbP;
                    VLDxtemp[ag][3][0][im][nm][rg][na] = V0[ag][3][0][im][nm][rg][na] * rTbP;
                    VLdx[ag][3][0][im][nm][rg][na] += V0[ag][3][0][im][nm][rg][na] * rTbP;

                    ///// Calculate the split between treatment completion and default.
                    temp  = VLDxtemp[ag][2][0][im][nm][rg][na]*LtTxParN[s][0]*(1-LtTxParN[s][1]); // tx completion
                    temp3 = VLDxtemp[ag][2][0][im][nm][rg][na]*LtTxParN[s][0]*LtTxParN[s][1]; // default

                    temp2 = VLDxtemp[ag][3][0][im][nm][rg][na]*LtTxParN[s][0]*(1-LtTxParN[s][1]); // tx completion
                    temp4 = VLDxtemp[ag][3][0][im][nm][rg][na]*LtTxParN[s][0]*LtTxParN[s][1]; // default

                    V1[ag][2][0][im][nm][rg][na]  -=  (temp+temp3);   //remove from latent slow
                    V1[ag][3][0][im][nm][rg][na]  -=  (temp2+temp4);  //remove from latent fast
                    //completion split between success and failure
                    V1[ag][1][1][im][nm][rg][na]  += (temp+temp2) * LtTxParN[s][2];  //exit to cure
                    V1[ag][2][1][im][nm][rg][na]  += (temp+temp2) * (1-LtTxParN[s][2]); //tx comp fail to latent slow
                    ///defaults are placed in tx naive because it is considered the same tb infection
                    V1[ag][2][0][im][nm][rg][na]  += (temp3+temp4); //latent tx default to latent slow

                    ///// how many treatment eligible remain

                    tx_eligible_pop += V1[ag][0][0][im][nm][rg][na] +
                                       V1[ag][1][0][im][nm][rg][na] +
                                       V1[ag][2][0][im][nm][rg][na] +
                                       V1[ag][3][0][im][nm][rg][na];

                    temp25 += VLdx[ag][2][0][im][nm][rg][na] + VLdx[ag][3][0][im][nm][rg][na];
                  } } } } }

          ///////////////////////////////////////////////////////////////////////////////
          ///////////////////////////////////////////////////////////////////////////////
          ///////////////////////////////////////////////////////////////////////////////
          ///////////////////////////////////////////////////////////////////////////////
          ///// RESET SEVERAL VARIABLES FOR TTT EXTRA SCREENING TO ZERO FOR BASECASE
          ///////////////////////////////////////////////////////////////////////////////
          ///////////////////////////////////////////////////////////////////////////////
          ///////////////////////////////////////////////////////////////////////////////
          ///////////////////////////////////////////////////////////////////////////////
          tempHR = 1;
          temp=temp2=temp3=temp4=temp5=temp6=temp7=temp8=temp9=temp10=0;
          ///// calculate the number of rows which are the number of populations screened
          // int rows = sizeof(ttt_samp_distN)/sizeof(ttt_samp_distN[0]);
          // int cols = 352;
          // ///// initialize variables
          // for (int l=0; l<rows; l++){
          //   rr_ltbi=1;
          //   for (int c=0; c<cols; c++){
          //     ttt_dist[l][c]=0;
          //   }
          // }
          ///////////////////////////////////////////////////////////////////////////////
          /////////////////// CHECK FOR ANY OF THE CASCADE COUNTERFACTUALS //////////////
          ///////////////////////////////////////////////////////////////////////////////
          ///// PROBABILITY OF LTBI TESTING ACCEPTANCE
          ttt_ltbi_acceptN=1;
          ///// PROBABILITY OF LTBI TREATMENT INITIATION
          if (ttt_ltbi_init==1){
            ttt_ltbi_initN=1;
          } else{
            ttt_ltbi_initN = LtTxParN[s][0];
          }
          ///// PROBABILITY OF LTBI TREATMENT COMPLETION
          if (ttt_ltbi_comp==1){
            ttt_ltbi_compN=1;
          } else{
            ttt_ltbi_compN=(1-LtTxParN[s][1]);
          }
          ///// PROBABILITY OF LTBI TREATMENT EFFECTIVENESS
          if (ttt_ltbi_eff==1){
            ttt_ltbi_effN=1;
          } else{
            ttt_ltbi_effN=LtTxParN[s][2];
          }
          // temp10 = 0; temp9=0; temp11=0;
          double tempSum = 0.0;
          ///////////////////////////////////////////////////////////////////////////////
          ///////////////////////////////////////////////////////////////////////////////
          ///////////////////////////////////////////////////////////////////////////////
          /////            START TTT EXTRA SCREENING INTERVENTION CODE              /////
          ///////////////////////////////////////////////////////////////////////////////
          ///////////////////////////////////////////////////////////////////////////////
          ///////////////////////////////////////////////////////////////////////////////
          ///// open a loop to initiate ttt when the month iterator falls within the
          ///// designated range
          ///////////////////////////////////////////////////////////////////////////////
          if (s == ttt_month.front()){
            ///////////////////////////////////////////////////////////////////////////////
            ///// CREATE OBJECTS TO HOLD THE TESTED AND DIAGNOSED FOR EACH POP
            ///// THESE DO NOT NEED TO BE SAVED BETWEEN MODEL STATE COMBOS
            ///////////////////////////////////////////////////////////////////////////////
            double ttt_test_all_vec[11][5][4][4][2][3] = {0};
            double ttt_diag_all_vec[11][5][4][4][2][3] = {0};

            rr_ltbi = ttt_ltbi[0];
            Rcpp::Rcout << "rr ltbi is " << rr_ltbi << "\n";
            for(int ag=0; ag<11; ag++) {
              for(int nm=0; nm<4; nm++) {
                for(int im=0; im<4; im++) {
                  for(int na=0; na<3; na++) {
                    if (na < 1) {
                      ni = 0;
                    } else{
                      ni = 1;
                    }
                    /// SET THE SENSITIVITY OF THE TEST BELOW
                    if (ttt_ltbi_acceptN == 1){
                      // ////////////// ALL US BORN  //////////////////////////
                      if(ni==0) {
                        rTbP_norm = 0.78;
                        rTbN_norm = (1 - 0.98);
                      }
                      ////////////// Young NUS (under 5)  /////////////////
                      if(ni > 0 && ag==0) {
                        rTbP_norm = 0.79;
                        rTbN_norm = (1 - 0.99);
                      }
                      //////////// NON US BORN  ///////////////////////////
                      if(ni > 0 && ag > 0) {
                        rTbP_norm = 0.79;
                        rTbN_norm = (1 - 0.99);
                      }
                    } else{
                      // ////////////// ALL US BORN  //////////////////////////
                      if(ni==0) {
                        rTbP_norm = LtDxPar_ltN[0][s];
                        rTbN_norm = LtDxPar_noltN[0][s];
                      }
                      ////////////// Young NUS (under 5)  /////////////////
                      if(ni > 0 && ag==0) {
                        rTbP_norm = LtDxPar_ltN[2][s];
                        rTbN_norm = LtDxPar_noltN[2][s];
                      }
                      //////////// NON US BORN  ///////////////////////////
                      if(ni > 0 && ag > 0) {
                        rTbP_norm = LtDxPar_ltN[3][s];
                        rTbN_norm = LtDxPar_noltN[3][s];
                      }
                    }
                    ///special care cascade (perfect)
                      if(ttt_ltbi_sens==1){
                        rTbP_norm=1;
                      }
                      if(ttt_ltbi_spec==1){
                        rTbN_norm=1;
                   }
                  /// define the fraction screened

                  ttt_frc_scrn = std::min(ttt_samp_distN[0][ag+(ni*11)+(nm*22)+(im*88)], 1.0);

                  for(int rg=0; rg<2; rg++) {
                    for(int tb=0; tb<5; tb++) {

                ///////////////////////////////////////////////////////////////////////////////
                /////                 CALCULATE THE NUMBER OF EXTRA TESTS                 /////
                ///////////////////////////////////////////////////////////////////////////////
                ///// second loop that iterates across all the populations being screened
                ///// need the VLdx to be the sum of all additional screening, but we can't sum
                ///// here because we need the population to change;
                ///////////////////////////////////////////////////////////////////////////////

                ///// test among the TB disease population first
                if (tb == 4){
                  ttt_test_all_vec[ag][tb][im][nm][rg][na] =  ttt_frc_scrn * V0[ag][tb][0][im][nm][rg][na] * ttt_ltbi_acceptN;
                ///// test among the LTBI population next
                } else if ((tb == 2) || (tb == 3)){
                  ttt_test_all_vec[ag][tb][im][nm][rg][na] =  ttt_frc_scrn * V0[ag][tb][0][im][nm][rg][na] * ttt_ltbi_acceptN * rr_ltbi;
                ///// test among the susceptible population last
                ///// adjust the susceptible population to have fewer susceptibles screened due to the elevated LTBI
                } else if (tb < 2){
                  double adj_susceptible_pop[11][2][4][4][2][3] = {0};
                  double frcStrata = V0[ag][tb][0][im][nm][rg][na] / (V0[ag][0][0][im][nm][rg][na] + V0[ag][1][0][im][nm][rg][na]);
                  adj_susceptible_pop[ag][tb][im][nm][rg][na] = V0[ag][tb][0][im][nm][rg][na] + ((1 - rr_ltbi)* (V0[ag][2][0][im][nm][rg][na] + V0[ag][3][0][im][nm][rg][na]) * frcStrata);
                  ttt_test_all_vec[ag][tb][im][nm][rg][na] =  ttt_frc_scrn * adj_susceptible_pop[ag][tb][im][nm][rg][na] * ttt_ltbi_acceptN;
                }
                ///// Update the test accumulator
                VLtest[ag][tb][0][im][nm][rg][na] += ttt_test_all_vec[ag][tb][im][nm][rg][na];
                /////////////////// CALCULATE THE NUMBER OF EXTRA DIAGNOSES
                ///// DO WE NEED TO ACCOUNT FOR HIGH RISK HERE?
                double testPerformance = 0;
                if ((tb == 2) || (tb == 3) || (tb == 4)){
                  testPerformance = rTbP_norm;
                } else {
                  testPerformance = rTbN_norm;
                }
                ttt_diag_all_vec[ag][tb][im][nm][rg][na] = std::max((std::min((ttt_test_all_vec[ag][tb][im][nm][rg][na] * testPerformance),
                                                                               ttt_test_all_vec[ag][tb][im][nm][rg][na])), 0.0);

                // Sum of total screened population for print below
                temp20 += ttt_test_all_vec[ag][tb][im][nm][rg][na];
                if (tb < 4){
                // Sum of total diagnosed population for print below
                temp22 += ttt_diag_all_vec[ag][tb][im][nm][rg][na];
                }

} // close TB loop
                ////////////////////////////////////////////////////////////////////////////
                /// REMOVE THE TREATED POPULATION
                /////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////
                /// remove those who test positive who are true LTBI negatives (false positives)
                /////////////////////////////////////////////////////////////////////////////
                /// Remove from the TB naive and PI states ////////////////////////
                V1[ag][0][0][im][nm][rg][na]  -= ttt_diag_all_vec[ag][0][im][nm][rg][na];
                V1[ag][1][0][im][nm][rg][na]  -= ttt_diag_all_vec[ag][1][im][nm][rg][na];
                ///// Update the number of diagnoses
                VLdx[ag][0][0][im][nm][rg][na] += ttt_diag_all_vec[ag][0][im][nm][rg][na];
                VLdx[ag][1][0][im][nm][rg][na] += ttt_diag_all_vec[ag][1][im][nm][rg][na];
                ///// Move to latent tx experienced /////////////////////////////////
                V1[ag][0][1][im][nm][rg][na]  += ttt_diag_all_vec[ag][0][im][nm][rg][na];
                V1[ag][1][1][im][nm][rg][na]  += ttt_diag_all_vec[ag][1][im][nm][rg][na];
                ///////////////////////////////////////////////////////////////////////////////
                ///// remove those who test positive who are true LTBI positive (true positives)
                ///////////////////////////////////////////////////////////////////////////////
                ///// Latent slow
                temp  = ttt_diag_all_vec[ag][2][im][nm][rg][na] * ttt_ltbi_initN * ttt_ltbi_compN; // tx completion
                temp3 = ttt_diag_all_vec[ag][2][im][nm][rg][na] * ttt_ltbi_initN * (1-ttt_ltbi_compN); // default
                ///// Update the number of diagnoses
                VLdx[ag][2][0][im][nm][rg][na] += ttt_diag_all_vec[ag][2][im][nm][rg][na];
                ///// Latent fast
                temp2 = ttt_diag_all_vec[ag][3][im][nm][rg][na]*ttt_ltbi_initN*ttt_ltbi_compN; // tx completion
                temp4 = ttt_diag_all_vec[ag][3][im][nm][rg][na]*ttt_ltbi_initN*(1-ttt_ltbi_compN); // default
                ///////////////////////////////////////////////////////////////////////////////
                ///// Update the number of diagnoses
                VLdx[ag][3][0][im][nm][rg][na] += ttt_diag_all_vec[ag][3][im][nm][rg][na];
                ///////////////////////////////////////////////////////////////////////////////
                V1[ag][2][0][im][nm][rg][na]  -=  (temp + temp3); //remove from latent slow
                V1[ag][3][0][im][nm][rg][na]  -=  (temp2 + temp4);  //remove from latent fast
                //completion split between success and failure
                V1[ag][1][1][im][nm][rg][na]  += (temp + temp2) * ttt_ltbi_effN;  //exit to cure
                V1[ag][2][1][im][nm][rg][na]  += (temp + temp2) * (1-ttt_ltbi_effN); //tx comp fail to latent slow
                ///defaults are placed in tx naive because it is considered the same tb infection
                V1[ag][2][0][im][nm][rg][na]  += (temp3 + temp4); //latent tx default to latent slow

                ///////////////////////////////////////////////////////////////////////////////
                ///// remove incidental tb disease cases found in the ttt intervention
                ///////////////////////////////////////////////////////////////////////////////
                V1 [ag][4][0][im][nm][rg][na]     -= ttt_diag_all_vec[ag][4][im][nm][rg][na];
                V1 [ag][5][0][im][nm][rg][na]     += ttt_diag_all_vec[ag][4][im][nm][rg][na];
                Vdx[ag][4][0][im][nm][rg][na]     += ttt_diag_all_vec[ag][4][im][nm][rg][na];

                // Sum of total screened population among latent for print below
                temp21 += ttt_test_all_vec[ag][2][im][nm][rg][na] + ttt_test_all_vec[ag][3][im][nm][rg][na];
                // Sum of total diagnosed population among latent for print below
                temp23 += ttt_diag_all_vec[ag][2][im][nm][rg][na] + ttt_diag_all_vec[ag][3][im][nm][rg][na];
                // Sum of total TB diagnosed population among latent for print below
                temp24 +=  VLdx[ag][2][0][im][nm][rg][na] + VLdx[ag][3][0][im][nm][rg][na];
} } } } }


            // this is for understanding the characteristics of the screened population
            // for(int ag=0; ag<11; ag++) {
            //   for(int im=0; im<4 ; im++) {
            //     for(int nm=0; nm<4; nm++) {
            //       for(int na=0; na<3; na++) {
            //         if (na < 1) {
            //           ni = 0;
            //         } else{
            //           ni = 1;
            //         }
            //         for (int tb = 0; tb < 5; tb++){
            //           for(int rg=0; rg<2; rg++) {
            //             ttt_frc_scrn = std::min(ttt_samp_distN[0][ag+(ni*11)+(nm*22)+(im*88)], 1.0);
            //           // VLtemp[ag][tb][0][im][nm][rg][na] = std::max((ttt_frc_scrn * std::max(V0[ag][tb][0][im][nm][rg][na] , 0.0) * ttt_ltbi_acceptN),0.0);
            //         }}}}}}

          } // close the screening month loop
/// print summary of screening statistics
          if (s == ttt_month.front()){
            Rcpp::Rcout<< "month = " << s << "\n";
            // Rcpp::Rcout<< "total population = " << temp19*1e6 << "\n";
            Rcpp::Rcout<< "total extra screening = " << temp20*1e6 << "\n";
            Rcpp::Rcout<< "total extra screening in latent pop = " << temp21 *1e6<< "\n";
            Rcpp::Rcout<< "total extra LTBI diagnoses = " << temp22*1e6 << "\n";
            Rcpp::Rcout<< "total extra LTBI diagnoses in latent pop = " << temp23*1e6 << "\n";
            Rcpp::Rcout<< "total vldx diagnoses = " << temp25*1e6 << "\n";
            Rcpp::Rcout<< "total extra vldx diagnoses = " << temp24*1e6 << "\n";
          }
          temp7=0;temp8=0;temp9=0;temp10=0; temp20=0; temp21=0; temp22=0; temp23=0; temp19=0; temp24=0;
          ///////////////////// TB DIAGNOSIS AND TX INITIATION  /////////////////////////
          temp=0;
          for(int ag=0; ag<11; ag++) {
            for(int lt=0; lt<2; lt++) {
              for(int im=0; im<4 ; im++) {
                for(int nm=0; nm<4; nm++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      if (lt == 0){
                        temp  = std::max(((V0[ag][4][lt][im][nm][rg][na] - VAtemp[ag][im][nm][rg][na]) * rDxtN[s][rg] / RRdxAge[ag]), 0.0);
                      } else {
                        temp  = std::max((V0[ag][4][lt][im][nm][rg][na] * rDxtN[s][rg] / RRdxAge[ag]), 0.0);
                      }
                      // if (temp == 0.0) { Rcpp::Rcout << "Not enough RegDiagAct in nat = " << na << " ag = " << ag << " rg = " << rg << " im = " << im << " nm = "<< nm << "\n";}

                      V1[ag][4][lt][im][nm][rg][na]  -= temp;
                      V1[ag][5][lt][im][nm][rg][na]  += temp;
                      Vdx[ag][4][lt][im][nm][rg][na] += temp;
                    } } } } } }

          for(int ag=0; ag<11; ag++) {
            for(int im=0; im<4 ; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    VAtemp[ag][im][nm][rg][na] = 0;
                  } } } } }
          //////////////////////// TB TREATMENT OUTCOMES /////////////////////////////
          temp=0;
          for(int ag=0; ag<11; ag++) {
            for(int lt=0; lt<2; lt++) {
              for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      ///// cure back to latent slow state
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
      /////////////////////////////////////////////////////////////////////////////
      /////                        FILL RESULTS TABLE                         /////
      /////////////////////////////////////////////////////////////////////////////
      //// start with those which can be scaled with the mid-month value     //////
      if(m == 6) {
        /////////////////////////////////// YEAR //////////////////////////////////
        Outputs[y][0]      = y+1950;
        ////////////////    COUNTS BY TOTAL, AGE, TB, RF, AND RG     //////////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      Outputs[y][1    ] += V0[ag][tb][lt][im][nm][rg][na];   // N_ALL
                      Outputs[y][2 +ag] += V0[ag][tb][lt][im][nm][rg][na];   // N_ by age (11)
                      Outputs[y][13+tb] += V0[ag][tb][lt][im][nm][rg][na];   // N_ by tb (6)
                      Outputs[y][19+im] += V0[ag][tb][lt][im][nm][rg][na];   // N_ by im (4)
                      Outputs[y][23+nm] += V0[ag][tb][lt][im][nm][rg][na];   // N_ by nm (4)
                      Outputs[y][27+rg] += V0[ag][tb][lt][im][nm][rg][na];   // N_ by rg (2)
                      Outputs[y][29+na] += V0[ag][tb][lt][im][nm][rg][na];   // N_ by na (3)
                    } } } } } } }
        ////////////////////    COUNTS BY NATIVITY AND AGE    ////////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                  for(int rg=0; rg<2; rg++) {
                    Outputs[y][32+ag] += V0[ag][tb][lt][im][nm][rg][0]; // N_ by age and US (11)
                    Outputs[y][43+ag] += V0[ag][tb][lt][im][nm][rg][1]+V0[ag][tb][lt][im][nm][rg][2];   // N_ by age and FB (11)
                  } } } } } }
        /////////////////    LATENT COUNTS BY NATIVITY, AGE, & LTBI    //////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  Outputs[y][54+ag] += (V0[ag][1][lt][im][nm][rg][0]*(1-pImmScen))+
                    V0[ag][2][lt][im][nm][rg][0]+V0[ag][3][lt][im][nm][rg][0];   // N_ by age and US (11) LATENT INFECTION
                  Outputs[y][65+ag] += ((V0[ag][1][lt][im][nm][rg][1]+V0[ag][1][lt][im][nm][rg][2])*(1-pImmScen))+
                    V0[ag][2][lt][im][nm][rg][1]+V0[ag][2][lt][im][nm][rg][2]+
                    V0[ag][3][lt][im][nm][rg][1]+V0[ag][3][lt][im][nm][rg][2]; // N_ by age and FB (11) LATENT INFECTION
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
                      Outputs[y][76+ag] += V0[ag][tb][lt][im][nm][rg][na];   // N_RF by age (11)
                    } } } } } } }
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
                  } } } } } }
        ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
        for(int i=131; i<133; i++) { Outputs[y][i] = Outputs[y][i]*12; }
        ///////////////////           FORCE OF INFECTION        //////////////////////
        for(int ag=0; ag<11; ag++) {
          Outputs[y][243] += VLkla[0][0][ag];
          Outputs[y][244] += VLkla[1][0][ag];
          Outputs[y][245] += VLkla[0][1][ag];
          Outputs[y][246] += VLkla[1][1][ag];
        }
        ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
        for(int i=243; i<247; i++) { Outputs[y][i] = Outputs[y][i]*12; }
        ///////////////////////          POPULATION          /////////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      if (na<1){
                        Outputs[y][296+im]  += V0[ag][tb][lt][im][nm][rg][na];
                        Outputs[y][304+rg]  += V0[ag][tb][lt][im][nm][rg][na];
                        Outputs[y][308+nm]  += V0[ag][tb][lt][im][nm][rg][na];
                      } else {
                        Outputs[y][300+im]  += V0[ag][tb][lt][im][nm][rg][na];
                        Outputs[y][306+rg]  += V0[ag][tb][lt][im][nm][rg][na];
                        Outputs[y][312+nm]  += V0[ag][tb][lt][im][nm][rg][na];

                      }
                    } } } } } } }
        ////////////////           POPULATION BY IM AND NM            ////////////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for(int nm=0; nm<4; nm++) {
                for(int im=0; im<4; im++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      Outputs[y][317+nm+(im*4)] += V0[ag][tb][lt][im][nm][rg][na];
                    } } } } } } }
        ////////////////        POPULATION BY AGE, IM AND NM          ////////////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      Outputs[y][333+nm+(im*4)+(ag*16)] += V0[ag][tb][lt][im][nm][rg][na];
                    } } } } } } }
        ////////////////           POPULATION BY AGE AND NM           ////////////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for(int nm=0; nm<4; nm++) {
                for(int im=0; im<4; im++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      Outputs[y][520+nm+(ag*4)] += V0[ag][tb][lt][im][nm][rg][na];
                    } } } } } } }
        ////////////////      NEW INFECTIONS + SUPER INFECTION        ////////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2 ; lt++) {
            for(int im=0; im<4 ; im++) {
              for(int nm=0; nm<4 ; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++){
                    Outputs[y][247+rg] += (V0[ag][0][lt][im][nm][rg][na]+V0[ag][1][lt][im][nm][rg][na]+V0[ag][2][lt][im][nm][rg][na])*VLkla[rg][na][ag]*NixTrans[s];
                    Outputs[y][249+na] += (V0[ag][0][lt][im][nm][rg][na]+V0[ag][1][lt][im][nm][rg][na]+V0[ag][2][lt][im][nm][rg][na])*VLkla[rg][na][ag]*NixTrans[s];
                  } } } } } }
        ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
        for(int i=247; i<252; i++) { Outputs[y][i] = Outputs[y][i]*12; }
        ////////////////   POPULATION BY NATIVITY, AGE, IM AND NM     ////////////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<5; tb++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    if (na<1){
                      if (ag <3){
                        Outputs[y][586+nm+(im*4)] += (V1[ag][tb][0][im][nm][rg][na] - V0[ag][4][0][im][nm][rg][na]*rDxtN[s][rg]/RRdxAge[ag]);
                      } if(2<ag && ag<7){
                        Outputs[y][602+nm+(im*4)] += (V1[ag][tb][0][im][nm][rg][na] - V0[ag][4][0][im][nm][rg][na]*rDxtN[s][rg]/RRdxAge[ag]);
                      } if (6<ag && ag<11){
                        Outputs[y][618+nm+(im*4)] += (V1[ag][tb][0][im][nm][rg][na] - V0[ag][4][0][im][nm][rg][na]*rDxtN[s][rg]/RRdxAge[ag]);
                      }
                    } else {
                      if (ag<3) {
                        Outputs[y][634+nm+(im*4)] += (V1[ag][tb][0][im][nm][rg][na] - V0[ag][4][0][im][nm][rg][na]*rDxtN[s][rg]/RRdxAge[ag]);
                      } if(2<ag && ag<7){
                        Outputs[y][650+nm+(im*4)] += (V1[ag][tb][0][im][nm][rg][na] - V0[ag][4][0][im][nm][rg][na]*rDxtN[s][rg]/RRdxAge[ag]);
                      } if (6<ag && ag<11){
                        Outputs[y][666+nm+(im*4)] += (V1[ag][tb][0][im][nm][rg][na] - V0[ag][4][0][im][nm][rg][na]*rDxtN[s][rg]/RRdxAge[ag]);
                      } } }
                } } } } }
        ////////////////      COUNTS BY LATENT TREATMENT STATUS     //////////////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    Outputs[y][682] += V0[ag][tb][0][im][nm][rg][na];   //TREATMENT NAIVE
                  } } } } } }
        ////////////////        LTBI TREATMENT NAIVE BY NATIVITY      ////////////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    Outputs[y][815+na] += V0[ag][tb][0][im][nm][rg][na];   //TREATMENT NAIVE
                  } } } } } }

        /////////////////    LATENT COUNTS BY NATIVITY & LATENT TX STATUS    //////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  Outputs[y][826 + lt]  += (V0[ag][1][lt][im][nm][rg][0]*(1-pImmScen))+
                    V0[ag][2][lt][im][nm][rg][0]+V0[ag][3][lt][im][nm][rg][0];   // N_ by age and US (11) LATENT INFECTION
                  Outputs[y][828 + lt] += ((V0[ag][1][lt][im][nm][rg][1]+V0[ag][1][lt][im][nm][rg][2])*(1-pImmScen))+
                    V0[ag][2][lt][im][nm][rg][1]+V0[ag][2][lt][im][nm][rg][2]+
                    V0[ag][3][lt][im][nm][rg][1]+V0[ag][3][lt][im][nm][rg][2]; // N_ by age and FB (11) LATENT INFECTION
                } } } } }

        /////////////////    LATENT COUNTS BY NATIVITY & LATENT TX STATUS    //////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int tb=0; tb<4; tb++) {
                    for(int na=0; na<3; na++) {
                      Outputs[y][830 + rg + (ag * 2)] +=  V0[ag][tb][lt][im][nm][rg][na];
                    } } } } } } }

        /////////////////    LATENT COUNTS BY NUSB without latent0  //////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  Outputs[y][857] +=
                    V0[ag][2][lt][im][nm][rg][1]+V0[ag][2][lt][im][nm][rg][2]+
                    V0[ag][3][lt][im][nm][rg][1]+V0[ag][3][lt][im][nm][rg][2]; // N_ by and FB (11) LATENT INFECTION
                } } } } }

        /////////////    Latent Tx Naive by Nativity and Age  /////////////////
        for(int ag=0; ag<11; ag++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                  if (na < 1){
                    Outputs[y][858 + ag] +=
                    V0[ag][0][0][im][nm][rg][0] +
                    V0[ag][1][0][im][nm][rg][0]+V0[ag][2][0][im][nm][rg][0]+
                    V0[ag][3][0][im][nm][rg][0]+V0[ag][4][0][im][nm][rg][0];
                  } else{
                    Outputs[y][869 + ag] +=
                      V0[ag][0][0][im][nm][rg][1]+V0[ag][0][0][im][nm][rg][2]+
                      V0[ag][1][0][im][nm][rg][1]+V0[ag][1][0][im][nm][rg][2]+
                      V0[ag][2][0][im][nm][rg][1]+V0[ag][2][0][im][nm][rg][2]+
                      V0[ag][3][0][im][nm][rg][1]+V0[ag][3][0][im][nm][rg][2] +
                      V0[ag][4][0][im][nm][rg][1]+V0[ag][4][0][im][nm][rg][2];
                  }
                } } } } }

        /////////////////    LATENT COUNTS BY NATIVITY & LATENT TX STATUS    //////////////////////
        // for(int ag=0; ag<11; ag++) {
        //   for(int lt=0; lt<2; lt++) {
        //     for(int im=0; im<4; im++) {
        //       for(int nm=0; nm<4; nm++) {
        //         for(int rg=0; rg<2; rg++) {
        //           Outputs[y][826 + lt]  += (V0[ag][1][lt][im][nm][rg][0])*(1-pImmScen)+
        //             V0[ag][2][lt][im][nm][rg][0]+V0[ag][3][lt][im][nm][rg][0];   // N_ by age and US (11) LATENT INFECTION
        //           Outputs[y][828 + lt] += (V0[ag][1][lt][im][nm][rg][1]+V0[ag][1][lt][im][nm][rg][2])*(1-pImmScen)+
        //             V0[ag][2][lt][im][nm][rg][1]+V0[ag][2][lt][im][nm][rg][2]+
        //             V0[ag][3][lt][im][nm][rg][1]+V0[ag][3][lt][im][nm][rg][2]; // N_ by age and FB (11) LATENT INFECTION
        //         } } } } }

        ///////////////////           POPS for FORCE OF INFECTION        //////////////////////
        // Outputs[y][858] = VNkl[0][0];
        // Outputs[y][859] = VNkl[1][0];
        // Outputs[y][860] = VNkl[0][1];
        // Outputs[y][861] = VNkl[1][1];

        ////////////////   POPULATION BY NATIVITY, AGE, IM AND NM     ////////////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<5; tb++) {
            for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    if (na < 1){
                        // if (tb == 4){
                        //   Outputs[y][880+nm+im*4+ag*16] += (V0[ag][4][0][im][nm][rg][na] - V0[ag][4][0][im][nm][rg][na]*rDxtN[s][rg]/RRdxAge[ag]);
                        // } else{
                          Outputs[y][880+nm+im*4+ag*16] += (V1[ag][tb][0][im][nm][rg][na]);
                          // }
                    } else {
                      // if (tb == 4){
                      //   Outputs[y][1056+nm+im*4+ag*16] += (V0[ag][tb][0][im][nm][rg][na] - V0[ag][4][0][im][nm][rg][na]*rDxtN[s][rg]/RRdxAge[ag]);
                      //  } else{
                         Outputs[y][1056+nm+im*4+ag*16] += (V1[ag][tb][0][im][nm][rg][na]);
                       // }
                    }
                } } } } } }
      } ////end of mid-year results bracket
      /////////////////////////////////////////////////////////////////////////////
      /////                       END MIDYEAR RESULTS                         /////
      /////////////////////////////////////////////////////////////////////////////
      ////// ADD IN THE RESULTS THAT ARE ACCUMULATED ACROSS THE YEAR OF THE ///////
      ////// MODEL RUN. THESE ARE TB SPECIFIC RESULTS SUCH AS LTBI TESTS,   ///////
      ////// DIAGNOSES, TX INITIATIONS, TX COMPLETIONS, TB DIAGNOSES AND    ///////
      ////// INFECTIONS, AND MORTALITY Outputs. ///////////////////////////////////
      if(m == res_m) {
        ///////////// TB MORTALITY COUNT BY AGE, NATIVITY///////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    if(na>0) {
                      ti = 11;
                    } else { ti = 0; }
                    Outputs[y][87+ag+ti]  += VMort[ag][4][lt][im][nm][rg][na];
                    Outputs[y][87+ag+ti]  += VMort[ag][5][lt][im][nm][rg][na];
                  } } } } } }
        ///////////////////////  RISK FACTOR MORTALITY BY AGE /////////////////////////
        ///this is horribly wrong, using it as a place holder
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      Outputs[y][109+ag]  += VMort[ag][tb][lt][im][nm][rg][na]/RRmuRFN[nm];
                    } } } } } } }
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
        ///////////////////////     TB TREATMENT OUTCOMES    /////////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    Outputs[y][133]  += VMort[ag][5 ][lt][im][nm][rg][na];// tx mort (not scaled)
                  } } } } } }
        ///////////////////////        TB NOTIFICATIONS      /////////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    Outputs[y][134   ] += Vdx[ag][4][lt][im][nm][rg][na];   // All dx (1)
                    Outputs[y][135+ag] += Vdx[ag][4][lt][im][nm][rg][na];   // dx by age (11)
                    Outputs[y][146+na] += Vdx[ag][4][lt][im][nm][rg][na];   // dx by na (3)
                    Outputs[y][149+rg] += Vdx[ag][4][lt][im][nm][rg][na];   // N_ by rg (2)
                  } } } } } }
        ///////////////////////    LTBI TREATMENT INITIATIONS   //////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++) {
              for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                  Outputs[y][151] += (VLdx[ag][2][0][im][nm][rg][na]+VLdx[ag][3][0][im][nm][rg][na])*LtTxParN[s][0] +
                    ((V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*rTbN*LtTxParN[s][0]);//all inits (((1- pop_frc)*rTbN) + (pop_frc*(1-(rTbP*rr_ltbi))))*LtTxParN[s][0]; //all init
                  if(na>0) {
                    Outputs[y][152] += (VLdx[ag][2][0][im][nm][rg][na]+VLdx[ag][3][0][im][nm][rg][na])*LtTxParN[s][0] +
                      ((V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*rTbN*LtTxParN[s][0]);} // FB inits
                  if(rg==1) {
                    Outputs[y][153] +=  (VLdx[ag][2][0][im][nm][rg][na]+VLdx[ag][3][0][im][nm][rg][na])*LtTxParN[s][0] +
                      ((V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*rTbN*LtTxParN[s][0]); } // high risk inits

                  Outputs[y][154] += (VLdx[ag][2][0][im][nm][rg][na]+VLdx[ag][3][0][im][nm][rg][na])*LtTxParN[s][0]; // inits with LTBI
                  // Calculate the number of LTBI Tests, Tx Initiations, and Tx Completions
                  if(na==0){
                    //// LTBI tests
                    Outputs[y][683+ag] +=  VLtest[ag][0][0][im][nm][rg][na]+VLtest[ag][1][0][im][nm][rg][na]+VLtest[ag][2][0][im][nm][rg][na]+VLtest[ag][3][0][im][nm][rg][na];
                    /// //// LTBI tests among true positive
                    Outputs[y][793+ag] +=  VLtest[ag][2][0][im][nm][rg][na]+VLtest[ag][3][0][im][nm][rg][na] + (VLtest[ag][1][0][im][nm][rg][na]*(1-pImmScen));
                    //// LTBI tests by progression risk
                    Outputs[y][818+im] +=  VLtest[ag][0][0][im][nm][rg][na]+VLtest[ag][1][0][im][nm][rg][na] + VLtest[ag][2][0][im][nm][rg][na] + VLtest[ag][3][0][im][nm][rg][na];

                    ///// LTBI treatment initiations
                    Outputs[y][705+ag] +=  (VLdx[ag][0][0][im][nm][rg][na] + VLdx[ag][1][0][im][nm][rg][na] + VLdx[ag][2][0][im][nm][rg][na] + VLdx[ag][3][0][im][nm][rg][na])*LtTxParN[s][0];
                    ///// LTBI treatment initiations among true positives
                    Outputs[y][1232+ag] += ((VLdx[ag][1][0][im][nm][rg][na]*(1-pImmScen)) + VLdx[ag][2][0][im][nm][rg][na] + VLdx[ag][3][0][im][nm][rg][na])*LtTxParN[s][0];

                    ///// LTBI treatment completions
                    Outputs[y][727+ag] +=  (VLdx[ag][0][0][im][nm][rg][na] + VLdx[ag][1][0][im][nm][rg][na] + VLdx[ag][2][0][im][nm][rg][na] + VLdx[ag][3][0][im][nm][rg][na])*LtTxParN[s][0]*(1-LtTxParN[s][1]) ;

                    ///// LTBI treatment completions among true positives
                    Outputs[y][1254+ag] +=  ((VLdx[ag][1][0][im][nm][rg][na]*(1-pImmScen)) + VLdx[ag][2][0][im][nm][rg][na] + VLdx[ag][3][0][im][nm][rg][na])*LtTxParN[s][0]*(1-LtTxParN[s][1]) ;
                  } else {
                    /// For the non-USB calculations, we need to allow for the addition of tests for the Intervention #1 scenario
                    /// in which every migrant is tested for LTBI prior to entry to the United States
                    Outputs[y][694+ag] += VLtest[ag][0][0][im][nm][rg][na] + VLtest[ag][1][0][im][nm][rg][na] +VLtest[ag][2][0][im][nm][rg][na] + VLtest[ag][3][0][im][nm][rg][na];
                    /////// LTBI tests among true positive
                    Outputs[y][804+ag] += VLtest[ag][2][0][im][nm][rg][na] + VLtest[ag][3][0][im][nm][rg][na]+(VLtest[ag][1][0][im][nm][rg][na]*(1-pImmScen));

                    Outputs[y][822+im] +=  VLtest[ag][0][0][im][nm][rg][na]+VLtest[ag][1][0][im][nm][rg][na]+VLtest[ag][2][0][im][nm][rg][na]+VLtest[ag][3][0][im][nm][rg][na];

                    ///// LTBI treatment initiations
                    Outputs[y][716+ag] += (VLdx[ag][0][0][im][nm][rg][na] + VLdx[ag][1][0][im][nm][rg][na] + VLdx[ag][2][0][im][nm][rg][na] + VLdx[ag][3][0][im][nm][rg][na])*LtTxParN[s][0];
                    ///// LTBI treatment initiations among true positives
                    Outputs[y][1243+ag] += ((VLdx[ag][1][0][im][nm][rg][na]*(1-pImmScen)) + VLdx[ag][2][0][im][nm][rg][na] + VLdx[ag][3][0][im][nm][rg][na])*LtTxParN[s][0];

                    ///// LTBI treatment completions
                    Outputs[y][738+ag] += (VLdx[ag][0][0][im][nm][rg][na] + VLdx[ag][1][0][im][nm][rg][na] + VLdx[ag][2][0][im][nm][rg][na] + VLdx[ag][3][0][im][nm][rg][na])*LtTxParN[s][0]*(1-LtTxParN[s][1]) ;
                    ///// LTBI treatment completions among true positives
                    Outputs[y][1265+ag] +=  ((VLdx[ag][1][0][im][nm][rg][na]*(1-pImmScen)) + VLdx[ag][2][0][im][nm][rg][na] + VLdx[ag][3][0][im][nm][rg][na])*LtTxParN[s][0]*(1-LtTxParN[s][1]) ;
                  }
                } } } }
          // Add in the additional tests for Intervention 1; If Int1 is not active, this adds a zero to every cell
          // come back and add a more refined way of doing this.
          Outputs[y][694+ag] += Int1TestN[y][ag];
          Outputs[y][716+ag] += Int1InitN[y][ag];
          Outputs[y][738+ag] += Int1TxN[y][ag];
        }
        ///////////////////      TB INCIDENCE, BY ALL VS RECENT      /////////////////
        // By recency (<2 years) == all immediate, 1-(1-rfast)^24 x all Lf
        // Break down from latent fast and slow
        for(int rg=0; rg<2; rg++)  {
          for(int ag=0; ag<11; ag++) {
            for(int im=0; im<4; im++) {
              temp4V[ag][im][rg] = (1-pow(1-(MrslowN[ag][im])-rRecov,24.0-(1.0/rDxtN[s][rg])))*(MrslowN[ag][im]);
            } }
          temp3V[rg] = (1-pow(1-rfast,24.0-(1.0/rDxtN[s][rg])))*rfast;
        }
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    temp2 = (Vinf[ag][0][lt][im][nm][rg][na] + Vinf[ag][1][lt][im][nm][rg][na])*((MpfastN[ag][im]*(1-pow(1-rfast,24.0)) + (1-MpfastN[ag][im])*(1-pow(1-MrslowN[ag][im],24.0))));// Progression from recent infection
                    temp =  Vinf[ag][2 ][lt][im][nm][rg][na]*(MrslowN[ag][im]) + V0[ag][3 ][lt][im][nm][rg][na]*rfast; // All progression
                    Outputs[y][155      ] += temp ;   // all incidence
                    Outputs[y][155+16   ] += temp2;   // all incidence, recent infection
                    Outputs[y][156+ag   ] += temp ;   // incidence by age
                    Outputs[y][156+ag+16] += temp2;   // incidence by age, recent infection
                    if(na<1) {
                      Outputs[y][167      ] += temp ;   //  incidence, US born
                      Outputs[y][167+16   ] += temp2;   //  incidence, US born, recent infection
                    } else {
                      Outputs[y][168      ] += temp ;   //  incidence, FB
                      Outputs[y][168+16   ] += temp2;   //  incidence, FB, recent infection
                    }
                    if(na==2) {
                      Outputs[y][169      ] += temp ;   //  incidence, FB2 born
                      Outputs[y][169+16   ] += temp2;    //  incidence, FB2 born, recent infection
                    }
                    if(rg==1) {
                      Outputs[y][170      ] += temp ;   //  incidence, HR
                      Outputs[y][170+16   ] += temp2;   //  incidence, HR, recent infection
                    }
                    // if(im>0) {
                    //   Outputs[y][171      ] += temp ;   //  incidence, HIV pos
                    //   Outputs[y][171+16   ] += temp2; } //  incidence, HIV pos, recent infection
                  } } } } } }
        /////////////////// TB NOTIFICATIONS (DEAD AT DIAGNOSIS) /////////////////////
        for(int nm=0; nm<4 ; nm++) {
          for(int ag=0; ag<11; ag++) {
            for(int lt=0; lt<2; lt++){
              for (int im=0; im<4; im++){
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++){
                    Outputs[y][187   ] += VMort[ag][4 ][lt][im][nm][rg][na];   // All dx (1)
                    Outputs[y][188+ag] += VMort[ag][4 ][lt][im][nm][rg][na];   // dx by age (11)
                    Outputs[y][199+na] += VMort[ag][4 ][lt][im][nm][rg][na];
                    // if(im>0) {
                    //   Outputs[y][228   ] += temp2;   // dx HIV pos (1)
                    // } else {
                    //   Outputs[y][229   ] += temp2; }  // dx HIV neg (1)
                    Outputs[y][202+rg] += VMort[ag][4 ][lt][im][nm][rg][na];   // N_ by rg (2)
                  } } } } } }
        ///////////////////      TB NOTIFICATIONS (US BORN)      /////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++){
            for (int im=0; im<4; im++){
              for (int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++) {
                  Outputs[y][204+ag] += Vdx[ag][4][lt][im][nm][rg][0];   // dx by age (11)
                } } } } }
        /////////////// TB NOTIFICATIONS (US BORN, DEAD AT DIAGNOSIS) ///////////////
        for(int nm=0; nm<4 ; nm++) {
          for (int im=0; im<4; im++){
            for(int ag=0; ag<11; ag++) {
              for(int lt=0; lt<2; lt++){
                for(int rg=0; rg<2; rg++) {
                  Outputs[y][215+ag] += VMort[ag][4 ][lt][im][nm][rg][0];   // dx by age (11)
                } } } } }
        /////////////////////  TOTAL MORTALITY BY AGE, HAVE TB   /////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++){
            for (int im=0; im<4; im++){
              for (int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++){
                    for(int tb=4; tb<6; tb++){
                      Outputs[y][226+ag]  += VMort[ag][tb][lt][im][nm][rg][na];
                    } } } } } } }
        ///////////////////        TB COUNTS (BY NATIVITY)       /////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++){
            for (int im=0; im<4; im++){
              for (int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++){
                    if(na<1) {
                      Outputs[y][237] += V0[ag][2][lt][im][nm][rg][na];
                      Outputs[y][238] += V0[ag][3][lt][im][nm][rg][na];
                      Outputs[y][239] += V0[ag][4][lt][im][nm][rg][na];
                    } else {
                      Outputs[y][240] += V0[ag][2][lt][im][nm][rg][na];
                      Outputs[y][241] += V0[ag][3][lt][im][nm][rg][na];
                      Outputs[y][242] += V0[ag][4][lt][im][nm][rg][na];
                    } } } } } } }
        ////////////////////// TB MORTALITY BY NATIVITY //////////////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++){
            for (int im=0; im<4; im++){
              for (int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++){
                    if(na > 0 ) {
                      ti = 1;
                    } else { ti = 0; }
                    // if(im > 2) { temp = muTbRF;
                    // } else { temp = 0; }
                    temp=0;
                    Outputs[y][252+ti]  += VMort[ag][4][lt][im][nm][rg][na];
                    Outputs[y][252+ti]  += VMort[ag][5][lt][im][nm][rg][na];
                  } } } } } }
        ///////////////////////    TOTAL MORTALITY BY AGE    /////////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      if (na < 1){
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
        ///////////////////////         TOTAL MORTALITY      /////////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      Outputs[y][316]  += VMort[ag][tb][lt][im][nm][rg][na];
                    } } } } } } }
        ////////////////             MORTALITY RATE BY AGE            ////////////////
        for(int ag=0; ag<11; ag++) {
          Outputs[y][509+ag] = (Outputs[y][120+ag]/12)/(Outputs[y][2+ag]);
        }
        ///////////  NEW INFECTIONS + SUPER INFECTION BY AGE AND NAT GROUP ///////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2 ; lt++) {
              for(int im=0; im<4 ; im++) {
                for(int nm=0; nm<4 ; nm++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++){
                      if (na <1){
                        Outputs[y][564+ag] += Vinf[ag][tb][lt][im][nm][rg][na];
                      } else {
                        Outputs[y][575+ag] += Vinf[ag][tb][lt][im][nm][rg][na];
                      }
                    } } } } } } }
        ////////////////        TB DIAGNOSES BY NATIVITY AND AGE       ///////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    if(na==0){
                      Outputs[y][749+ag] += Vdx[ag][4][lt][im][nm][rg][na]; //how to handle TB tests among non-TB persons
                      Outputs[y][771+ag] += Vdx[ag][4][lt][im][nm][rg][na]*TxVec[0]; //check with TxVec[5]
                    } else {
                      Outputs[y][760+ag] += Vdx[ag][4][lt][im][nm][rg][na];
                      Outputs[y][782+ag] += Vdx[ag][4][lt][im][nm][rg][na]*TxVec[0];
                    }
                  } } } } } }
        ///////////////////        TB COUNTS (BY NATIVITY)       /////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++){
            for (int im=0; im<4; im++){
              for (int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++){
                    if(na<1) {
                      Outputs[y][852] += V0[ag][0][lt][im][nm][rg][na];
                      Outputs[y][853] += V0[ag][1][lt][im][nm][rg][na];
                    } else {
                      Outputs[y][854] += V0[ag][0][lt][im][nm][rg][na];
                      Outputs[y][855] += V0[ag][1][lt][im][nm][rg][na];
                    } } } } } } }
        //// LTBI TREATMENT INITIATIONS with LATENT ZERO INCLUDED ///////////////////////////
        for(int ag=0; ag<11; ag++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++) {
              for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                  Outputs[y][856] += ((VLdx[ag][1][0][im][nm][rg][na]*(1-pImmScen)) + VLdx[ag][2][0][im][nm][rg][na]+VLdx[ag][3][0][im][nm][rg][na])*LtTxParN[s][0]; // inits with LTBI
                } } } } }
        ////////////////      RESET THE CUMULATIVE ANNUAL COUNTS      ////////////////
        for(int ag=0; ag<11; ag++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++) {
              for(int nm=0; nm<4; nm++) {
                for(int im=0; im<4; im++) {
                  for(int lt=0; lt<2 ; lt++) {
                    for(int tb=0; tb<6; tb++) {
                      VMort [ag][tb][lt][im][nm][rg][na] = 0;
                      Vdx   [ag][tb][lt][im][nm][rg][na] = 0;
                      VLdx  [ag][tb][lt][im][nm][rg][na] = 0;
                      VLtest[ag][tb][lt][im][nm][rg][na] = 0;
                      VLtemp[ag][tb][lt][im][nm][rg][na] = 0;
                      VLDxtemp[ag][tb][lt][im][nm][rg][na] = 0;
                      Vinf[ag][tb][lt][im][nm][rg][na]  = 0;
                    } } } } } } }
      } ////end of end-year results bracket
      /////////////////////////////////////////////////////////////////////////////
      /////                   END END OF YEAR RESULTS                         /////
      /////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////
      ////////////////////// REBALANCE THE POPULATION /////////////////////////
      /////////////////////////////////////////////////////////////////////////
      ///// Define the current distribution of persons across the risk groups
      ///// [nm] and [im] at this timestep
      if (reblnc==1){
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for (int im=0; im<4; im++){
                for (int nm=0; nm<4; nm++){
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++){
                      // V2[ag][tb][0][im][nm][rg][na]=0;
                      for (int m2=0; m2<4; m2++){
                        for (int p2=0; p2<4; p2++){
                          temp = V1[ag][tb][lt][p2][m2][rg][na]*(trans_mat_tot_agesN[m2+p2*4][(16*(ag+1))-(16-(nm+im*4))]);
                          V1[ag][tb][lt][im][nm][rg][na] += temp;

                          V1[ag][tb][lt][p2][m2][rg][na] -= temp;
                        } } } } } } } } }
      } //end of rebalancing loop
      for(int ag=0; ag<11; ag++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++){
            for(int rg=0; rg<2; rg++){
              for (int na=0; na<3; na++){
                for(int lt=0; lt<2; lt++){
                  for(int tb=0; tb<6; tb++) {
                    V0[ag][tb][lt][im][nm][rg][na] = V1[ag][tb][lt][im][nm][rg][na];
                  } } } } } } }

    } //// end of month loop!//////////////////////////////////////////////////////////
  } //// end of year loop!///////////////////////////////////////////////////////////

  for (int i=0; i<4; i++){
    for (int j=0; j<4; j++){
      temp_mat[i][j]=0; } }
  mat_sum=0;
  for(int nm=0; nm<4; nm++){
    for(int ag=0; ag<11; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++){
          for(int im=0; im<4; im++) {
            for(int rg=0; rg<2; rg++){
              for(int na=0; na<3; na++){
                temp_mat[nm][im]  += V1[ag][tb][lt][im][nm][rg][na];
              } } } } } }
  }
  for(int nm=0; nm<4; nm++){
    for(int im=0; im<4; im++){
      mat_sum+=temp_mat[nm][im];
    } }
  for(int nm=0; nm<4; nm++){
    for(int im=0; im<4; im++){

      dist_mat(nm,im) = temp_mat[nm][im]/mat_sum;
    } }
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
      Rcpp::Named("V1")= CheckV,
      Rcpp::Named("V") = CheckV0,
      Rcpp::Named("dist_mat") = dist_mat

    );

}
