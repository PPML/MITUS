#include <Rcpp.h>
#include <math.h>
#include <vector>
using namespace Rcpp;
//'@name fin2_cSim
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
//'@param Birthst Births over time
//'@param HrEntEx Matrix of Entry and Exit rates into the High Risk population
//'@param ImmNon Immigration with no TB
//'@param ImmLat Immigration with Latent TB
//'@param ImmAct Immigration with Active TB
//'@param ImmFst Immigration with Fast Progressing TB
//'@param net_mig_usb net internal migration usb
//'@param net_mig_nusb net internal migration nusb
//'@param mubt background mortality over time
//'@param RelInf beta
//'@param RelInfRg beta based off of risk group
//'@param Vmix 1-sigma
//'@param rEmmigFB rate of emmigration in non-US born population
//'@param TxVec vector of parameters for TB Tx
//'@param TunTxMort Tuning parameter for mortality on TB Tx
//'@param rDeft rate of default from TB treatment over time
//'@param rLtScrt rate of latent screening over time
//'@param ttt_samp_dist probabilities of screening for ttt intervention
//'@param ttt_ag which age groups to apply ttt
//'@param ttt_na which nativity groups to apply ttt
//'@param ttt_ltbi how much to increase ltbi
//'@param ttt_month when to apply ttt interventions
//'@param ttt_pop_frc population size to apply the ltbi prev to
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
//'@param dist_gen general distribution across tb progression and mort
//'@param trans_mat_tot_ages
//'@return Outputs a list of outputs
//[[Rcpp::export]]

Rcpp::List fin2_cSim(
    int                 nYrs,
    int                 nRes,
    Rcpp::NumericMatrix rDxt,
    std::vector<double> TxQualt,
    Rcpp::NumericMatrix       InitPop,
    Rcpp::NumericMatrix       Mpfast,
    std::vector<double> ExogInf,
    Rcpp::NumericMatrix       MpfastPI,
    Rcpp::NumericMatrix       Mrslow,
    std::vector<double> rrSlowFB,
    double              rfast,
    double              RRcurDef,
    double              rSlfCur,
    double              p_HR,
    Rcpp::NumericMatrix       vTMort,
    std::vector<double> RRmuRF,
    std::vector<double> RRmuHR,
    std::vector<double> Birthst,
    Rcpp::NumericMatrix       HrEntEx,
    Rcpp::NumericMatrix       ImmNon,
    Rcpp::NumericMatrix       ImmLat,
    Rcpp::NumericMatrix       ImmAct,
    Rcpp::NumericMatrix       ImmFst,
    std::vector<double> net_mig_usb,
    std::vector<double> net_mig_nusb,
    Rcpp::NumericMatrix       mubt,
    std::vector<double> RelInf,
    std::vector<double> RelInfRg,
    std::vector<double> Vmix,
    std::vector<double> rEmmigFB,
    std::vector<double> TxVec,
    double              TunTxMort,
    std::vector<double> rDeft,
    std::vector<double> rLtScrt,
    Rcpp::NumericMatrix LtTxPar,
    Rcpp::NumericMatrix ttt_samp_dist,
    std::vector<double> ttt_ag,
    std::vector<double> ttt_na,
    std::vector<double> ttt_month,
    double              ttt_pop_frc,
    double              ttt_ltbi,
    Rcpp::NumericMatrix LtDxPar_lt,
    Rcpp::NumericMatrix LtDxPar_nolt,
    std::vector<double> RRdxAge,
    double              rRecov,
    double              pImmScen,
    std::vector<double>  EarlyTrend,
    std::vector<double> pReTx,
    Rcpp::NumericMatrix ag_den,
    std::vector<double> NixTrans,
    Rcpp::NumericMatrix       dist_gen,
    Rcpp::NumericMatrix       trans_mat_tot_ages
) {
  ////////////////////////////////////////////////////////////////////////////////
  ////////    BELOW IS A LIST OF THE VARIABLES CREATED INTERNALLY IN MODEL   /////
  ////////////////////////////////////////////////////////////////////////////////
  int           ti; int n2;
  int           s;
  double        InitPopN[InitPop.nrow()][InitPop.ncol()];
  double        InitPopZ[InitPop.nrow()][InitPop.ncol()];
  double        MpfastN[Mpfast.nrow()][Mpfast.ncol()];
  double        MpslowN[Mpfast.nrow()][Mpfast.ncol()];
  double        MpfastPIN[MpfastPI.nrow()][MpfastPI.ncol()];
  double        MpslowPIN[MpfastPI.nrow()][MpfastPI.ncol()];
  double        MrslowN[Mrslow.nrow()][Mrslow.ncol()];
  double        vTMortN[vTMort.nrow()][vTMort.ncol()];
  double        ag_denN[ag_den.nrow()][ag_den.ncol()];

  // double        RRs_muN[RRs_mu.nrow()][RRs_mu.ncol()];
  // double        vIsxtoIsyN[vIsxtoIsy.nrow()][vIsxtoIsy.ncol()];
  // double        vNmxtoNmyN[vNmxtoNmy.nrow()][vNmxtoNmy.ncol()];
  // double        vrgxtorgyN[vrgxtorgy.nrow()][vrgxtorgy.ncol()];
  double        HrEntExN[HrEntEx.nrow()][HrEntEx.ncol()];
  double        ImmNonN[ImmNon.nrow()][ImmNon.ncol()];
  double        ImmLatN[ImmLat.nrow()][ImmLat.ncol()];
  double        ImmFstN[ImmFst.nrow()][ImmFst.ncol()];
  double        ImmActN[ImmAct.nrow()][ImmAct.ncol()];
  double        mubtN[mubt.nrow()][mubt.ncol()];
  //  double        rIntvInitN[rIntvInit.nrow()][rIntvInit.ncol()];
  double        rDxtN[rDxt.nrow()][rDxt.ncol()];
  double        ttt_samp_distN[ttt_samp_dist.nrow()][ttt_samp_dist.ncol()];
  double        LtTxParN[LtTxPar.nrow()][LtTxPar.ncol()];
  double        LtDxPar_ltN[LtDxPar_lt.nrow()][LtDxPar_lt.ncol()];
  double        LtDxPar_noltN[LtDxPar_nolt.nrow()][LtDxPar_nolt.ncol()];
  double        TxVecZ[6];
  double        temp;
  double        temp2;
  double        temp3;
  double        temp4;
  double        indextemp;
  double        temp4V[11][5];
  double        rTbP;
  double        rTbN;
  double        pop_frc;
  double        rr_ltbi;
  double        Outputs[nYrs][nRes];
  double   V0[11][6][2][4][4][2][3];
  double   V1[11][6][2][4][4][2][3];
  double   VMort[11][6][2][4][4][2][3];
  double        Vdx[11][6][2][4][4][2][3];
  double        VLdx[11][6][2][4][4][2][3];
  double        VNkl[2][2];  ///HIGH AND LOW RISK, NATIVITY
  double        VGjkl[2][2]; ///HIGH AND LOW RISK, NATIVITY
  double        Vjaf[4];     ///BY NUMBER OF MIXING GROUPS
  double        VLjkl[2][2];  ///HIGH AND LOW RISK, NATIVITY
  int           N;
  double        dist_genN[dist_gen.nrow()][dist_gen.ncol()];
  double        temp_vec[4];
  double        temp_mat[4][4];
  double   temp_mat2[4][4];
  double trans_mat_tot_agesN[trans_mat_tot_ages.nrow()][trans_mat_tot_ages.ncol()];
  // double burnvec[12672];
  // double        pop_t;
  // double        sse;
  // long double        rowsum[16];
  double        mat_sum;
  double temp_vec2[4];
  int reblnc; int tb_dyn;
  Rcpp::NumericMatrix Outputs2(nYrs,nRes);
  Rcpp::NumericMatrix dist_mat(4,4);
  double RRmuRFN[4];
  double mort_dist[4];
  // double func_dist[16];
  // long double temp_vec3[11];
  // long double temp_vec4[11];
  // long double temp_vec5[11];

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
  // for(int i=0; i<RRs_mu.nrow(); i++) {
  //   for(int j=0; j<RRs_mu.ncol(); j++) {
  //     RRs_muN[i][j] = RRs_mu(i,j);
  //   } }
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

  for(int i=0; i<mubt.nrow(); i++) {
    for(int j=0; j<mubt.ncol(); j++) {
      mubtN[i][j] = mubt(i,j);

    } }
  for(int i=0; i<rDxt.nrow(); i++) {
    for(int j=0; j<rDxt.ncol(); j++) {
      rDxtN[i][j] = rDxt(i,j);
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
        for(int tb=0; tb<6; tb++) {
          for(int rg=0; rg<2; rg++) {
            for(int lt=0; lt<2; lt++){
              for(int na=0; na<3; na++){
                V0[ag][tb][lt][im][nm][rg][na]    = 0;
                V1[ag][tb][lt][im][nm][rg][na]    = 0;
                VMort[ag][tb][lt][im][nm][rg][na] = 0;
                Vdx[ag][tb][lt][im][nm][rg][na]   = 0;
                VLdx[ag][tb][lt][im][nm][rg][na]  = 0;
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
    mort_dist[i]=0;
    temp_vec2[i]=0;
    temp_vec[i]=0;

    RRmuRFN[i]=RRmuRF[i];
  }
  for(int ag=0; ag<11; ag++) {
    for(int im=0; im<4; im++) {
      temp4V[ag][im] = (1-pow(1-(MrslowN[ag][im])-rRecov,24.0))*(MrslowN[ag][im]);
    } }

  for(int i=0; i<trans_mat_tot_ages.nrow(); i++) {
    for(int j=0; j<trans_mat_tot_ages.ncol(); j++) {
      trans_mat_tot_agesN[i][j] = trans_mat_tot_ages(i,j);
    } }

  // for(int i=0; i<16; i++) {
  //   dist(i)=0;
  //   func_dist[i]=0;
  //   for(int j=0; j<16; j++) {
  //     did_goN[i][j] = 0;
  //   } }
  //
  // for(int i=0; i<4; i++) {
  //   for(int j=0; j<4; j++) {
  //     temp_mat[i][j] = 0;
  //
  //     dist_orig[i][j]= 0;
  //   } }
  //
  for(int i=0; i<4; i++) {
    for(int j=0; j<4; j++) {
      temp_mat2[i][j] = 0;
      //     trans_mat[i][j] = 0;
      //     trans_mat_tot[i][j] = 0;
      //     trans_mat_fin(i,j)=0;
    } }
  //
  // for(int i=0; i<16; i++) {
  //   dist_i_v[i]=0;
  //   temp_vec[i]=0;
  // }
  N=30; indextemp=0;
  reblnc=1;
  tb_dyn=1;
  temp=0; n2=0;
  rr_ltbi=1;
  pop_frc=0;
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
    // for(int i=0; i < 11; i++)
    // {  temp_vec3[i]=0;
    //   temp_vec4[i]=0;
    //   temp_vec5[i]=0;
    // }
  }

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
  for(int m=0; m<1201; m++) {
    // Rcpp::Rcout << m << "\n";
    /////////////////////////////////START BURN IN//////////////////////////////////
    ////////////////////////////////////BIRTHS//////////////////////////////////////
    ///////////USE DISTRIBUTION TO POPULATE THE MODEL ACROSS RISK GROUPS////////////
    for(int im=0; im<4; im++) {
      for(int nm=0; nm<4; nm++){
        V1[0][0][0][im][nm][0][0]  += Birthst[0]*dist_genN[nm][im]*(1-p_HR);
        V1[0][0][0][im][nm][1][0]  += Birthst[0]*dist_genN[nm][im]*(p_HR);
      } }
    // //////////////////////////////////IMMIGRATION///////////////////////////////////
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

          V1[ag][3][0][im][nm][0][1]   += ImmFstN[0][ag]*dist_genN[nm][im]*(1-p_HR);   // LATENT FAST, low risk
          V1[ag][3][0][im][nm][1][1]   += ImmFstN[0][ag]*dist_genN[nm][im]*(p_HR);   // LATENT FAST, high risk

          V1[ag][4][0][im][nm][0][1]   += ImmActN[0][ag]*dist_genN[nm][im]*(1-p_HR);   //ACTIVE TB, low risk
          V1[ag][4][0][im][nm][1][1]   += ImmActN[0][ag]*dist_genN[nm][im]*(p_HR);   //ACTIVE TB, high risk
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
    // for (int i=0; i<4; i++){
    //   temp_vec2[i]=0; }
    // mat_sum=0;
    // ////make a count of # of ppl in each mortality group
    // for(int ag=0; ag<11; ag++) {
    //   for(int tb=0; tb<6; tb++) {
    //     for(int im=0; im<4; im++) {
    //       for(int nm=0; nm<4; nm++){
    //         for(int rg=0; rg<2; rg++){
    //           for(int na=0; na<3; na++){
    //             temp_vec2[nm] += V1[ag][tb][0][im][nm][rg][na];
    //           } } } } } }
    // ////create a population total at this time point
    // for(int nm=0; nm<4; nm++){
    //   mat_sum+=temp_vec2[nm];}
    // ///calculate the mortality
    // for(int nm=0; nm<4; nm++){
    //   mort_dist[nm] = temp_vec2[nm]/mat_sum; }
    // // Rcpp::Rcout << "mort dist is" << mort_dist[nm] << "@m= "<< m<< "for age = "<<ag<<"\n";}
    // temp=0;
    // for(int nm=0; nm<4; nm++){
    //   temp+=(RRmuRF[nm]*mort_dist[nm]);}
    // mat_sum=0;
    // for(int nm=0; nm<4; nm++){
    //   RRmuRFN[nm]=RRmuRF[nm]/temp;
    //   // Rcpp::Rcout << "RRmuRF is" << RRmuRFN[nm] << "@m= "<< m<< "\n";
    // }
    temp=0;
    for(int ag=0; ag<11; ag++) {
      for(int im=0; im<4; im++) {
        for(int nm=0; nm<4; nm++) {
          for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++){
              for(int tb=0; tb<5; tb++) {
                if ((ag<9) | ((RRmuRFN[nm]*RRmuHR[rg]) < 5)){
                  temp = ((RRmuRFN[nm]*RRmuHR[rg])*mubtN[0][ag]);
                } else {
                  temp =  (mubtN[0][ag]*5);
                }
                V1[ag][tb][0][im][nm][rg][na]  -= (V0[ag][tb][0][im][nm][rg][na]*(temp+vTMortN[ag][tb]));

              }//close the tb loop

              ////////////////          MORTALITY WITH TB TREATMENT         ////////////////////
              V1[ag][5 ][0][im][nm][rg][na]  -= V0[ag][5 ][0][im][nm][rg][na]*(temp+vTMortN[ag][5 ]*pow(1.0-TxVecZ[1],TunTxMort)); //check the mortality in param
            } } } } }
    /////////////////////////////////////AGING///////////////////////////////////////
    for(int ag=0; ag<10; ag++) {
      for(int tb=0; tb<6; tb++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++) {
                // /////          IF AGE > 4, IT TAKES 120 MONTHS TO LEAVE AGE GROUP          /////
                // if(ag>0) {
                //   temp2 = 120;
                //   /////          IF AGE < 4, IT TAKES 60 MONTHS TO LEAVE AGE GROUP           /////
                // } else {
                //   temp2 = 60;
                // }
                temp = V0[ag][tb][0][im][nm][rg][na]/ag_denN[0][ag];
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
              temp  = V0[ag][tb][0][im][nm][0][na]*HrEntExN[ag][0]; //entry
              temp2 = V0[ag][tb][0][im][nm][1][na]*HrEntExN[ag][1];  //exit
              //THESE CODES WERE UPDATED, BUT REMAIN ALMOST THE SAME
              V1[ag][tb][0][im][nm][0][na]  += temp2-temp;
              V1[ag][tb][0][im][nm][1][na]  += temp-temp2;
            } } } } }
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
          for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
              // LOW RISK US BORN
              VNkl[0][0]  += V0[ag][tb][0][im][nm][0][0];
              ///////// HIGH RISK US BORN
              ///////// ranges under 1 every time step
              VNkl[1][0]  += V0[ag][tb][0][im][nm][1][0];
              ///////// LOW RISK NON US BORN
              ///////// ranges from .5 -10 people between time steps
              VNkl[0][1]  += V0[ag][tb][0][im][nm][0][1] + V0[ag][tb][0][im][nm][0][2];
              ///////// also ranges from .5 -10 people between time steps
              ///////// HIGH RISK NON US BORN
              VNkl[1][1]  += V0[ag][tb][0][im][nm][1][1] + V0[ag][tb][0][im][nm][1][2];
            } } } }

      // for (int i=0; i<2; i++){
      //   for (int j=0; j<2; j++){
      //     Rcpp::Rcout <<"sums are " <<  VNkl[i][j] << "for i= " << i << " and j = " <<  j << "\n";
      //   }
      // }
      // Step 2  (active TB)
      // Number of active cases* relative infectiousness
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
              /////////  LOW RISK US BORN
              VGjkl[0][0]  +=  V0[ag][tb][0][im][nm][0][0]                               *RelInf[tb];
              ///////// HIGH RISK US BORN
              VGjkl[1][0]  +=  V0[ag][tb][0][im][nm][1][0]                               *RelInf[tb];
              ///////// LOW RISK NON US BORN
              VGjkl[0][1]  += (V0[ag][tb][0][im][nm][0][1] + V0[ag][tb][0][im][nm][0][2])*RelInf[tb];
              /////////  HIGH RISK NON US BORN
              VGjkl[1][1]  += (V0[ag][tb][0][im][nm][1][1] + V0[ag][tb][0][im][nm][1][2])*RelInf[tb];
            } } } }

      // for (int i=0; i<2; i++){
      //   for (int j=0; j<2; j++){
      //     Rcpp::Rcout <<"VGs are " <<  VGjkl[i][j] << "for i= " << i << " and j = " <<  j << "\n";
      //   }
      // }
      // Step 2 (treated TB)
      // No contribution to force of infection

      // Step 3
      Vjaf[0]  =( RelInfRg[0]*VGjkl[0][0] +
        RelInfRg[1]*VGjkl[1][0]*Vmix[0]+
        RelInfRg[2]*VGjkl[0][1]*Vmix[1]+
        RelInfRg[3]*VGjkl[1][1]*Vmix[1]*Vmix[0]) /
          (RelInfRg[0]*VNkl[0][0] +
            RelInfRg[1]*VNkl[1][0]*Vmix[0]+
            RelInfRg[2]*VNkl[0][1]*Vmix[1]+
            RelInfRg[3]*VNkl[1][1]*Vmix[1]*Vmix[0] + 1e-12);

      Vjaf[1] = ((RelInfRg[1]*VGjkl[0][1]) + ( RelInfRg[3]*VGjkl[1][1]*Vmix[0])) /
        ((RelInfRg[1]*VNkl[0][1])  +  (RelInfRg[3]*VNkl[1][1]*Vmix[0]) + 1e-12);

      Vjaf[2] = ((RelInfRg[2]*VGjkl[1][0]) +  (RelInfRg[3]*VGjkl[1][1]*Vmix[1])) /
        ((RelInfRg[2]*VNkl[1][1] ) +  (RelInfRg[3]*VNkl[1][1]*Vmix[1]) + 1e-12);

      Vjaf[3] = VGjkl[1][1] / (VNkl[1][1] + 1e-12);


      // for (int i=0; i<4; i++){
      //  Rcpp::Rcout << "Vjaf is "<<  Vjaf[i] << "\n";
      // }

      // Step 4
      /// LOW RISK US BORN
      VLjkl[0 ][0 ]  = (RelInfRg[0]*Vjaf[0]);
      ///////// HIGH RISK US BORN
      VLjkl[1 ][0 ]  = (RelInfRg[1]*Vjaf[1]*(1-Vmix[0]) + RelInfRg[0]*Vjaf[0]*Vmix[0]);
      ///////// LOW RISK NON US BORN
      VLjkl[0 ][1 ]  = ((RelInfRg[2]*Vjaf[2]*(1-Vmix[1]) + RelInfRg[0]*Vjaf[0]*Vmix[1])) + ExogInf[0];
      ///////// HIGH RISK NON US BORN
      ///check the use of RelInfRg here as beta, might need to be a combo param but unclear check the old param file
      VLjkl[1 ][1 ]  = ((RelInfRg[3]*Vjaf[3]*(1-Vmix[0])*(1-Vmix[1]) + RelInfRg[2]*Vjaf[2]*Vmix[0]*(1-Vmix[1]) + RelInfRg[1]*Vjaf[1]*Vmix[1]*(1-Vmix[0]) + RelInfRg[0]*Vjaf[0]*Vmix[0]*Vmix[1])) + ExogInf[0];

      // for (int i=0; i<2; i++){
      //   for (int j=0; j<2; j++){
      //     Rcpp::Rcout <<"VLs are " <<  VLjkl[i][j] << "for i= " << i << " and j = " <<  j << "\n";
      //   }
      // }
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

                ///////////////////////////////   SUCEPTIBLE  /////™////////////////////////////
                temp = V0[ag][0][0][im][nm][rg][na]*(VLjkl[rg][n2])*EarlyTrend[m];
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
                temp = V0[ag][1][0][im][nm][rg][na]*(VLjkl[rg][n2]);
                V1[ag][1][0][im][nm][rg][na] -= temp;
                V1[ag][2][0][im][nm][rg][na] += temp*MpslowPIN[ag][im];
                V1[ag][3][0][im][nm][rg][na] += temp*MpfastPIN[ag][im];
                ///////////////////////////////////////////////////////////////////////////////

                /////////////////////////////// SUPER-INFECTION LS ////////////////////////////
                temp = V0[ag][2][0][im][nm][rg][na]*(VLjkl[rg][n2]);
                V1[ag][2][0][im][nm][rg][na]  -= temp;
                V1[ag][2][0][im][nm][rg][na]  += temp*MpslowPIN[ag][im];
                V1[ag][3][0][im][nm][rg][na]  += temp*MpfastPIN[ag][im];

              } } } } }

      // temp=0;
      // for(int ag=0; ag<11; ag++) {
      //   for(int im=0; im<4; im++) {
      //     for(int nm=0; nm<4; nm++){
      //       for(int rg=0; rg<2; rg++){
      //         for (int na=0; na<3; na++){
      //           temp+=V1[ag][4][0][im][nm][rg][0];
      //
      //         } } } } }
      // Rcpp::Rcout << "at m= "<< m<< "number of active US: "<< temp << "\n";

      ///////////////////////////         BREAK DOWN      /////////////////////////////
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

      ////////////////////TB DIAGNOSIS AND TX INITIATION PUBLIC ///////////////////////
      ///////////////////////for all age groups, living cond///////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++) {
                temp  = V0[ag][4 ][0][im][nm][rg][na]*rDxtN[0][rg]/RRdxAge[ag]/EarlyTrend[m];
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

    } //end of tb dynamics loop

    ///////                 REBALANCE THE POPULATION
    ////////           NOW FINALLY UPDATE THE DISTRIBUTION           ///////////
    if (reblnc==1){
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for (int im=0; im<4; im++){
            for (int nm=0; nm<4; nm++){
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++){
                  // V2[ag][tb][0][im][nm][rg][na]=0;
                  for (int m2=0; m2<4; m2++){
                    for (int p2=0; p2<4; p2++){
                      temp = V1[ag][tb][0][p2][m2][rg][na]*(trans_mat_tot_agesN[m2+p2*4][(16*(ag+1))-(16-(nm+im*4))]);
                      V1[ag][tb][0][im][nm][rg][na] += temp;
                      V1[ag][tb][0][p2][m2][rg][na] -= temp;

                    } } } } } } } }
    } //end of rebalancing loop
    /////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////RESET POPULATION SIZE/////////////////////////////

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


  } ///////////////////////////END BURN IN///////////////////////////////////////

  // for(int ag=0; ag<11; ag++) {
  //
  //   temp=0;
  //   for (int i=0; i<4; i++){
  //     for (int j=0; j<4; j++){
  //
  //       temp_mat[i][j]=0; } }
  //   mat_sum=0;
  //   ////make a count of # of ppl in each mortality group
  //   // for(int ag=0; ag<11; ag++) {
  //   for(int tb=0; tb<6; tb++) {
  //     for(int na=0; na<3; na++){
  //
  //       for(int im=0; im<4; im++) {
  //         for(int nm=0; nm<4; nm++){
  //           for(int rg=0; rg<2; rg++){
  //             if (V1[ag][tb][0][im][nm][rg][na]<0){
  //               Rcout<< "V1 is neg @ ag " << ag << "\n";
  //             }
  //             temp_mat[nm][im] += V1[ag][tb][0][im][nm][rg][na];
  //             temp+=V1[ag][tb][0][im][nm][1][na];
  //           } } } } }
  //   ////create a population total at this time point
  //   for(int nm=0; nm<4; nm++){
  //     for(int im=0; im<4; im++){
  //
  //       mat_sum+=temp_mat[nm][im];
  //       // Rcpp::Rcout << "mat sum is" << mat_sum << "at ag = "<< ag<<"\n";
  //     }}
  //
  //   ///calculate the mortality
  //   for(int nm=0; nm<4; nm++){
  //     for(int im=0; im<4; im++){
  //       temp_mat2[nm][im] = temp_mat[nm][im]/mat_sum;
  //       Rcpp::Rcout << "mort dist is" << temp_mat2[nm][im] << "at ag = "<< ag <<  " mort = "<< nm << " prog= "<<im <<"\n";}
  //   }
  //   Rcpp::Rcout<<"proportion of high risk at ag = " << ag << "is " << temp/mat_sum<< "\n";
  //   }


  //       for (int i=0; i<4; i++){
  //         for (int j=0; j<4; j++){
  //
  //       Rcpp::Rcout << "dist diff after reblnc is" << dist_mat(i,j) << "@m= "<< m<< "\n";} }
  ///////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////
  // temp=0;
  Rcpp::NumericVector  CheckV0(12672);
  for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
      for(int lt=0; lt<2; lt++){
        for(int im=0; im<4; im++){
          for(int nm=0; nm<4; nm++){
            for(int rg=0; rg<2; rg++) {
              for(int na=0; na<3; na++) {
                //                 if (std::any_of(V1[ag][tb][lt][im][nm][rg][na]<0)){
                // // Rcpp::Rcout << "After burn in pop is negative at ag = " << ag << " tb = "<< tb << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
                //  Rcpp::Rcout << "After burn in pop is negative /n" ;
                //
                //                 }
                CheckV0(ag+tb*11+lt*66+im*132+nm*528+rg*2112+na*4224) = V1[ag][tb][lt][im][nm][rg][na];
              } } } } } } }
  // Rcpp::Rcout << "pop is" << temp << "\n";
  /////////////////////////////////////////////////////////////////////////////
  ///////////////////////////BEGIN MODEL///////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  ///////////////////CREATE TIME LOOPS FOR MODEL RUN///////////////////////////
  ////////////////////////////////YEAR LOOP////////////////////////////////////
  for(int y=0; y<nYrs; y++) {
    /////////////////////////////////MONTH LOOP////////////////////////////////////
    for(int m=0; m<12; m++) {
      /////////////////////CREATE A COUNTER OF MONTHS SINCE START////////////////////
      s = y*12+m;
      // Rcpp::Rcout << s <<"\n";
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
      TxVecZ[3] = TxVec[0]*(1-TxVecZ[1])*(1-pReTx[s]) + rDeft[s]*(1-TxVecZ[1]*RRcurDef);
      ///////// RATE OF TREATMENT EXIT TO RE TREATMENT////////////////////////////////
      TxVecZ[4] = TxVec[0]*(1-TxVecZ[1])*(pReTx[s]);
      ////////////////////// P(TREATMENT COMPLETION) /////////////////////////////////
      TxVecZ[5] = TxVec[0]*(1-(1.0-TxVecZ[1])*pReTx[s]);
      // /////////////////////////////////////////////////////////////////////////////////
      // /////////////////////////////////////BIRTHS//////////////////////////////////////
      // ///////////////////////////////////////////////////////////////////////////////
      for (int im=0; im<4; im++) {
        for (int nm=0; nm<4; nm++) {
          /////LOW RISK GROUP BIRTHS//////////////////////////////////////////////////////
          V1[0][0][0][im][nm][0][0]  += Birthst[s]*dist_genN[nm][im]*(1-p_HR);
          /////HIGH RISK GROUP BIRTHS//////////////////////////////////////////////////////
          V1[0][0][0][im][nm][1][0]  += Birthst[s]*dist_genN[nm][im]*(p_HR);
        } }
      // for(int ag=0; ag<11; ag++) {
      //   for(int tb=0; tb<6; tb++) {
      //     for(int lt=0; lt<2; lt++){
      //       for(int im=0; im<4; im++){
      //         for(int nm=0; nm<4; nm++){
      //           for(int rg=0; rg<2; rg++) {
      //             for(int na=0; na<3; na++) {
      //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
      //                 Rcpp::Rcout << "After Births pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
      //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
      //
      //               }
      //             } } } } } } }

      /////////////////////////////// IMMIGRATION ///////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {

            V1[ag][0][0][im][nm][0][1]   += ImmNonN[s][ag]*(1-p_HR)*dist_genN[nm][im];  // NO TB, low risk
            V1[ag][0][0][im][nm][1][1]   += ImmNonN[s][ag]*(p_HR)*dist_genN[nm][im];    // NO TB, high risk

            V1[ag][2][0][im][nm][0][1]   += ImmLatN[s][ag]*(1-p_HR)*dist_genN[nm][im]; // LATENT SLOW TB, low risk
            V1[ag][2][0][im][nm][1][1]   += ImmLatN[s][ag]*(p_HR)*dist_genN[nm][im];   // LATENT SLOW TB, high risk

            V1[ag][3][0][im][nm][0][1]   += ImmFstN[s][ag]*(1-p_HR)*dist_genN[nm][im];   // LATENT FAST, low risk
            V1[ag][3][0][im][nm][1][1]   += ImmFstN[s][ag]*(p_HR)*dist_genN[nm][im];   // LATENT FAST, high risk

            V1[ag][4][0][im][nm][0][1]   += ImmActN[s][ag]*(1-p_HR)*dist_genN[nm][im];   //ACTIVE TB, low risk
            V1[ag][4][0][im][nm][1][1]   += ImmActN[s][ag]*(p_HR)*dist_genN[nm][im];   //ACTIVE TB, high risk
          } } }
      // for(int ag=0; ag<11; ag++) {
      //   for(int tb=0; tb<6; tb++) {
      //     for(int lt=0; lt<2; lt++){
      //       for(int im=0; im<4; im++){
      //         for(int nm=0; nm<4; nm++){
      //           for(int rg=0; rg<2; rg++) {
      //             for(int na=0; na<3; na++) {
      //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
      //                 Rcpp::Rcout << "after immigration pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
      //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
      //
      //               }
      //             } } } } } } }
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
      // for(int ag=0; ag<11; ag++) {
      //   for(int tb=0; tb<6; tb++) {
      //     for(int lt=0; lt<2; lt++){
      //       for(int im=0; im<4; im++){
      //         for(int nm=0; nm<4; nm++){
      //           for(int rg=0; rg<2; rg++) {
      //             for(int na=0; na<3; na++) {
      //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
      //                 Rcpp::Rcout << "after emigration pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
      //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
      //
      //               }
      //             } } } } } } }
      /////////////////////////////////  NET INTERNAL MIGRATION  ///////////////////////////////////
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++){
            for(int im=0; im<4; im++){
              for(int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++) {
                  V1[ag][tb][lt][im][nm][rg][0]  += V0[ag][tb][lt][im][nm][rg][0]*net_mig_usb[ag];      // US
                  V1[ag][tb][lt][im][nm][rg][1]  += V0[ag][tb][lt][im][nm][rg][1]*net_mig_nusb[ag];      // FB1
                  V1[ag][tb][lt][im][nm][rg][2]  += V0[ag][tb][lt][im][nm][rg][2]*net_mig_nusb[ag];      // FB2
                } } } } } }
      /////////////////////////////////  MORTALITY  ///////////////////////////////////

      // for (int i=0; i<4; i++){
      //   temp_vec2[i]=0; }
      // mat_sum=0;
      // ///create a vector of the counts of # ppl in each mort group
      //   for(int nm=0; nm<4; nm++){
      //     for(int tb=0; tb<6; tb++) {
      //       for(int lt=0; lt<2; lt++){
      //         for(int im=0; im<4; im++) {
      //           for(int rg=0; rg<2; rg++){
      //             for(int na=0; na<3; na++){
      //               temp_vec2[nm] += V1[ag][tb][lt][im][nm][rg][na];
      //             } } } } } }
      // ////create a population total
      // for(int nm=0; nm<4; nm++){
      //   mat_sum+=temp_vec2[nm];
      // }
      // ///mortality distribution at mortality time step
      // for(int nm=0; nm<4; nm++){
      //   mort_dist[nm] = temp_vec2[nm]/mat_sum; }
      // //reset the temp variable
      // temp=0;
      // ///
      // for(int nm=0; nm<4; nm++){
      //   temp+=(RRmuRF[nm]*mort_dist[nm]);}
      // for(int nm=0; nm<4; nm++){
      //   RRmuRFN[nm]=RRmuRF[nm]/temp;
      // }
      temp=0;
      for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++){
          for(int im=0; im<4 ; im++) {
            for(int nm=0; nm<4; nm++) {
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                  for(int tb=0; tb<4; tb++) {
                    if ((ag<9) | ((RRmuRFN[nm]*RRmuHR[rg]) < 5)){
                      temp = ((RRmuRFN[nm]*RRmuHR[rg])*mubtN[s][ag]);}
                    else {
                      temp =  (mubtN[s][ag]*5);
                    }

                    VMort[ag][tb ][lt][im][nm][rg][na]  = V0[ag][tb][lt][im][nm][rg][na]*temp;
                  }//close the tb loop
                  ////////////////////////      ACTIVE TB         /////////////////////////////////
                  VMort[ag][4 ][lt][im][nm][rg][na]  = V0[ag][4 ][lt][im][nm][rg][na]*(temp+vTMortN[ag][4 ] );

                  ////////////////////////    TB TREATMENT        /// //////////////////////////////
                  VMort[ag][5 ][lt][im][nm][rg][na]  = V0[ag][5 ][lt][im][nm][rg][na]*(temp+vTMortN[ag][5 ]*pow(1.0-TxVecZ[1],TunTxMort));

                  ///////////// UPDATE THE PRIMARY VECTOR BY REMOVING MORTALITY /////////////////
                  for(int tb=0; tb<6; tb++) {
                    V1[ag][tb][lt][im][nm][rg][na]  -= VMort[ag][tb][lt][im][nm][rg][na];
                    // temp += VMort[ag][tb][lt][im][nm][rg][na];
                    // Rcpp::Rcout << "total mortality at time" << s << "is" << temp << "\n";
                  }

                } } } } } }

      // for(int ag=0; ag<11; ag++) {
      //   for(int tb=0; tb<6; tb++) {
      //     for(int lt=0; lt<2; lt++){
      //       for(int im=0; im<4; im++){
      //         for(int nm=0; nm<4; nm++){
      //           for(int rg=0; rg<2; rg++) {
      //             for(int na=0; na<3; na++) {
      //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
      //                 Rcpp::Rcout << "after mort pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
      //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
      //
      //               }
      //             } } } } } } }
      /////////////////////////////////////AGING///////////////////////////////////////
      for(int ag=0; ag<10; ag++) {
        /////          IF AGE > 4, IT TAKES 120 MONTHS TO LEAVE AGE GROUP          /////
        // if(ag>0) {
        //   temp2 = 120;
        //   /////          IF AGE < 4, IT TAKES 60 MONTHS TO LEAVE AGE GROUP           /////
        // } else {
        //   temp2 = 60;
        // }
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
      // for(int ag=0; ag<11; ag++) {
      //   for(int tb=0; tb<6; tb++) {
      //     for(int lt=0; lt<2; lt++){
      //       for(int im=0; im<4; im++){
      //         for(int nm=0; nm<4; nm++){
      //           for(int rg=0; rg<2; rg++) {
      //             for(int na=0; na<3; na++) {
      //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
      //                 Rcpp::Rcout << "after aging pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
      //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
      //
      //               }
      //             } } } } } } }
      ///////////////////////// NEW FB -> ESTABLISHED FB ///////////////////////////////
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
      // for(int ag=0; ag<11; ag++) {
      //   for(int tb=0; tb<6; tb++) {
      //     for(int lt=0; lt<2; lt++){
      //       for(int im=0; im<4; im++){
      //         for(int nm=0; nm<4; nm++){
      //           for(int rg=0; rg<2; rg++) {
      //             for(int na=0; na<3; na++) {
      //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
      //                 Rcpp::Rcout << "after fb trans pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
      //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
      //
      //               }
      //             } } } } } } }
      //////////////////////////// HIGH-RISK ENTRY/EXIT ////////////////////////////////
      //is there ever a way for high risk exit to be greater than high risk entry
      for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int na=0; na<3; na++) {
                  temp  = V0[ag][tb][lt][im][nm][0][na]*HrEntExN[ag][0];
                  temp2 = V0[ag][tb][lt][im][nm][1][na]*HrEntExN[ag][1];
                  V1[ag][tb][lt][im][nm][0][na]  += (temp2-temp);
                  V1[ag][tb][lt][im][nm][1][na]  += (temp-temp2);
                } } } } } }
      // for(int ag=0; ag<11; ag++) {
      //   for(int tb=0; tb<6; tb++) {
      //     for(int lt=0; lt<2; lt++){
      //       for(int im=0; im<4; im++){
      //         for(int nm=0; nm<4; nm++){
      //           for(int rg=0; rg<2; rg++) {
      //             for(int na=0; na<3; na++) {
      //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
      //                 Rcpp::Rcout << "after hrentex pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
      //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
      //
      //               }
      //             } } } } } } }
      // for(int ag=0; ag<11; ag++) {
      // for(int tb=0; tb<6; tb++) {
      //   for(int im=0; im<4; im++) {
      //     for(int lt=0; lt<2; lt++){
      //       for(int na=0; na<3; na++){
      //         for(int nm=0; nm<4; nm++){
      //           for(int rg=0; rg<2; rg++){
      //             if (V1[ag][tb][lt][im][nm][rg][na]<0){
      //               Rcout<< "V1 is neg @ na " << na << "ag = " << ag << "tb="<< tb << "rg = " << rg << "nm = "<<nm<< "im = "<< im<< "\n";
      //             } } } } } } } }

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

                  // LOW RISK US BORN
                  VNkl[0][0]  += V0[ag][tb][lt][im][nm][0][0];
                  ///////// HIGH RISK US BORN
                  ///////// ranges under 1 every time step
                  VNkl[1][0]  += V0[ag][tb][lt][im][nm][1][0];
                  ///////// LOW RISK NON US BORN
                  ///////// ranges from .5 -10 people between time steps
                  VNkl[0][1]  += V0[ag][tb][lt][im][nm][0][1] + V0[ag][tb][lt][im][nm][0][2];
                  ///////// also ranges from .5 -10 people between time steps
                  ///////// HIGH RISK NON US BORN
                  VNkl[1][1]  += V0[ag][tb][lt][im][nm][1][1] + V0[ag][tb][lt][im][nm][1][2];
                } } } } }

        // for (int i=0; i<2; i++){
        //   for (int j=0; j<2; j++){
        //     Rcpp::Rcout <<"sums are " <<  VNkl[i][j] << "for i= " << i << " and j = " <<  j << "\n";
        //   }
        // }
        // Step 2  (active TB)
        // Number of active cases* relative infectiousness
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int lt=0; lt<2; lt++) {

                  /////////  LOW RISK US BORN
                  VGjkl[0][0]  +=  V0[ag][tb][lt][im][nm][0][0]                                *RelInf[tb];
                  ///////// HIGH RISK US BORN
                  VGjkl[1][0]  +=  V0[ag][tb][lt][im][nm][1][0]                                *RelInf[tb];
                  ///////// LOW RISK NON US BORN
                  VGjkl[0][1]  += (V0[ag][tb][lt][im][nm][0][1] + V0[ag][tb][lt][im][nm][0][2])*RelInf[tb];
                  /////////  HIGH RISK NON US BORN
                  VGjkl[1][1]  += (V0[ag][tb][lt][im][nm][1][1] + V0[ag][tb][lt][im][nm][1][2])*RelInf[tb];
                } } } } }

        // for (int i=0; i<2; i++){
        //   for (int j=0; j<2; j++){
        //     Rcpp::Rcout <<"VGs are " <<  VGjkl[i][j] << "for i= " << i << " and j = " <<  j << "\n";
        //   }
        // }
        // Step 2 (treated TB)
        // No contribution to force of infection

        // Step 3
        Vjaf[0]  =( RelInfRg[0]*VGjkl[0][0] +
          RelInfRg[1]*VGjkl[1][0]*Vmix[0]+
          RelInfRg[2]*VGjkl[0][1]*Vmix[1]+
          RelInfRg[3]*VGjkl[1][1]*Vmix[1]*Vmix[0]) /
            ((RelInfRg[0]*VNkl[0][0] +
              RelInfRg[1]*VNkl[1][0]*Vmix[0]+
              RelInfRg[2]*VNkl[0][1]*Vmix[1]+
              RelInfRg[3]*VNkl[1][1]*Vmix[1]*Vmix[0]) + 1e-12);

        Vjaf[1] = ((RelInfRg[1]*VGjkl[0][1]) + ( RelInfRg[3]*VGjkl[1][1]*Vmix[0])) /
          ((RelInfRg[1]*VNkl[0][1])  +  (RelInfRg[3]*VNkl[1][1]*Vmix[0]) + 1e-12);

        Vjaf[2] = ((RelInfRg[2]*VGjkl[1][0]) +  (RelInfRg[3]*VGjkl[1][1]*Vmix[1])) /
          ((RelInfRg[2]*VNkl[1][1] ) +  (RelInfRg[3]*VNkl[1][1]*Vmix[1]) + 1e-12);

        Vjaf[3] = VGjkl[1][1] / (VNkl[1][1] + 1e-12);


        // for (int i=0; i<4; i++){
        //   Rcpp::Rcout << "Vjaf is "<<  Vjaf[i] << "\n";
        // }

        // Step 4
        /// LOW RISK US BORN
        VLjkl[0 ][0 ]  = (RelInfRg[0]*Vjaf[0]);
        ///////// HIGH RISK US BORN
        VLjkl[1 ][0 ]  = (RelInfRg[1]*Vjaf[1]*(1-Vmix[0]) + RelInfRg[0]*Vjaf[0]*Vmix[0]);
        ///////// LOW RISK NON US BORN
        VLjkl[0 ][1 ]  = ((RelInfRg[2]*Vjaf[2]*(1-Vmix[1]) + RelInfRg[0]*Vjaf[0]*Vmix[1])) + ExogInf[s];
        ///////// HIGH RISK NON US BORN
        ///check the use of RelInfRg here as beta, might need to be a combo param but unclear check the old param file
        VLjkl[1 ][1 ]  = ((RelInfRg[3]*Vjaf[3]*(1-Vmix[0])*(1-Vmix[1]) + RelInfRg[2]*Vjaf[2]*Vmix[0]*(1-Vmix[1]) + RelInfRg[1]*Vjaf[1]*Vmix[1]*(1-Vmix[0]) + RelInfRg[0]*Vjaf[0]*Vmix[0]*Vmix[1]) ) + ExogInf[s];


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
                    temp = V0[ag][0][lt][im][nm][rg][na]*(VLjkl[rg][n2])*NixTrans[s];
                    //////////////////////////// REMOVE FROM SUSCEPTIBLE //////////////////////////
                    V1[ag][0][lt][im][nm][rg][na]  -= temp;
                    //////////////////////////////// LATENT TB SLOW ///////////////////////////////
                    V1[ag][2][lt][im][nm][rg][na]  += temp*MpslowN[ag][im];
                    //////////////////////////////// LATENT TB FAST ///////////////////////////////
                    V1[ag][3][lt][im][nm][rg][na]  += temp*MpfastN[ag][im];
                    ///////////////////////////////////////////////////////////////////////////////

                    /////////////////////////////// SUPER-INFECTION SP ////////////////////////////
                    temp = V0[ag][1][lt][im][nm][rg][na]*(VLjkl[rg][n2])*NixTrans[s];
                    V1[ag][1][lt][im][nm][rg][na] -= temp;
                    V1[ag][2][lt][im][nm][rg][na] += temp*MpslowPIN[ag][im];
                    V1[ag][3][lt][im][nm][rg][na] += temp*MpfastPIN[ag][im];
                    ///////////////////////////////////////////////////////////////////////////////

                    /////////////////////////////// SUPER-INFECTION LS ////////////////////////////
                    temp = V0[ag][2][lt][im][nm][rg][na]*(VLjkl[rg][n2])*NixTrans[s];
                    V1[ag][2][lt][im][nm][rg][na]  -= temp;
                    V1[ag][2][lt][im][nm][rg][na]  += temp*MpslowPIN[ag][im];
                    V1[ag][3][lt][im][nm][rg][na]  += temp*MpfastPIN[ag][im];

                  } } } } } }

        // for(int ag=0; ag<11; ag++) {
        //   for(int tb=0; tb<6; tb++) {
        //     for(int lt=0; lt<2; lt++){
        //       for(int im=0; im<4; im++){
        //         for(int nm=0; nm<4; nm++){
        //           for(int rg=0; rg<2; rg++) {
        //             for(int na=0; na<3; na++) {
        //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
        //                 Rcpp::Rcout << "after infection pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
        //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
        //
        //               }
        //             } } } } } } }
        ///////////////////////////////////BREAKDOWN///////////////////////////////////
        ///////////////////////for all age groups, risk groups/////////////////////////
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
        // for(int ag=0; ag<11; ag++) {
        //   for(int tb=0; tb<6; tb++) {
        //     for(int lt=0; lt<2; lt++){
        //       for(int im=0; im<4; im++){
        //         for(int nm=0; nm<4; nm++){
        //           for(int rg=0; rg<2; rg++) {
        //             for(int na=0; na<3; na++) {
        //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
        //                 Rcpp::Rcout << "after breakdown pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
        //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
        //
        //               }
        //             } } } } } } }
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
        // for(int ag=0; ag<11; ag++) {
        //   for(int tb=0; tb<6; tb++) {
        //     for(int lt=0; lt<2; lt++){
        //       for(int im=0; im<4; im++){
        //         for(int nm=0; nm<4; nm++){
        //           for(int rg=0; rg<2; rg++) {
        //             for(int na=0; na<3; na++) {
        //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
        //                 Rcpp::Rcout << "after lsts pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
        //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
        //
        //               }
        //             } } } } } } }
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
        // for(int ag=0; ag<11; ag++) {
        //   for(int tb=0; tb<6; tb++) {
        //     for(int lt=0; lt<2; lt++){
        //       for(int im=0; im<4; im++){
        //         for(int nm=0; nm<4; nm++){
        //           for(int rg=0; rg<2; rg++) {
        //             for(int na=0; na<3; na++) {
        //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
        //                 Rcpp::Rcout << "after slfcur pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
        //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
        //
        //               }
        //             } } } } } } }

        /// LTBI SCREENING AND TLTBI INITIATION /// only for no previous TB or LTBI tx
        int agi = ttt_ag.size(); int nai = ttt_na.size(); int si = ttt_month.size();
        rr_ltbi=1; pop_frc=0;
        for(int ag=0; ag<11; ag++) {
          for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
            ////////////// US BORN, LOW RISK  //////////////////
            if( rg==0 & na==0) {
              rTbP = rLtScrt[s]*LtDxPar_ltN[0][s];
              rTbN = rLtScrt[s]*LtDxPar_noltN[0][s];
            }
            ////////////// US BORN, HIGH RISK  /////////////////
            if(rg==1 & na==0) {
              rTbP = rLtScrt[s]*LtDxPar_ltN[1][s];
              rTbN = rLtScrt[s]*LtDxPar_noltN[1][s];
            }
              ////////////// Young NUS (under 5)  /////////////////
              if(rg==0 & na > 0 & ag==0) {
                rTbP = rLtScrt[s]*LtDxPar_ltN[2][s];
                rTbN = rLtScrt[s]*LtDxPar_noltN[2][s];
              }
              //////////// NON US BORN  ////////////////
              if(rg==0 & na > 0 & ag > 0) {
                rTbP = rLtScrt[s]*LtDxPar_ltN[3][s];
                rTbN = rLtScrt[s]*LtDxPar_noltN[3][s];
              }
              ////////////// NON US BORN, HIGH RISK  /////////////////
              if(rg==1 & na >0) {
                rTbP = rLtScrt[s]*LtDxPar_ltN[4][s];
                rTbN = rLtScrt[s]*LtDxPar_noltN[4][s];
              }
              for (int i=0; i<agi; i++){
              for (int j=0; j<nai; j++){
                for (int k=0; k<si; k++){
                  if (ag==ttt_ag[i] & na==ttt_na[j] & s==ttt_month[si]){

                    ////////////// US BORN, LOW RISK  //////////////////
                    if(rg==0 & na==0) {
                      rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[0][s];
                      rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[0][s];
                    }
                    ////////////// US BORN, HIGH RISK  /////////////////
                    if(rg==1 & na==0) {
                      rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[1][s];
                      rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[1][s];
                    }
                    ////////////// Young NUS (under 5)  /////////////////
                    if(rg==0 & na > 0 & ag==0) {
                      rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[2][s];
                      rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[2][s];
                    }
                    //////////// NON US BORN  ////////////////
                    if(rg==0 & na > 0 & ag > 0) {
                      rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[3][s];
                      rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[3][s];
                    }
                    ////////////// NON US BORN, HIGH RISK  /////////////////
                    if(rg==1 & na >0) {
                      rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[4][s];
                      rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[4][s];
                    }
                    rr_ltbi=ttt_ltbi;
                    pop_frc=ttt_pop_frc;
                  }
                } } }

                  ////////////// Dont have LTBI
                temp= V0[ag][0][0][im][nm][rg][na]*(((1- pop_frc)*rTbN) +
                  (pop_frc*(1-(rTbP*rr_ltbi))));
                temp2= (V0[ag][1][0][im][nm][rg][na]*((1- pop_frc)*rTbN )+
                 ( pop_frc*(1-(rTbP*rr_ltbi))));
                  V1[ag][0][0][im][nm][rg][na]  -= temp;
                  V1[ag][1][0][im][nm][rg][na]  -= temp2;
                  ///////moving to latent tx experienced as in last model -- is this correct?
                  V1[ag][0][1][im][nm][rg][na]  += temp;
                  V1[ag][1][1][im][nm][rg][na]  += temp2;
                } } } } }
        // for(int ag=0; ag<11; ag++) {
        //   for(int tb=0; tb<6; tb++) {
        //     for(int lt=0; lt<2; lt++){
        //       for(int im=0; im<4; im++){
        //         for(int nm=0; nm<4; nm++){
        //           for(int rg=0; rg<2; rg++) {
        //             for(int na=0; na<3; na++) {
        //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
        //                 Rcpp::Rcout << "after screen pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
        //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
        //
        //               }
        //             } } } } } } }
        /// TLTBI: TX COMPLETION + DEFAULT /// only need to consider tx naive compartment
        rr_ltbi=1; pop_frc=0;
              for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                  ////////////// US BORN, LOW RISK  //////////////////
                  if( rg==0 & na==0) {
                    rTbP = rLtScrt[s]*LtDxPar_ltN[0][s];
                    rTbN = rLtScrt[s]*LtDxPar_noltN[0][s];
                  }
                  ////////////// US BORN, HIGH RISK  /////////////////
                  if(rg==1 & na==0) {
                    rTbP = rLtScrt[s]*LtDxPar_ltN[1][s];
                    rTbN = rLtScrt[s]*LtDxPar_noltN[1][s];
                  }
                  for(int ag=0; ag<11; ag++) {

                  ////////////// Young NUS (under 5)  /////////////////
                  if(rg==0 & na > 0 & ag==0) {
                    rTbP = rLtScrt[s]*LtDxPar_ltN[2][s];
                    rTbN = rLtScrt[s]*LtDxPar_noltN[2][s];
                  }
                  //////////// NON US BORN  ////////////////
                  if(rg==0 & na > 0 & ag > 0) {
                    rTbP = rLtScrt[s]*LtDxPar_ltN[3][s];
                    rTbN = rLtScrt[s]*LtDxPar_noltN[3][s];
                  }
                  ////////////// NON US BORN, HIGH RISK  /////////////////
                  if(rg==1 & na >0) {
                    rTbP = rLtScrt[s]*LtDxPar_ltN[4][s];
                    rTbN = rLtScrt[s]*LtDxPar_noltN[4][s];
                  }
                  for(int im=0; im<4 ; im++) {
                    for(int nm=0; nm<4; nm++) {
                      for (int i=0; i<agi; i++){
                        for (int j=0; j<nai; j++){
                          for (int k=0; k<si; k++){
                            if (ag==ttt_ag[i] & na==ttt_na[j] & s==ttt_month[si]){

                              ////////////// US BORN, LOW RISK  //////////////////
                              if(rg==0 & na==0) {
                                rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[0][s];
                                rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[0][s];
                              }
                              ////////////// US BORN, HIGH RISK  /////////////////
                              if(rg==1 & na==0) {
                                rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[1][s];
                                rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[1][s];
                              }
                              ////////////// Young NUS (under 5)  /////////////////
                              if(rg==0 & na > 0 & ag==0) {
                                rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[2][s];
                                rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[2][s];
                              }
                              //////////// NON US BORN  ////////////////
                              if(rg==0 & na > 0 & ag > 0) {
                                rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[3][s];
                                rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[3][s];
                              }
                              ////////////// NON US BORN, HIGH RISK  /////////////////
                              if(rg==1 & na >0) {
                                rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[4][s];
                                rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[4][s];
                              }
                              rr_ltbi=ttt_ltbi;
                              pop_frc=ttt_pop_frc;
                            }
                          } } }
                  //N(latent)*r(posLTBIscreen)*p(TLTBI Initiation)*p(TLTBI completion)
                  temp  = V0[ag][2][0][im][nm][rg][na]*(((1-pop_frc)*rTbP)+
                          (pop_frc*(rTbP*rr_ltbi)))
                          *LtTxParN[s][0]*(1-LtTxParN[s][1]); // tx completion
                        temp2  = V0[ag][3][0][im][nm][rg][na]*(((1-pop_frc)*rTbP)+
                        (pop_frc*(rTbP*rr_ltbi)))
                          *LtTxParN[s][0]*(1-LtTxParN[s][1]); // tx completion

                        temp3 =V0[ag][2][0][im][nm][rg][na]*(((1-pop_frc)*rTbP)+
                        (pop_frc*(rTbP*rr_ltbi)))
                          *LtTxParN[s][0]*LtTxParN[s][1]; // default
                        temp4 = V0[ag][3][0][im][nm][rg][na]*(((1-pop_frc)*rTbP)+
                        (pop_frc*(rTbP*rr_ltbi)))
                          *LtTxParN[s][0]*LtTxParN[s][1]; // default
                  V1[ag][2][0][im][nm][rg][na]  -=  (temp+temp3); //remove from latent slow
                  V1[ag][3][0][im][nm][rg][na]  -=  (temp2+temp4);  //remove from latent fast

                  V1[ag][1][1][im][nm][rg][na]   += (temp+temp2)*LtTxParN[s][2]; //*EffLt0[s]; //exit to cure
                  V1[ag][2][1][im][nm][rg][na]   += (temp+temp2)*(1-LtTxParN[s][2]); //*(1-EffLt0[s]) //tx comp fail to latent slow

                  ///defaults are placed in tx naive because it is considered the same tb infection
                  V1[ag][2][0][im][nm][rg][na]  += (temp3+temp4); //latent tx default to latent slow
                } } } } }
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
        // for(int ag=0; ag<11; ag++) {
        //   for(int tb=0; tb<6; tb++) {
        //     for(int lt=0; lt<2; lt++){
        //       for(int im=0; im<4; im++){
        //         for(int nm=0; nm<4; nm++){
        //           for(int rg=0; rg<2; rg++) {
        //             for(int na=0; na<3; na++) {
        //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
        //                 Rcpp::Rcout << "after tbdx pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
        //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
        //
        //               }
        //             } } } } } } }
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
        // for(int ag=0; ag<11; ag++) {
        //   for(int tb=0; tb<6; tb++) {
        //     for(int lt=0; lt<2; lt++){
        //       for(int im=0; im<4; im++){
        //         for(int nm=0; nm<4; nm++){
        //           for(int rg=0; rg<2; rg++) {
        //             for(int na=0; na<3; na++) {
        //               if (V1[ag][tb][lt][im][nm][rg][na]<0){
        //                 Rcpp::Rcout << "after txout pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
        //                 Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
        //
        //               }
        //             } } } } } } }
      }//end of TB loop
      //
      // // temp=0;
      // //   temp2=0;
      // //   for(int ag=0; ag<11; ag++) {
      // //       for(int im=0; im<4; im++) {
      // //         for(int nm=0; nm<4; nm++){
      // //           for(int rg=0; rg<2; rg++){
      // //             for (int na=0; na<3; na++){
      // //               for(int lt=0; lt<2; lt++){
      // //                 for(int tb=0; tb<6; tb++) {
      // //
      // //                 temp2+=V1[ag][tb][lt][im][nm][rg][na];
      // //                   if (V1[ag][tb][0][im][nm][rg][na] <0 ){
      // //                     Rcpp::Rcout << "pop is negative at ag " << ag <<"tb=" << tb << "na " << na << "im " << im << "nm "<<nm << "& rg" << rg<< "\n";
      // //                     Rcpp::Rcout << "pop is "<< V1[ag][tb][0][im][nm][rg][na]<< "at ag " << ag <<"tb=" << tb << "na " << na << "im " << im << "nm "<<nm << "& rg" << rg<< "\n";
      // //                   } }
      // //
      // //               temp+=V1[ag][0][lt][im][nm][rg][na];
      // //
      // //               }
      // //
      // //
      // //           } } } } }
      // //
      // //   Rcpp::Rcout << "at m= "<< m<< "population total: "<< temp2 << "\n";
      // //
      // //   Rcpp::Rcout << "at m= "<< m<< "number of susceptibles: "<< temp << "\n";

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
        ///////////// TB MORTALITY COUNT BY AGE, NATIVITY///////////////
        ///////////// This output will likely be updated but it is not calib'd
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    if(na>0) {
                      ti = 11;
                    } else { ti = 0; }
                    // if(im==3) {
                    //   temp = muTbRF;
                    // } else { temp = 0;
                    // }
                    temp=0;
                    Outputs[y][87+ag+ti]  += V0[ag][4 ][lt][im][nm][rg][na]*(vTMortN[ag][4 ]+temp);
                    Outputs[y][87+ag+ti]  += V0[ag][5 ][lt][im][nm][rg][na]*(vTMortN[ag][5 ]+temp)*pow(1.0-TxVecZ[1],TunTxMort); //last term is approx .6
                  } } } } } }
        ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
        for(int i=87; i<109; i++) { Outputs[y][i] = Outputs[y][i]*12; }

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
        for(int i=134; i<151; i++) { Outputs[y][i] = Outputs[y][i]*12; } //yes these are updated
        /// TLTBI INITS ///
        int agi = ttt_ag.size(); int nai = ttt_na.size(); int si = ttt_month.size();
        rr_ltbi=1; pop_frc=0;
        for(int ag=0; ag<11; ag++) {
        for(int im=0; im<4; im++) {
          for(int nm=0; nm<4; nm++) {
        for(int rg=0; rg<2; rg++) {
          for(int na=0; na<3; na++) {
            ////////////// US BORN, LOW RISK  //////////////////
            if( rg==0 & na==0) {
              rTbP = rLtScrt[s]*LtDxPar_ltN[0][s];
              rTbN = rLtScrt[s]*LtDxPar_noltN[0][s];
            }
            ////////////// US BORN, HIGH RISK  /////////////////
            if(rg==1 & na==0) {
              rTbP = rLtScrt[s]*LtDxPar_ltN[1][s];
              rTbN = rLtScrt[s]*LtDxPar_noltN[1][s];
            }
              ////////////// Young NUS (under 5)  /////////////////
              if(rg==0 & na > 0 & ag==0) {
                rTbP = rLtScrt[s]*LtDxPar_ltN[2][s];
                rTbN = rLtScrt[s]*LtDxPar_noltN[2][s];
              }
              //////////// NON US BORN  ////////////////
              if(rg==0 & na > 0 & ag > 0) {
                rTbP = rLtScrt[s]*LtDxPar_ltN[3][s];
                rTbN = rLtScrt[s]*LtDxPar_noltN[3][s];
              }
              ////////////// NON US BORN, HIGH RISK  /////////////////
              if(rg==1 & na >0) {
                rTbP = rLtScrt[s]*LtDxPar_ltN[4][s];
                rTbN = rLtScrt[s]*LtDxPar_noltN[4][s];
              }
              for (int i=0; i<agi; i++){
                for (int j=0; j<nai; j++){
                  for (int k=0; k<si; k++){
                    if (ag==ttt_ag[i] & na==ttt_na[j] & s==ttt_month[si]){

                      ////////////// US BORN, LOW RISK  //////////////////
                      if(rg==0 & na==0) {
                        rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[0][s];
                        rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[0][s];
                      }
                      ////////////// US BORN, HIGH RISK  /////////////////
                      if(rg==1 & na==0) {
                        rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[1][s];
                        rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[1][s];
                      }
                      ////////////// Young NUS (under 5)  /////////////////
                      if(rg==0 & na > 0 & ag==0) {
                        rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[2][s];
                        rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[2][s];
                      }
                      //////////// NON US BORN  ////////////////
                      if(rg==0 & na > 0 & ag > 0) {
                        rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[3][s];
                        rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[3][s];
                      }
                      ////////////// NON US BORN, HIGH RISK  /////////////////
                      if(rg==1 & na >0) {
                        rTbP = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_ltN[4][s];
                        rTbN = (rLtScrt[s]+ttt_samp_distN[nm][im])*LtDxPar_noltN[4][s];
                      }
                      rr_ltbi=ttt_ltbi;
                      pop_frc=ttt_pop_frc;
                    }
                  } } }

                  Outputs[y][151] +=(V0[ag][2][0][im][nm][rg][na]+V0[ag][3][0][im][nm][rg][na])*(((1-pop_frc)*rTbP)+(pop_frc*(rTbP*rr_ltbi)))*LtTxParN[s][0] +
                    (V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*(((1- pop_frc)*rTbN) + (pop_frc*(1-(rTbP*rr_ltbi)))) *LtTxParN[s][0]; //all init
                  if(na>0) {
                    Outputs[y][152] += (V0[ag][2][0][im][nm][rg][na]+V0[ag][3][0][im][nm][rg][na])*(((1-pop_frc)*rTbP)+(pop_frc*(rTbP*rr_ltbi)))*LtTxParN[s][0] +
                      (V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*(((1- pop_frc)*rTbN) + (pop_frc*(1-(rTbP*rr_ltbi)))) *LtTxParN[s][0]; } // FB inits
                  if(rg==1) {
                    Outputs[y][153] +=  (V0[ag][2][0][im][nm][rg][na]+V0[ag][3][0][im][nm][rg][na])*(((1-pop_frc)*rTbP)+(pop_frc*(rTbP*rr_ltbi)))*LtTxParN[s][0] +
                      (V0[ag][1 ][0 ][im][nm][rg][na]+V0[ag][0 ][0 ][im][nm][rg][na])*(((1- pop_frc)*rTbN) + (pop_frc*(1-(rTbP*rr_ltbi)))) *LtTxParN[s][0]; } // high risk inits

                  Outputs[y][154] += (V0[ag][2][0][im][nm][rg][na]+V0[ag][3][0][im][nm][rg][na])*(((1-pop_frc)*rTbP)+(pop_frc*(rTbP*rr_ltbi)))*LtTxParN[s][0]; // inits with LTBI
                } } } } }
        for(int i=151; i<155; i++) { Outputs[y][i] = Outputs[y][i]*12; } // annualize


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
                    temp =  V0[ag][2 ][lt][im][nm][rg][na]*(MrslowN[ag][im]) + V0[ag][3 ][lt][im][nm][rg][na]*rfast; // All progression
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
        for(int i=155; i<187; i++) { Outputs[y][i] = Outputs[y][i]*12; } // annualize

        // NOTIFICATIONS, dead at diagnosis
        for(int nm=0; nm<4 ; nm++) {
          // if(nm > 2) {
          //   temp = muTbRF;
          // } else {
          //   temp = 0;
          // }
          temp = 0;
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
            // if(im > 2) { temp = muTbRF;  } else { temp = 0;  }
            temp=0;
            for(int ag=0; ag<11; ag++) {
              for(int lt=0; lt<2; lt++){
                for(int rg=0; rg<2; rg++) {
                  temp2 = V0[ag][4 ][lt][im][nm][rg][0]*(vTMortN[ag][4]+temp); //vTMort[ag][tb]
                  Outputs[y][215+ag] += temp2;   // dx by age (11)
                } } } } }
        for(int i=215; i<226; i++) { Outputs[y][i] = Outputs[y][i]*12; }

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
        for(int i=247; i<252; i++) { Outputs[y][i] = Outputs[y][i]*12; }

        ////////////////////// TB MORTALITY BY NATIVITY //////////////////////////////
        for(int ag=0; ag<11; ag++){
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

        ////total mortality
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

        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for(int nm=0; nm<4; nm++) {
                for(int im=0; im<4; im++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      Outputs[y][317+nm+(im*4)] += V1[ag][tb][lt][im][nm][rg][na];
                    } } } } } } }

        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {

                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      Outputs[y][333+nm+(im*4)+(ag*16)] += VMort[ag][tb][lt][im][nm][rg][na];
                    } } } } } } }

        for(int i=333; i<509; i++) { Outputs[y][i] = Outputs[y][i]*12; }

        for(int ag=0; ag<11; ag++) {
          Outputs[y][509+ag] = (Outputs[y][120+ag]/12)/(Outputs[y][2+ag]);
        }

        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++) {
              for(int nm=0; nm<4; nm++) {
                for(int im=0; im<4; im++) {
                  for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                      Outputs[y][520+nm+(ag*4)] += V1[ag][tb][lt][im][nm][rg][na];
                    } } } } } } }
        ///  NEW INFECTIONS + SUPER INFECTION BY AGE AND NAT GROUP
        for(int ag=0; ag<11; ag++) {
          for(int lt=0; lt<2 ; lt++) {
            for(int im=0; im<4 ; im++) {
              for(int nm=0; nm<4 ; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++){
                    if (na <1){
                      Outputs[y][564+ag] += (V0[ag][0][lt][im][nm][rg][na]+V0[ag][1][lt][im][nm][rg][na]+V0[ag][2][lt][im][nm][rg][na])*VLjkl[rg][na]*NixTrans[s];
                    } else {
                      Outputs[y][575+ag] += (V0[ag][0][lt][im][nm][rg][na]+V0[ag][1][lt][im][nm][rg][na]+V0[ag][2][lt][im][nm][rg][na])*VLjkl[rg][na]*NixTrans[s];}

                  }
                } } } } }
        ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
        for(int i=564; i<586; i++) { Outputs[y][i] = Outputs[y][i]*12; }

        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    if (na<1){
                      if (ag <3){
                        Outputs[y][586+nm+(im*4)] += V1[ag][tb][0][im][nm][rg][na];
                      } if(2<ag & ag<7){
                        Outputs[y][602+nm+(im*4)] += V1[ag][tb][0][im][nm][rg][na];
                      } if (6<ag & ag<11){
                        Outputs[y][618+nm+(im*4)] += V1[ag][tb][0][im][nm][rg][na];
                      }
                    } else {
                      if (ag<3) {
                        Outputs[y][634+nm+(im*4)] += V1[ag][tb][0][im][nm][rg][na];
                      } if(2<ag & ag<7){
                        Outputs[y][650+nm+(im*4)] += V1[ag][tb][0][im][nm][rg][na];
                      } if (6<ag & ag<11){
                        Outputs[y][666+nm+(im*4)] += V1[ag][tb][0][im][nm][rg][na];
                      } } }
                } } } } }
        ////////////////    COUNTS BY LATENT TREATMENT STATUS     //////////////////
        for(int ag=0; ag<11; ag++) {
          for(int tb=0; tb<6; tb++) {
            for(int im=0; im<4; im++) {
              for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                  for(int na=0; na<3; na++) {
                    Outputs[y][682] += V1[ag][tb][0][im][nm][rg][na];   //TREATMENT NAIVE
                  } } } } } }
      } ////end of mid-year results bracket
      ///////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////END MIDYEAR RESULTS//////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////
      ////////////////////////// REBALANCE THE POPULATION //////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////

      ////// need to define the current distribution of persons across the RG at this timestep
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
                  }
                } } } } }
      }

      // for(int ag=0; ag<11; ag++) {
      //   for(int tb=0; tb<6; tb++) {
      //     for(int lt=0; lt<2; lt++){
      //       for(int im=0; im<4; im++){
      //         for(int nm=0; nm<4; nm++){
      //           for(int rg=0; rg<2; rg++) {
      //             for(int na=0; na<3; na++) {
      //               if (std::any_of(V1[ag][tb][lt][im][nm][rg][na]<0)){
      //                 //Rcpp::Rcout << "after rblnc pop is negative at ag = " << ag << " tb = "<< tb << "lt = "<< lt << " im = " << im << " nm = " << nm << " rg = " << rg << " na = " << na << "/n";
      //              //   Rcpp::Rcout << "V1 is = "<<  V1[ag][tb][lt][im][nm][rg][na] << "\n";
      //              Rcpp::Rcout << "after rblnc pop is negative /n";
      //               }
      //             } } } } } } }


      //     //
      //     //    ///////////////////////////////////////////////////////////////////////////////////
      //     //    ///////////                       UPDATE V0 as V1                       ///////////
      // //     //    ///////////////////////////////////////////////////////////////////////////////////
      // for(int ag=0; ag<11; ag++) {
      //   for(int tb=0; tb<6; tb++) {
      //     for(int lt=0; lt<2; lt++){
      //       for (int im=0; im<4; im++){
      //         for (int nm=0; nm<4; nm++){
      //           for(int rg=0; rg<2; rg++) {
      //             for(int na=0; na<3; na++){
      //               // if ((reblnc == 1) & ((m==5)|(m==11))){
      //               //   V0[ag][tb][lt][im][nm][rg][na] = V2[ag][tb][lt][im][nm][rg][na];
      //               //   V1[ag][tb][lt][im][nm][rg][na] = V2[ag][tb][lt][im][nm][rg][na];
      //               // } else {
      //                 V0[ag][tb][lt][im][nm][rg][na] = V1[ag][tb][lt][im][nm][rg][na];
      //               // }
      //             } } } } } } }
      // for(int ag=0; ag<11; ag++) {
      //
      // for (int i=0; i<4; i++){
      //   for (int j=0; j<4; j++){
      //     temp_mat[i][j]=0; } }
      // mat_sum=0;
      // for(int nm=0; nm<4; nm++){
      //     for(int tb=0; tb<6; tb++) {
      //       for(int lt=0; lt<2; lt++){
      //
      //         for(int im=0; im<4; im++) {
      //           for(int rg=0; rg<2; rg++){
      //             for(int na=0; na<3; na++){
      //               temp_mat[nm][im]  += V1[ag][tb][lt][im][nm][rg][na];
      //             } } } } }
      // }
      // for(int nm=0; nm<4; nm++){
      //   for(int im=0; im<4; im++){
      //     mat_sum+=temp_mat[nm][im];
      //   } }
      // for(int nm=0; nm<4; nm++){
      //   for(int im=0; im<4; im++){
      //
      //     temp_mat2[nm][im] = temp_mat[nm][im]/mat_sum;
      //   } }
      // for(int nm=0; nm<4; nm++){
      //   for(int im=0; im<4; im++){
      // Rcpp::Rcout <<"at s "<<s <<"dist is = "<< temp_mat2[nm][im] << "at nm =" << nm <<" at im =" << im << "& ag = " <<ag<< "\n";
      // } }
      // }


    } //// end of month loop!//////////////////////////////////////////////////////////
  } //// end of year loop!///////////////////////////////////////////////////////////
  // for(int ag=0; ag<11; ag++) {
  //
  //   for (int i=0; i<4; i++){
  //     for (int j=0; j<4; j++){
  //
  //       temp_mat[i][j]=0; } }
  //   mat_sum=0;
  //   for(int nm=0; nm<4; nm++){
  //     for(int tb=0; tb<6; tb++) {
  //       for(int lt=0; lt<2; lt++){
  //         for(int im=0; im<4; im++) {
  //           for(int rg=0; rg<2; rg++){
  //             for(int na=0; na<3; na++){
  //               temp_mat[nm][im]  += V1[ag][tb][lt][im][nm][rg][na];
  //             } } } } }
  //   }
  //   for(int nm=0; nm<4; nm++){
  //     for(int im=0; im<4; im++){
  //       mat_sum+=temp_mat[nm][im];
  //     } }
  //   for(int nm=0; nm<4; nm++){
  //     for(int im=0; im<4; im++){
  //
  //       temp_mat[nm][im] = temp_mat[nm][im]/mat_sum;
  //       Rcpp::Rcout << temp_mat[nm][im] << "at nm = " << nm << "& im "<< im << "at ag "<< ag<< "\n";
  //
  //     } } }
  //

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

  // for (int i=0; i <12672; i++){
  // if (CheckV(i) <0){
  //   for(int ag=0; ag<11; ag++) {
  //     for(int tb=0; tb<6; tb++) {
  //       for(int lt=0; lt<2; lt++){
  //         for(int im=0; im<4; im++){
  //           for(int nm=0; nm<4; nm++){
  //             for(int rg=0; rg<2; rg++) {
  //               for(int na=0; na<3; na++) {
  //   Rcout <<"population is negative at ag = "<< ag << " tb = " << tb <<
  //     " lt = " << lt << " im = " << im << " nm = " << nm << " rg = " << rg << " & na = " << na << "\n";
  //               } } } } } } }
  // // } else {Rcout << "no negatives \n" ;
  // //        }
  // } }
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
