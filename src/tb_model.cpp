using namespace Rcpp;

//[[Rcpp::export]]

List cSim(
          int                 nYrs,
          int                 nRes,
          NumericMatrix       rDxt,
          std::vector<double> TxQualt,
          NumericMatrix       InitPop,
          NumericMatrix       Mpfast,
          NumericMatrix       ExogInf,
          NumericMatrix       MpfastPI,
          NumericMatrix       Mrslow,
          std::vector<double> rrSlowFB,
          double              rfast,
          double              RRcurDef,
          std::vector<double> VrSlfCur,
          NumericMatrix       vTMort,
          NumericMatrix       vNmMort, //vector of non-TB mortality
          double              muTbNm   //factor for comorbidity btw TB and non-TB,
          std::vector<double> Birthst,
          NumericMatrix       ImmNon,
          NumericMatrix       ImmLat,
          NumericMatrix       ImmAct,
          NumericMatrix       ImmFst,
          std::vector<double> TxExpAge,
          NumericMatrix       mubt,
          NumericMatrix       RelInf,
          std::vector<double> RelInfRg,
          std::vector<double> Vmix,
          NumericMatrix       vIsxtoIsy, //matrix for transitions within the Is dimension
          NumericMatrix       vNmxtoNmy, //matrix for transitions within the Nm dimension
          NumericMatrix       vrgxtorgy, //matrix for transitions within the rg dimension
          double              rArtDef,
          NumericMatrix       rRFt,//rate of risk factor (population) of interest over time
          std::vector<double> rEmmigFB,
          NumericMatrix       rIntvInit,//rate of intervention for RF of interest over time
          NumericMatrix       TxMat,
          double              TunTxMort,
          std::vector<double> rDeft,
          std::vector<double> rDeftH,
          std::vector<double> LtTxPar,
          NumericMatrix       LtDxPar,
          std::vector<double> RelInfHivt,
          std::vector<double> RRdxAge,
          double              rRecov,
          double              pImmScen,
          std::vector<double> EarlyTrend,
          NumericMatrix       EffLtX,
          double              EffLt,
          std::vector<double> dLtt,
          std::vector<double> NixTrans
          ) {
    
    ////////////////////////////////////////////////////////////////////////////////
    ////////    BELOW IS A LIST OF THE VARIABLES CREATED INTERNALLY IN MODEL   /////
    ////////////////////////////////////////////////////////////////////////////////
    int           ti;
    int           ti2;
    int           ti3;
    int           ti4;
    int           ti5;
    int           s;
    int           extrV[10];
    int           h2;
    int           r2;
    double        ExogInfN[ExogInf.nrow()][ExogInf.ncol()];
    double        InitPopN[InitPop.nrow()][InitPop.ncol()];
    double        InitPopZ[InitPop.nrow()][InitPop.ncol()];
    double        MpfastN[Mpfast.nrow()][Mpfast.ncol()];
    double        MpslowN[Mpfast.nrow()][Mpfast.ncol()];
    double        MpfastPIN[MpfastPI.nrow()][MpfastPI.ncol()];
    double        MpslowPIN[MpfastPI.nrow()][MpfastPI.ncol()];
    double        MrslowN[Mrslow.nrow()][Mrslow.ncol()];
    double        MpSmPosN[Mpfast.nrow()][Mpfast.ncol()];
    double        vTMortN[vTMort.nrow()][vTMort.ncol()];
    double        vNmMortN[vNmMort.nrow()][vNmMort.ncol()];
   // double        vIsxtoIsyN[vIsxtoIsy.nrow()][vIsxtoIsy.ncol()];
   // double        vNmxtoNmyN[vNmxtoNmy.nrow()][vNmxtoNmy.ncol()];
   // double        vrgxtorgyN[vrgxtorgy.nrow()][vrgxtorgy.ncol()];
    double        ImmNonN[ImmNon.nrow()][ImmNon.ncol()];
    double        ImmLatN[ImmLat.nrow()][ImmLat.ncol()];
    double        ImmFstN[ImmFst.nrow()][ImmFst.ncol()];
    double        ImmActN[ImmAct.nrow()][ImmAct.ncol()];
    double        mubtN[mubt.nrow()][mubt.ncol()];
    double        RelInfN[RelInf.nrow()][RelInf.ncol()];
  //  double        rIntvInitN[rIntvInit.nrow()][rIntvInit.ncol()];
    double        rDxtN[rDxt.nrow()][rDxt.ncol()];
    double        TxMatN[TxMat.nrow()][TxMat.ncol()];
    double        LtDxParN[LtDxPar.nrow()][LtDxPar.ncol()];
    double        EffLtXN[EffLtX.nrow()][EffLtX.ncol()];
    double        TxMatZ[23][10][2];
    double        temp;
    double        temp2;
    double        temp3;
    double        temp4V[11][5];
    double        rTbP;
    double        rTbN;
    double        p2ndL;
    double        Outputs[nYrs][nRes];
    double        V0[11][11][5][3][5][4];
    double        V1[11][11][5][3][5][4];
    double        VMort[11][11][5][3][5][4];
    double        Vdx[11][11][5][3][5][4];
    double        VNkl[3][2];
    double        VGjkl[3][2]; ///removed drug resistance dimension
    double        Vjaf[6];     ///removed drug resistance dimension
    double        VLjk[3][2];  ///removed drug resistance dimension
    NumericMatrix Outputs2(nYrs,nRes);
    
    ///////////////////////////////////////////////////////////////////////////////
    ///////                            INITIALIZE                             /////
    ///////////////////////////////////////////////////////////////////////////////
    for(int i=0; i<InitPop.nrow(); i++) {
        for(int j=0; j<InitPop.ncol(); j++) {
            InitPopN[i][j] = InitPop(i,j);
        } }
    for(int i=0; i<ExogInf.nrow(); i++) {
        for(int j=0; j<ExogInf.ncol(); j++) {
            ExogInfN[i][j] = ExogInf(i,j);
        } }
    for(int i=0; i<Mpfast.nrow(); i++) {
        for(int j=0; j<Mpfast.ncol(); j++) {
            MpfastN[i][j]   = Mpfast(i,j);
            MpfastPIN[i][j]   = MpfastPI(i,j);
            MpslowN[i][j]   = 1-Mpfast(i,j);
            MpslowPIN[i][j]   = 1-MpfastPI(i,j);
            MrslowN[i][j]   = Mrslow(i,j);
        } }
    for(int i=0; i<LtDxPar.nrow(); i++) {
        for(int j=0; j<LtDxPar.ncol(); j++) {
            LtDxParN[i][j] = LtDxPar(i,j);
        } }
    for(int i=0; i<EffLtX.nrow(); i++) {
        for(int j=0; j<EffLtX.ncol(); j++) {
            EffLtXN[i][j] = EffLtX(i,j);
        } }
    for(int i=0; i<vTMort.nrow(); i++) {
        for(int j=0; j<vTMort.ncol(); j++) {
            vTMortN[i][j] = vTMort(i,j);
        } }
///this was HIV mortality temporarily replaced with RF of interest
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
    for(int i=0; i<mubt.nrow(); i++) {
        for(int j=0; j<mubt.ncol(); j++) {
            mubtN[i][j] = mubt(i,j);
        } }
    for(int i=0; i<RelInf.nrow(); i++) {
        for(int j=0; j<RelInf.ncol(); j++) {
            RelInfN[i][j] = RelInf(i,j);
        } }
//////was the rate of HIV and has been replaced as rate of risk factor of interest////
    for(int i=0; i<rRFt.nrow(); i++) {
        for(int j=0; j<rRFt.ncol(); j++) {
            rRFtN[i][j][0] = rRFt(i,j);
            rRFtN[i][j][1] = rRFt(i,j)*RFHrPar;
        } }
    for(int i=0; i<rIntvInit.nrow(); i++) {
        for(int j=0; j<rIntvInit.ncol(); j++) {
            rIntvInitN[i][j] = rIntvInit(i,j);
        } }
    for(int i=0; i<TxMat.nrow(); i++) {
        for(int j=0; j<TxMat.ncol(); j++) {
            TxMatN[i][j] = TxMat(i,j);
        } }
    for(int i=0; i<rDxt.nrow(); i++) {
        for(int j=0; j<rDxt.ncol(); j++) {
            rDxtN[i][j] = rDxt(i,j);
        } }
    for(int i=0; i<HrEntEx.nrow(); i++) {
        for(int j=0; j<HrEntEx.ncol(); j++) {
            HrEntExN[i][j][0] = HrEntEx(i,j);
            HrEntExN[i][j][1] = HrEntEx(i,j)*HivHrPar;
        } }
    for(int i=0; i<22; i++) {
        for(int j=0; j<2; j++) {
            TxMatZ[i][j] = 0.0;
        } } }
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
for(int i=0; i<3; i++) {   //high risk low risk foreign us born (expand)
    for(int j=0; j<2; j++) { //????
        VNkl[i][j] = 0;
} }
for(int j=0; j<3; j++) {     //high risk low risk foreign us born (expand)
    for(int k=0; k<2; k++) { //????
        VGjkl[j][k] = 0;
        VLjkl[j][k] = 0;
} } }
for(int i=0; i<5; i++) { //drug resistance?
    for(int j=0; j<6; j++) { ///????
        Vjaf[i][j] = 0;
} }
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
for(int j=0; j<2; j++) {
    if(j==0) {
        temp = rDeft[0];
    } else {
        temp = rDeftH[0];
    }
    //////// TREATMENT EFFICACY UPDATED FOR TREATMENT QUALITY //////////////////////
    TxMatZ[1][j] = TxMatN[1][j]*TxQualt[0];
    ///////// RATE OF TREATMENT EXIT TO CURE (LS) //////////////////////////////////
    TxMatZ[7][j] = TxMatN[0][j]*TxMatZ[1][j] + temp*TxMatZ[1][j]*RRcurDef;
    //////// RATE OF TREATMENT EXIT TO FAILURE (IN/IP) /////////////////////////////
    TxMatZ[8][j] = TxMatN[0][j]*(1-TxMatZ[1][j]) + temp*(1-TxMatZ[1][j]*RRcurDef);
} }
////////////////////////////////////////////////////////////////////////////////
//////                             StatList                                /////
//////                             BURN IN                                 /////
//////                           Populate model                            /////
////////////////////////////////////////////////////////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int is=0; is<4; is++) {
        for(int nm=0; nm<4; nm++){
            for(int rg=0; rg<3; rg++){
////////////////////        UNINFECTED/SUSCEPTIBLE POP /////////////////////////
                V0[ag][0][0][0][0][0][0] = InitPopN[ag][0]*0.40*(1-p_HR); //US born
                V0[ag][0][0][0][0][0][2] = InitPopN[ag][1]*0.40;          //new non-US born
/////////////////////////   LATENT SLOW INFECTED POP  //////////////////////////
                V0[ag][2][0][0][0][0][0] = InitPopN[ag][0]*0.60*(1-p_HR);
                V0[ag][2][0][0][0][0][2] = InitPopN[ag][1]*0.60;
            }
//////create a 2nd array with same dimensions as V0 & populate w/ same values//
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int is=0; is<4; is++) {
            for(int nm=0; nm<4; nm++){
                for(int rg=0; rg<2; rg++){
                    for (int na=0; na<3; na++){
                        V1[ag][tb][0][is][nm][rg][na]  = V1[ag][tb][0][is][nm][rg][na];
} } } } } }
////////////////////////RUN THE MODEL FOR 3000 MONTHS /////////////////////////
for(int m=0; m<3001; m++) {
/////////////////////////////////START BURN IN//////////////////////////////////
////////////////////////////////////BIRTHS//////////////////////////////////////
///////////USE DISTRIBUTION TO POPULATE THE MODEL ACROSS RISK GROUPS////////////
for (i=0, i<4; i++){
    for (j=0;i<4; j++){
        for (k=0;k<4;k++){
            V0[0][0][0][i][j][k][0]  += Birthst[0]*dist[i,j]*rg[k] //check this
} } }
//////////////////////////////////IMMIGRATION///////////////////////////////////
////////  Assume 0% RF prev, 0% NONTBMORT Assume All Pansens   ////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        V1[ag][0][0][0][0][0][1]  += ImmNonN[0][ag];      // NO TB
        V1[ag][2][0][0][0][0][1]  += ImmLatN[0][ag];      // LATENT SLOW TB
        V1[ag][3][0][im][nm][rg][1]    += (TBImm[ag][0][s][1] + TBImm[ag][1][s][1]);   // LATENT FAST
        V1[ag][4][0][im][nm][rg][1]    += (TBImm[ag][0][s][0] + TBImm[ag][1][s][0]);   //ACTIVE TB
} }
///////////////////////////////EMMIGRATION//////////////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        V1[ag][tb][0][0][0][0][1]  -= V0[ag][tb][0][0][0][1]*rEmmigFB[0];   // FB1
        V1[ag][tb][0][0][0][0][2]  -= V0[ag][tb][0][0][0][2]*rEmmigFB[1];   // FB2
} }
////////////////////////////////MORTALITY//////////////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int rg=0; rg<2; rg++) {
        for(int tb=0; tb<6; tb++) {
            V1[ag][tb][0][0][0][rg][0]  -= V0[ag][tb][0][0][rg][0]*(mubtN[0][ag]*RRmuHR[rg]+vTMortN[ag][tb]);
}
////////////////          MORTALITY WITH TB TREATMENT         ////////////////////
            V1[ag][5 ][0][0][0][rg][0]  -= V0[ag][5 ][0][0][0]rg[0]*(mubtN[0][ag]*RRmuHR[rg]+vTMortN[ag][5 ]*pow(1.0-TxMatZ[1][0][0],TunTxMort)); //check the mortality in param
} }
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
        temp = V0[ag][tb][0][0][0][1]/24;
        V1[ag][tb][0][0][0][0][1]  -= temp;
        V1[ag][tb][0][0][0][0][2]  += temp;
} }
//////////////////////////// HIGH-RISK ENTRY/EXIT ////////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        temp  = V0[ag][tb][0][0][0][0][0]*HrEntExN[ag][0][0];
        temp2 = V0[ag][tb][0][0][0][0][1]*HrEntExN[ag][1][0];
//THESE CODES WERE UPDATED, BUT REMAIN ALMOST THE SAME
        V1[ag][tb][0][0][0][0][0]  += temp2-temp;
        V1[ag][tb][0][0][0][0][1]  += temp-temp2;
} }
////////////////////////////  TRANSMISSION RISK  ////////////////////////////////
for(int i=0; i<3; i++) {
    VNkl[i][0] = 0;
    VGjkl[0][i][0] = 0;   // set to zero
}
// Step 1
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++){
///////// RISK FACTOR FREE & LOW RISK US BORN
            VNkl[0][0][0]  += V0[ag][tb][lt][0][0][0][0];
///////// RISK FACTOR FREE & HIGH RISK US BORN
            VNkl[0][1][0]  += V0[ag][tb][lt][0][0][0][0];
///////// RISK FACTOR FREE & LOW RISK NON US BORN
            VNkl[1][0][0]  += V0[ag][tb][lt][0][0][0][1] + V0[ag][tb][lt][0][0][0][2];
///////// RISK FACTOR FREE & HIGH RISK NON US BORN
            VNkl[1][1][0]  += V0[ag][tb][lt][0][0][1][1] + V0[ag][tb][lt][0][0][1][2];
} } }
// Step 2  (active TB)
for(int ag=0; ag<11; ag++) {
///////// RISK FACTOR FREE & LOW RISK US BORN
        VGjkl[0][0][0][0]  += V0[ag][4][0][0][0][0][0]*RelInfN[tb][0];
///////// RISK FACTOR FREE & HIGH RISK US BORN
        VGjkl[0][0][1][0]  += V0[ag][4][0][0][0][1][1]*RelInfN[tb][0];
///////// RISK FACTOR FREE & LOW RISK NON US BORN
        VGjkl[0][1][0][0]  += (V0[ag][4][0][0][0][0][1] + V0[ag][4][0][0][0][0][2])*RelInfN[tb][0];
///////// RISK FACTOR FREE & HIGH RISK NON US BORN
        VGjkl[0][1][1][0]  += (V0[ag][4][0][0][0][1][1] + V0[ag][4][0][0][0][1][2])*RelInfN[tb][0];

}
// Step 2 (treated TB)
// No contribution to force of infection
                            
// Step 3
    Vjaf[0][0]  = (RelInfRg[0]*VGjkl[0][0][0]+RelInfRg[1]*VGjkl[0][1][0]*Vmix[1]+RelInfRg[2]*VGjkl[0][2][0]*Vmix[2]) /(RelInfRg[0]*VNkl[0][0]    +RelInfRg[1]*VNkl[1][0]*Vmix[1]    +RelInfRg[2]*VNkl[2][0]*Vmix[2]+1e-12);
    Vjaf[0][1]  = VGjkl[0][2][0] / (VNkl[2][0]+1e-12);
    Vjaf[0][2]  = VGjkl[0][1][0] / (VNkl[1][0]+1e-12);
// Step 4 CREATE LAMBDA FORCE OF INFECTION
/// VLjkl[dr][na][rg][hv]
    VLjkl[0 ][0 ][0 ][0 ]  = RelInfRg[0]*Vjaf[0][0];
///////// RISK FACTOR FREE & HIGH RISK US BORN
    VLjkl[0 ][0 ][1 ][0 ]  = RelInfRg[1]*((1-Vmix[1])*Vjaf[0][2]+Vmix[1]*Vjaf[0][0]);
///////// RISK FACTOR FREE & LOW RISK NON US BORN
    VLjkl[0 ][1 ][0 ][0 ]  = RelInfRg[2]*((1-Vmix[2])*Vjaf[0][1]+Vmix[2]*Vjaf[0][0])+ExogInfN[0][0];
///////// RISK FACTOR FREE & HIGH RISK NON US BORN
    VLjkl[0 ][1 ][1 ][0 ]  = RelInfRg[2]*((1-Vmix[2])*Vjaf[0][1]+Vmix[2]*Vjaf[0][0])+ExogInfN[0][0];
///////////////////////////////    INFECTION   /////////////////////////////////
for(int rg=0; rg<2; rg++) {
    for (int na=0; na<3; na++){
        if(na<2) {
            r2 = rg;
        } else {
            r2 = 2;
} }
for(int ag=0; ag<11; ag++) {
///////////////////////////////   SUCEPTIBLE  /////////////////////////////////
    temp = V0[ag][0][0][0][rg][na]*VLjkl[0][na][r2][0]*EarlyTrend[m];
//////////////////////////// REMOVE FROM SUSCEPTIBLE //////////////////////////
    V1[ag][0][0][0][0][rg][na]  -= temp;
//////////////////////////////// LATENT TB SLOW ///////////////////////////////
    V1[ag][2][0][0][0][rg][na]  += temp*MpslowN[ag][0];
//////////////////////////////// LATENT TB FAST ///////////////////////////////
    V1[ag][3][0][0][0][rg][na]  += temp*MpfastN[ag][0];
///////////////////////////////////////////////////////////////////////////////
                                    
/////////////////////////////// SUPER-INFECTION SP ////////////////////////////
    temp = V0[ag][1][0][0][0][rg][na]*VLjkl[0][na][r2][0];
    V1[ag][1][0][0][0][rg][na] -= temp;
    V1[ag][2][0][0][0][rg][na] += temp*MpslowPIN[ag][0];
    V1[ag][3][0][0][0][rg][na] += temp*MpfastPIN[ag][0];
///////////////////////////////////////////////////////////////////////////////
                                    
/////////////////////////////// SUPER-INFECTION LS ////////////////////////////
    temp = V0[ag][2][0][0][0][rg][na]*VLjkl[0][na][r2][0];
    V1[ag][2][0][0][0][rg][na]  -= temp;
    V1[ag][2][0][0][0][rg][na]  += temp*MpslowPIN[ag][0];
    V1[ag][3][0][0][0][rg][na]  += temp*MpfastPIN[ag][0];
}
////////////////////////////////////////////////////////////////////////////////
                            
///////////////////////////////////BREAKDOWN///////////////////////////////////
///////////////////////for all age groups, risk groups/////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int rg=0; rg<2; rg++) {
        for (int na=0; na<3; na++){
        temp  = V0[ag][2][0][0][0][rg][na]*MrslowN[ag][0]*rrSlowFB[rg];
        temp2 = V0[ag][3][0][0][0][rg][na]*rfast;
        V1[ag][2][0][0][0][rg][na]  -= temp;
        V1[ag][3][0][0][0[rg][na]  -= temp2;
        V1[ag][4][0][0][0][rg][na]  += (temp+temp2);
} } }
///////////////////////////     LATENT SLOW TO SAFE /////////////////////////////
///////////////////////for all age groups, risk groups///////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int rg=0; rg<2; rg++) {
        temp  = V0[ag][2][0][0][0][rg][na]*rRecov;
        V1[ag][2][0][0][0][rg][na]  -= temp;
        V1[ag][1][0][0][0][rg][na]  += temp;
} }
                    
////////////////////////////////// SELF CURE/////////////////////////////////////
/////////////////for all age groups, risk groups, only TB 4 /////////////////////
for(int ag=0; ag<11; ag++) {
    for(int rg=0; rg<2; rg++) {
        for (int na=0; na<3; na++){
        temp = V0[ag][4][0][0][0][rg][na]*VrSlfCur[0];
               V1[ag][4][0][0][0][rg][na]  -= temp;
               V1[ag][2 ][0][0][0][rg][na]  += temp;
} } }
////////////////////TB DIAGNOSIS AND TX INITIATION PUBLIC ///////////////////////
///////////////////////for all age groups, living cond///////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int rg=0; rg<2; rg++) {
       if(rg!=1) {
            ti = 0;
     } else {
            ti = 1;
     }
    temp  = V0[ag][4 ][0][0][0][0][rg]*rDxtN[0][ti]/RRdxAge[ag]/EarlyTrend[m];
            V1[ag][4 ][0][0][0][0][rg]  -= temp;
            V1[ag][5 ][0][0][0][0][rg]  += temp;
} }
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
            temp  = V0[ag][5][0][0][0][rg][na]*TxMatZ[7][0][ti2];
            V1[ag][5][0][0][0][rg][na]  -= temp;
            V1[ag][2][0][0][0][rg][na]  += temp;
//////////FAILURES(INCLUDING TREATMENT DEFAULT)//////////////////////////////////
            temp  = V0[ag][5][0][0][0][rg][na]*TxMatZ[8][0][ti2];
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
            InitPopZ[ag][0]  += V1[ag][tb][0][0][0][0][0]+V1[ag][tb][0][0][0][0][1];
            InitPopZ[ag][1]  += V1[ag][tb][0][0][0][0][2]+V1[ag][tb][0][0][0][0][3];
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
    V0[ag][tb][0][0][0][rg][na]  = V1[ag][tb][0][0][0][rg][na];
} } }
}  /////////////////////  END OF MONTH LOOP  //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////////////////////////////END BURN IN///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
NumericVector  CheckV0(14784);
    for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; lt++){
                for(int im=0; im<4; im++){
                    for(int nm=0; nm<4; nm++){
                        for(int rg=0; rg<2; rg++) {
                            for(int na=0; na<2; na++){
CheckV0(ag+tb*11+lt*66+tx*132+im*528+nm*2112+rg*4928+na*12672) = V1[ag][tb][lt][im][nm][rg][na];
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
                                ///// NEED TO UPDATE THE TXMATZ IN THE PARAM FILE /////////////////////////////
                                ///// TxMatZ: 0=completion rate, 1 = tx success, 2:6 = AR probabilities, //////
                                ///// 7 = rate of exit to cure 8 = adj factor for AR with tx completion, //////
                                ///// 9 = adj factor for AR with default, 10:14 = exit rates to active TB;/////
                                ///// 15:19 = exit rates to retx; 20 = sum of AR for active TB, ///////////////
                                ///// 21 = sum of AR for retx /////////////////////////////////////////////////
                                ///////////////////////////////////////////////////////////////////////////////
                                for(int j=0; j<10; j++) {
                                    for(int k=0; k<2; k++) {
                                        if(k==0) { temp = rDeft[s];
                                        } else { temp = rDeftH[s]; } ///UPDATE rDeftH/////////////////////////////
                                        //////// TREATMENT EFFICACY UPDATED FOR TREATMENT QUALITY //////////////////////
                                        TxMatZ[1][j][k]            = TxMatN[1][j]*TxQualt[s];
                                        ///////// RATE OF TREATMENT EXIT TO CURE (LS) //////////////////////////////////
                                        TxMatZ[7][j][k]            = TxMatN[0][j]*TxMatZ[1][j][k] + temp*TxMatZ[1][j][k]*RRcurDef;
                                        ///REMOVE TXMATZ[8-9] AS THEY ARE AR PROBABILITIES
                                        for(int i=0; i<5; i++) {  // Rates of exit to active disease from completion, potentially with AR
                                            TxMatZ[10+i][j][k]       = (TxMatN[0][j]*TxMatZ[8][j][k]*(1-pReTx[s]) + temp*TxMatZ[9][j][k])*TxMatN[2+i][j];
                                            TxMatZ[15+i][j][k]       = TxMatN[0][j]*TxMatN[2+i][j]*TxMatZ[8][j][k] * pReTx[s]; } // rate_complete * p(AR|complete,fail,0) * Adj_factor * p(retx|complete)
                                        TxMatZ[20][j][k]           = TxMatZ[10][j][k]+TxMatZ[11][j][k]+TxMatZ[12][j][k]+TxMatZ[13][j][k]+TxMatZ[14][j][k]; // Sum total of AR exit rates to ACTIVE TB
                                        TxMatZ[21][j][k]           = TxMatZ[15][j][k]+TxMatZ[16][j][k]+TxMatZ[17][j][k]+TxMatZ[18][j][k]+TxMatZ[19][j][k]; // Sum total of AR exit rates to RETX
                                        TxMatZ[10+extrV[j]][j][k]  = (TxMatN[0][j]*(1.0-TxMatZ[1][j][k])*(1-pReTx[s]) + temp*(1.0-TxMatZ[1][j][k]*RRcurDef)) - TxMatZ[20][j][k]; // exit to active TB, no AR
                                        TxMatZ[15+extrV[j]][j][k]  = TxMatN[0][j]*(1.0-TxMatZ[1][j][k])*pReTx[s] - TxMatZ[21][j][k]; // exit to RETX, no AR
                                        TxMatZ[22][j][k]           = TxMatN[0][j]*(1-(1.0-TxMatZ[1][j][k])*pReTx[s]); // p(tx completion)
                                    } }
                                /////////////////////////////////////BIRTHS//////////////////////////////////////
                                /////ALL BIRTHS ENTER IN 0-4 AGE, UNINF&SUSC, PANSENSITIVE, TREAT.NAIVE, HIV-NEG
                                /////LOW RISK GROUP BIRTHS//////////////////////////////////////////////////////
                                V1[0][0][0][0][0][0][0]  += Birthst[s]*(1-p_HR);
                                /////HIGH RISK GROUP BIRTHS//////////////////////////////////////////////////////
                                V1[0][0][0][0][0][1][0]  += Birthst[s]*p_HR;
                                ///////////////////////////////// IMMIGRATION ///////////////////////////////////
                                for(int ag=0; ag<11; ag++) {
                                    V1[ag][0][0][0 ][0 ][0 ][1]    += ImmNonN[s][ag];      // NO TB
                                    V1[ag][2][0][0 ][0 ][0 ][1]    += ImmLatN[s][ag];      // LATENT TB
                                    V1[ag][3][0][im][nm][rg][1]    += TBImm[ag][0][s][1];   // tx naive, LAT FAST
                                    V1[ag][3][1][im][nm][rg][1]    += TBImm[ag][1][s][1];   // tx exp, LAT FAST
                                    V1[ag][4][0][im][nm][rg][1]    += TBImm[ag][0][s][0];   // tx naive, active
                                    V1[ag][4][1][im][nm][rg][1]    += TBImm[ag][1][s][0];   // tx exp, active
                                    }
                                /////////////////////////////////  IMMIGRATION ///////////////////////////////////
                                for(int ag=0; ag<11; ag++) {
                                    for(int tb=0; tb<6; tb++) {
                                        for(int lt=0; lt<2; lt++){
                                            for(int im=0; im<4; im++){
                                                for(int nm=0; nm<4; nm++){
                                                    for(int rg=0; rg<2; rg++) {
                                                        V1[ag][tb][lt][im][nm][rg][1]  -= V0[ag][tb][lt][im][nm][rg][1]*rEmmigFB[0];      // FB1
                                                        V1[ag][tb][lt][im][nm][rg][2]  -= V0[ag][tb][lt][im][nm][rg][2]*rEmmigFB[1];      // FB2
                                                    } } } } } } }
                            
/////////////////////////////////  EMIGRATION ///////////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2; lt++){
            for(int im=0; im<5 ; im++) {
                if(im==3) {
                    temp = muTbRF;
                } else { temp = 0;
                }
for(int rg=0; rg<2; rg++) {
for(int na=0; na<3; na++) {
////////////////////////UNINFECTED, SUSCEPTIBLE//////////////////////////////////
    VMort[ag][0 ][lt][im][nm][rg][na]  = V0[ag][0 ][lt][im][nm][rg][na]*
    (mubtN[s][ag]*RRmuHR[rg][na]+vNmMortN[ag][nm] );
////////////////////////UNINFECTED, PART. IMMUNE/////////////////////////////////
                                                    VMort[ag][1 ][lt][im][nm][rg][na]  = V0[ag][1 ][lt][im][nm][rg][na]*
                                                    (mubtN[s][ag]*RRmuHR[rg][na]+vNmMortN[ag][nm] );
////////////////////////    LATENT TB SLOW      /////////////////////////////////
                                                    VMort[ag][2 ][lt][im][nm][rg][na]  = V0[ag][2 ][lt][im][nm][rg][na]*
                                                    (mubtN[s][ag]*RRmuHR[rg][na]+vNmMortN[ag][nm] );
////////////////////////    LATENT TB FAST      /////////////////////////////////
                                                    VMort[ag][3 ][lt][im][nm][rg][na]  = V0[ag][3 ][lt][im][nm][rg][na]*
                                                    (mubtN[s][ag]*RRmuHR[rg][na]+vNmMortN[ag][nm] );
////////////////////////      ACTIVE TB         /////////////////////////////////
                                                    VMort[ag][4 ][lt][im][nm][rg][na]  = V0[ag][4 ][lt][im][nm][rg][na]*
                                                    (mubtN[s][ag]*RRmuHR[rg][na]+vNmMortN[ag][nm]+vTMortN[ag][4 ]+temp );
////////////////////////    TB TREATMENT        /////////////////////////////////
                                                    VMort[ag][5 ][lt][im][nm][rg][na]  = V0[ag][5 ][lt][im][nm][rg][na]*
                                                    (mubtN[s][ag]*RRmuHR[rg][na]+vNmMortN[ag][nm]+(vTMortN[ag][5 ]+temp)*pow(1.0-TxMatZ[1][dr*2+0][0],TunTxMort));
                                                    
for(int tb=0; tb<6; tb++) {
    V1[ag][tb][lt][im][nm][rg][na]  -= VMort[ag][tb][lt][im][nm][rg][na];
}
} } } } }
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
                            temp = V0[ag][tb][lt][im][nm][rg][na]/temp2;
                                V1[ag  ][tb][lt][im][nm][rg][na]  -= temp;
                                    V1[ag+1][tb][lt][im][nm][rg][na]  += temp;
} } } } } } }
///////////////////////// NEW FB -> ESTABLISHED FB ///////////////////////////////
///////////////////////// TWO YEARS FOR TRANSITION ///////////////////////////////
for(int ag=0; ag<11; ag++)
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++){
            for (int im=0; im<4; im++){
                for (int nm=0; nm<4; nm++){
                    for(int rg=0; rg<2; rg++) {
                        temp = V0[ag][tb][lt][im][nm][rg][1]/24;
                        V1[ag][tb][lt][im][nm][rg][1]  -= temp;
                        V1[ag][tb][lt][im][nm][rg][2]  += temp;
                    } } } } }
//////////////////////////// HIGH-RISK ENTRY/EXIT ////////////////////////////////
///////////////////REVIEW AND UPDATE PROPERLY; JUST EDITED INDEXES////////////////
for(int ag=0; ag<11; ag++)
    for(int tb=0; tb<6; tb++) {
        for(int tx=0; tx<2; tx++){
            temp  = V0[ag][tb][dr][tx][0][0]*HrEntExN[ag][0][0];
            temp2 = V0[ag][tb][dr][tx][0][1]*HrEntExN[ag][1][0];
            V1[ag][tb][dr][tx][0][0]  += temp2-temp;
            V1[ag][tb][dr][tx][0][1]  += temp-temp2;
            for(int im=1; im<5 ; im++) {
                temp  = V0[ag][tb][dr][tx][im][0]*HrEntExN[ag][0][1];
                temp2 = V0[ag][tb][dr][tx][im][1]*HrEntExN[ag][1][0]; // note this stays the same as HIV neg ([0] suffix)
                V1[ag][tb][dr][tx][im][0]  += temp2-temp;
                V1[ag][tb][dr][tx][im][1]  += temp-temp2;
} } } }
////////////////////////////  TRANSMISSION RISK  ////////////////////////////////
for(int i=0; i<4; i++) {
    for(int j=0; j<2; j++) {
        VNkl[i][j] = 0; // set to zero
        VGjkl[i][j] = 0; // set to zero //removed dr
} }
////////////////////////////          Step 1         ////////////////////////////
for(int ag=0; ag<11; ag++) {
 	for(int tb=0; tb<6; tb++) {
    	for(int lt=0; lt<2; tx++) {
            for(int nm=0; nm<4; nm++){
////////////////////////////    LOW RISK, US BORN    ////////////////////////////
                VNkl[0][0]  += V0[ag][tb][lt][0][nm][0][0];
////////////////////////////   HIGH RISK, US BORN    ////////////////////////////
                VNkl[1][0]  += V0[ag][tb][lt][0][nm][1][0];
////////////////////////////  LOW RISK, NON US BORN  ////////////////////////////
                VNkl[2][0]  += V0[ag][tb][lt][0][nm][0][1] + V0[ag][tb][lt][0][nm][0][2];
//////////////////////////// HIGH RISK, NON US BORN  ////////////////////////////
                VNkl[3][0]  += V0[ag][tb][lt][0][nm][1][1] + V0[ag][tb][lt][0][nm][1][2];
                for(int im=1; im<4 ; nm++) {
////////////////////////////    LOW RISK, US BORN    ////////////////////////////
                VNkl[0][1]  += V0[ag][tb][lt][im][nm][0][0];
////////////////////////////   HIGH RISK, US BORN    ////////////////////////////
                VNkl[1][1]  += V0[ag][tb][lt][im][nm][1][0];
////////////////////////////  LOW RISK, NON US BORN  ////////////////////////////
                VNkl[2][1]  += V0[ag][tb][lt][im][nm][0][1] + V0[ag][tb][lt][im][nm][0][2];
//////////////////////////// HIGH RISK, NON US BORN  ////////////////////////////
                VNkl[3][1]  += V0[ag][tb][lt][im][nm][1][1] + V0[ag][tb][lt][im][nm][1][2];
                }
} } } }
/////////////////////////// /   Step 2  (ACTIVE TB)   ////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=4; tb<6; tb++) {
        for(int lt=0; lt<2 ; lt++) {
            for(int nm=0; nm<4; nm++){
////////////////////////////    LOW RISK, US BORN    ////////////////////////////
                VGjkl[0][0]  += V0[ag][tb][lt][0][nm][0][0]*RelInfN[tb];
////////////////////////////   HIGH RISK, US BORN    ////////////////////////////
                VGjkl[1][0]  += V0[ag][tb][lt][0][nm][1][0]*RelInfN[tb];
////////////////////////////  LOW RISK, NON US BORN  ////////////////////////////
                VGjkl[2][0]  += (V0[ag][tb][lt][0][nm][0][1] + V0[ag][lt][0][nm][0][2])*RelInfN[tb];
//////////////////////////// HIGH RISK, NON US BORN  ////////////////////////////
                VGjkl[3][0]  += (V0[ag][tb][lt][0][nm][1][1] + V0[ag][tb][lt][0][nm][1][2])*RelInfN[tb];
                    for(int im=1; im<5 ; im++) {
////////////////////////////    LOW RISK, US BORN    ////////////////////////////
                VGjkl[0][1]  += V0[ag][tb][dr][tx][im][0]*RelInfN[tb];
////////////////////////////   HIGH RISK, US BORN    ////////////////////////////
                VGjkl[1][1]  += V0[ag][tb][dr][tx][im][0]*RelInfN[tb];
////////////////////////////  LOW RISK, NON US BORN  ////////////////////////////
                VGjkl[2][1]  += (V0[ag][tb][lt][im][nm][0][1] + V0[ag][tb][lt][im][nm][0][2]])*RelInfN[tb];
//////////////////////////// HIGH RISK, NON US BORN  ////////////////////////////
                VGjkl[3][1]  += (V0[ag][tb][lt][im][nm][1][1] + V0[ag][tb][lt][im][nm][1][2])*RelInfN[tb];
                    }
} } } }
////////////////////////////   Step 2 (TREATED TB)   ////////////////////////////
////////////////////  No contribution to force of infection  ////////////////////
                            
////////////////////////////          Step 3         ////////////////////////////
 for(int dr=0; dr<5; dr++) {
                                Vjaf[dr][0]  = (RelInfRg[0]*(VGjkl[dr][0][0]+RelInfHivt[s]*Vmix[0]*VGjkl[dr][0][1])+
                                                RelInfRg[1]*(VGjkl[dr][1][0]+Vmix[0]*VGjkl[dr][1][1])*Vmix[1]+
                                                RelInfRg[2]*(VGjkl[dr][2][0]+RelInfHivt[s]*Vmix[0]*VGjkl[dr][2][1])*Vmix[2]) /
                                (RelInfRg[0]*(VNkl[0][0]+RelInfHivt[s]*Vmix[0]*VGjkl[dr][0][1])+
                                 RelInfRg[1]*(VNkl[1][0]+Vmix[0]*VNkl[1][1])*Vmix[1]+
                                 RelInfRg[2]*(VNkl[2][0]+RelInfHivt[s]*Vmix[0]*VNkl[2][1])*Vmix[2]+1e-12);
                                Vjaf[dr][1]  = (Vmix[0]*VGjkl[dr][2][1]*RelInfHivt[s]+VGjkl[dr][2][0]) / (Vmix[0]*VNkl[2][1]*RelInfHivt[s]+VNkl[2][0]+1e-12);
                                Vjaf[dr][2]  = (Vmix[0]*VGjkl[dr][1][1]+VGjkl[dr][1][0]) / (Vmix[0]*VNkl[1][1]+VNkl[1][0]+1e-12);
                                Vjaf[dr][3]  = (RelInfRg[0]*VGjkl[dr][0][1]+Vmix[1]*RelInfRg[1]*VGjkl[dr][1][1]+Vmix[2]*RelInfRg[2]*VGjkl[dr][2][1]) /
                                (RelInfRg[0]*VNkl[0][1]+Vmix[1]*RelInfRg[1]*VNkl[1][1]+Vmix[2]*RelInfRg[2]*VNkl[2][1]+1e-12);
                                Vjaf[dr][4]  = VGjkl[dr][2][1]/(VNkl[2][1]+1e-12);
                                Vjaf[dr][5]  = VGjkl[dr][1][1]/(VNkl[1][1]+1e-12);  }
                            // Step 4
                            for(int dr=0; dr<5; dr++) {
                                VLjkl[dr][0][0]  = RelInfRg[0]*Vjaf[dr][0];
                                VLjkl[dr][0][1]  = RelInfRg[0]*RelInfHivt[s]*((1-Vmix[0])*Vjaf[dr][3]+Vmix[0]*Vjaf[dr][0]);
                                VLjkl[dr][1][0]  = RelInfRg[1]*((1-Vmix[1])*Vjaf[dr][2]+Vmix[1]*Vjaf[dr][0]);
                                VLjkl[dr][1][1]  = RelInfRg[1]*((1-Vmix[0])*((1-Vmix[1])*Vjaf[dr][5]+Vmix[1]*Vjaf[dr][3])+
                                                                Vmix[0] *((1-Vmix[1])*Vjaf[dr][2]+Vmix[1]*Vjaf[dr][0]));
                                VLjkl[dr][2][0]  = RelInfRg[2]*((1-Vmix[2])*Vjaf[dr][1]+Vmix[2]*Vjaf[dr][0])+ExogInfN[s][dr];
                                VLjkl[dr][2][1]  = RelInfRg[2]*RelInfHivt[s]*((1-Vmix[0])*((1-Vmix[2])*Vjaf[dr][4]+Vmix[2]*Vjaf[dr][3])+
                                                                              Vmix[0] *((1-Vmix[2])*Vjaf[dr][1]+Vmix[2]*Vjaf[dr][0]))+ExogInfN[s][dr];  }
                            
///////////////////////////////    INFECTION   /////////////////////////////////
for(int im=0; im<5 ; im++) {
                                for(int rg=0; rg<4; rg++) {
                                    if(im>0) { h2 = 1;
                                    } else { h2 = 0; }
                                    if(rg<3) { r2 = rg;
                                    } else { r2 = 2; }
                                    for(int ag=0; ag<11; ag++) {
                                        for(int dr=0; dr<5; dr++) {
                                            for(int tx=0; tx<3 ; tx++) {
                                                for(int d2=0; d2<5; d2++) {
                                                    ///////////////////////////////   SUCEPTIBLE  /////////////////////////////////
                                                    temp = V0[ag][0][dr][tx][im][rg]*VLjkl[d2][r2][h2]*NixTrans[s];  // Su
                                                    V1[ag][0][dr][tx][im][rg]  -= temp;
                                                    V1[ag][2][d2][tx][im][rg]  += temp*MpslowN[ag][im];
                                                    V1[ag][3][d2][tx][im][rg]  += temp*MpfastN[ag][im];
                                                    ///////////////////////////////   SUCEPTIBLE, PI  /////////////////////////////////
                                                    temp = V0[ag][1][dr][tx][im][rg]*VLjkl[d2][r2][h2]*NixTrans[s];  // Sp
                                                    V1[ag][1][dr][tx][im][rg]  -= temp;
                                                    V1[ag][2][d2][tx][im][rg]  += temp*MpslowPIN[ag][im];
                                                    V1[ag][3][d2][tx][im][rg]  += temp*MpfastPIN[ag][im];
                                                    /////////////////SUPER INFECTION LATENT TB SLOW ///////////////////////////////
                                                    temp = V0[ag][2][dr][tx][im][rg]*VLjkl[d2][r2][h2]*NixTrans[s];  // Ls
                                                    V1[ag][2][dr][tx][im][rg]  -= temp;
                                                    V1[ag][2][dr][tx][im][rg]  += temp*MpslowPIN[ag][im]/2;
                                                    V1[ag][2][d2][tx][im][rg]  += temp*MpslowPIN[ag][im]/2;
                                                    V1[ag][3][d2][tx][im][rg]  += temp*MpfastPIN[ag][im];
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
                                                    temp = V0[ag][6 ][dr][tx][im][rg]*MrslowN[ag][im]*rrSlowFB[rg]*(1-EffLtXN[s][dr]);
                                                    V1[ag][6][dr][tx][im][rg]  -= temp;
                                                    V1[ag][4][dr][tx][im][rg]  += temp*(1-MpSmPosN[ag][im]);
                                                    V1[ag][5][dr][tx][im][rg]  += temp*MpSmPosN[ag][im];
                                                } } } } }
                                
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
                                                        temp  = V0[ag][4 ][lt][im][nm][rg][na]*VrSlfCur[im];
                                                        V1[ag][4 ][lt][im][nm][rg][na]  -= temp;
                                                        V1[ag][2 ][lt][im][nm][rg][na]  += temp;
                                                    } } } } } }
                                
                                /// LTBI SCREENING AND TLTBI INITIATION /// only for no previous TB or LTBI tx
                                for(int im=0; im<5; im++) {
                                    for(int rg=0; rg<2; rg++) {
                                        for(int na=0; na<3; na++) {
                                            //////////////  NO INCREASE TB REACT RISK; US BORN, LOW RISK  //////////////////
                                            if(im==0 & rg==0 & na==0) {
                                                rTbP = rLtScrt[s]*LtDxParN[0][0]; ///check these dimensions LtDxParN
                                                rTbN = rLtScrt[s]*LtDxParN[0][1];
                                            }
                                            //////////////  NO INCREASE TB REACT RISK; US BORN, HIGH RISK  /////////////////
                                            if(im==0 & rg==1 & na==0 ) {
                                                rTbP = rLtScrt[s]*LtDxParN[1][0];
                                                rTbN = rLtScrt[s]*LtDxParN[1][1];
                                            }
                                            ////////////  NO INCREASE TB REACT RISK; NON US BORN, LOW RISK  ////////////////
                                            if(im==0 & rg==0 & na>0) {
                                                rTbP = rLtScrt[s]*LtDxParN[2][0];
                                                rTbN = rLtScrt[s]*LtDxParN[2][1];
                                            }
                                            ////////////  NO INCREASE TB REACT RISK; NON US BORN, HIGH RISK  ///////////////
                                            if(im==0 & rg==1 & na>0 ) {
                                                rTbP = rLtScrt[s]*LtDxParN[2][0];
                                                rTbN = rLtScrt[s]*LtDxParN[2][1];
                                            }
                                            //////////////   INCREASE TB REACT RISK; US BORN, LOW RISK   //////////////////
                                            if(im>0  & rg==0 & na==0 ) {
                                                rTbP = rLtScrt[s]*LtDxParN[3][0];
                                                rTbN = rLtScrt[s]*LtDxParN[3][1];
                                            }
                                            //////////////   INCREASE TB REACT RISK; US BORN, HIGH RISK  //////////////////
                                            if(im>0 & rg==1 & na==0 ) {
                                                rTbP = rLtScrt[s]*LtDxParN[4][0];
                                                rTbN = rLtScrt[s]*LtDxParN[4][1];
                                            }
                                            ////////////   INCREASE TB REACT RISK; NON US BORN, LOW RISK    ///////////////
                                            if(im>0 & rg==0 & na>0 ) {
                                                rTbP = rLtScrt[s]*LtDxParN[5][0];
                                                rTbN = rLtScrt[s]*LtDxParN[5][1];
                                            }
////////////   INCREASE TB REACT RISK; NON US BORN, HIGH RISK   ///////////////
if(im>0 & rg==1 & na>0 ) {
    rTbP = rLtScrt[s]*LtDxParN[5][0];
    rTbN = rLtScrt[s]*LtDxParN[5][1];
}
for(int ag=0; ag<11; ag++) {
    // Have LTBI
    temp  = V0[ag][2][0][im][nm][rg][na]*rTbP;
    temp2 = V0[ag][3][0][im][nm][rg][na]*rTbP;
    V1[ag][2][0][im][nm][rg][na]  -= temp;
    V1[ag][3][0][im][nm][rg][na]  -= temp2;
    V1[ag][6][0][im][nm][rg][na]  += temp+temp2;
// Dont have LTBI
                                                temp  = V0[ag][0][0][im][nm][rg][na]*rTbN;
                                                temp2 = V0[ag][1][0][im][nm][rg][na]*rTbN;
                                                V1[ag][0][0][im][nm][rg][na]  -= temp;
                                                V1[ag][1][0][im][nm][rg][na]  -= temp2;
                                                V1[ag][0][0][im][nm][rg][na]  += temp;
                                                V1[ag][1][0][im][nm][rg][na]  += temp2;
                                            } } } }
                                
                                /// TLTBI: TX COMPLETION + DEFAULT /// only need to consider tx naive compartment
                                for(int ag=0; ag<11; ag++) {
                                    for(int im=0; im<4 ; im++) {
                                        for(int nm=0; nm<4; nm++) {
                                            for(int rg=0; rg<2; rg++) {
                                                for(int na=0; na<3; na++) {
                                                    temp  = V0[ag][6][0][im][nm][rg][na]*dLtt[s]; // tx completion
                                                    temp2 = V0[ag][6][0][im][nm][rg][na]*LtTxPar[1]; // default
                                                    V1[ag][6 ][0][im][nm][rg][na]  -= temp+temp2;
                                                    V1[ag][1 ][1][im][nm][rg][na]  += temp*EffLtXN[s]*EffLt;
                                                    V1[ag][2 ][1][im][nm][rg][na]  += temp*(1-EffLtXN[s]*EffLt);
                                                    V1[ag][2 ][0][im][nm][rg][na]  += temp2;
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
                        temp  = V0[ag][4][dr][tx][im][rg]*rDxtN[s][ti  ]/RRdxAge[ag];
                        V1[ag][4][lt][im][nm][rg][na]      -= temp;
                        V1[ag][5][lt][im][nm][rg][na]      += temp *(pTbTx);
                        Vdx[ag][4][lt][im][nm][rg][na]   = temp;
} } } } } }
/////////////////////////// TB REACT PROGRESSION  ////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++) {
            for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                        temp  = V0[ag][tb][lt][1][nm][rg][na]*vImxtoImyN[ag][0];
                        temp2 = V0[ag][tb][lt][2][nm][rg][na]*vImxtoImyN[ag][1];
                        temp3 = V0[ag][tb][lt][3][nm][rg][na]*vImxtoImyN[ag][2];
                                                        V1[ag][tb][lt][1][nm][rg][na]  -= temp;
                                                        V1[ag][tb][lt][2][nm][rg][na] += temp;
                                                        V1[ag][tb][lt][2][nm][rg][na] -= temp2;
                                                        V1[ag][tb][lt][3][nm][rg][na]  += temp2;
                                                        V1[ag][tb][lt][3][nm][rg][na] -= temp3;
                                                        V1[ag][tb][lt][4][nm][rg][na]  += temp3;
                                                    } } } } } }
/////////////////////////// NON TB MORTALITY PROGRESSION  ////////////////////
                                for(int ag=0; ag<11; ag++) {
                                    for(int tb=0; tb<6; tb++) {
                                        for(int lt=0; lt<2; lt++) {
                                            for(int im=0; im<4; im++) {
                                                for(int rg=0; rg<2; rg++) {
                                                    for(int na=0; na<3; na++) {
                                                        temp  = V0[ag][tb][lt][im][1][rg][na]*vNmxtoNmyN[ag][0];
                                                        temp2 = V0[ag][tb][lt][im][2][rg][na]*vNmxtoNmyN[ag][1];
                                                        temp3 = V0[ag][tb][lt][im][3][rg][na]*vNmxtoNmyN[ag][2];
                                                        V1[ag][tb][lt][im][1][rg][na]  -= temp;
                                                        V1[ag][tb][lt][im][2][rg][na] += temp;
                                                        V1[ag][tb][lt][im][2][rg][na] -= temp2;
                                                        V1[ag][tb][lt][im][3][rg][na]  += temp2;
                                                        V1[ag][tb][lt][im][3][rg][na] -= temp3;
                                                        V1[ag][tb][lt][im][4][rg][na]  += temp3;
                                                    } } } } } }
//////////////////// RISK FACTOR OF INTEREST INCIDENCE //////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++) {
            for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                    if(rg!=1) {
                        temp = V0[ag][tb][lt][0][0][rg][na]*rRFtN[s][ag][0];
                    } else {
                        temp = V0[ag][tb][lt][0][0][rg][na]*rRFtN[s][ag][1];
                    }
                    V1[ag][tb][lt][0][0][rg][na]  -= temp;
                    V1[ag][tb][lt][1][1][rg][na]  += temp;
} } } } }
//////////////////////// INTERVENTION ENROLLMENT /////////////////////////////
// rate differs by RISK category, rg group
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++) {
            for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                    for(int rf2=0; rf2<4; rf2++) {
                        if( tb>5 | rf2==1 ) {
                            ti = 1;
                        } else { ti = 0;
                        }   // ti = 0 for slow rate, 1 for higher rate (advanced HIV and/or diag TB)
                        if( rg==1 ) {
                            ti += 2; } // 3rd and 4th col of rIntvInit for HR group (will need to update)
                        temp = V0[ag][tb][lt][1+rf2][1+rf2][rg][na]*rIntvInitN[s][ti];
                        V1[ag][tb][lt][1+rf2][1+rf2][rg][na]  -= temp;
                        V1[ag][tb][lt][2+rf2][2+rf2][rg][na]  += temp;
} } } } } }
////////////////////////  INTERVENTION DEFAULT  /////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++) {
            for(int rg=0; rg<2; rg++) {
                for(int na=0; na<3; na++) {
                    if(rg!=1) { temp3=rIntvDef; } else { temp3=rIntvDef*2; }
                    temp  = V0[ag][tb][dr][tx][2 ][rg]*temp3;
                    temp2 = V0[ag][tb][dr][tx][4 ][rg]*temp3;
                    V1[ag][tb][dr][tx][2 ][rg]  -= temp;
                    V1[ag][tb][dr][tx][1 ][rg]  += temp;
                    V1[ag][tb][dr][tx][4 ][rg]  -= temp2;
                    V1[ag][tb][dr][tx][3 ][rg]  += temp2;
} } } }
                                    
//////////////////////// TB TREATMENT OUTCOMES /////////////////////////////
for(int rg=0; rg<2; rg++) {
    if(rg!=1) { ti2 = 0; }
    else { ti2 = 1; }
    for(int ag=0; ag<11; ag++) {
        for(int lt=0; lt<2; lt++) {
            for(int nm=0; nm<4; nm++) {
                for(int na=0; na<3; na++) {
// Cures back to Ls state
                    for(int i=0; i<2; i++) {
                        ti3 = dr*2+i;
                        ti4 = 7+i*2;
                        ti5 = 8+i*2;
                        if(i==0) { ti = 7; }
                        else { ti = 9; }
                        temp  = V0[ag][ti  ][dr][tx][im][rg]*TxMatZ[7][ti3][ti2];
                        temp2 = V0[ag][ti+1][dr][tx][im][rg]*TxMatZ[7][ti3][ti2];
                        V1[ag][ti4][dr][tx][im][rg]  -= temp;
                        V1[ag][ti5][dr][tx][im][rg]  -= temp2;
                        V1[ag][2  ][dr][2 ][im][rg]  += temp + temp2;
                                                            // Failures stay in naive (if naive) as will still be one episode of tb
                                                            for(int d2=0; d2<5; d2++) {  if(idAR[d2][ti3]==1) {
                                                                // To active disease state, with AR distribution
                                                                temp  = V0[ag][ti  ][dr][tx][im][rg]*TxMatZ[10+d2][ti3][ti2];  // sm neg
                                                                temp2 = V0[ag][ti+1][dr][tx][im][rg]*TxMatZ[10+d2][ti3][ti2];  // sm pos
                                                                V1[ag][ti4][dr][tx][im][rg]  -= temp;
                                                                V1[ag][4  ][d2][tx][im][rg]  += temp;
                                                                V1[ag][ti5][dr][tx][im][rg]  -= temp2;
                                                                V1[ag][5  ][d2][tx][im][rg]  += temp2;
                                                                // To retreatment, with AR distribution
                                                                temp  = V0[ag][ti  ][dr][tx][im][rg]*TxMatZ[15+d2][ti3][ti2];   // sm neg
                                                                temp2 = V0[ag][ti+1][dr][tx][im][rg]*TxMatZ[15+d2][ti3][ti2];   // sm pos
                                                                V1[ag][ti4][dr][tx][im][rg]  -= temp;
                                                                V1[ag][9  ][d2][tx][im][rg]  += temp;
                                                                V1[ag][ti5][dr][tx][im][rg]  -= temp2;
                                                                V1[ag][10 ][d2][tx][im][rg]  += temp2;    } } }
                                                    } } } } }
///////////////////////////////////////////////////////////////////////////////
/////////////////////////    FILL RESULTS TABLE    ////////////////////////////
///////////////////////////////////////////////////////////////////////////////
/////////////////////////     MID-YEAR RESULTS     ////////////////////////////
if(m==6) {
/////////////////////////////////// YEAR //////////////////////////////////////
Outputs[y][0]      = y+1950;  // Year
////////////////    COUNTS BY TOTAL, AGE, TB, RF, AND RG     //////////////////
for(int ag=0; ag<11; ag++) {
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
} } } } } }
////////////////////    COUNTS BY NATIVITY AND AGE    ////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++) {
            for(int im=0; im<4; im++) {
                for(int nm=0; nm<4; nm++) {
                    for(int rg=0; rg<2; rg++) {
                        Outputs[y][33+ag] += V1[ag][tb][lt][im][nm][rg][0] // N_ by age and US (11)
                        Outputs[y][44+ag] += V1[ag][tb][lt][im][nm][rg][1]+V1[ag][tb][lt][im][nm][rg][2];   // N_ by age and FB (11)
} } } } }
/////////////////    COUNTS BY NATIVITY, AGE, & LTBI    //////////////////////
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2; lt++) {
        for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                    Outputs[y][55+ag] += (V1[ag][1][lt][im][nm][rg][0])*(1-pImmScen)+
                        V1[ag][2][lt][im][nm][rg][0]+V1[ag][3][lt][im][nm][rg][0];   // N_ by age and US (11) LATENT INFECTION
                        V1[ag][1][lt][im][nm][rg][1]+V1[1][tb][lt][im][nm][rg][2]
                    Outputs[y][66+ag] += (V1[ag][1][lt][im][nm][rg][1]+V1[1][tb][lt][im][nm][rg][2])*(1-pImmScen)+
                        V1[ag][2][lt][im][nm][rg][1]+V1[2][tb][lt][im][nm][rg][2]
                        V1[ag][3][lt][im][nm][rg][1]+V1[3][tb][lt][im][nm][rg][2]; // N_ by age and FB (11) LATENT INFECTION
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
} } } } } }
///////////// TB MORTALITY COUNT BY AGE, RISK FACTOR OF INTEREST///////////////
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2; lt++) {
        for(int im=0; im<4; im++) {
            for(int nm=0; nm<4; nm++) {
                for(int rg=0; rg<2; rg++) {
                    for(int na=0; na<3; na++) {
                        if(im>0) { ti = 11; } else { ti = 0; }
                        if(im==3) { temp = muTbH; } else { temp = 0; }
                        Outputs[y][88+ag+ti]  += V0[ag][4 ][lt][im][nm][rg][na]*(vTMortN[ag][4 ]+temp);
                        Outputs[y][88+ag+ti]  += V0[ag][5 ][lt][im][nm][rg][na]*(vTMortN[ag][5 ]+temp)*pow(1.0-TxMatZ[1][dr*2+0][0],TunTxMort); ///check txmatz
} } } } }
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
                            Outputs[y][110+ag]  += V0[ag][tb][lt][im][nm][rg][na]*vHMortN[ag][im];
                                                        } } } } }
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
                                                                                if(rg!=1) { temp = rDeft[s]; } else { temp = rDeftH[s]; }
                                                                                Outputs[y][132]  += V0[ag][5 ][lt][im][nm][rg][na]*TxMatZ[22][dr*2+0][0] // tx completion  check this!!!
                                                                                Outputs[y][133]  += V0[ag][5 ][lt][im][nm][rg][na] // tx discontinuation
                                                                                Outputs[y][134]  += VMort[ag][5 ][lt][im][nm][rg][na]// tx mort
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
                                        if(tx<2) {
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
                                                        for(int im=0; im<5 ; im++) { for(int rg=0; rg<2; rg++) {
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
                                                                
                                                                    Outputs[y][179] += (V0[ag][3 ][0 ][im][nm][rg][na]+V0[ag][2 ][0 ][im][nm][rg][na])*rTbP // inits with LTBI
                                                            } } }
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
                        temp = V0[ag][2 ][lt][im][nm][rg][na]*MrslowN[ag][im] + V0[ag][3 ][lt][im][nm][rg][na]**rfast; // All progression
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
                
                    } } } } }
    
for(int i=180; i<214; i++) { Outputs[y][i] = Outputs[y][i]*12; } // annualize
                                                        
// NOTIFICATIONS, dead at diagnosis
for(int nm=0; nm<4 ; im++) {
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
                        temp2 = V0[ag][4 ][lt][im][nm][rg][na]*(vTMortN[ag][tb]+temp)
                        Outputs[y][214   ] += temp2;   // All dx (1)
                        Outputs[y][215+ag] += temp2;   // dx by age (11)
                        if(tx<1) {
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
    for(int ag=0; ag<11; ag++)
        for(int lt=0; lt<2; lt++){
            for (int im=0; im<4; im++){
                 for (int nm=0; nm<4; nm++){
                     for(int rg=0; rg<2; rg++) {
                         for(int na=0; na<3; na++){
                            Outputs[y][232+ag] += Vdx[ag][4 ][lt][im][nm][rg][na];   // dx by age (11)
} } } } }
   for(int i=232; i<243; i++) { Outputs[y][i] = Outputs[y][i]*12; }
// NOTIFICATIONS US, dead at diagnosis
for(int nm=0; nm<4 ; nm++) {
if(im > 2) { temp = muTbRF;  } else { temp = 0;  }
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2; lt++){
        for (int im=0; im<4; im++){
        for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++){
                temp2 = V0[ag][4 ][lt][im][nm][rg][na]*(vTMortN[ag][tb]+temp);
                Outputs[y][243+ag] += temp2;   // dx by age (11)
} } } } } }
for(int i=243; i<254; i++) { Outputs[y][i] = Outputs[y][i]*12; }
                                                        
// NOTIFICATIONS, dead at diagnosis  HIV_NEGATIVE
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2; lt++){
        for(int rg=0; rg<2; rg++) {
            for(int na=0; na<3; na++){
                temp2 = V0[ag][4 ][lt][0][0][rg][na]*(vTMortN[ag][tb]+temp);
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
} } } } } }
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
        for(int ag=0; ag<11; ag++)
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
                                        Outputs[y][296] += V1[ag][5][lt][im][nm][rg][na]; }
                                } } } } }
                // Force of infection
                for(int dr=0; dr<5; dr++) {
                                                            Outputs[y][297] += VLjkl[dr][0][0];
                                                            Outputs[y][298] += VLjkl[dr][1][0];
                                                            Outputs[y][299] += VLjkl[dr][2][0];
                                                            Outputs[y][300] += VLjkl[dr][0][1];
                                                            Outputs[y][301] += VLjkl[dr][1][1];
                                                            Outputs[y][302] += VLjkl[dr][2][1];      }
                                                        ////////////     CREATE YEARLY VALUES FROM THE MONTH ESTIMATE     ////////////
                                                        for(int i=297; i<303; i++) { Outputs[y][i] = Outputs[y][i]*12; }
                                                        
///  NEW INFECTIONS + SUPER INFECTIOPN ///
for(int im=0; im<4 ; im++) {
    for(int rg=0; rg<2; rg++) {
        for(int na=0; na<3; na++){
            if(im>0) {
                h2 = 1;
            } else {
                h2 = 0;
            }
if(rg<3) { r2 = rg; } else { r2 = 2; }
for(int ag=0; ag<11; ag++) {
    for(int dr=0; dr<5; dr++) {
        for(int tx=0; tx<3 ; tx++) {
            for(int d2=0; d2<5; d2++) {
                Outputs[y][303+rg] += (V0[ag][0][lt][im][nm][rg][na]+V0[ag][1][lt][im][nm][rg][na]+V0[ag][2][lt][im][nm][rg][na])*VLjkl[d2][r2][h2]*NixTrans[s];
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
                        Outputs[y][307+ti]  += V0[ag][5][lt][im][nm][rg][na]*(vTMortN[ag][5 ]+temp)*pow(1.0-TxMatZ[1][dr*2+0][0],TunTxMort); //check the TxMatZ portion
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
NumericVector  CheckV0(14784);
for(int ag=0; ag<11; ag++) {
    for(int tb=0; tb<6; tb++) {
        for(int lt=0; lt<2; lt++){
            for(int im=0; im<4; im++){
                for(int nm=0; nm<4; nm++){
                    for(int rg=0; rg<2; rg++) {
                        for(int na=0; na<2; na++){
                            CheckV0(ag+tb*11+lt*66+tx*132+im*528+nm*2112+rg*4928+na*12672) = V1[ag][tb][lt][im][nm][rg][na];
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
);      }
                                    
                                    

