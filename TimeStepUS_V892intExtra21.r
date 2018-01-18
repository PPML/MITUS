library(Rcpp)
############

StatNam2 <- c("0-4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p")
StatNam2 <- paste(rep(c("Su","Sp","Ls","Lf","In","Ip","Tl","Fn","Fp","Mn","Mp"),each=length(StatNam2)),rep(StatNam2,11),sep="_")
StatNam2 <- paste(rep(1:5,each=length(StatNam2)),rep(StatNam2,5),sep="_")
StatNam2 <- paste(rep(c("NT","TL","TX"),each=length(StatNam2)),rep(StatNam2,2),sep="_")
StatNam2 <- paste(rep(c("N0","H1","T1","H2","T2"),each=length(StatNam2)),rep(StatNam2,4),sep="_")
StatNam2 <- paste(rep(c("LR","HR","F1","F2"),each=length(StatNam2)),rep(StatNam2,4),sep="_")
length(StatNam2) # 36300

StatList <- noquote(list(c("0_4",paste(0:8*10+5,1:9*10+4,sep="_"),"95p"),
            c("Su","Sp","Ls","Lf","In","Ip","Tl","Fn","Fp","Mn","Mp"),
            1:5,c("NT","TL","TX"),c("N0","H1","T1","H2","T2"),c("LR","HR","F1","F2")))

########################################### dim(rArtInit)

cppFunction('
  List cSim(
    int                 nYrs,
    int                 nRes,
    NumericMatrix       rDxt,
    std::vector<double> TxQualt,
    NumericMatrix       InitPop,
    double              p_HR,
    NumericMatrix       Mpfast,
    NumericMatrix       ExogInf,
    NumericMatrix       MpfastPI,
    double              pimmed,
    NumericMatrix       MpSmPos,
    double              HivHrPar,
    std::vector<double> RRmuHR,
    NumericMatrix       Mrslow,
    std::vector<double> rrSlowFB,
    double              rfast,
    double              RRcurDef,
    double              fAR,
    std::vector<double> VrSlfCur,
    std::vector<double> VrSmConv,
    NumericMatrix       vTMort,  
    NumericMatrix       vHMort,
    double              muTbH,
    std::vector<double> Birthst,
    NumericMatrix       ImmNon,
    NumericMatrix       ImmLat,
    NumericMatrix       ImmAct,
    NumericMatrix       ImmFst,
    NumericMatrix       DrN,
    NumericMatrix       DrE,
    std::vector<double> TxExpAge,
    double              p_Imm_SP,
    NumericMatrix       mubt,
    NumericMatrix       RelInf,
    std::vector<double> RelInfRg,
    std::vector<double> Vmix,
    NumericMatrix       vHxtoHy,
    double              rArtDef,
    NumericMatrix       rHIVt,
    std::vector<double> rEmmigFB,
    NumericMatrix       rArtInit,
    std::vector<double> pDstt,
    NumericMatrix       TxMat,
    double              TunTxMort,
    std::vector<double> rDeft,
    std::vector<double> rDeftH,
    std::vector<double> pReTx,
    std::vector<double> LtTxPar,
    NumericMatrix       LtDxPar,
    std::vector<double> rLtScrt,
    NumericMatrix       HrEntEx,
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
 
  /////////////////////////////////////////////////////////////////////////////
  //// CREATED INTERNALLY 
    int           idAR[5][10];
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
    double        MpimmedNp[Mpfast.nrow()][Mpfast.ncol()];
    double        MpimmedNn[Mpfast.nrow()][Mpfast.ncol()];
    double        MpslowN[Mpfast.nrow()][Mpfast.ncol()];
    double        MpfastPIN[MpfastPI.nrow()][MpfastPI.ncol()];
    double        MpimmedPINp[MpfastPI.nrow()][MpfastPI.ncol()];
    double        MpimmedPINn[MpfastPI.nrow()][MpfastPI.ncol()];
    double        MpslowPIN[MpfastPI.nrow()][MpfastPI.ncol()];
    double        MrslowN[Mrslow.nrow()][Mrslow.ncol()];
    double        MpSmPosN[Mpfast.nrow()][Mpfast.ncol()];
    double        vTMortN[vTMort.nrow()][vTMort.ncol()];
    double        vHMortN[vHMort.nrow()][vHMort.ncol()];
    double        vHxtoHyN[vHxtoHy.nrow()][vHxtoHy.ncol()];
    double        ImmNonN[ImmNon.nrow()][ImmNon.ncol()];
    double        ImmLatN[ImmLat.nrow()][ImmLat.ncol()];
    double        ImmFstN[ImmFst.nrow()][ImmFst.ncol()];
    double        ImmActN[ImmAct.nrow()][ImmAct.ncol()];
    double        DrImm[11][5][2][DrN.nrow()][2];
    double        DrNN[DrN.nrow()][DrN.ncol()];
    double        DrEN[DrE.nrow()][DrE.ncol()];
    double        mubtN[mubt.nrow()][mubt.ncol()];
    double        RelInfN[RelInf.nrow()][RelInf.ncol()];
    double        rHIVtN[rHIVt.nrow()][rHIVt.ncol()][2];
    double        rArtInitN[rArtInit.nrow()][rArtInit.ncol()];
    double        rDxtN[rDxt.nrow()][rDxt.ncol()];
    double        TxMatN[TxMat.nrow()][TxMat.ncol()];
    double        HrEntExN[HrEntEx.nrow()][HrEntEx.ncol()][2];
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
    double        VGjkl[5][3][2]; 
    double        Vjaf[5][6];
    double        VLjkl[5][3][2]; 
    NumericMatrix Outputs2(nYrs,nRes);

  //// INITIALIZE 
    for(int i=0; i<5; i++) { for(int j=0; j<10; j++) { idAR[i][j] = 0; } }
    idAR[0][0] = 1;  idAR[0][1] = 1;
    idAR[1][0] = 1;  idAR[1][1] = 1;  idAR[1][2] = 1;  idAR[1][3] = 1;
    idAR[2][0] = 1;  idAR[2][1] = 1;  idAR[2][4] = 1;  idAR[2][5] = 1;
    idAR[3][0] = 1;  idAR[3][1] = 1;  idAR[3][2] = 1;  idAR[3][3] = 1; idAR[3][4] = 1;  idAR[3][5] = 1; idAR[3][6] = 1;  idAR[3][7] = 1;
    idAR[4][7] = 1;  idAR[4][8] = 1;  idAR[4][9] = 1;

    for(int i=0; i<InitPop.nrow(); i++) { for(int j=0; j<InitPop.ncol(); j++) { InitPopN[i][j] = InitPop(i,j); } }
    for(int i=0; i<ExogInf.nrow(); i++) { for(int j=0; j<ExogInf.ncol(); j++) { ExogInfN[i][j] = ExogInf(i,j); } }
    for(int i=0; i<Mpfast.nrow(); i++) { for(int j=0; j<Mpfast.ncol(); j++) { 
      MpfastN[i][j]   = Mpfast(i,j)*(1-pimmed);               MpfastPIN[i][j]   = MpfastPI(i,j)*(1-pimmed);
      MpimmedNp[i][j] = Mpfast(i,j)*pimmed*MpSmPos(i,j);      MpimmedPINp[i][j] = MpfastPI(i,j)*pimmed*MpSmPos(i,j);
      MpimmedNn[i][j] = Mpfast(i,j)*pimmed*(1-MpSmPos(i,j));  MpimmedPINn[i][j] = MpfastPI(i,j)*pimmed*(1-MpSmPos(i,j));
      MpslowN[i][j]   = 1-Mpfast(i,j);                        MpslowPIN[i][j]   = 1-MpfastPI(i,j);
      MrslowN[i][j]   = Mrslow(i,j);                          MpSmPosN[i][j]    = MpSmPos(i,j);  } }
    for(int i=0; i<LtDxPar.nrow(); i++) { for(int j=0; j<LtDxPar.ncol(); j++) { LtDxParN[i][j] = LtDxPar(i,j); } }
    for(int i=0; i<EffLtX.nrow(); i++) { for(int j=0; j<EffLtX.ncol(); j++) { EffLtXN[i][j] = EffLtX(i,j); } }
    for(int i=0; i<vTMort.nrow(); i++) { for(int j=0; j<vTMort.ncol(); j++) { vTMortN[i][j] = vTMort(i,j); } }
    for(int i=0; i<vHMort.nrow(); i++) { for(int j=0; j<vHMort.ncol(); j++) { vHMortN[i][j] = vHMort(i,j); } }
    for(int i=0; i<vHxtoHy.nrow(); i++) { for(int j=0; j<vHxtoHy.ncol(); j++) { vHxtoHyN[i][j] = vHxtoHy(i,j); } }
    for(int i=0; i<ImmNon.nrow(); i++) { for(int j=0; j<ImmNon.ncol(); j++) { ImmNonN[i][j] = ImmNon(i,j); } }
    for(int i=0; i<ImmLat.nrow(); i++) { for(int j=0; j<ImmLat.ncol(); j++) { ImmLatN[i][j] = ImmLat(i,j); } }
    for(int i=0; i<ImmFst.nrow(); i++) { for(int j=0; j<ImmFst.ncol(); j++) { ImmFstN[i][j] = ImmFst(i,j); } }
    for(int i=0; i<ImmAct.nrow(); i++) { for(int j=0; j<ImmAct.ncol(); j++) { ImmActN[i][j] = ImmAct(i,j); } }
    for(int i=0; i<DrN.nrow(); i++) { for(int j=0; j<DrN.ncol(); j++) { DrNN[i][j] = DrN(i,j); } }
    for(int i=0; i<DrE.nrow(); i++) { for(int j=0; j<DrE.ncol(); j++) { DrEN[i][j] = DrE(i,j); } }
    for(int ag=0; ag<11; ag++) {  for(int dr=0; dr<5; dr++) { for(int s=0; s<DrN.nrow(); s++) {
        		DrImm[ag][dr][0][s][0] = (1-TxExpAge[ag])*DrNN[s][dr]*ImmActN[s][ag];
      			DrImm[ag][dr][1][s][0] = TxExpAge[ag]    *DrEN[s][dr]*ImmActN[s][ag];
        		DrImm[ag][dr][0][s][1] = (1-TxExpAge[ag])*DrNN[s][dr]*ImmFstN[s][ag];
      			DrImm[ag][dr][1][s][1] = TxExpAge[ag]    *DrEN[s][dr]*ImmFstN[s][ag]; }	}	}
    for(int i=0; i<mubt.nrow(); i++) { for(int j=0; j<mubt.ncol(); j++) { mubtN[i][j] = mubt(i,j); } }
    for(int i=0; i<RelInf.nrow(); i++) { for(int j=0; j<RelInf.ncol(); j++) { RelInfN[i][j] = RelInf(i,j); } }
    for(int i=0; i<rHIVt.nrow(); i++) { for(int j=0; j<rHIVt.ncol(); j++) { rHIVtN[i][j][0] = rHIVt(i,j); rHIVtN[i][j][1] = rHIVt(i,j)*HivHrPar; } }
    for(int i=0; i<rArtInit.nrow(); i++) { for(int j=0; j<rArtInit.ncol(); j++) { rArtInitN[i][j] = rArtInit(i,j); } }
    for(int i=0; i<TxMat.nrow(); i++) { for(int j=0; j<TxMat.ncol(); j++) { TxMatN[i][j] = TxMat(i,j); } }
    for(int i=0; i<rDxt.nrow(); i++) { for(int j=0; j<rDxt.ncol(); j++) { rDxtN[i][j] = rDxt(i,j); } }
    for(int i=0; i<HrEntEx.nrow(); i++) { for(int j=0; j<HrEntEx.ncol(); j++) { HrEntExN[i][j][0] = HrEntEx(i,j); HrEntExN[i][j][1] = HrEntEx(i,j)*HivHrPar; } }
    for(int i=0; i<22; i++) { for(int j=0; j<10; j++) { for(int k=0; k<2; k++) { TxMatZ[i][j][k] = 0.0; } } }
    for(int i=0; i<nYrs; i++) {  for(int j=0; j<nRes; j++) {  Outputs[i][j] = 0;  }  }
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      			V0[ag][tb][dr][tx][hv][rg]    = 0; 	V1[ag][tb][dr][tx][hv][rg]  = 0; 
    				VMort[ag][tb][dr][tx][hv][rg] = 0;  Vdx[ag][tb][dr][tx][hv][rg] = 0;   }	}	}	} } }
    for(int i=0; i<3; i++) { for(int j=0; j<2; j++) { VNkl[i][j] = 0; } }
    for(int i=0; i<5; i++) { for(int j=0; j<3; j++) { for(int k=0; k<2; k++) { VGjkl[i][j][k] = 0; VLjkl[i][j][k] = 0; } } }
    for(int i=0; i<5; i++) { for(int j=0; j<6; j++) { Vjaf[i][j] = 0; } }
    for(int i=0; i<5; i++) { extrV[i*2] = extrV[i*2+1] = i; }
    for(int ag=0; ag<11; ag++) { for(int hv=0; hv<5; hv++) {  temp4V[ag][hv] = (1-pow(1-MrslowN[ag][hv]-rRecov,24.0))*MrslowN[ag][hv]; } }

  /////////////////////////////////////////////////////////////////////////////
 
  //// UPDATING TREATMENT PARAMETERS   ****  THIS DIFFERENT TO MAIN MODEL DUE TO SIMPLIFIED OUTCOMES ***
    for(int j=0; j<10; j++) { for(int k=0; k<2; k++) {
      if(k==0) { temp = rDeft[0]; } else { temp = rDeftH[0]; }
      TxMatZ[1][j][k] = TxMatN[1][j]*TxQualt[0];            // TxEff updated for txqual
      TxMatZ[7][j][k] = TxMatN[0][j]*TxMatZ[1][j][k] + temp*TxMatZ[1][j][k]*RRcurDef; // Rate of treatment exit to Ls (cure)
      TxMatZ[8][j][k] = TxMatN[0][j]*(1-TxMatZ[1][j][k]) + temp*(1-TxMatZ[1][j][k]*RRcurDef); // Rate of treatment exit to In/Ip (failure)  ** THIS CHANGED COMPARED TO MAIN MODEL ***
      } }

  /////////////////////////////////////////////////////////////////////////////// StatList
  //// BURN IN  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// Populate model
    for(int ag=0; ag<11; ag++) {
      V0[ag][0][0][0][0][0] = InitPopN[ag][0]*0.40*(1-p_HR); V0[ag][0][0][0][0][1] = InitPopN[ag][0]*0.40*p_HR; V0[ag][0][0][0][0][3] = InitPopN[ag][1]*0.40;
      V0[ag][2][0][0][0][0] = InitPopN[ag][0]*0.60*(1-p_HR); V0[ag][2][0][0][0][1] = InitPopN[ag][0]*0.60*p_HR; V0[ag][2][0][0][0][3] = InitPopN[ag][1]*0.60;  }

    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) {  for(int rg=0; rg<4; rg++) {
        V1[ag][tb][0][0][0][rg]  = V0[ag][tb][0][0][0][rg];  } } }

  for(int m=0; m<3001; m++) {  // ~~~ start burn in ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  /// BIRTHS /// 
     V1[0][0][0][0][0][0]  += Birthst[0]*(1-p_HR); 
     V1[0][0][0][0][0][1]  += Birthst[0]*p_HR;

  /// IMMIGRATION /// Assume 0% HIV prev, Assume 0% prior tx, Assume All Pansens
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { 
      V1[ag][0][0][0][0][2]  += ImmNonN[0][ag];      // NO TB
      V1[ag][2][0][0][0][2]  += ImmLatN[0][ag];      // LATENT TB
      V1[ag][3][0][0][0][2]  += (DrImm[ag][0][0][0][1]+DrImm[ag][0][1][0][1]);  // LATENT FAST
      V1[ag][4][0][0][0][2]  += (DrImm[ag][0][0][0][0]+DrImm[ag][0][1][0][0])*(1-p_Imm_SP);
      V1[ag][5][0][0][0][2]  += (DrImm[ag][0][0][0][0]+DrImm[ag][0][1][0][0])*p_Imm_SP;  } }
 
  /// EMIGRATION /// 
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { 
      V1[ag][tb][0][0][0][2]  -= V0[ag][tb][0][0][0][2]*rEmmigFB[0];        // FB1
      V1[ag][tb][0][0][0][3]  -= V0[ag][tb][0][0][0][3]*rEmmigFB[1];  } }   // FB2

  /// MORTALITY ///
    for(int ag=0; ag<11; ag++) { for(int rg=0; rg<4; rg++) {
      for(int tb=0; tb<7; tb++) {
      V1[ag][tb][0][0][0][rg]  -= V0[ag][tb][0][0][0][rg]*(mubtN[0][ag]*RRmuHR[rg]+vTMortN[ag][tb]); }
      V1[ag][7 ][0][0][0][rg]  -= V0[ag][7 ][0][0][0][rg]*(mubtN[0][ag]*RRmuHR[rg]+vTMortN[ag][7 ]*pow(1.0-TxMatZ[1][0][0],TunTxMort));
      V1[ag][8 ][0][0][0][rg]  -= V0[ag][8 ][0][0][0][rg]*(mubtN[0][ag]*RRmuHR[rg]+vTMortN[ag][8 ]*pow(1.0-TxMatZ[1][0][0],TunTxMort));
      V1[ag][9 ][0][0][0][rg]  -= V0[ag][9 ][0][0][0][rg]*(mubtN[0][ag]*RRmuHR[rg]+vTMortN[ag][9 ]*pow(1.0-TxMatZ[1][1][0],TunTxMort));
      V1[ag][10][0][0][0][rg]  -= V0[ag][10][0][0][0][rg]*(mubtN[0][ag]*RRmuHR[rg]+vTMortN[ag][10]*pow(1.0-TxMatZ[1][1][0],TunTxMort));   } }

  /// AGING ///
    for(int ag=0; ag<10; ag++) { for(int tb=0; tb<11; tb++) { for(int rg=0; rg<4; rg++) { 
      if(ag>0) { temp2 = 120; } else { temp2 = 60; }
      temp = V0[ag][tb][0][0][0][rg]/temp2;
      V1[ag  ][tb][0][0][0][rg]  -= temp; 
      V1[ag+1][tb][0][0][0][rg]  += temp;  } } }

  /// NEW FB -> ESTABLISHED FB ///
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { 
      temp = V0[ag][tb][0][0][0][2]/24;
      V1[ag][tb][0][0][0][2]  -= temp; 
      V1[ag][tb][0][0][0][3]  += temp;   } }

  /// HIGH-RISK ENTRY/EXIT ///
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) {  
      temp  = V0[ag][tb][0][0][0][0]*HrEntExN[ag][0][0];
      temp2 = V0[ag][tb][0][0][0][1]*HrEntExN[ag][1][0];
      V1[ag][tb][0][0][0][0]  += temp2-temp; 
      V1[ag][tb][0][0][0][1]  += temp-temp2;  } } 

  /// TRANSMISSION RISK ///  
    for(int i=0; i<3; i++) { VNkl[i][0] = 0; VGjkl[0][i][0] = 0;  } // set to zero
    // Step 1
      for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { 
      VNkl[0][0]  += V0[ag][tb][0][0][0][0];
      VNkl[1][0]  += V0[ag][tb][0][0][0][1];
      VNkl[2][0]  += V0[ag][tb][0][0][0][2] + V0[ag][tb][0][0][0][3]; } }
    // Step 2  (active TB)
    for(int ag=0; ag<11; ag++) { for(int tb=4; tb<6; tb++) { 
      VGjkl[0][0][0]  += V0[ag][tb][0][0][0][0]*RelInfN[tb][0];
      VGjkl[0][1][0]  += V0[ag][tb][0][0][0][1]*RelInfN[tb][0];
      VGjkl[0][2][0]  += (V0[ag][tb][0][0][0][2] + V0[ag][tb][0][0][0][3])*RelInfN[tb][0];  } }
    // Step 2 (treated TB)
      // No contribution to force of infection 

    // Step 3
      Vjaf[0][0]  = (RelInfRg[0]*VGjkl[0][0][0]+RelInfRg[1]*VGjkl[0][1][0]*Vmix[1]+RelInfRg[2]*VGjkl[0][2][0]*Vmix[2]) /
                    (RelInfRg[0]*VNkl[0][0]    +RelInfRg[1]*VNkl[1][0]*Vmix[1]    +RelInfRg[2]*VNkl[2][0]*Vmix[2]+1e-12);
      Vjaf[0][1]  = VGjkl[0][2][0] / (VNkl[2][0]+1e-12);
      Vjaf[0][2]  = VGjkl[0][1][0] / (VNkl[1][0]+1e-12);
    // Step 4
      VLjkl[0][0][0]  = RelInfRg[0]*Vjaf[0][0];
      VLjkl[0][1][0]  = RelInfRg[1]*((1-Vmix[1])*Vjaf[0][2]+Vmix[1]*Vjaf[0][0]);
      VLjkl[0][2][0]  = RelInfRg[2]*((1-Vmix[2])*Vjaf[0][1]+Vmix[2]*Vjaf[0][0])+ExogInfN[0][0];

  ///  INFECTION /// StatList
      for(int rg=0; rg<4; rg++) {
        if(rg<3) { r2 = rg; } else { r2 = 2; }
        for(int ag=0; ag<11; ag++) { 
        temp = V0[ag][0][0][0][0][rg]*VLjkl[0][r2][0]*EarlyTrend[m];  // Su
        V1[ag][0][0][0][0][rg]  -= temp; 
        V1[ag][2][0][0][0][rg]  += temp*MpslowN[ag][0]; 
        V1[ag][3][0][0][0][rg]  += temp*MpfastN[ag][0];
        V1[ag][4][0][0][0][rg]  += temp*MpimmedNn[ag][0];
        V1[ag][5][0][0][0][rg]  += temp*MpimmedNp[ag][0];
      // SUPER-INFECTION SP
        temp = V0[ag][1 ][0][0][0][rg]*VLjkl[0][r2][0];  // Sp
        V1[ag][1][0][0][0][rg]  -= temp; 
        V1[ag][2][0][0][0][rg]  += temp*MpslowPIN[ag][0]; 
        V1[ag][3][0][0][0][rg]  += temp*MpfastPIN[ag][0]; 
        V1[ag][4][0][0][0][rg]  += temp*MpimmedPINn[ag][0]; 
        V1[ag][5][0][0][0][rg]  += temp*MpimmedPINp[ag][0];
        // SUPER-INFECTION LS
        temp = V0[ag][2 ][0][0][0][rg]*VLjkl[0][r2][0];  // Ls
        V1[ag][2][0][0][0][rg]  -= temp; 
        V1[ag][2][0][0][0][rg]  += temp*MpslowPIN[ag][0]; 
        V1[ag][3][0][0][0][rg]  += temp*MpfastPIN[ag][0]; 
        V1[ag][4][0][0][0][rg]  += temp*MpimmedPINn[ag][0]; 
        V1[ag][5][0][0][0][rg]  += temp*MpimmedPINp[ag][0]; } } 
 
  /// BREAKDOWN ///  
  for(int ag=0; ag<11; ag++) {   for(int rg=0; rg<4; rg++) { 
      temp  = V0[ag][2][0][0][0][rg]*MrslowN[ag][0]*rrSlowFB[rg];
      temp2 = V0[ag][3][0][0][0][rg]*rfast;
      V1[ag][2][0][0][0][rg]  -= temp; 
      V1[ag][3][0][0][0][rg]  -= temp2; 
      V1[ag][4][0][0][0][rg]  += (temp+temp2)*(1-MpSmPosN[ag][0]); 
      V1[ag][5][0][0][0][rg]  += (temp+temp2)*MpSmPosN[ag][0];    } }

  /// LATENT SLOW TO SAFE ///  
  for(int ag=0; ag<11; ag++) {   for(int rg=0; rg<4; rg++) { 
      temp  = V0[ag][2][0][0][0][rg]*rRecov;
      V1[ag][2][0][0][0][rg]  -= temp;
      V1[ag][1][0][0][0][rg]  += temp;   } }

  /// SMEAR CONVERSION ///
    for(int ag=0; ag<11; ag++) {  for(int rg=0; rg<4; rg++) { 
      temp = V0[ag][4][0][0][0][rg]*VrSmConv[0];
      V1[ag][4][0][0][0][rg]  -= temp; 
      V1[ag][5][0][0][0][rg]  += temp;  } }

  /// SELF-CURE /// 
    for(int ag=0; ag<11; ag++) {  for(int tb=4; tb<6; tb++) {  for(int rg=0; rg<4; rg++) { 
      temp = V0[ag][tb][0][0][0][rg]*VrSlfCur[0];
      V1[ag][tb][0][0][0][rg]  -= temp; 
      V1[ag][2 ][0][0][0][rg]  += temp;    } } }

  /// TB DIAGNOSIS AND TX INITIATION ///  
    for(int ag=0; ag<11; ag++) {  for(int rg=0; rg<4; rg++) {
      if(rg!=1) { ti = 0; } else { ti = 1;  }
      temp  = V0[ag][4 ][0][0][0][rg]*rDxtN[0][ti  ]/RRdxAge[ag]/EarlyTrend[m]; 
      temp2 = V0[ag][5 ][0][0][0][rg]*rDxtN[0][ti+2]/RRdxAge[ag]/EarlyTrend[m]; 
      V1[ag][4 ][0][0][0][rg]  -= temp;
      V1[ag][5 ][0][0][0][rg]  -= temp2;
      V1[ag][7 ][0][0][0][rg]  += temp; 
      V1[ag][8 ][0][0][0][rg]  += temp2;  } }

  /// TB TREATMENT OUTCOMES ///  
    for(int ag=0; ag<11; ag++) {  for(int rg=0; rg<4; rg++) {
      if(rg!=1) { ti2 = 0; } else { ti2 = 1; } 
    // Cures
      temp  = V0[ag][7][0][0][0][rg]*TxMatZ[7][0][ti2]; // sm neg
      temp2 = V0[ag][8][0][0][0][rg]*TxMatZ[7][0][ti2];  // sm pos
      V1[ag][7][0][0][0][rg]  -= temp;
      V1[ag][8][0][0][0][rg]  -= temp2;
      V1[ag][2][0][0][0][rg]  += temp + temp2;  
    // Failures 
      temp  = V0[ag][7][0][0][0][rg]*TxMatZ[8][0][ti2]; // sm neg
      temp2 = V0[ag][8][0][0][0][rg]*TxMatZ[8][0][ti2];  // sm pos
      V1[ag][7][0][0][0][rg]  -= temp;
      V1[ag][8][0][0][0][rg]  -= temp2;
      V1[ag][4][0][0][0][rg]  += temp;
      V1[ag][5][0][0][0][rg]  += temp2;  } }

  /// RESET POP SIZE ///
    for(int ag=0; ag<11; ag++) { for(int i=0; i<2; i++) { InitPopZ[ag][i] = 0; } }
     for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) {
      InitPopZ[ag][0]  += V1[ag][tb][0][0][0][0]+V1[ag][tb][0][0][0][1];
      InitPopZ[ag][1]  += V1[ag][tb][0][0][0][2]+V1[ag][tb][0][0][0][3];    } }
    for(int ag=0; ag<11; ag++) { for(int i=0; i<2; i++) { // factor for pop size reset
      InitPopZ[ag][i] = InitPopN[ag][i]/(InitPopZ[ag][i]+1e-12); } }
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) {   // reset pop to InitPop
      V1[ag][tb][0][0][0][0]  = V1[ag][tb][0][0][0][0]*InitPopZ[ag][0];
      V1[ag][tb][0][0][0][1]  = V1[ag][tb][0][0][0][1]*InitPopZ[ag][0];
      V1[ag][tb][0][0][0][2]  = V1[ag][tb][0][0][0][2]*InitPopZ[ag][1];
      V1[ag][tb][0][0][0][3]  = V1[ag][tb][0][0][0][3]*InitPopZ[ag][1];

      for(int rg=0; rg<4; rg++) {
        V0[ag][tb][0][0][0][rg]  = V1[ag][tb][0][0][0][rg];  } } }
    }  //  ~~~ end burn in ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ///////////////////////////////////////////////////////////////////////////////
   NumericVector        CheckV0(36300); 
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
     CheckV0(ag+tb*11+dr*121+tx*605+hv*1815+rg*9075) = V1[ag][tb][dr][tx][hv][rg];
      } } } } } } 

///////////////////////////////////////////////////////////////////////////////

  //// YEAR LOOP  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(int y=0; y<nYrs; y++) {
  //// MONTH LOOP  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  for(int m=0; m<12; m++) {  
    s = y*12+m; // counter of months since start
            
  //// UPDATING TREATMENT PARAMETERS
    // TxMatZ: 0=completion rate, 1 = tx success, 2:6 = AR probabilities, 7 = rate of exit to cure
    //         8 = adj factor for AR with tx completion, 9 = adj factor for AR with default
    //         10:14 = exit rates to active TB; 15:19 = exit rates to retx; 20 = sum of AR for active TB, 21 = sum of AR for retx

  for(int j=0; j<10; j++) { for(int k=0; k<2; k++) {
    if(k==0) { temp = rDeft[s]; } else { temp = rDeftH[s]; }
    TxMatZ[1][j][k]            = TxMatN[1][j]*TxQualt[s];            // TxEff updated for txqual
    TxMatZ[7][j][k]            = TxMatN[0][j]*TxMatZ[1][j][k] + temp*TxMatZ[1][j][k]*RRcurDef; // Rate of treatment exit to Ls (cure)
    TxMatZ[8][j][k]            = (1-exp(-(1.0-TxMatZ[1][j][k]         )*fAR)) / fAR; // Adjustment factor for AR probabilities COMPLETION = p(AR | complete) / pAR0 = p(fail)*p(AR|complete and fail) / pAR0
    TxMatZ[9][j][k]            = (1-exp(-(1.0-TxMatZ[1][j][k]*RRcurDef)*fAR)) / fAR; // Adjustment factor for AR probabilities DEFAULT
    for(int i=0; i<5; i++) {  // Rates of exit to active disease from completion, potentially with AR
      TxMatZ[10+i][j][k]       = (TxMatN[0][j]*TxMatZ[8][j][k]*(1-pReTx[s]) + temp*TxMatZ[9][j][k])*TxMatN[2+i][j];
      TxMatZ[15+i][j][k]       = TxMatN[0][j]*TxMatN[2+i][j]*TxMatZ[8][j][k] * pReTx[s]; } // rate_complete * p(AR|complete,fail,0) * Adj_factor * p(retx|complete)
    TxMatZ[20][j][k]           = TxMatZ[10][j][k]+TxMatZ[11][j][k]+TxMatZ[12][j][k]+TxMatZ[13][j][k]+TxMatZ[14][j][k]; // Sum total of AR exit rates to ACTIVE TB
    TxMatZ[21][j][k]           = TxMatZ[15][j][k]+TxMatZ[16][j][k]+TxMatZ[17][j][k]+TxMatZ[18][j][k]+TxMatZ[19][j][k]; // Sum total of AR exit rates to RETX
    TxMatZ[10+extrV[j]][j][k]  = (TxMatN[0][j]*(1.0-TxMatZ[1][j][k])*(1-pReTx[s]) + temp*(1.0-TxMatZ[1][j][k]*RRcurDef)) - TxMatZ[20][j][k]; // exit to active TB, no AR
    TxMatZ[15+extrV[j]][j][k]  = TxMatN[0][j]*(1.0-TxMatZ[1][j][k])*pReTx[s] - TxMatZ[21][j][k]; // exit to RETX, no AR
    TxMatZ[22][j][k]           = TxMatN[0][j]*(1-(1.0-TxMatZ[1][j][k])*pReTx[s]); // p(tx completion)
  } }

  /// BIRTHS ///
     V1[0][0][0][0][0][0]  += Birthst[s]*(1-p_HR);  
     V1[0][0][0][0][0][1]  += Birthst[s]*p_HR;

  /// IMMIGRATION /// 
    for(int ag=0; ag<11; ag++) { 
      V1[ag][0][0][0][0][2]  += ImmNonN[s][ag];      // NO TB
      V1[ag][2][0][0][0][2]  += ImmLatN[s][ag];      // LATENT TB
      for(int dr=0; dr<5; dr++) {       // ACTIVE TB + LAT FAST
        V1[ag][3][dr][0][0][2]  += DrImm[ag][dr][0][s][1];   // tx naive, LAT FAST
        V1[ag][3][dr][2][0][2]  += DrImm[ag][dr][1][s][1];   // tx exp, LAT FAST
        V1[ag][4][dr][0][0][2]  += DrImm[ag][dr][0][s][0]*(1-p_Imm_SP);   // tx naive, sn
        V1[ag][5][dr][0][0][2]  += DrImm[ag][dr][0][s][0]*p_Imm_SP;       // tx naive, sp
        V1[ag][4][dr][2][0][2]  += DrImm[ag][dr][1][s][0]*(1-p_Imm_SP);   // tx exp, sn
        V1[ag][5][dr][2][0][2]  += DrImm[ag][dr][1][s][0]*p_Imm_SP;  } }  // tx exp, sp

  /// EMIGRATION /// 
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { 
      V1[ag][tb][dr][tx][hv][2]  -= V0[ag][tb][dr][tx][hv][2]*rEmmigFB[0];      // FB1
      V1[ag][tb][dr][tx][hv][3]  -= V0[ag][tb][dr][tx][hv][3]*rEmmigFB[1];      // FB2
      } } } } }

  /// MORTALITY /// 
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { 
    if(hv==3) { temp = muTbH;  } else { temp = 0;  }
    for(int rg=0; rg<4; rg++) {  
      VMort[ag][0 ][dr][tx][hv][rg]  = V0[ag][0 ][dr][tx][hv][rg]*(mubtN[s][ag]*RRmuHR[rg]+vHMortN[ag][hv] );  
      VMort[ag][1 ][dr][tx][hv][rg]  = V0[ag][1 ][dr][tx][hv][rg]*(mubtN[s][ag]*RRmuHR[rg]+vHMortN[ag][hv] );  
      VMort[ag][2 ][dr][tx][hv][rg]  = V0[ag][2 ][dr][tx][hv][rg]*(mubtN[s][ag]*RRmuHR[rg]+vHMortN[ag][hv] );  
      VMort[ag][3 ][dr][tx][hv][rg]  = V0[ag][3 ][dr][tx][hv][rg]*(mubtN[s][ag]*RRmuHR[rg]+vHMortN[ag][hv] );  
      VMort[ag][4 ][dr][tx][hv][rg]  = V0[ag][4 ][dr][tx][hv][rg]*(mubtN[s][ag]*RRmuHR[rg]+vHMortN[ag][hv]+vTMortN[ag][4 ]+temp );  
      VMort[ag][5 ][dr][tx][hv][rg]  = V0[ag][5 ][dr][tx][hv][rg]*(mubtN[s][ag]*RRmuHR[rg]+vHMortN[ag][hv]+vTMortN[ag][5 ]+temp );  
      VMort[ag][6 ][dr][tx][hv][rg]  = V0[ag][6 ][dr][tx][hv][rg]*(mubtN[s][ag]*RRmuHR[rg]+vHMortN[ag][hv] );  
      VMort[ag][7 ][dr][tx][hv][rg]  = V0[ag][7 ][dr][tx][hv][rg]*(mubtN[s][ag]*RRmuHR[rg]+vHMortN[ag][hv]+(vTMortN[ag][7 ]+temp)*pow(1.0-TxMatZ[1][dr*2+0][0],TunTxMort));
      VMort[ag][8 ][dr][tx][hv][rg]  = V0[ag][8 ][dr][tx][hv][rg]*(mubtN[s][ag]*RRmuHR[rg]+vHMortN[ag][hv]+(vTMortN[ag][8 ]+temp)*pow(1.0-TxMatZ[1][dr*2+0][0],TunTxMort));
      VMort[ag][9 ][dr][tx][hv][rg]  = V0[ag][9 ][dr][tx][hv][rg]*(mubtN[s][ag]*RRmuHR[rg]+vHMortN[ag][hv]+(vTMortN[ag][9 ]+temp)*pow(1.0-TxMatZ[1][dr*2+1][0],TunTxMort));
      VMort[ag][10][dr][tx][hv][rg]  = V0[ag][10][dr][tx][hv][rg]*(mubtN[s][ag]*RRmuHR[rg]+vHMortN[ag][hv]+(vTMortN[ag][10]+temp)*pow(1.0-TxMatZ[1][dr*2+1][0],TunTxMort));

    for(int tb=0; tb<11; tb++) {
      V1[ag][tb][dr][tx][hv][rg]  -= VMort[ag][tb][dr][tx][hv][rg]; }
      } } } } }

  /// AGING ///
    for(int ag=0; ag<10; ag++) {
      if(ag>0) { temp2 = 120; } else { temp2 = 60; }
    for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      temp = V0[ag][tb][dr][tx][hv][rg]/temp2;
      V1[ag  ][tb][dr][tx][hv][rg]  -= temp; 
      V1[ag+1][tb][dr][tx][hv][rg]  += temp; 
      } } } } } }

  /// NEW FB -> ESTABLISHED FB ///
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { 
      temp = V0[ag][tb][dr][tx][hv][2]/24;
      V1[ag][tb][dr][tx][hv][2]  -= temp; 
      V1[ag][tb][dr][tx][hv][3]  += temp; 
      } } } } }

  /// HIGH-RISK ENTRY/EXIT ///
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { 
      temp  = V0[ag][tb][dr][tx][0][0]*HrEntExN[ag][0][0];
      temp2 = V0[ag][tb][dr][tx][0][1]*HrEntExN[ag][1][0];
      V1[ag][tb][dr][tx][0][0]  += temp2-temp; 
      V1[ag][tb][dr][tx][0][1]  += temp-temp2;
      for(int hv=1; hv<5 ; hv++) { 
      temp  = V0[ag][tb][dr][tx][hv][0]*HrEntExN[ag][0][1];
      temp2 = V0[ag][tb][dr][tx][hv][1]*HrEntExN[ag][1][0]; // note this stays the same as HIV neg ([0] suffix)
      V1[ag][tb][dr][tx][hv][0]  += temp2-temp; 
      V1[ag][tb][dr][tx][hv][1]  += temp-temp2;
      } } } } }

  /// TRANSMISSION RISK ///  
    for(int i=0; i<3; i++) { for(int j=0; j<2; j++) { VNkl[i][j] = 0; } } // set to zero
    for(int i=0; i<3; i++) { for(int j=0; j<2; j++) { for(int k=0; k<5; k++) { VGjkl[k][i][j] = 0; } } } // set to zero
    // Step 1
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { 
      VNkl[0][0]  += V0[ag][tb][dr][tx][0][0];
      VNkl[1][0]  += V0[ag][tb][dr][tx][0][1];
      VNkl[2][0]  += V0[ag][tb][dr][tx][0][2] + V0[ag][tb][dr][tx][0][3];
    for(int hv=1; hv<5 ; hv++) { 
      VNkl[0][1]  += V0[ag][tb][dr][tx][hv][0];
      VNkl[1][1]  += V0[ag][tb][dr][tx][hv][1];
      VNkl[2][1]  += V0[ag][tb][dr][tx][hv][2] + V0[ag][tb][dr][tx][hv][3];  }
      } } } }
    // Step 2 (active TB)
    for(int ag=0; ag<11; ag++) { for(int tb=4; tb<6; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) {
      VGjkl[dr][0][0]  += V0[ag][tb][dr][tx][0][0]*RelInfN[tb][dr];
      VGjkl[dr][1][0]  += V0[ag][tb][dr][tx][0][1]*RelInfN[tb][dr];
      VGjkl[dr][2][0]  += (V0[ag][tb][dr][tx][0][2] + V0[ag][tb][dr][tx][0][3])*RelInfN[tb][dr];
    for(int hv=1; hv<5 ; hv++) { 
      VGjkl[dr][0][1]  += V0[ag][tb][dr][tx][hv][0]*RelInfN[tb][dr];
      VGjkl[dr][1][1]  += V0[ag][tb][dr][tx][hv][1]*RelInfN[tb][dr];
      VGjkl[dr][2][1]  += (V0[ag][tb][dr][tx][hv][2] + V0[ag][tb][dr][tx][hv][3])*RelInfN[tb][dr];  }
      } } } }
    // Step 2 (treated TB)
      // No contribution to force of infection 

    // Step 3
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

  ///  INFECTION ///
    for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
    if(hv>0) { h2 = 1; } else { h2 = 0; }
    if(rg<3) { r2 = rg; } else { r2 = 2; }
    for(int ag=0; ag<11; ag++) {  for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { 
    for(int d2=0; d2<5; d2++) {
      temp = V0[ag][0][dr][tx][hv][rg]*VLjkl[d2][r2][h2]*NixTrans[s];  // Su
      V1[ag][0][dr][tx][hv][rg]  -= temp; 
      V1[ag][2][d2][tx][hv][rg]  += temp*MpslowN[ag][hv]; 
      V1[ag][3][d2][tx][hv][rg]  += temp*MpfastN[ag][hv];
      V1[ag][4][d2][tx][hv][rg]  += temp*MpimmedNn[ag][hv];
      V1[ag][5][d2][tx][hv][rg]  += temp*MpimmedNp[ag][hv];
      temp = V0[ag][1][dr][tx][hv][rg]*VLjkl[d2][r2][h2]*NixTrans[s];  // Sp
      V1[ag][1][dr][tx][hv][rg]  -= temp; 
      V1[ag][2][d2][tx][hv][rg]  += temp*MpslowPIN[ag][hv]; 
      V1[ag][3][d2][tx][hv][rg]  += temp*MpfastPIN[ag][hv];
      V1[ag][4][d2][tx][hv][rg]  += temp*MpimmedPINn[ag][hv];
      V1[ag][5][d2][tx][hv][rg]  += temp*MpimmedPINp[ag][hv];
      // SUPER-INFECTION LATENT SLOW
      temp = V0[ag][2][dr][tx][hv][rg]*VLjkl[d2][r2][h2]*NixTrans[s];  // Ls
      V1[ag][2][dr][tx][hv][rg]  -= temp; 
      V1[ag][2][dr][tx][hv][rg]  += temp*MpslowPIN[ag][hv]/2; 
      V1[ag][2][d2][tx][hv][rg]  += temp*MpslowPIN[ag][hv]/2; 
      V1[ag][3][d2][tx][hv][rg]  += temp*MpfastPIN[ag][hv];
      V1[ag][4][d2][tx][hv][rg]  += temp*MpimmedPINn[ag][hv];
      V1[ag][5][d2][tx][hv][rg]  += temp*MpimmedPINp[ag][hv];
      } } } } } }

  /// BREAKDOWN ///  
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      temp  = V0[ag][2][dr][tx][hv][rg]*MrslowN[ag][hv]*rrSlowFB[rg];  // Ls
      temp2 = V0[ag][3][dr][tx][hv][rg]*rfast;  // Lf
      V1[ag][2][dr][tx][hv][rg]  -= temp;
      V1[ag][3][dr][tx][hv][rg]  -= temp2;
      V1[ag][4][dr][tx][hv][rg]  += (temp+temp2)*(1-MpSmPosN[ag][hv]);
      V1[ag][5][dr][tx][hv][rg]  += (temp+temp2)*MpSmPosN[ag][hv];
      // Tl progression if INH resistant
      temp = V0[ag][6 ][dr][tx][hv][rg]*MrslowN[ag][hv]*rrSlowFB[rg]*(1-EffLtXN[s][dr]);
      V1[ag][6][dr][tx][hv][rg]  -= temp;
      V1[ag][4][dr][tx][hv][rg]  += temp*(1-MpSmPosN[ag][hv]);
      V1[ag][5][dr][tx][hv][rg]  += temp*MpSmPosN[ag][hv];
      } } } } }

  /// LATENT SLOW TO SAFE ///  
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      temp  = V0[ag][2][dr][tx][hv][rg]*rRecov;  // Ls
      V1[ag][2][dr][tx][hv][rg]  -= temp;
      V1[ag][1][dr][tx][hv][rg]  += temp;
      } } } } }

  /// SMEAR CONVERSION /// 
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      temp = V0[ag][4 ][dr][tx][hv][rg]*VrSmConv[hv]; 
      V1[ag][4 ][dr][tx][hv][rg]  -= temp;
      V1[ag][5 ][dr][tx][hv][rg]  += temp;
      } } } } }

  /// SELF CURE /// 
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      temp  = V0[ag][4 ][dr][tx][hv][rg]*VrSlfCur[hv]; 
      temp2 = V0[ag][5 ][dr][tx][hv][rg]*VrSlfCur[hv];
      V1[ag][4 ][dr][tx][hv][rg]  -= temp;
      V1[ag][5 ][dr][tx][hv][rg]  -= temp2;
      V1[ag][2 ][dr][tx][hv][rg]  += temp+temp2;
      } } } } }

  /// LTBI SCREENING AND TLTBI INITIATION /// only for no previous TB or LTBI tx
    for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      if(hv==0 & rg==0 ) { rTbP = rLtScrt[s]*LtDxParN[0][0]; rTbN = rLtScrt[s]*LtDxParN[0][1]; } 
      if(hv==0 & rg==1 ) { rTbP = rLtScrt[s]*LtDxParN[1][0]; rTbN = rLtScrt[s]*LtDxParN[1][1]; } 
      if(hv==0 & rg>1  ) { rTbP = rLtScrt[s]*LtDxParN[2][0]; rTbN = rLtScrt[s]*LtDxParN[2][1]; } 
      if(hv>0  & rg==0 ) { rTbP = rLtScrt[s]*LtDxParN[3][0]; rTbN = rLtScrt[s]*LtDxParN[3][1]; } 
      if(hv>0  & rg==1 ) { rTbP = rLtScrt[s]*LtDxParN[4][0]; rTbN = rLtScrt[s]*LtDxParN[4][1]; } 
      if(hv>0  & rg>1  ) { rTbP = rLtScrt[s]*LtDxParN[5][0]; rTbN = rLtScrt[s]*LtDxParN[5][1]; } 
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
      // Have LTBI
      temp  = V0[ag][2][dr][0][hv][rg]*rTbP; 
      temp2 = V0[ag][3][dr][0][hv][rg]*rTbP; 
      V1[ag][2][dr][0][hv][rg]  -= temp;
      V1[ag][3][dr][0][hv][rg]  -= temp2;
      V1[ag][6][dr][0][hv][rg]  += temp+temp2; 
      // Dont have LTBI
      temp  = V0[ag][0][dr][0][hv][rg]*rTbN; 
      temp2 = V0[ag][1][dr][0][hv][rg]*rTbN; 
      V1[ag][0][dr][0][hv][rg]  -= temp;
      V1[ag][1][dr][0][hv][rg]  -= temp2;
      V1[ag][0][dr][1][hv][rg]  += temp; 
      V1[ag][1][dr][1][hv][rg]  += temp2; 
    } } } }

  /// TLTBI: TX COMPLETION + DEFAULT /// only need to consider tx naive compartment
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
    for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      temp  = V0[ag][6][dr][0][hv][rg]*dLtt[s]; // tx completion
      temp2 = V0[ag][6][dr][0][hv][rg]*LtTxPar[1]; // default
      V1[ag][6 ][dr][0][hv][rg]  -= temp+temp2;
      V1[ag][1 ][dr][1][hv][rg]  += temp*EffLtXN[s][dr]*EffLt;
      V1[ag][2 ][dr][1][hv][rg]  += temp*(1-EffLtXN[s][dr]*EffLt);
      V1[ag][2 ][dr][0][hv][rg]  += temp2;
      } } } }

  /// TB DIAGNOSIS AND TX INITIATION ///  
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      if(dr==0) { p2ndL = 0; } else { p2ndL = pDstt[s];  }
      if(rg!=1) { ti = 0; } else { ti = 1;  }
      temp  = V0[ag][4 ][dr][tx][hv][rg]*rDxtN[s][ti  ]/RRdxAge[ag]; 
      temp2 = V0[ag][5 ][dr][tx][hv][rg]*rDxtN[s][ti+2]/RRdxAge[ag]; 
      V1[ag][4 ][dr][tx][hv][rg]  -= temp;
      V1[ag][5 ][dr][tx][hv][rg]  -= temp2;
      V1[ag][7 ][dr][tx][hv][rg]  += temp *(1-p2ndL); 
      V1[ag][8 ][dr][tx][hv][rg]  += temp2*(1-p2ndL); 
      V1[ag][9 ][dr][tx][hv][rg]  += temp *p2ndL; 
      V1[ag][10][dr][tx][hv][rg]  += temp2*p2ndL; 
      Vdx[ag][4][dr][tx][hv][rg]   = temp;
      Vdx[ag][5][dr][tx][hv][rg]   = temp2;
      } } } } }

  /// HIV PROGRESSION ///
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int rg=0; rg<4; rg++) {
      temp  = V0[ag][tb][dr][tx][1 ][rg]*vHxtoHyN[ag][0]; 
      temp2 = V0[ag][tb][dr][tx][2 ][rg]*vHxtoHyN[ag][1]; 
      V1[ag][tb][dr][tx][1 ][rg]  -= temp;
      V1[ag][tb][dr][tx][3 ][rg]  += temp;
      V1[ag][tb][dr][tx][2 ][rg]  -= temp2;
      V1[ag][tb][dr][tx][4 ][rg]  += temp2;
      } } } } } 

  /// HIV INCIDENCE /// 
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) {  for(int rg=0; rg<4; rg++) {
      if(rg!=1) { 
        temp = V0[ag][tb][dr][tx][0 ][rg]*rHIVtN[s][ag][0]; 
      } else {
        temp = V0[ag][tb][dr][tx][0 ][rg]*rHIVtN[s][ag][1]; }
      V1[ag][tb][dr][tx][0 ][rg]  -= temp;
      V1[ag][tb][dr][tx][1 ][rg]  += temp;
      } } } } }

  /// ART ENROLLMENT ///
    // rate differs by HIV category, risk group, tb treatment 
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int rg=0; rg<4; rg++) { for(int hv2=0; hv2<2; hv2++) {
      if( tb>6 | hv2==1 ) { ti = 1; } else { ti = 0; }   // ti = 0 for slow rate, 1 for higher rate (advanced HIV and/or diag TB)
      if( rg==1 ) { ti += 2; } // 3rd and 4th col of rArtInit for HR group
      temp = V0[ag][tb][dr][tx][1+hv2*2][rg]*rArtInitN[s][ti];
      V1[ag][tb][dr][tx][1+hv2*2][rg]  -= temp;
      V1[ag][tb][dr][tx][2+hv2*2][rg]  += temp;
      } } } } } }

  /// ART DEFAULT ///   
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int rg=0; rg<4; rg++) {
      if(rg!=1) { temp3=rArtDef; } else { temp3=rArtDef*2; }
      temp  = V0[ag][tb][dr][tx][2 ][rg]*temp3; 
      temp2 = V0[ag][tb][dr][tx][4 ][rg]*temp3; 
      V1[ag][tb][dr][tx][2 ][rg]  -= temp;
      V1[ag][tb][dr][tx][1 ][rg]  += temp;
      V1[ag][tb][dr][tx][4 ][rg]  -= temp2;
      V1[ag][tb][dr][tx][3 ][rg]  += temp2;
      } } } } }

  /// TB TREATMENT OUTCOMES ///  
    for(int rg=0; rg<4; rg++) {
      if(rg!=1) { ti2 = 0; } else { ti2 = 1; }
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { 
  // Cures back to Ls state, treatment experienced subdivision
      for(int i=0; i<2; i++) {
      ti3 = dr*2+i; ti4 = 7+i*2; ti5 = 8+i*2;
      if(i==0) { ti = 7; } else { ti = 9; }
      temp  = V0[ag][ti  ][dr][tx][hv][rg]*TxMatZ[7][ti3][ti2]; 
      temp2 = V0[ag][ti+1][dr][tx][hv][rg]*TxMatZ[7][ti3][ti2]; 
      V1[ag][ti4][dr][tx][hv][rg]  -= temp;
      V1[ag][ti5][dr][tx][hv][rg]  -= temp2;
      V1[ag][2  ][dr][2 ][hv][rg]  += temp + temp2;  
  // Failures stay in naive (if naive) as will still be one episode of tb
      for(int d2=0; d2<5; d2++) {  if(idAR[d2][ti3]==1) {
    // To active disease state, with AR distribution
      temp  = V0[ag][ti  ][dr][tx][hv][rg]*TxMatZ[10+d2][ti3][ti2];  // sm neg
      temp2 = V0[ag][ti+1][dr][tx][hv][rg]*TxMatZ[10+d2][ti3][ti2];  // sm pos
      V1[ag][ti4][dr][tx][hv][rg]  -= temp;
      V1[ag][4  ][d2][tx][hv][rg]  += temp;
      V1[ag][ti5][dr][tx][hv][rg]  -= temp2;
      V1[ag][5  ][d2][tx][hv][rg]  += temp2;
    // To retreatment, with AR distribution
      temp  = V0[ag][ti  ][dr][tx][hv][rg]*TxMatZ[15+d2][ti3][ti2];   // sm neg
      temp2 = V0[ag][ti+1][dr][tx][hv][rg]*TxMatZ[15+d2][ti3][ti2];   // sm pos
      V1[ag][ti4][dr][tx][hv][rg]  -= temp;
      V1[ag][9  ][d2][tx][hv][rg]  += temp;
      V1[ag][ti5][dr][tx][hv][rg]  -= temp2;
      V1[ag][10 ][d2][tx][hv][rg]  += temp2;    } } }
      } } } } }

  //// FILL RESULTS TABLE   ~~~~~~ 

    // MID-YEAR RESULTS
    if(m==6) {
    // YEAR
      Outputs[y][0]      = y+1950;  // Year
    // COUNTS BY TOTAL, AGE, TB, HIV, AND RG
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      Outputs[y][1    ] += V1[ag][tb][dr][tx][hv][rg];   // N_ALL
      Outputs[y][2 +ag] += V1[ag][tb][dr][tx][hv][rg];   // N_ by age (11)
      Outputs[y][13+tb] += V1[ag][tb][dr][tx][hv][rg];   // N_ by tb (11)
      Outputs[y][24+hv] += V1[ag][tb][dr][tx][hv][rg];   // N_ by hv (5)
      Outputs[y][29+rg] += V1[ag][tb][dr][tx][hv][rg];   // N_ by rg (4)
      } } } } } }

    // US AND FB COUNTS BY AGE
    for(int ag=0; ag<11; ag++) { 
    	for(int tb=0; tb<11; tb++) { 
    		for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { 
    	for(int hv=0; hv<5 ; hv++) { 
      Outputs[y][33+ag] += V1[ag][tb][dr][tx][hv][0]+V1[ag][tb][dr][tx][hv][1];   // N_ by age and US (11)
      Outputs[y][44+ag] += V1[ag][tb][dr][tx][hv][2]+V1[ag][tb][dr][tx][hv][3];   // N_ by age and FB (11)
      } } } } }

    // US AND FB COUNTS BY AGE, LTBI
    for(int ag=0; ag<11; ag++) {  for(int dr=0; dr<5; dr++) { for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { 
      Outputs[y][55+ag] += (V1[ag][1][dr][tx][hv][0]+V1[ag][1][dr][tx][hv][1])*(1-pImmScen)+
                            V1[ag][2][dr][tx][hv][0]+V1[ag][2][dr][tx][hv][1]+
                            V1[ag][3][dr][tx][hv][0]+V1[ag][3][dr][tx][hv][1]+
                            V1[ag][6][dr][tx][hv][0]+V1[ag][6][dr][tx][hv][1];   // N_ by age and US (11) LATENT INFECTION
      Outputs[y][66+ag] += (V1[ag][1][dr][tx][hv][2]+V1[ag][1][dr][tx][hv][3])*(1-pImmScen)+
                            V1[ag][2][dr][tx][hv][2]+V1[ag][2][dr][tx][hv][3]+
                            V1[ag][3][dr][tx][hv][2]+V1[ag][3][dr][tx][hv][3]+
                            V1[ag][6][dr][tx][hv][2]+V1[ag][6][dr][tx][hv][3];   // N_ by age and FB (11) LATENT INFECTION
      } } } } 

    // HIV COUNT BY AGE
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=1; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      Outputs[y][77+ag] += V1[ag][tb][dr][tx][hv][rg];   // N_HIV by age (11)
      } } } } } }

    // TB MORTALITY BY AGE AND HIV
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      if(hv>0) { ti = 11; } else { ti = 0; }
      if(hv==3) { temp = muTbH; } else { temp = 0; }
      Outputs[y][88+ag+ti]  += V0[ag][4 ][dr][tx][hv][rg]*(vTMortN[ag][4 ]+temp);  
      Outputs[y][88+ag+ti]  += V0[ag][5 ][dr][tx][hv][rg]*(vTMortN[ag][5 ]+temp);  
      Outputs[y][88+ag+ti]  += V0[ag][7 ][dr][tx][hv][rg]*(vTMortN[ag][7 ]+temp)*pow(1.0-TxMatZ[1][dr*2+0][0],TunTxMort);
      Outputs[y][88+ag+ti]  += V0[ag][8 ][dr][tx][hv][rg]*(vTMortN[ag][8 ]+temp)*pow(1.0-TxMatZ[1][dr*2+0][0],TunTxMort);
      Outputs[y][88+ag+ti]  += V0[ag][9 ][dr][tx][hv][rg]*(vTMortN[ag][9 ]+temp)*pow(1.0-TxMatZ[1][dr*2+1][0],TunTxMort);
      Outputs[y][88+ag+ti]  += V0[ag][10][dr][tx][hv][rg]*(vTMortN[ag][10]+temp)*pow(1.0-TxMatZ[1][dr*2+1][0],TunTxMort);
      } } } } }
    for(int i=88; i<110; i++) { Outputs[y][i] = Outputs[y][i]*12; }

    // HIV MORTALITY BY AGE 
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      Outputs[y][110+ag]  += V0[ag][tb][dr][tx][hv][rg]*vHMortN[ag][hv];
      } } } } } }
    for(int i=110; i<121; i++) { Outputs[y][i] = Outputs[y][i]*12;  }

    // TOTAL MORTALITY BY AGE
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      Outputs[y][121+ag]  += VMort[ag][tb][dr][tx][hv][rg];
      } } } } } }
    for(int i=121; i<132; i++) { Outputs[y][i] = Outputs[y][i]*12; }

    // TB TREATMENT OUTCOMES
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      if(rg!=1) { temp = rDeft[s]; } else { temp = rDeftH[s]; }
      Outputs[y][132]  += (V0[ag][7 ][dr][tx][hv][rg]+V0[ag][8 ][dr][tx][hv][rg])*TxMatZ[22][dr*2+0][0]+
                          (V0[ag][9 ][dr][tx][hv][rg]+V0[ag][10][dr][tx][hv][rg])*TxMatZ[22][dr*2+1][0];  // tx completion
      Outputs[y][133]  += (V0[ag][7 ][dr][tx][hv][rg]+V0[ag][8 ][dr][tx][hv][rg]+
                           V0[ag][9 ][dr][tx][hv][rg]+V0[ag][10][dr][tx][hv][rg])*temp;  // tx discontinuation
      Outputs[y][134]  += VMort[ag][7 ][dr][tx][hv][rg]+VMort[ag][8 ][dr][tx][hv][rg]+
                          VMort[ag][9 ][dr][tx][hv][rg]+VMort[ag][10][dr][tx][hv][rg];  // tx mort
      } } } } }
    for(int i=132; i<135; i++) { Outputs[y][i] = Outputs[y][i]*12; }

    // NOTIFICATIONS 
    for(int ag=0; ag<11; ag++) { for(int tb=4; tb<6; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      Outputs[y][135   ] += Vdx[ag][tb][dr][tx][hv][rg];   // All dx (1)
      Outputs[y][136+ag] += Vdx[ag][tb][dr][tx][hv][rg];   // dx by age (11)
      if(tx<2) {
      Outputs[y][147   ] += Vdx[ag][tb][dr][tx][hv][rg];   // dx by tx history - naive (2)
      } else {
      Outputs[y][148   ] += Vdx[ag][tb][dr][tx][hv][rg]; } // dx by tx history - experienced (2)
      if(hv>0) {
      Outputs[y][149   ] += Vdx[ag][tb][dr][tx][hv][rg];   // dx HIV pos (1)
      } else {
      Outputs[y][150   ] += Vdx[ag][tb][dr][tx][hv][rg]; }  // dx HIV neg (1)
      Outputs[y][151+rg] += Vdx[ag][tb][dr][tx][hv][rg];   // N_ by rg (4)
      } } } } } }
    for(int i=135; i<155; i++) { Outputs[y][i] = Outputs[y][i]*12; }

   // DR DISTRIBUTION IN NOTIFICATIONS 
    for(int ag=0; ag<11; ag++) { for(int tb=4; tb<6; tb++) {  
    for(int dr=0; dr<5; dr++) {for(int hv=0; hv<5 ; hv++) {
      Outputs[y][155+dr] += Vdx[ag][tb][dr][0][hv][0]+Vdx[ag][tb][dr][0][hv][1]+Vdx[ag][tb][dr][1][hv][0]+Vdx[ag][tb][dr][1][hv][1];   // US N dx by dr (5) 
      Outputs[y][160+dr] += Vdx[ag][tb][dr][2][hv][0]+Vdx[ag][tb][dr][2][hv][1];   // US E dx by dr (5)
      Outputs[y][165+dr] += Vdx[ag][tb][dr][0][hv][2]+Vdx[ag][tb][dr][0][hv][3]+Vdx[ag][tb][dr][1][hv][2]+Vdx[ag][tb][dr][1][hv][3];   // FB N dx by dr (5)
      Outputs[y][170+dr] += Vdx[ag][tb][dr][2][hv][2]+Vdx[ag][tb][dr][2][hv][3];   // FB E dx by dr (5)
      } } } }
    for(int i=155; i<175; i++) { Outputs[y][i] = Outputs[y][i]*12; }

  /// TLTBI INITS /// 
    for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      if(hv==0 & rg==0 ) { rTbP = rLtScrt[s]*LtDxParN[0][0]; rTbN = rLtScrt[s]*LtDxParN[0][1]; } 
      if(hv==0 & rg==1 ) { rTbP = rLtScrt[s]*LtDxParN[1][0]; rTbN = rLtScrt[s]*LtDxParN[1][1]; } 
      if(hv==0 & rg>1  ) { rTbP = rLtScrt[s]*LtDxParN[2][0]; rTbN = rLtScrt[s]*LtDxParN[2][1]; } 
      if(hv>0  & rg==0 ) { rTbP = rLtScrt[s]*LtDxParN[3][0]; rTbN = rLtScrt[s]*LtDxParN[3][1]; } 
      if(hv>0  & rg==1 ) { rTbP = rLtScrt[s]*LtDxParN[4][0]; rTbN = rLtScrt[s]*LtDxParN[4][1]; } 
      if(hv>0  & rg>1  ) { rTbP = rLtScrt[s]*LtDxParN[5][0]; rTbN = rLtScrt[s]*LtDxParN[5][1]; } 
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
      Outputs[y][175] += (V0[ag][3][dr][0][hv][rg]+V0[ag][2][dr][0][hv][rg])*rTbP + (V0[ag][1][dr][0][hv][rg]+V0[ag][0][dr][0][hv][rg])*rTbN;  // all inits
      if(rg>1) {
      Outputs[y][176] += (V0[ag][3][dr][0][hv][rg]+V0[ag][2][dr][0][hv][rg])*rTbP + (V0[ag][1][dr][0][hv][rg]+V0[ag][0][dr][0][hv][rg])*rTbN; } // FB inits
      if(rg==1) {
      Outputs[y][177] += (V0[ag][3][dr][0][hv][rg]+V0[ag][2][dr][0][hv][rg])*rTbP + (V0[ag][1][dr][0][hv][rg]+V0[ag][0][dr][0][hv][rg])*rTbN; } // FB inits
      if(hv>0) {
      Outputs[y][178] += (V0[ag][3][dr][0][hv][rg]+V0[ag][2][dr][0][hv][rg])*rTbP + (V0[ag][1][dr][0][hv][rg]+V0[ag][0][dr][0][hv][rg])*rTbN; } // FB inits
      Outputs[y][179] += (V0[ag][3][dr][0][hv][rg]+V0[ag][2][dr][0][hv][rg])*rTbP;  // inits with LTBI
    } } } }
    for(int i=175; i<180; i++) { Outputs[y][i] = Outputs[y][i]*12; } // annualize

  /// TB INCIDENCE, BY ALL VS RECENT  ///
  // By recency (<2 years) == all immediate, 1-(1-rfast)^24 x all Lf
  //  INFECTION, IMMEDIATE PROGRESSION 
    for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
    if(hv>0) { h2 = 1; } else { h2 = 0; }
    if(rg<3) { r2 = rg; } else { r2 = 2; }
    for(int ag=0; ag<11; ag++) {  for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { 
    for(int d2=0; d2<5; d2++) {
      temp  = V0[ag][0][dr][tx][hv][rg]*VLjkl[d2][r2][h2]*(MpimmedNn[ag][hv]  +MpimmedNp[ag][hv]  ) +  // immediate transition Su
              V0[ag][1][dr][tx][hv][rg]*VLjkl[d2][r2][h2]*(MpimmedPINn[ag][hv]+MpimmedPINp[ag][hv]) +  // immediate transition Sp
              V0[ag][2][dr][tx][hv][rg]*VLjkl[d2][r2][h2]*(MpimmedPINn[ag][hv]+MpimmedPINp[ag][hv]);   // immediate transition Ls
    Outputs[y][180]    += temp;  // all incidence
    Outputs[y][181+ag] += temp;  // incidence by age
    if(rg<2) {
    Outputs[y][192] += temp;      //  incidence, US born
    } else { 
    Outputs[y][193] += temp;      //  incidence, FB
    if(rg==3) {
    Outputs[y][194] += temp; } }  //  incidence, FB2
    if(rg==1) {
    Outputs[y][195] += temp; }    //  incidence, HR
    if(hv>0) {
    Outputs[y][196] += temp; }    //  incidence, HIV pos
      } } } } } }
    for(int i=180; i<197; i++) { Outputs[y][i+17] = Outputs[y][i]; } // copy these over to recent infection lines

// BREAKDOWN from Ls, Lf 
    temp3 = (1-pow(1-rfast,24.0))*rfast;
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      temp2 = V0[ag][2][dr][tx][hv][rg]*temp4V[ag][hv] + V0[ag][3][dr][tx][hv][rg]*temp3;  // Progression from recent infection
      temp  = V0[ag][2][dr][tx][hv][rg]*MrslowN[ag][hv] + V0[ag][3][dr][tx][hv][rg]*rfast;  // All progression
      if( (dr==1) | (dr>2) ) { 
      temp = temp + V0[ag][6][dr][tx][hv][rg]*MrslowN[ag][hv]; } // Progression from Tl for INH resistant

    Outputs[y][180      ] += temp ;   // all incidence
    Outputs[y][180+17   ] += temp2;   // all incidence, recent infection
    Outputs[y][181+ag   ] += temp ;   // incidence by age
    Outputs[y][181+ag+17] += temp2;   // incidence by age, recent infection
    if(rg<2) {
    Outputs[y][192      ] += temp ;   //  incidence, US born
    Outputs[y][192+17   ] += temp2;   //  incidence, US born, recent infection
    } else {
    Outputs[y][193      ] += temp ;   //  incidence, FB
    Outputs[y][193+17   ] += temp2;  //  incidence, FB, recent infection
    if(rg==3) {
    Outputs[y][194      ] += temp ;   //  incidence, FB2 born
    Outputs[y][194+17   ] += temp2; } } //  incidence, FB2 born, recent infection
    if(rg==1) {
    Outputs[y][195      ] += temp ;   //  incidence, HR
    Outputs[y][195+17   ] += temp2; } //  incidence, HR, recent infection
    if(hv>0) {
    Outputs[y][196      ] += temp ;   //  incidence, HIV pos
    Outputs[y][196+17   ] += temp2; } //  incidence, HIV pos, recent infection
    } } } } }

    for(int i=180; i<214; i++) { Outputs[y][i] = Outputs[y][i]*12; } // annualize

    // NOTIFICATIONS, dead at diagnosis
    for(int hv=0; hv<5 ; hv++) { 
    if(hv==3) { temp = muTbH;  } else { temp = 0;  }
    for(int ag=0; ag<11; ag++) { for(int tb=4; tb<6; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int rg=0; rg<4; rg++) {
      temp2 = V0[ag][tb][dr][tx][hv][rg]*(vTMortN[ag][tb]+temp);
      Outputs[y][214   ] += temp2;   // All dx (1)
      Outputs[y][215+ag] += temp2;   // dx by age (11)
      if(tx<2) {
      Outputs[y][226   ] += temp2;   // dx by tx history - naive (2)
      } else {
      Outputs[y][227   ] += temp2; } // dx by tx history - experienced (2)
      if(hv>0) {
      Outputs[y][228   ] += temp2;   // dx HIV pos (1)
      } else {
      Outputs[y][229   ] += temp2; }  // dx HIV neg (1)
      Outputs[y][230+rg] += temp2;   // N_ by rg (4)
      } } } } } }
    for(int i=214; i<234; i++) { Outputs[y][i] = Outputs[y][i]*12; }

    // NOTIFICATIONS US
    for(int ag=0; ag<11; ag++) { for(int tb=4; tb<6; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<2; rg++) {
      Outputs[y][234+ag] += Vdx[ag][tb][dr][tx][hv][rg];   // dx by age (11)
      } } } } } }
    for(int i=234; i<245; i++) { Outputs[y][i] = Outputs[y][i]*12; }

    // NOTIFICATIONS US, dead at diagnosis
    for(int hv=0; hv<5 ; hv++) { 
    if(hv==3) { temp = muTbH;  } else { temp = 0;  }
    for(int ag=0; ag<11; ag++) { for(int tb=4; tb<6; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int rg=0; rg<2; rg++) {
      temp2 = V0[ag][tb][dr][tx][hv][rg]*(vTMortN[ag][tb]+temp);
      Outputs[y][245+ag] += temp2;   // dx by age (11)
      } } } } } }
    for(int i=245; i<256; i++) { Outputs[y][i] = Outputs[y][i]*12; }

    // NOTIFICATIONS, dead at diagnosis  HIV_NEGATIVE
    for(int ag=0; ag<11; ag++) { for(int tb=4; tb<6; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int rg=0; rg<4; rg++) {
      temp2 = V0[ag][tb][dr][tx][0 ][rg]*(vTMortN[ag][tb]+temp);
      Outputs[y][256+ag] += temp2;   // dx by age (11)
      } } } } }
    for(int i=256; i<267; i++) { Outputs[y][i] = Outputs[y][i]*12; }

    // TOTAL MORTALITY BY AGE, HAVE HIV
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=1; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      Outputs[y][267+ag]  += VMort[ag][tb][dr][tx][hv][rg];
      } } } } } }
    for(int i=267; i<278; i++) { Outputs[y][i] = Outputs[y][i]*12; }

    // TOTAL MORTALITY BY AGE, HAVE TB
    for(int ag=0; ag<11; ag++) {  for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
    for(int tb=4; tb<6; tb++) {
      Outputs[y][278+ag]  += VMort[ag][tb][dr][tx][hv][rg]; }
    for(int tb=7; tb<11; tb++) {
      Outputs[y][278+ag]  += VMort[ag][tb][dr][tx][hv][rg]; }
      } } } } }
    for(int i=278; i<289; i++) { Outputs[y][i] = Outputs[y][i]*12; }

  // COUNTS TB BY US/FB 
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
      if(rg<2) {
      Outputs[y][289] += V1[ag][2][dr][tx][hv][rg];
      Outputs[y][290] += V1[ag][3][dr][tx][hv][rg];
      Outputs[y][291] += V1[ag][4][dr][tx][hv][rg];
      Outputs[y][292] += V1[ag][5][dr][tx][hv][rg];
      } else {
      Outputs[y][293] += V1[ag][2][dr][tx][hv][rg];
      Outputs[y][294] += V1[ag][3][dr][tx][hv][rg];
      Outputs[y][295] += V1[ag][4][dr][tx][hv][rg];
      Outputs[y][296] += V1[ag][5][dr][tx][hv][rg]; }
      } } } } }
  // Force of infection
    for(int dr=0; dr<5; dr++) { 
      Outputs[y][297] += VLjkl[dr][0][0];
      Outputs[y][298] += VLjkl[dr][1][0];
      Outputs[y][299] += VLjkl[dr][2][0];
      Outputs[y][300] += VLjkl[dr][0][1];
      Outputs[y][301] += VLjkl[dr][1][1];
      Outputs[y][302] += VLjkl[dr][2][1];      }
    for(int i=297; i<303; i++) { Outputs[y][i] = Outputs[y][i]*12; }

///  NEW INFECTIONS + SUPER INFECTIOPN ///
    for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
            if(hv>0) { h2 = 1; } else { h2 = 0; }
            if(rg<3) { r2 = rg; } else { r2 = 2; }
            for(int ag=0; ag<11; ag++) {  for(int dr=0; dr<5; dr++) { 
            for(int tx=0; tx<3 ; tx++) { 
            for(int d2=0; d2<5; d2++) {
            Outputs[y][303+rg] += (V0[ag][0][dr][tx][hv][rg]+V0[ag][1][dr][tx][hv][rg]+V0[ag][2][dr][tx][hv][rg])*VLjkl[d2][r2][h2]*NixTrans[s];
            } } } } } }
    for(int i=303; i<307; i++) { Outputs[y][i] = Outputs[y][i]*12; }

    // TB MORTALITY BY NATIVITY
    for(int ag=0; ag<11; ag++) { for(int dr=0; dr<5; dr++) { 
            for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
            if(rg>1) { ti = 1; } else { ti = 0; }
            if(hv==3) { temp = muTbH; } else { temp = 0; }
            Outputs[y][307+ti]  += V0[ag][4 ][dr][tx][hv][rg]*(vTMortN[ag][4 ]+temp);  
            Outputs[y][307+ti]  += V0[ag][5 ][dr][tx][hv][rg]*(vTMortN[ag][5 ]+temp);  
            Outputs[y][307+ti]  += V0[ag][7 ][dr][tx][hv][rg]*(vTMortN[ag][7 ]+temp)*pow(1.0-TxMatZ[1][dr*2+0][0],TunTxMort);
            Outputs[y][307+ti]  += V0[ag][8 ][dr][tx][hv][rg]*(vTMortN[ag][8 ]+temp)*pow(1.0-TxMatZ[1][dr*2+0][0],TunTxMort);
            Outputs[y][307+ti]  += V0[ag][9 ][dr][tx][hv][rg]*(vTMortN[ag][9 ]+temp)*pow(1.0-TxMatZ[1][dr*2+1][0],TunTxMort);
            Outputs[y][307+ti]  += V0[ag][10][dr][tx][hv][rg]*(vTMortN[ag][10]+temp)*pow(1.0-TxMatZ[1][dr*2+1][0],TunTxMort);
            } } } } }
            for(int i=307; i<309; i++) { Outputs[y][i] = Outputs[y][i]*12; }

 //  LTBI 
    for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) { for(int ag=0; ag<11; ag++) { 
    for(int dr=0; dr<5; dr++) { for(int tx=0; tx<3 ; tx++) { 
            temp  =  (V1[ag][1][dr][tx][hv][rg])*(1-pImmScen)+V1[ag][2][dr][tx][hv][rg]+V1[ag][3][dr][tx][hv][rg]+V1[ag][6][dr][tx][hv][rg];
            Outputs[y][309]    += temp;  // all LTBI
            Outputs[y][310+ag] += temp;  // LTBI by age
            if(rg<2) {
            Outputs[y][321] += temp;      //  LTBI, US born
            } else { 
            Outputs[y][322] += temp;      //  LTBI, FB
            if(rg==3) {
            Outputs[y][323] += temp; } }  //  LTBI, FB2
            if(rg==1) {
            Outputs[y][324] += temp; }    //  LTBI, HR
            if(hv>0) {
            Outputs[y][325] += temp; }    //  LTBI, HIV pos
            } } } } }

            } // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ end mid-year stock results

     //// UPDATE V0 as V1
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
     V0[ag][tb][dr][tx][hv][rg] = V1[ag][tb][dr][tx][hv][rg];
      } } } } } }           
            } //// end of month loop! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            } //// end of year loop! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Cumulative incidence plots
    NumericMatrix CumInc(480,8);
      for(int aa=0; aa<480; aa++) { for(int bb=0; bb<8; bb++) { CumInc(aa,bb) = 0; } }
    double        Vzz1[11][4]; // Su[0 to 0], Ls[2 to 1], Lf[3 to 2], TB[4/5 to 3]
    double        Vzz0[11][4]; 

/// 25-35 adult uninfected USB   ############################################################################
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb] = 0; } }
  Vzz0[3 ][0] = Vzz1[3 ][0] = 1;  // assign starting position

  for(int ss=0; ss<480; ss++) {  // start month loop

///  INFECTION /// StatList
       if(ss==0) { for(int ag=0; ag<11; ag++) { 
        temp          = Vzz0[ag][0];
        Vzz1[ag][0]  -= temp; 
        Vzz1[ag][1]  += temp*MpslowN[ag][0]; 
        Vzz1[ag][2]  += temp*MpfastN[ag][0];
        Vzz1[ag][3]  += temp*(MpimmedNn[ag][0]+MpimmedNp[ag][0]); }  }

/// BREAKDOWN ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*MrslowN[ag][0]*rrSlowFB[0];
      temp2 = Vzz0[ag][2]*rfast;
      Vzz1[ag][1]  -= temp; 
      Vzz1[ag][2]  -= temp2; 
      Vzz1[ag][3]  += temp+temp2;   } 

  /// LATENT SLOW TO SAFE ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*rRecov;
      Vzz1[ag][1]  -= temp;  } 

  /// AGING ///
    for(int ag=0; ag<10; ag++) { for(int tb=0; tb<4; tb++) {
      if(ag>0) { temp2 = 120; } else { temp2 = 60; }
      temp = Vzz0[ag][tb]/temp2;
      Vzz1[ag  ][tb]  -= temp; 
      Vzz1[ag+1][tb]  += temp;  } }
  
  /// add up results ///
  for(int ag=0; ag<11; ag++) {
    CumInc(ss,0) +=  Vzz1[ag][3]; }

  /// reset vectors ///
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb]; } }

  } // end month loop

/// 0-4 adult uninfected USB   ############################################################################
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb] = 0; } }
  Vzz0[0 ][0] = Vzz1[0 ][0] = 1;  // assign starting position

  for(int ss=0; ss<480; ss++) {  // start month loop

///  INFECTION /// StatList
       if(ss==0) { for(int ag=0; ag<11; ag++) { 
        temp          = Vzz0[ag][0];
        Vzz1[ag][0]  -= temp; 
        Vzz1[ag][1]  += temp*MpslowN[ag][0]; 
        Vzz1[ag][2]  += temp*MpfastN[ag][0];
        Vzz1[ag][3]  += temp*(MpimmedNn[ag][0]+MpimmedNp[ag][0]); }  }

/// BREAKDOWN ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*MrslowN[ag][0]*rrSlowFB[0];
      temp2 = Vzz0[ag][2]*rfast;
      Vzz1[ag][1]  -= temp; 
      Vzz1[ag][2]  -= temp2; 
      Vzz1[ag][3]  += temp+temp2;   } 

  /// LATENT SLOW TO SAFE ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*rRecov;
      Vzz1[ag][1]  -= temp;  } 

  /// AGING ///
    for(int ag=0; ag<10; ag++) { for(int tb=0; tb<4; tb++) {
      if(ag>0) { temp2 = 120; } else { temp2 = 60; }
      temp = Vzz0[ag][tb]/temp2;
      Vzz1[ag  ][tb]  -= temp; 
      Vzz1[ag+1][tb]  += temp;  } }
  
  /// add up results ///
  for(int ag=0; ag<11; ag++) {
    CumInc(ss,1) +=  Vzz1[ag][3]; }

  /// reset vectors ///
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb]; } }

  } // end month loop

/// 5-14 adult uninfected USB   ############################################################################
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb] = 0; } }
  Vzz0[1 ][0] = Vzz1[1 ][0] = 1;  // assign starting position

  for(int ss=0; ss<480; ss++) {  // start month loop

///  INFECTION /// StatList
       if(ss==0) { for(int ag=0; ag<11; ag++) { 
        temp          = Vzz0[ag][0];
        Vzz1[ag][0]  -= temp; 
        Vzz1[ag][1]  += temp*MpslowN[ag][0]; 
        Vzz1[ag][2]  += temp*MpfastN[ag][0];
        Vzz1[ag][3]  += temp*(MpimmedNn[ag][0]+MpimmedNp[ag][0]); }  }

/// BREAKDOWN ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*MrslowN[ag][0]*rrSlowFB[0];
      temp2 = Vzz0[ag][2]*rfast;
      Vzz1[ag][1]  -= temp; 
      Vzz1[ag][2]  -= temp2; 
      Vzz1[ag][3]  += temp+temp2;   } 

  /// LATENT SLOW TO SAFE ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*rRecov;
      Vzz1[ag][1]  -= temp;  } 

  /// AGING ///
    for(int ag=0; ag<10; ag++) { for(int tb=0; tb<4; tb++) {
      if(ag>0) { temp2 = 120; } else { temp2 = 60; }
      temp = Vzz0[ag][tb]/temp2;
      Vzz1[ag  ][tb]  -= temp; 
      Vzz1[ag+1][tb]  += temp;  } }
  
  /// add up results ///
  for(int ag=0; ag<11; ag++) {
    CumInc(ss,2) +=  Vzz1[ag][3]; }

  /// reset vectors ///
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb]; } }

  } // end month loop

/// 75-85 adult uninfected USB   ############################################################################
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb] = 0; } }
  Vzz0[8 ][0] = Vzz1[8 ][0] = 1;  // assign starting position

  for(int ss=0; ss<480; ss++) {  // start month loop

///  INFECTION /// StatList
       if(ss==0) { for(int ag=0; ag<11; ag++) { 
        temp          = Vzz0[ag][0];
        Vzz1[ag][0]  -= temp; 
        Vzz1[ag][1]  += temp*MpslowN[ag][0]; 
        Vzz1[ag][2]  += temp*MpfastN[ag][0];
        Vzz1[ag][3]  += temp*(MpimmedNn[ag][0]+MpimmedNp[ag][0]); }  }

/// BREAKDOWN ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*MrslowN[ag][0]*rrSlowFB[0];
      temp2 = Vzz0[ag][2]*rfast;
      Vzz1[ag][1]  -= temp; 
      Vzz1[ag][2]  -= temp2; 
      Vzz1[ag][3]  += temp+temp2;   } 

  /// LATENT SLOW TO SAFE ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*rRecov;
      Vzz1[ag][1]  -= temp;  } 

  /// AGING ///
    for(int ag=0; ag<10; ag++) { for(int tb=0; tb<4; tb++) {
      if(ag>0) { temp2 = 120; } else { temp2 = 60; }
      temp = Vzz0[ag][tb]/temp2;
      Vzz1[ag  ][tb]  -= temp; 
      Vzz1[ag+1][tb]  += temp;  } }
  
  /// add up results ///
  for(int ag=0; ag<11; ag++) {
    CumInc(ss,3) +=  Vzz1[ag][3]; }

  /// reset vectors ///
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb]; } }

  } // end month loop

/// 25-35 adult uninfected FB   ############################################################################
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb] = 0; } }
  Vzz0[3 ][0] = Vzz1[3 ][0] = 1;  // assign starting position

  for(int ss=0; ss<480; ss++) {  // start month loop

///  INFECTION /// StatList
       if(ss==0) { for(int ag=0; ag<11; ag++) { 
        temp          = Vzz0[ag][0];
        Vzz1[ag][0]  -= temp; 
        Vzz1[ag][1]  += temp*MpslowN[ag][0]; 
        Vzz1[ag][2]  += temp*MpfastN[ag][0];
        Vzz1[ag][3]  += temp*(MpimmedNn[ag][0]+MpimmedNp[ag][0]); }  }

/// BREAKDOWN ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*MrslowN[ag][0]*rrSlowFB[2];
      temp2 = Vzz0[ag][2]*rfast;
      Vzz1[ag][1]  -= temp; 
      Vzz1[ag][2]  -= temp2; 
      Vzz1[ag][3]  += temp+temp2;   } 

  /// LATENT SLOW TO SAFE ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*rRecov;
      Vzz1[ag][1]  -= temp;  } 

  /// AGING ///
    for(int ag=0; ag<10; ag++) { for(int tb=0; tb<4; tb++) {
      if(ag>0) { temp2 = 120; } else { temp2 = 60; }
      temp = Vzz0[ag][tb]/temp2;
      Vzz1[ag  ][tb]  -= temp; 
      Vzz1[ag+1][tb]  += temp;  } }
  
  /// add up results ///
  for(int ag=0; ag<11; ag++) {
    CumInc(ss,4) +=  Vzz1[ag][3]; }

  /// reset vectors ///
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb]; } }

  } // end month loop

/// 25-35 adult LTBI USB   ############################################################################
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb] = 0; } }
  Vzz0[3 ][1] = Vzz1[3 ][1] = 1;  // assign starting position

  for(int ss=0; ss<480; ss++) {  // start month loop

///  INFECTION /// StatList
       if(ss==0) { for(int ag=0; ag<11; ag++) { 
        temp          = Vzz0[ag][1];
        Vzz1[ag][1]  -= temp; 
        Vzz1[ag][1]  += temp*MpslowPIN[ag][0]; 
        Vzz1[ag][2]  += temp*MpfastPIN[ag][0];
        Vzz1[ag][3]  += temp*(MpimmedPINn[ag][0]+MpimmedPINp[ag][0]); }  }

/// BREAKDOWN ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*MrslowN[ag][0]*rrSlowFB[0];
      temp2 = Vzz0[ag][2]*rfast;
      Vzz1[ag][1]  -= temp; 
      Vzz1[ag][2]  -= temp2; 
      Vzz1[ag][3]  += temp+temp2;   } 

  /// LATENT SLOW TO SAFE ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*rRecov;
      Vzz1[ag][1]  -= temp;  } 

  /// AGING ///
    for(int ag=0; ag<10; ag++) { for(int tb=0; tb<4; tb++) {
      if(ag>0) { temp2 = 120; } else { temp2 = 60; }
      temp = Vzz0[ag][tb]/temp2;
      Vzz1[ag  ][tb]  -= temp; 
      Vzz1[ag+1][tb]  += temp;  } }
  
  /// add up results ///
  for(int ag=0; ag<11; ag++) {
    CumInc(ss,5) +=  Vzz1[ag][3]; }

  /// reset vectors ///
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb]; } }

  } // end month loop

/// 25-35 adult uninfected USB ART2  ############################################################################
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb] = 0; } }
  Vzz0[3 ][0] = Vzz1[3 ][0] = 1;  // assign starting position

  for(int ss=0; ss<480; ss++) {  // start month loop

///  INFECTION /// StatList
       if(ss==0) { for(int ag=0; ag<11; ag++) { 
        temp          = Vzz0[ag][0];
        Vzz1[ag][0]  -= temp; 
        Vzz1[ag][1]  += temp*MpslowN[ag][4]; 
        Vzz1[ag][2]  += temp*MpfastN[ag][4];
        Vzz1[ag][3]  += temp*(MpimmedNn[ag][4]+MpimmedNp[ag][4]); }  }

/// BREAKDOWN ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*MrslowN[ag][4]*rrSlowFB[0];
      temp2 = Vzz0[ag][2]*rfast;
      Vzz1[ag][1]  -= temp; 
      Vzz1[ag][2]  -= temp2; 
      Vzz1[ag][3]  += temp+temp2;   } 

  /// LATENT SLOW TO SAFE ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*rRecov;
      Vzz1[ag][1]  -= temp;  } 

  /// AGING ///
    for(int ag=0; ag<10; ag++) { for(int tb=0; tb<4; tb++) {
      if(ag>0) { temp2 = 120; } else { temp2 = 60; }
      temp = Vzz0[ag][tb]/temp2;
      Vzz1[ag  ][tb]  -= temp; 
      Vzz1[ag+1][tb]  += temp;  } }
  
  /// add up results ///
  for(int ag=0; ag<11; ag++) {
    CumInc(ss,6) +=  Vzz1[ag][3]; }

  /// reset vectors ///
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb]; } }

  } // end month loop

/// 25-35 adult uninfected USB ADVANCED HIV  ############################################################################
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb] = 0; } }
  Vzz0[3 ][0] = Vzz1[3 ][0] = 1;  // assign starting position

  for(int ss=0; ss<480; ss++) {  // start month loop

///  INFECTION /// StatList
       if(ss==0) { for(int ag=0; ag<11; ag++) { 
        temp          = Vzz0[ag][0];
        Vzz1[ag][0]  -= temp; 
        Vzz1[ag][1]  += temp*MpslowN[ag][3]; 
        Vzz1[ag][2]  += temp*MpfastN[ag][3];
        Vzz1[ag][3]  += temp*(MpimmedNn[ag][3]+MpimmedNp[ag][3]); }  }

/// BREAKDOWN ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*MrslowN[ag][3]*rrSlowFB[0];
      temp2 = Vzz0[ag][2]*rfast;
      Vzz1[ag][1]  -= temp; 
      Vzz1[ag][2]  -= temp2; 
      Vzz1[ag][3]  += temp+temp2;   } 

  /// LATENT SLOW TO SAFE ///  
  for(int ag=0; ag<11; ag++) {
      temp  = Vzz0[ag][1]*rRecov;
      Vzz1[ag][1]  -= temp;  } 

  /// AGING ///
    for(int ag=0; ag<10; ag++) { for(int tb=0; tb<4; tb++) {
      if(ag>0) { temp2 = 120; } else { temp2 = 60; }
      temp = Vzz0[ag][tb]/temp2;
      Vzz1[ag  ][tb]  -= temp; 
      Vzz1[ag+1][tb]  += temp;  } }
  
  /// add up results ///
  for(int ag=0; ag<11; ag++) {
    CumInc(ss,7) +=  Vzz1[ag][3]; }

  /// reset vectors ///
  for(int ag=0; ag<11; ag++) { for(int tb=0; tb<4; tb++) { Vzz0[ag][tb] = Vzz1[ag][tb]; } }

  } // end month loop

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 


    
   NumericVector        CheckV(36300); 
    for(int ag=0; ag<11; ag++) { for(int tb=0; tb<11; tb++) { for(int dr=0; dr<5; dr++) { 
    for(int tx=0; tx<3 ; tx++) { for(int hv=0; hv<5 ; hv++) { for(int rg=0; rg<4; rg++) {
     CheckV(ag+tb*11+dr*121+tx*605+hv*1815+rg*9075) = V1[ag][tb][dr][tx][hv][rg];
      } } } } } }  

   //// RETURN STUFF
     for(int i=0; i<nYrs; i++) {  for(int j=0; j<nRes; j++) {  Outputs2(i,j)  = Outputs[i][j];  }  }

    return Rcpp::List::create(
      Rcpp::Named("Outputs") = Outputs2,
      Rcpp::Named("V1") = CheckV,
      Rcpp::Named("V0") = CheckV0,
      Rcpp::Named("CumInc") = CumInc
            );      }')


