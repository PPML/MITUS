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
                Outputs[y][110+ag]  += V0[ag][tb][lt][im][nm][rg][na]*vRFMort[nm];
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
