// for (int intv=0; intv<ttt_sampling_dist.nrow(); intv++){
// /// LTBI SCREENING AND TLTBI INITIATION /// only for no previous TB or LTBI tx
// int agi = ttt_ag.size(); int nai = ttt_na.size(); int si = ttt_month.size();
// // if (s==1){
// //   Rcpp::Rcout<<agi<<" "<<nai<<" "<<si<<"\n";
// // }
// for(int rg=0; rg<2; rg++) {
//   for(int na=0; na<3; na++) {
//     for(int im=0; im<4; im++) {
//       for(int nm=0; nm<4; nm++) {
//         for(int ag=0; ag<11; ag++) {
//
//           rr_ltbi=1; pop_frc=0;
//           for(int i=0; i<4; i++) {
//             for(int j=0; j<4; j++) {
//               ttt_dist[i][j] =0;
//             } }
//           if(std::find(std::begin(ttt_month), std::end(ttt_month), s) != std::end(ttt_month)){
//             for (int i=0; i<agi; i++){
//               for (int j=0; j<nai; j++){
//                 if (ag==ttt_ag[i] & na==ttt_na[j]){
//                   rr_ltbi=ttt_ltbi;
//                   ttt_dist[nm][im]=ttt_samp_distN[nm][im];
//                 } } } }
//           ////////////// US BORN, LOW RISK  //////////////////
//           if(rg==0 & na==0) {
//             rTbP = rLtScrt[s]*LtDxPar_ltN[0][s];
//             rTbP_norm=(ttt_dist[nm][im]*LtDxPar_ltN[0][s]*rr_ltbi)/(((rLtScrt[s]+ttt_dist[nm][im])*LtDxPar_ltN[0][s]*rr_ltbi)+(1-(rLtScrt[s]+ttt_dist[nm][im])*LtDxPar_ltN[0][s]));
//             rTbN = rLtScrt[s]*LtDxPar_noltN[0][s];
//           }
//           ////////////// US BORN, HIGH RISK  /////////////////
//           if(rg==1 & na==0) {
//             rTbP = rLtScrt[s]*LtDxPar_ltN[1][s];
//             rTbP_norm=(ttt_dist[nm][im]*LtDxPar_ltN[1][s]*rr_ltbi)/(((rLtScrt[s]+ttt_dist[nm][im])*LtDxPar_ltN[1][s]*rr_ltbi)+(1-(rLtScrt[s]+ttt_dist[nm][im])*LtDxPar_ltN[1][s]));
//             rTbN = rLtScrt[s]*LtDxPar_noltN[1][s];
//           }
//           ////////////// Young NUS (under 5)  /////////////////
//           if(rg==0 & na > 0 & ag==0) {
//             rTbP = rLtScrt[s]*LtDxPar_ltN[2][s];
//             rTbP_norm=(ttt_dist[nm][im]*LtDxPar_ltN[2][s]*rr_ltbi)/(((rLtScrt[s]+ttt_dist[nm][im])*LtDxPar_ltN[2][s]*rr_ltbi)+(1-(rLtScrt[s]+ttt_dist[nm][im])*LtDxPar_ltN[2][s]));
//             rTbN = rLtScrt[s]*LtDxPar_noltN[2][s];
//           }
//
//           //////////// NON US BORN  ////////////////
//           if(rg==0 & na > 0 & ag > 0) {
//             rTbP = rLtScrt[s]*LtDxPar_ltN[3][s];
//             rTbP_norm=(ttt_dist[nm][im]*LtDxPar_ltN[3][s]*rr_ltbi)/(((rLtScrt[s]+ttt_dist[nm][im])*LtDxPar_ltN[3][s]*rr_ltbi)+(1-(rLtScrt[s]+ttt_dist[nm][im])*LtDxPar_ltN[3][s]));
//             rTbN = rLtScrt[s]*LtDxPar_noltN[3][s];
//           }
//           ////////////// NON US BORN, HIGH RISK  /////////////////
//           if(rg==1 & na >0) {
//             rTbP = rLtScrt[s]*LtDxPar_ltN[4][s];
//             rTbP_norm=(ttt_dist[nm][im]*LtDxPar_ltN[4][s]*rr_ltbi)/(((rLtScrt[s]+ttt_dist[nm][im])*LtDxPar_ltN[4][s]*rr_ltbi)+(1-(rLtScrt[s]+ttt_dist[nm][im])*LtDxPar_ltN[4][s]));
//             rTbN = rLtScrt[s]*LtDxPar_noltN[4][s];
//           }
//
//           // if (s==858){
//           //   Rcpp::Rcout << "rTbPnorm = " << rTbP << " @ ag = "<< ag << " & na = "<< na << " & rg = " << rg <<" \n";
//           //   if ((V0[ag][2][0][im][nm][rg][na]+V0[ag][3][0][im][nm][rg][na])*ttt_dist[nm][im]*rTbP !=0){
//           //     Rcpp::Rcout << "mod additional = " << (V0[ag][2][0][im][nm][rg][na]+V0[ag][3][0][im][nm][rg][na])*ttt_dist[nm][im]*rTbP << " @ ag = "<< ag << " & na = "<< na << " & rg = " << rg <<" \n";
//           //   }}
//           ////////////// True Status -- LTBI Negative
//           temp5= V0[ag][0][0][im][nm][rg][na]*rTbN;
//           temp6= V0[ag][1][0][im][nm][rg][na]*rTbN;
//           ///remove from the TB naive and PI states
//           V1[ag][0][0][im][nm][rg][na]  -= temp5;
//           V1[ag][1][0][im][nm][rg][na]  -= temp6;
//           ///////moving to latent tx experienced as in last model -- is this correct?
//           V1[ag][0][1][im][nm][rg][na]  += temp5;
//           V1[ag][1][1][im][nm][rg][na]  += temp6;
//           ////True Status -- LTBI Positive
//           base_diag=V0[ag][2][0][im][nm][rg][na]*rTbP;
//           temp  =(base_diag + (V0[ag][2][0][im][nm][rg][na]*rTbP_norm))*LtTxParN[s][0]*(1-LtTxParN[s][1]);// tx completion
//           temp3 =(base_diag + (V0[ag][2][0][im][nm][rg][na]*rTbP_norm))*LtTxParN[s][0]*LtTxParN[s][1]; // default
//
//           base_diag=V0[ag][3][0][im][nm][rg][na]*rTbP;
//           temp2  = (base_diag + (V0[ag][3][0][im][nm][rg][na]*rTbP_norm))*LtTxParN[s][0]*(1-LtTxParN[s][1]);// tx completion
//           temp4  = (base_diag + (V0[ag][3][0][im][nm][rg][na]*rTbP_norm))*LtTxParN[s][0]*LtTxParN[s][1]; // default
//
//           V1[ag][2][0][im][nm][rg][na]  -=  (temp+temp3); //remove from latent slow
//           V1[ag][3][0][im][nm][rg][na]  -=  (temp2+temp4);  //remove from latent fast
//           //completion split between success and failure
//           V1[ag][1][1][im][nm][rg][na]  += (temp+temp2)*LtTxParN[s][2]; //*EffLt0[s]; //exit to cure
//           V1[ag][2][1][im][nm][rg][na]  += (temp+temp2)*(1-LtTxParN[s][2]); //*(1-EffLt0[s]) //tx comp fail to latent slow
//           ///defaults are placed in tx naive because it is considered the same tb infection
//           V1[ag][2][0][im][nm][rg][na]  += (temp3+temp4); //latent tx default to latent slow
//         } } } } }
// } //end intv
