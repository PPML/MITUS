//
//  mixing.cpp
//  
//
//  Created by Nicole Swartwood on 2/20/18.
//

#include <Rcpp.h>

using namespace Rcpp;

//[[Rcpp::export]]


double        VNkl[2][2];  ///HIGH AND LOW RISK, NATIVITY
double        VGjkl[2][2]; ///HIGH AND LOW RISK, NATIVITY
double        Vjaf[4];     ///BY NUMBER OF MIXING GROUPS
double        VLjk[2][2];  ///HIGH AND LOW RISK, NATIVITY

////////INITIALIZE
for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
            VNkl[i][j] = 0;
} }
for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
        VGjkl[i][j] = 0;
        VLjkl[i][j] = 0;
    } } }
for(int i=0; i<4; i++) { ///effective contact rates
        Vjaf[i]= 0;
    } }


///BURN IN

for(int i=0; i<2; i++) {
    for(int j=0; j<2; j++) {
        VNkl [i][j] = 0;
        VGjkl[i][j] = 0;   // set to zero
    }
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
        VGjkl[0][0] += V0[ag][4][0][0][0][0][0]*RelInfN[tb][0];
        ///////// HIGH RISK US BORN
        VGjkl[1][0]  += V0[ag][4][0][0][0][1][1]*RelInfN[tb][0];
        ///////// LOW RISK NON US BORN
        VGjkl[0][1]  += (V0[ag][4][0][0][0][0][1] + V0[ag][4][0][0][0][0][2])*RelInfN[tb][0];
        /////////  HIGH RISK NON US BORN
        VGjkl[1][1]  += (V0[ag][4][0][0][0][1][1] + V0[ag][4][0][0][0][1][2])*RelInfN[tb][0];
        
    }
    // Step 2 (treated TB)
    // No contribution to force of infection
    
    // Step 3
    Vjaf[0]  = (RelInf[0]*VGjkl[0][0]         +       //LOW RISK US BORN
                RelInf[1]*VGjkl[1][0]*Vmix[0] +       //HIGH RISK US BORN
                RelInf[2]*VGjkl[0][1]*Vmix[1])+       //LOW RISK FOREIGN BORN
                RelInf[3]*VGjkl[1][1]*Vmix[0]*Vmix[1] //HIGH RISK FOREIGN BORN
            /  (RelInf[0]*VNkl[0][0]                    +
                RelInf[1]*VNkl[1][0]*Vmix[0]            +
                RelInf[2]*VNkl[0][1]*Vmix[1]            +
                RelInf[2]*VNkl[0][1]*Vmix[0]*Vmix[1]    +
                1e-12);

    Vjaf[1]  = (Vmix[0]*VGjkl[1][1] + Vmix[1]*VGjkl[1][0]) / (Vmix[0]*VNkl[1][1] +VNkl[1][0]+1e-12);
    
    Vjaf[2]  = (Vmix[1]*VGjkl[1][1] + Vmix[0]*VGjkl[0][1]) / (Vmix[1]*VNkl[1][1] + VNkl[0][1]+1e-12);
    
    Vjaf[3]  = VGjkl[1][1] / (VNkl[1][1]+1e-12);

    // Step 4
    /// LOW RISK US BORN
    VLjkl[0 ][0 ]  = RelInf[0]*Vjaf[0];
    ///////// HIGH RISK US BORN
    VLjkl[1 ][0 ]  = RelInf[1]*((1-Vmix[0])*Vjaf[1]+Vmix[0]*Vjaf[0]);
    ///////// LOW RISK NON US BORN
    VLjkl[0 ][1 ]  = RelInf[2]*((1-Vmix[1])*Vjaf[2]+Vmix[1]*Vjaf[0])+ExogInfN[0];
    ///////// HIGH RISK NON US BORN
    VLjkl[1 ][1 ]  = RelInf[3]*((1-Vmix[0])*(1-Vmix[1])*Vjaf[2]+Vmix[0]*Vmix[1]*Vjaf[0])+ExogInfN[0];

    ///real model
    
    ////////////////////////////  TRANSMISSION RISK  ////////////////////////////////
    for(int i=0; i<2; i++) {
        for(int j=0; j<2; j++) {
            VNkl[i][j] = 0; // set to zero
            VGjkl[i][j] = 0; // set to zero
        } }
////////////////////////////          Step 1         ////////////////////////////
    for(int ag=0; ag<11; ag++) {
        for(int tb=0; tb<6; tb++) {
            for(int lt=0; lt<2; tx++) {
                for(int im=0; im<4; im++){
                    for(int nm=0; nm<4; nm++){
////////////////////////////    LOW RISK, US BORN    ////////////////////////////
                    VNkl[0][0]  += V0[ag][tb][lt][im][nm][0][0];
////////////////////////////   HIGH RISK, US BORN    ////////////////////////////
                    VNkl[1][0]  += V0[ag][tb][lt][im][nm][1][0];
////////////////////////////  LOW RISK, NON US BORN  ////////////////////////////
                    VNkl[2][0]  += V0[ag][tb][lt][im][nm][0][1] + V0[ag][tb][lt][0][nm][0][2];
//////////////////////////// HIGH RISK, NON US BORN  ////////////////////////////
                    VNkl[3][0]  += V0[ag][tb][lt][im][nm][1][1] + V0[ag][tb][lt][0][nm][1][2];
} } } } }
/////////////////////////// /   Step 2  (ACTIVE TB)   ////////////////////////////
for(int ag=0; ag<11; ag++) {
    for(int lt=0; lt<2 ; lt++) {
        for(int im=0; im<4 ; im++) {
            for(int nm=0; nm<4; nm++){
////////////////////////////    LOW RISK, US BORN    ////////////////////////////
                VGjkl[0][0]  +=  V0[ag][4][lt][im][nm][0][0]*RelInfN[tb];
////////////////////////////   HIGH RISK, US BORN    ////////////////////////////
                VGjkl[1][0]  +=  V0[ag][4][lt][im][nm][1][0]*RelInfN[tb];
////////////////////////////  LOW RISK, NON US BORN  ////////////////////////////
                VGjkl[0][1]  += (V0[ag][4][lt][im][nm][0][1] + V0[ag][4][lt][im][nm][0][2])*RelInfN[tb];
//////////////////////////// HIGH RISK, NON US BORN  ////////////////////////////
                VGjkl[1][1]  += (V0[ag][4][lt][im][nm][1][1] + V0[ag][4][lt][im][nm][1][2])*RelInfN[tb];
} } } }
////////////////////////////   Step 2 (TREATED TB)   ////////////////////////////
////////////////////  No contribution to force of infection  ////////////////////
    
////////////////////////////          Step 3         ////////////////////////////
    Vjaf[0]  = (RelInf[0]*VGjkl[0][0]         +       //LOW RISK US BORN
                RelInf[1]*VGjkl[1][0]*Vmix[0] +       //HIGH RISK US BORN
                RelInf[2]*VGjkl[0][1]*Vmix[1])+       //LOW RISK FOREIGN BORN
                RelInf[3]*VGjkl[1][1]*Vmix[0]*Vmix[1] //HIGH RISK FOREIGN BORN
            /  (RelInf[0]*VNkl[0][0]                    +
                RelInf[1]*VNkl[1][0]*Vmix[0]            +
                RelInf[2]*VNkl[0][1]*Vmix[1]            +
                RelInf[2]*VNkl[0][1]*Vmix[0]*Vmix[1]    +
                1e-12);
    
    Vjaf[1]  = (Vmix[0]*VGjkl[1][1] + Vmix[1]*VGjkl[1][0]) / (Vmix[0]*VNkl[1][1] +VNkl[1][0]+1e-12);
    
    Vjaf[2]  = (Vmix[1]*VGjkl[1][1] + Vmix[0]*VGjkl[0][1]) / (Vmix[1]*VNkl[1][1] + VNkl[0][1]+1e-12);
    
    Vjaf[3]  = VGjkl[1][1] / (VNkl[1][1]+1e-12);

    // Step 4
    //CREATE LAMBDA FORCE OF INFECTION
    /// LOW RISK US BORN
    VLjkl[0 ][0 ]  = RelInf[0]*Vjaf[0];
    ///////// HIGH RISK US BORN
    VLjkl[1 ][0 ]  = RelInf[1]*((1-Vmix[0])*Vjaf[1]+Vmix[0]*Vjaf[0]);
    ///////// LOW RISK NON US BORN
    VLjkl[0 ][1 ]  = RelInf[2]*((1-Vmix[1])*Vjaf[2]+Vmix[1]*Vjaf[0])+ExogInfN[0];
    ///////// HIGH RISK NON US BORN
    VLjkl[1 ][1 ]  = RelInf[3]*((1-Vmix[0])*(1-Vmix[1])*Vjaf[2]+Vmix[0]*Vmix[1]*Vjaf[0])+ExogInfN[0];
