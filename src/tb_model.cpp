using namespace Rcpp;

//[[Rcpp::export]]  

List cSim(
    int                 nYrs,     //number of years to run the model
    int                 nRes,     //number of results
    NumericMatrix       rDxt,     /rate of tb diagnosis over time
    std::vector<double> TxQualt,  //tb treatment quality
    NumericMatrix       InitPop,  //initial population
    NumericMatrix       Mpfast,   //matrix of probabilities of fast tb progression
    NumericMatrix       ExogInf,  //
    NumericMatrix       MpfastPI, //matrix of probabilities of fast tb progression in partially immune
    double              pimmed,   //probability of immediate progression to active tb upon infection
    NumericMatrix       MpSmPos,  //matrix of probabilities of smear positive tb
    NumericMatrix       Mrslow,   //matrix of rates of slow tb progression to active disease
    std::vector<double> rrSlowFB, //vector of rate ratios for modification of rate of slow progression among non-US population
    double              rfast,    //rate of fast progression 
    double              RRcurDef, //rate ratio for modification of rate of cure among those who default tb treatment
    std::vector<double> VrSlfCur, //vector of the rates of self-cure of tb 
    std::vector<double> VrSmConv, //vector of the rates of smear conversion
    NumericMatrix       vTMort,   //vector of tb mortality
    NumericMatrix       vNmMort,  //vector of non-TB mortality
    double              muTbNm    //factor for comorbidity btw TB and non-TB
    std::vector<double> Birthst,  //vector of number of US births by month
    NumericMatrix       ImmNon,   //
    NumericMatrix       ImmLat,   //
    NumericMatrix       ImmAct,   // 
    NumericMatrix       ImmFst,   // 
    std::vector<double> TxExpAge, //
    double              p_Imm_SP, //probability of immigrants smear positivity (might be unneeded)
    NumericMatrix       mubt,     //matrix of background mortality
    NumericMatrix       RelInf,   //relative infectiousness
    std::vector<double> RelInfRg, //relative infectiousness by risk group (needs modification) 
    std::vector<double> Vmix,     //vector of mixing parameters
    NumericMatrix       vIsxtoIsy, //matrix for transitions within the Is dimension
    NumericMatrix       vNmxtoNmy, //matrix for transitions within the Nm dimension
    NumericMatrix       vLcxtoLcy, //matrix for transitions within the Lc dimension   
    double              rIntvDef,  //rate of default of the intervention for RF of interest 
    NumericMatrix       rRFt,      //rate of risk factor (population) of interest over time
    std::vector<double> rEmmigFB,  //rate of emmigration among non-US born population
    NumericMatrix       rIntvInit, //rate of intervention inititation for RF of interest over time
    NumericMatrix       TxMat,     //Matrix of various tb treatment characteristics (likely to become vector/double)
    double              TunTxMort, //tuning parameter for mortality when treatment is considered
    std::vector<double> rDeft,     //rate of tb treatment default over time
    std::vector<double> rDeftH,    //rate of tb treatment default 
    std::vector<double> LtTxPar,   //latent tb treatment 
    NumericMatrix       LtDxPar,   //latent diagnosis
    std::vector<double> RelInfHivt,//relative infectiousness by HIV over time (needs modification)
    std::vector<double> RRdxAge,   //rate ratio of tb diagnosis by age group
    double              rRecov,    //rate of tb recovery
    double              pImmScen,  //
    std::vector<double> EarlyTrend, 
    NumericMatrix       EffLtX,    
    double              EffLt,
    std::vector<double> dLtt,
    std::vector<double> NixTrans
    ) {  
 
