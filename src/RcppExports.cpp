// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cSim_noTB
Rcpp::List cSim_noTB(int nYrs, int nRes, Rcpp::NumericMatrix rDxt, std::vector<double> TxQualt, Rcpp::NumericMatrix InitPop, Rcpp::NumericMatrix Mpfast, std::vector<double> ExogInf, Rcpp::NumericMatrix MpfastPI, Rcpp::NumericMatrix Mrslow, std::vector<double> rrSlowFB, double rfast, double RRcurDef, double rSlfCur, double p_HR, Rcpp::NumericMatrix vTMort, std::vector<double> RRmuRF, std::vector<double> RRmuHR, double muTbRF, std::vector<double> Birthst, Rcpp::NumericMatrix HrEntEx, Rcpp::NumericMatrix ImmNon, Rcpp::NumericMatrix ImmLat, Rcpp::NumericMatrix ImmAct, Rcpp::NumericMatrix ImmFst, Rcpp::NumericMatrix mubt, std::vector<double> RelInf, std::vector<double> RelInfRg, std::vector<double> Vmix, std::vector<double> rEmmigFB, std::vector<double> TxVec, double TunTxMort, std::vector<double> rDeft, std::vector<double> rLtScrt, std::vector<double> LtTxPar, Rcpp::NumericMatrix LtDxPar, std::vector<double> RRdxAge, double rRecov, double pImmScen, std::vector<double> EarlyTrend, std::vector<double> pReTx, std::vector<double> NixTrans, Rcpp::NumericMatrix dist_gen, Rcpp::NumericMatrix trans_mat_tot_ages);
RcppExport SEXP _MITUS_cSim_noTB(SEXP nYrsSEXP, SEXP nResSEXP, SEXP rDxtSEXP, SEXP TxQualtSEXP, SEXP InitPopSEXP, SEXP MpfastSEXP, SEXP ExogInfSEXP, SEXP MpfastPISEXP, SEXP MrslowSEXP, SEXP rrSlowFBSEXP, SEXP rfastSEXP, SEXP RRcurDefSEXP, SEXP rSlfCurSEXP, SEXP p_HRSEXP, SEXP vTMortSEXP, SEXP RRmuRFSEXP, SEXP RRmuHRSEXP, SEXP muTbRFSEXP, SEXP BirthstSEXP, SEXP HrEntExSEXP, SEXP ImmNonSEXP, SEXP ImmLatSEXP, SEXP ImmActSEXP, SEXP ImmFstSEXP, SEXP mubtSEXP, SEXP RelInfSEXP, SEXP RelInfRgSEXP, SEXP VmixSEXP, SEXP rEmmigFBSEXP, SEXP TxVecSEXP, SEXP TunTxMortSEXP, SEXP rDeftSEXP, SEXP rLtScrtSEXP, SEXP LtTxParSEXP, SEXP LtDxParSEXP, SEXP RRdxAgeSEXP, SEXP rRecovSEXP, SEXP pImmScenSEXP, SEXP EarlyTrendSEXP, SEXP pReTxSEXP, SEXP NixTransSEXP, SEXP dist_genSEXP, SEXP trans_mat_tot_agesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nYrs(nYrsSEXP);
    Rcpp::traits::input_parameter< int >::type nRes(nResSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type rDxt(rDxtSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type TxQualt(TxQualtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type InitPop(InitPopSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Mpfast(MpfastSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type ExogInf(ExogInfSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type MpfastPI(MpfastPISEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Mrslow(MrslowSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rrSlowFB(rrSlowFBSEXP);
    Rcpp::traits::input_parameter< double >::type rfast(rfastSEXP);
    Rcpp::traits::input_parameter< double >::type RRcurDef(RRcurDefSEXP);
    Rcpp::traits::input_parameter< double >::type rSlfCur(rSlfCurSEXP);
    Rcpp::traits::input_parameter< double >::type p_HR(p_HRSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type vTMort(vTMortSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RRmuRF(RRmuRFSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RRmuHR(RRmuHRSEXP);
    Rcpp::traits::input_parameter< double >::type muTbRF(muTbRFSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type Birthst(BirthstSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type HrEntEx(HrEntExSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ImmNon(ImmNonSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ImmLat(ImmLatSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ImmAct(ImmActSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ImmFst(ImmFstSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mubt(mubtSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RelInf(RelInfSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RelInfRg(RelInfRgSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type Vmix(VmixSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rEmmigFB(rEmmigFBSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type TxVec(TxVecSEXP);
    Rcpp::traits::input_parameter< double >::type TunTxMort(TunTxMortSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rDeft(rDeftSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rLtScrt(rLtScrtSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type LtTxPar(LtTxParSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type LtDxPar(LtDxParSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RRdxAge(RRdxAgeSEXP);
    Rcpp::traits::input_parameter< double >::type rRecov(rRecovSEXP);
    Rcpp::traits::input_parameter< double >::type pImmScen(pImmScenSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type EarlyTrend(EarlyTrendSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type pReTx(pReTxSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type NixTrans(NixTransSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dist_gen(dist_genSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type trans_mat_tot_ages(trans_mat_tot_agesSEXP);
    rcpp_result_gen = Rcpp::wrap(cSim_noTB(nYrs, nRes, rDxt, TxQualt, InitPop, Mpfast, ExogInf, MpfastPI, Mrslow, rrSlowFB, rfast, RRcurDef, rSlfCur, p_HR, vTMort, RRmuRF, RRmuHR, muTbRF, Birthst, HrEntEx, ImmNon, ImmLat, ImmAct, ImmFst, mubt, RelInf, RelInfRg, Vmix, rEmmigFB, TxVec, TunTxMort, rDeft, rLtScrt, LtTxPar, LtDxPar, RRdxAge, rRecov, pImmScen, EarlyTrend, pReTx, NixTrans, dist_gen, trans_mat_tot_ages));
    return rcpp_result_gen;
END_RCPP
}
// reblncd
Rcpp::NumericMatrix reblncd(Rcpp::NumericMatrix mubt, Rcpp::NumericMatrix can_go, double RRmuHR, std::vector<double> RRmuRF, std::vector<double> HRdist, std::vector<double> dist_gen_v, std::vector<double> adj_fact);
RcppExport SEXP _MITUS_reblncd(SEXP mubtSEXP, SEXP can_goSEXP, SEXP RRmuHRSEXP, SEXP RRmuRFSEXP, SEXP HRdistSEXP, SEXP dist_gen_vSEXP, SEXP adj_factSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mubt(mubtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type can_go(can_goSEXP);
    Rcpp::traits::input_parameter< double >::type RRmuHR(RRmuHRSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RRmuRF(RRmuRFSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type HRdist(HRdistSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type dist_gen_v(dist_gen_vSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type adj_fact(adj_factSEXP);
    rcpp_result_gen = Rcpp::wrap(reblncd(mubt, can_go, RRmuHR, RRmuRF, HRdist, dist_gen_v, adj_fact));
    return rcpp_result_gen;
END_RCPP
}
// cSim
Rcpp::List cSim(int nYrs, int nRes, Rcpp::NumericMatrix rDxt, std::vector<double> TxQualt, Rcpp::NumericMatrix InitPop, Rcpp::NumericMatrix Mpfast, std::vector<double> ExogInf, Rcpp::NumericMatrix MpfastPI, Rcpp::NumericMatrix Mrslow, std::vector<double> rrSlowFB, double rfast, double RRcurDef, double rSlfCur, double p_HR, Rcpp::NumericMatrix vTMort, std::vector<double> RRmuRF, std::vector<double> RRmuHR, double muTbRF, std::vector<double> Birthst, Rcpp::NumericMatrix HrEntEx, Rcpp::NumericMatrix ImmNon, Rcpp::NumericMatrix ImmLat, Rcpp::NumericMatrix ImmAct, Rcpp::NumericMatrix ImmFst, Rcpp::NumericMatrix mubt, std::vector<double> RelInf, std::vector<double> RelInfRg, std::vector<double> Vmix, std::vector<double> rEmmigFB, std::vector<double> TxVec, double TunTxMort, std::vector<double> rDeft, std::vector<double> rLtScrt, std::vector<double> LtTxPar, Rcpp::NumericMatrix LtDxPar, std::vector<double> RRdxAge, double rRecov, double pImmScen, std::vector<double> EarlyTrend, std::vector<double> pReTx, std::vector<double> NixTrans, Rcpp::NumericMatrix dist_gen, Rcpp::NumericMatrix trans_mat_tot_ages);
RcppExport SEXP _MITUS_cSim(SEXP nYrsSEXP, SEXP nResSEXP, SEXP rDxtSEXP, SEXP TxQualtSEXP, SEXP InitPopSEXP, SEXP MpfastSEXP, SEXP ExogInfSEXP, SEXP MpfastPISEXP, SEXP MrslowSEXP, SEXP rrSlowFBSEXP, SEXP rfastSEXP, SEXP RRcurDefSEXP, SEXP rSlfCurSEXP, SEXP p_HRSEXP, SEXP vTMortSEXP, SEXP RRmuRFSEXP, SEXP RRmuHRSEXP, SEXP muTbRFSEXP, SEXP BirthstSEXP, SEXP HrEntExSEXP, SEXP ImmNonSEXP, SEXP ImmLatSEXP, SEXP ImmActSEXP, SEXP ImmFstSEXP, SEXP mubtSEXP, SEXP RelInfSEXP, SEXP RelInfRgSEXP, SEXP VmixSEXP, SEXP rEmmigFBSEXP, SEXP TxVecSEXP, SEXP TunTxMortSEXP, SEXP rDeftSEXP, SEXP rLtScrtSEXP, SEXP LtTxParSEXP, SEXP LtDxParSEXP, SEXP RRdxAgeSEXP, SEXP rRecovSEXP, SEXP pImmScenSEXP, SEXP EarlyTrendSEXP, SEXP pReTxSEXP, SEXP NixTransSEXP, SEXP dist_genSEXP, SEXP trans_mat_tot_agesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nYrs(nYrsSEXP);
    Rcpp::traits::input_parameter< int >::type nRes(nResSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type rDxt(rDxtSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type TxQualt(TxQualtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type InitPop(InitPopSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Mpfast(MpfastSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type ExogInf(ExogInfSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type MpfastPI(MpfastPISEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Mrslow(MrslowSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rrSlowFB(rrSlowFBSEXP);
    Rcpp::traits::input_parameter< double >::type rfast(rfastSEXP);
    Rcpp::traits::input_parameter< double >::type RRcurDef(RRcurDefSEXP);
    Rcpp::traits::input_parameter< double >::type rSlfCur(rSlfCurSEXP);
    Rcpp::traits::input_parameter< double >::type p_HR(p_HRSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type vTMort(vTMortSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RRmuRF(RRmuRFSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RRmuHR(RRmuHRSEXP);
    Rcpp::traits::input_parameter< double >::type muTbRF(muTbRFSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type Birthst(BirthstSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type HrEntEx(HrEntExSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ImmNon(ImmNonSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ImmLat(ImmLatSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ImmAct(ImmActSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ImmFst(ImmFstSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mubt(mubtSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RelInf(RelInfSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RelInfRg(RelInfRgSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type Vmix(VmixSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rEmmigFB(rEmmigFBSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type TxVec(TxVecSEXP);
    Rcpp::traits::input_parameter< double >::type TunTxMort(TunTxMortSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rDeft(rDeftSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rLtScrt(rLtScrtSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type LtTxPar(LtTxParSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type LtDxPar(LtDxParSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RRdxAge(RRdxAgeSEXP);
    Rcpp::traits::input_parameter< double >::type rRecov(rRecovSEXP);
    Rcpp::traits::input_parameter< double >::type pImmScen(pImmScenSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type EarlyTrend(EarlyTrendSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type pReTx(pReTxSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type NixTrans(NixTransSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dist_gen(dist_genSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type trans_mat_tot_ages(trans_mat_tot_agesSEXP);
    rcpp_result_gen = Rcpp::wrap(cSim(nYrs, nRes, rDxt, TxQualt, InitPop, Mpfast, ExogInf, MpfastPI, Mrslow, rrSlowFB, rfast, RRcurDef, rSlfCur, p_HR, vTMort, RRmuRF, RRmuHR, muTbRF, Birthst, HrEntEx, ImmNon, ImmLat, ImmAct, ImmFst, mubt, RelInf, RelInfRg, Vmix, rEmmigFB, TxVec, TunTxMort, rDeft, rLtScrt, LtTxPar, LtDxPar, RRdxAge, rRecov, pImmScen, EarlyTrend, pReTx, NixTrans, dist_gen, trans_mat_tot_ages));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MITUS_cSim_noTB", (DL_FUNC) &_MITUS_cSim_noTB, 43},
    {"_MITUS_reblncd", (DL_FUNC) &_MITUS_reblncd, 7},
    {"_MITUS_cSim", (DL_FUNC) &_MITUS_cSim, 43},
    {NULL, NULL, 0}
};

RcppExport void R_init_MITUS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
