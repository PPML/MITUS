// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// national_cSim
Rcpp::List national_cSim(std::vector<int> setup_pars, Rcpp::NumericMatrix rDxt, std::vector<double> TxQualt, Rcpp::NumericMatrix InitPop, Rcpp::NumericMatrix Mpfast, std::vector<double> ExogInf, Rcpp::NumericMatrix MpfastPI, Rcpp::NumericMatrix Mrslow, std::vector<double> rrSlowFB, double rfast, double RRcurDef, double rSlfCur, double p_HR, Rcpp::NumericMatrix vTMort, std::vector<double> RRmuRF, std::vector<double> RRmuHR, std::vector<double> RRmuTBPand, std::vector<double> Birthst, Rcpp::NumericMatrix HrEntEx, Rcpp::NumericMatrix ImmNon, Rcpp::NumericMatrix ImmLat, Rcpp::NumericMatrix ImmAct, Rcpp::NumericMatrix ImmFst, Rcpp::NumericMatrix SpImmNon, std::vector<double> net_mig_usb, std::vector<double> net_mig_nusb, Rcpp::NumericMatrix mubt, std::vector<double> RelInf, std::vector<double> RelInfRg, std::vector<double> RRcrAG, std::vector<double> Vmix, std::vector<double> rEmmigFB, std::vector<double> TxVec, double TunTxMort, std::vector<double> rDeft, Rcpp::NumericMatrix rLtScrt, Rcpp::NumericMatrix ttt_samp_dist, std::vector<int> ttt_month, double ttt_pop_scrn, std::vector<double> ttt_ltbi, double ttt_ltbi_accept, double ttt_ltbi_init, double ttt_ltbi_comp, double ttt_ltbi_eff, std::vector<double> ttt_ltbi_sens, std::vector<double> ttt_ltbi_spec, Rcpp::NumericMatrix LtTxPar, Rcpp::NumericMatrix LtDxPar_lt, Rcpp::NumericMatrix LtDxPar_nolt, double rrTestLrNoTb, double rrTestHr, Rcpp::NumericMatrix Int1Test, Rcpp::NumericMatrix Int1Init, Rcpp::NumericMatrix Int1Tx, std::vector<double> RRdxAge, double rRecov, double pImmScen, std::vector<double> EarlyTrend, std::vector<double> pReTx, Rcpp::NumericMatrix ag_den, std::vector<double> NixTrans, std::vector<double> NixTb, Rcpp::NumericMatrix dist_gen, Rcpp::NumericMatrix trans_mat_tot_ages);
RcppExport SEXP _MITUS_national_cSim(SEXP setup_parsSEXP, SEXP rDxtSEXP, SEXP TxQualtSEXP, SEXP InitPopSEXP, SEXP MpfastSEXP, SEXP ExogInfSEXP, SEXP MpfastPISEXP, SEXP MrslowSEXP, SEXP rrSlowFBSEXP, SEXP rfastSEXP, SEXP RRcurDefSEXP, SEXP rSlfCurSEXP, SEXP p_HRSEXP, SEXP vTMortSEXP, SEXP RRmuRFSEXP, SEXP RRmuHRSEXP, SEXP RRmuTBPandSEXP, SEXP BirthstSEXP, SEXP HrEntExSEXP, SEXP ImmNonSEXP, SEXP ImmLatSEXP, SEXP ImmActSEXP, SEXP ImmFstSEXP, SEXP SpImmNonSEXP, SEXP net_mig_usbSEXP, SEXP net_mig_nusbSEXP, SEXP mubtSEXP, SEXP RelInfSEXP, SEXP RelInfRgSEXP, SEXP RRcrAGSEXP, SEXP VmixSEXP, SEXP rEmmigFBSEXP, SEXP TxVecSEXP, SEXP TunTxMortSEXP, SEXP rDeftSEXP, SEXP rLtScrtSEXP, SEXP ttt_samp_distSEXP, SEXP ttt_monthSEXP, SEXP ttt_pop_scrnSEXP, SEXP ttt_ltbiSEXP, SEXP ttt_ltbi_acceptSEXP, SEXP ttt_ltbi_initSEXP, SEXP ttt_ltbi_compSEXP, SEXP ttt_ltbi_effSEXP, SEXP ttt_ltbi_sensSEXP, SEXP ttt_ltbi_specSEXP, SEXP LtTxParSEXP, SEXP LtDxPar_ltSEXP, SEXP LtDxPar_noltSEXP, SEXP rrTestLrNoTbSEXP, SEXP rrTestHrSEXP, SEXP Int1TestSEXP, SEXP Int1InitSEXP, SEXP Int1TxSEXP, SEXP RRdxAgeSEXP, SEXP rRecovSEXP, SEXP pImmScenSEXP, SEXP EarlyTrendSEXP, SEXP pReTxSEXP, SEXP ag_denSEXP, SEXP NixTransSEXP, SEXP NixTbSEXP, SEXP dist_genSEXP, SEXP trans_mat_tot_agesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int> >::type setup_pars(setup_parsSEXP);
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
    Rcpp::traits::input_parameter< std::vector<double> >::type RRmuTBPand(RRmuTBPandSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type Birthst(BirthstSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type HrEntEx(HrEntExSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ImmNon(ImmNonSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ImmLat(ImmLatSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ImmAct(ImmActSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ImmFst(ImmFstSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type SpImmNon(SpImmNonSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type net_mig_usb(net_mig_usbSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type net_mig_nusb(net_mig_nusbSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type mubt(mubtSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RelInf(RelInfSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RelInfRg(RelInfRgSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RRcrAG(RRcrAGSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type Vmix(VmixSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rEmmigFB(rEmmigFBSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type TxVec(TxVecSEXP);
    Rcpp::traits::input_parameter< double >::type TunTxMort(TunTxMortSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type rDeft(rDeftSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type rLtScrt(rLtScrtSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ttt_samp_dist(ttt_samp_distSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type ttt_month(ttt_monthSEXP);
    Rcpp::traits::input_parameter< double >::type ttt_pop_scrn(ttt_pop_scrnSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type ttt_ltbi(ttt_ltbiSEXP);
    Rcpp::traits::input_parameter< double >::type ttt_ltbi_accept(ttt_ltbi_acceptSEXP);
    Rcpp::traits::input_parameter< double >::type ttt_ltbi_init(ttt_ltbi_initSEXP);
    Rcpp::traits::input_parameter< double >::type ttt_ltbi_comp(ttt_ltbi_compSEXP);
    Rcpp::traits::input_parameter< double >::type ttt_ltbi_eff(ttt_ltbi_effSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type ttt_ltbi_sens(ttt_ltbi_sensSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type ttt_ltbi_spec(ttt_ltbi_specSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type LtTxPar(LtTxParSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type LtDxPar_lt(LtDxPar_ltSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type LtDxPar_nolt(LtDxPar_noltSEXP);
    Rcpp::traits::input_parameter< double >::type rrTestLrNoTb(rrTestLrNoTbSEXP);
    Rcpp::traits::input_parameter< double >::type rrTestHr(rrTestHrSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Int1Test(Int1TestSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Int1Init(Int1InitSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Int1Tx(Int1TxSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type RRdxAge(RRdxAgeSEXP);
    Rcpp::traits::input_parameter< double >::type rRecov(rRecovSEXP);
    Rcpp::traits::input_parameter< double >::type pImmScen(pImmScenSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type EarlyTrend(EarlyTrendSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type pReTx(pReTxSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ag_den(ag_denSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type NixTrans(NixTransSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type NixTb(NixTbSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type dist_gen(dist_genSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type trans_mat_tot_ages(trans_mat_tot_agesSEXP);
    rcpp_result_gen = Rcpp::wrap(national_cSim(setup_pars, rDxt, TxQualt, InitPop, Mpfast, ExogInf, MpfastPI, Mrslow, rrSlowFB, rfast, RRcurDef, rSlfCur, p_HR, vTMort, RRmuRF, RRmuHR, RRmuTBPand, Birthst, HrEntEx, ImmNon, ImmLat, ImmAct, ImmFst, SpImmNon, net_mig_usb, net_mig_nusb, mubt, RelInf, RelInfRg, RRcrAG, Vmix, rEmmigFB, TxVec, TunTxMort, rDeft, rLtScrt, ttt_samp_dist, ttt_month, ttt_pop_scrn, ttt_ltbi, ttt_ltbi_accept, ttt_ltbi_init, ttt_ltbi_comp, ttt_ltbi_eff, ttt_ltbi_sens, ttt_ltbi_spec, LtTxPar, LtDxPar_lt, LtDxPar_nolt, rrTestLrNoTb, rrTestHr, Int1Test, Int1Init, Int1Tx, RRdxAge, rRecov, pImmScen, EarlyTrend, pReTx, ag_den, NixTrans, NixTb, dist_gen, trans_mat_tot_ages));
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

static const R_CallMethodDef CallEntries[] = {
    {"_MITUS_national_cSim", (DL_FUNC) &_MITUS_national_cSim, 64},
    {"_MITUS_reblncd", (DL_FUNC) &_MITUS_reblncd, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_MITUS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
