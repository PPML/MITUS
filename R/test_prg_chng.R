# defprg<-def_prgchng(Par[10,])
# print(defprg)
#
# #double screening coverage
# doublescreen<-defprg
# doublescreen[2]<-doublescreen[2]*2
# identical(fin_param(P,doublescreen),fin_param_init(P,prg_chng=doublescreen))
#
# #raise ltbi frc to .9
# ltbi_frc<-defprg
# ltbi_frc[4]<-.9
# identical(fin_param(P,ltbi_frc),fin_param_init(P,prg_chng=ltbi_frc))
#
# #raise ltbi complete frc to .9
# ltbi_cfrc<-defprg
# ltbi_cfrc[5]<-.9
# identical(fin_param(P,ltbi_cfrc),fin_param_init(P,prg_chng=ltbi_cfrc))
#
# #raise ltbi effec frc to .95
# ltbi_efrc<-defprg
# ltbi_efrc[6]<-.95
# identical(fin_param(P,ltbi_efrc),fin_param_init(P,prg_chng=ltbi_efrc))
#
# #reduce tb_tim2tx_frc to .8
# timtx<-defprg
# timtx[7]<-80
# identical(fin_param(P,timtx),fin_param_init(P,prg_chng=timtx))
#
# #tb_txdef_frc is doubled
# txdef<-defprg
# txdef[8]<-defprg[8]*2
# identical(fin_param(P,txdef),fin_param_init(P,prg_chng=txdef))
#
