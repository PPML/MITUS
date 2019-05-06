# yo4<-yo3<-yo2<-yo<-matrix(NA,11,4)
#
# for (ag in 1:ncol(prms$mubt)){
#   for (nm in 1:length(prms$RRmuRF)){
#     yo[ag,nm]<-(1-(prms$mubt[1,ag]*prms$RRmuRF[nm]))
#     yo2[ag,nm]<-mort_v[nm]*(1-(prms$mubt[1,ag]*prms$RRmuRF[nm]))
#     yo3[ag,nm]<-yo[ag,nm]
#
#   }}
# for (ag in 1:ncol(prms$mubt)){
#   for (nm in 1:length(prms$RRmuRF)){
#
#     yo3[ag,nm]<-yo[ag,nm]/rowSums(yo)[ag]
#     yo4[ag,nm]<-yo2[ag,nm]/rowSums(yo2)[ag]
#
#   }}
