redo_load<-function(loc){
model_load(loc)
View(Opt)
sm <<- readRDS(paste0("~/MITUS/inst/", loc, "/sm_resTab_2019-04-23.rds"))
nr<-10
mResTabfb <- mResTabus <- mResTab <- array(NA,dim=c(9,32,7,11))

for(l in 1:9){
  for (i in 1:32){
    for (j in 1:7){
      for (k in 1:11){
        mResTab[l,i,j,k]<-mean(na.omit(sm[["ResTab"]][1:nr+((l-1)*10),i,j,k]))
        mResTabus[l,i,j,k]<-mean(na.omit(sm[["ResTabus"]][1:nr+((l-1)*10),i,j,k]))
        mResTabfb[l,i,j,k]<-mean(na.omit(sm[["ResTabfb"]][1:nr+((l-1)*10),i,j,k]))
      } } }}

smResTabC <- list(mResTab,mResTabus,mResTabfb)
saveRDS(smResTabC, file=paste0("~/MITUS/inst/", loc, "/sm_resTab_", Sys.Date(), ".rds"))


bg <<- readRDS(paste0("~/MITUS/inst/", loc, "/bg_resTab_2019-04-23.rds"))

mResTabfb <- mResTabus <- mResTab <- array(NA,dim=c(9,32,7,4))

for(l in 1:9){
  for (i in 1:32){
    for (j in 1:7){
      for (k in 1:11){
        mResTab[l,i,j,k]<-mean(na.omit(bg[["ResTab"]][1:nr+((l-1)*10),i,j,k]))
        mResTabus[l,i,j,k]<-mean(na.omit(bg[["ResTabus"]][1:nr+((l-1)*10),i,j,k]))
        mResTabfb[l,i,j,k]<-mean(na.omit(bg[["ResTabfb"]][1:nr+((l-1)*10),i,j,k]))
      } } }}

bgResTabC <- list(mResTab,mResTabus,mResTabfb)

saveRDS(bgResTabC, file=paste0("~/MITUS/inst/", loc, "/bgm_resTab_", Sys.Date(), ".rds"))


}
