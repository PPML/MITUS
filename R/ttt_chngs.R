#'This function applies
#'@name ttt_chngs
#'@param
#'@param
#'@return
ttt_chngs <- function(UI_ttt_chngs,P) {
  if (ttt_chng==TRUE){
    #which nativity
    #which age
    if (UI_ttt_chngs$age==1){
    } else if (UI_ttt_chngs$age==2){
    } else if (UI_ttt_chngs$age==3){
    } else if (UI_ttt_chngs$age==4){
    }
    #how many people
    for (ag in 1:11){
      for (na in 1:3){
    tot_pop<-sum()
      }}
    pop_frc<-UI_ttt_chngs$scrn_pop/tot_pop
    #fraction screened
    #intervention years
    #change in LTBI prevalence
    #change in mortality
    #change in LTBI progression
    #will need to construct a new sub distribution
    #example for HIV
    pars1 = c(0.5,1,0.8); #  parameters for HIV dist
    dist1 <- matrix(NA,4,4)
    for(i in 1:4) {
      for(j in 1:4) {
        dist1[i,j] <- pmvnorm(lower = cuts[c(i,j)],
                              upper = cuts[c(i,j)+1],
                              mean  = pars1[c(1,1)],
                              sigma = matrix(pars1[2]*c(1,pars1[3],pars1[3],1),2,2) )[[1]]
      }
    }
    dist1 <- dist1/sum(dist1)
  }
}
