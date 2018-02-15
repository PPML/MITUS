trans_react = matrix(NA, nrow=4, ncol=4)
trans_mort  = matrix(NA, nrow=4, ncol=4)
trans_tot   = matrix(NA, nrow=16, ncol=16)
for (i in 1:4){
  for (j in 1:4){
     for (k in 1:16){
	for (l in 1:16){
   trans_tot[k,l]<-trans_react[i,j]*trans_mort
  }

}}}
trans_tot
