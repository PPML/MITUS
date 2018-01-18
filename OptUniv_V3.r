
## Load Par

 UnivOptim <- function(parr) {
   print(posterior(parr))
   parr2 <- parr
   for(fi in 1:ncol(StartVal)) {
     posterior1 <- function(gg) { parr2[fi] <- gg; posterior(parr2) }
    res <- optimize(posterior1,parr[fi]*c(0.1,10),tol = .Machine$double.eps^0.25/10)
    if(res$objective < posterior(parr2) ) {  parr2[fi] <- res$minimum }
     print(paste(fi,"of",ncol(StartVal),"complete, objective = ",posterior(parr2)))   }
     return(list(par=parr2,value=res$objective)) }
  


