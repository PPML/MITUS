######
library(mnormt)
library(numDeriv)
library(lhs)


#  ASSESS DENSITY OF PROPOSAL DISTRIBUTION
dpropnew <- function(samps,mu,vcv,logg=F) {
  if(dim(as.data.frame(samps))[2]==1) { samps <- t(as.data.frame(samps)) }
  ##
  lpz_y <- dmnorm(samps,mu,vcv,log=T)
  ##
  lpy_y <- dmnorm(samps,rep(0,length(mu)),diag(length(mu)),log=T)
  ##
  lpx_invf_y <- rep(0,nrow(samps))
  for(samp_i in 1:nrow(samps)) {
    # norm2unif
    Par3 <- Par2 <- pnorm(samps[samp_i,],0,1)
    # unif2true
    Par3[idZ0] <- qbeta( Par2[idZ0], shape1  = ParamInitZ[idZ0,6], shape2 = ParamInitZ[idZ0,7])
    Par3[idZ1] <- qgamma(Par2[idZ1], shape   = ParamInitZ[idZ1,6], rate   = ParamInitZ[idZ1,7])
    Par3[idZ2] <- qnorm( Par2[idZ2], mean    = ParamInitZ[idZ2,6], sd     = ParamInitZ[idZ2,7])
    lpx_invf_y[samp_i] <- lPrior2(samps[samp_i,],Par3)
    if(lpx_invf_y[samp_i]%in%c(-Inf,Inf)) { lpx_invf_y[samp_i] <- NA }  }
  ldpn <- lpz_y - lpy_y + lpx_invf_y
  ret <- if(logg==T) { ldpn } else { exp(ldpn) }
  return(ret)  }


Imis2 <- function (B = 1000, B.re = 10000, number_k = 100, seed_val=1,SaveNam="xx",n_cores=2,OptInit) {
  # lprior = log prior for matrix with rows = samples, cols = pars (or a vector)
  # llikelihood = log likelihood for matrix with rows = samples,cols = pars
  # seed_val = seed
  # VCVinit = VCV of prior
  D = 1
  B0               = B * 10
  Sig2_global      = diag(length(OptInit))
  set.seed(seed_val+10^5)
  zd               = sample.prior(B0-1)
  for(ij in 1:length(OptInit)) zd[,ij] <- zd[,ij]+OptInit[ij]
  X_all            = X_k = rbind(OptInit,zd)
  stat_all         = matrix(NA,number_k,7); colnames(stat_all) <- c("Stage","UniquePoint","MaxWeight","ESS","MaxPost","NoDropped","NoSamp")
  center_all       = lprior_all = llike_all = lgaus1_all = NULL
  sigma_all        = list()
 option.opt = 0
  for (k in 1:number_k) { # k=1
    ptm.like       = proc.time()
    lgaus1_all     = c(lgaus1_all, dpropnew(X_k,OptInit,diag(length(OptInit)),logg=T))
    lprior_all     = c(lprior_all, lprior(X_k))
    llike_all      = c(llike_all, llikelihood(X_k,n_cores))
    ptm.use        = (proc.time() - ptm.like)[3]
    if (k == 1) { print(paste(B0, "likelihoods are evaluated in", round(ptm.use/60,1), "minutes")); flush.console() }
    if (k == 1) {  lenvelop_all = lgaus1_all } else { lenvelop_all = log(apply(rbind(exp(lgaus1_all)*B0/B,gaussian_all),2,sum)/(B0/B+D+(k-2))) }
    lWeights       = lprior_all + llike_all - lenvelop_all
    ### DROP THE WORST ##########
    # which to drop
    lWeightsZ = lWeights-max(lWeights,na.rm=T)
    lWeightsZ[is.na(lWeightsZ)] = -Inf
    Drop      = lWeightsZ < (-10000)
    lWeightsZ
    if(sum(Drop==F)>1e5)  Drop[rank(-lWeightsZ)>1e5] = T
    # drop them...
    X_all        = X_all[Drop==F,]
    lgaus1_all   = lgaus1_all[Drop==F]
    lprior_all   = lprior_all[Drop==F]
    llike_all    = llike_all[Drop==F]
    lWeights     = lWeights[Drop==F]
    lenvelop_all = lenvelop_all[Drop==F]
    if (k > 1) {
    gaussian_all = gaussian_all[,Drop==F] }
    ################## ########
    Weights        = exp(lWeights-max(lWeights))
    Weights        = Weights/sum(Weights)
    stat_all[k, 1] = k
    stat_all[k, 2] = sum(1-(1-Weights)^B.re)
    stat_all[k, 3] = max(Weights)
    stat_all[k, 4] = 1/sum(Weights^2)
    stat_all[k, 5] = max(lprior_all + llike_all,na.rm=T)[1]
    stat_all[k, 6] = sum(Drop)
    stat_all[k, 7] = length(Weights)

    print(round(stat_all[1:k,], 3)); flush.console()

    ### TEMP SAVE
    nonzero        = which(Weights>0)
    if(length(nonzero)>1) {
    which_X        = sample(nonzero,B.re,replace = T,prob=Weights[nonzero])
    resample_X     = X_all[which_X, ]
    resamp_lprior  = lprior_all[which_X]
    resamp_llike   = llike_all[which_X]
    resamp_lenvelop = lenvelop_all[which_X]
    imis_res_tmp   = list(stat            = stat_all,
                          resample        = resample_X,
                          resamp_lprior   = resamp_lprior,
                          resamp_llike    = resamp_llike,
                          resamp_lenvelop = resamp_lenvelop)
    save(imis_res_tmp,file=paste("imis_result_temp_",SaveNam,".rData",sep=""))
    }

    ###
      important = which(Weights == max(Weights))
      if (length(important) > 1) { important = important[1] }
        X_imp              = X_all[important, ]
        center_all         = rbind(center_all, X_imp)
        distance_all       = mahalanobis(X_all, X_imp, diag(diag(Sig2_global)))
        Bz                 = min(B,length(distance_all))
        which_var          = sort(distance_all,index=T)$ix[1:Bz]
        Sig2               = cov.wt(X_all[which_var, ],wt=Weights[which_var]+1/length(Weights),center=X_imp)$cov
        sigma_all[[D+k-1]] = Sig2
        jr                 = tryCatch({ X_k = rmnorm(B, X_imp, Sig2) }, error = function(e) NA)
        while(is.na(jr)[1]) {
        gf                 = 0
        jr                 = tryCatch({ X_k = rmnorm(B,X_imp,Sig2+diag(nrow(Sig2))/10^(12-gf)) }, error = function(e) NA)
        sigma_all[[D+k-1]] = Sig2+diag(nrow(Sig2))/10^(12-gf)
        gf                 = gf+1  }
        X_all              = rbind(X_all, X_k)

      ###
    if (k == 1) {
      gaussian_all = matrix(NA,D,dim(X_all)[1])
      for (i in 1:D) { gaussian_all[i, ] = dpropnew(X_all,center_all[i,],sigma_all[[i]],logg=F) }
     }
    if (k > 1) {
        gaussian_new = matrix(0, D + k - 1, dim(X_all)[1])
        gaussian_new[1:(D+k-2),1:(dim(X_all)[1]-B)] = gaussian_all
        gaussian_new[D+k-1,] = dpropnew(X_all,X_imp,sigma_all[[D+k-1]],logg=F)
        for (j in 1:(D+k-2)) { gaussian_new[j,(dim(X_all)[1]-B+1):dim(X_all)[1]] = dpropnew(X_k, center_all[j,],sigma_all[[j]],logg=F) }
       gaussian_all = gaussian_new
    }
    if(stat_all[k,3]>(1-exp(-1))*B.re) {  break }
  }

   return(list(stat            = stat_all,
               resample        = resample_X,
               resamp_lprior   = resamp_lprior,
               resamp_llike    = resamp_llike,
               resamp_lenvelop = resamp_lenvelop))
  }

