##########################################################
####      Posterior/Bootstrap Summary Function        ####
##########################################################


##### function to calculate mean/median/CI/sd of a vector of posterior/bootstrap samples

postresults <- function(posteriorsamp){
  toreturn <- vector()
  toreturn["mean"]   <- mean(posteriorsamp, na.rm=TRUE)
  toreturn["sd"]     <- sd(posteriorsamp, na.rm=TRUE)
  toreturn["lower"]  <- quantile(posteriorsamp, probs=0.025, na.rm=TRUE)
  toreturn["median"]  <- quantile(posteriorsamp, probs=0.5, na.rm=TRUE)
  toreturn["upper"] <- quantile(posteriorsamp, probs=0.975, na.rm=TRUE)
  return(toreturn)
}


##########################################################
####             Estimate TE for BKMR                 ####
##########################################################


TE.bkmr <- function(a, astar, fit.y.TE, X.predict, sel, seed){
  toreturn <- list()
  
  set.seed(seed)
  newz  <- rbind(a,astar)
  TE.mat <- SamplePred(fit.y.TE, Znew = newz, Xnew = X.predict, sel=sel) 

  Ya     <- TE.mat[,"znew1"]
  Yastar <- TE.mat[,"znew2"]
  
  toreturn$TE.samp     <- as.vector(Ya - Yastar)
  toreturn$Ya.samp     <- as.vector(Ya) 
  toreturn$Yastar.samp <- as.vector(Yastar) 
  
  TE.sum <- postresults(toreturn$TE)
  toreturn$est  <- TE.sum[c("mean","median","lower","upper")]
  return(toreturn)
}



##########################################################
####           Estimate CDE for BKMR                  ####
##########################################################


CDE.bkmr <- function(a, astar, m.quant, fit.y, X.predict=rep(0,ncol(fit.y$X)), sel, seed){
  toreturn <- list()
  m <- fit.y$Z[,ncol(fit.y$Z)]  ### okay as long as m is the LAST variable in Zm
  Z <- fit.y$Z[,-ncol(fit.y$Z)]
  
  toreturn$est <- matrix(NA, nrow=length(m.quant), ncol=4, dimnames=list(paste0("CDE",m.quant*100), c("mean","median","lower","upper")))
  for(i in seq_along(m.quant)){
    print(paste(i, "out of", length(m.quant)))
    mnew <- quantile(m, probs=m.quant[i])
    
    set.seed(seed)
    newZm  <- rbind(c(a,mnew),c(astar,mnew))
    CDE.mat <- SamplePred(fit.y, Znew = newZm, Xnew = X.predict, sel=sel) 
    
    Yam     <- CDE.mat[,"znew1"]
    Yastarm <- CDE.mat[,"znew2"]
    
    CDE <- as.vector(Yam - Yastarm)
    toreturn[[paste0("CDE",m.quant[i]*100,".samp")]] <- CDE
    toreturn$est[paste0("CDE",m.quant[i]*100),] <- postresults(CDE)[c("mean","median","lower","upper")]
  }
  
  return(toreturn)
}



  
##########################################################
####           YaMastar counterfactual                ####
##########################################################


YaMastar.SamplePred <- function(a, astar, fit.m, fit.y, X.predict.M, X.predict.Y, sel, seed, K){
  start.time <- proc.time()
 
  set.seed(seed)
  EM.samp <- SamplePred(fit.m, Znew = astar, Xnew = X.predict.M, sel=sel) 
  Mastar     <- as.vector(EM.samp)
  
  sigma.samp  <- sqrt(fit.m$sigsq.eps[sel])
  random.samp <- matrix(rnorm(length(sel)*K),nrow=length(sel),ncol=K)
  
  Mastar.samp  <- Mastar + sigma.samp*random.samp
  
  YaMastar.samp.mat     <- matrix(NA,nrow=length(sel),ncol=K)
  for(j in 1:length(sel)){
    Mastar.j <-  Mastar.samp[j,]
    aMastar.j <- cbind(matrix(a, nrow=K, ncol=length(a), byrow=TRUE), Mastar.j)
    row.seed <- j + 10000
    set.seed(row.seed)
    YaMastar.j <- SamplePred(fit.y, Znew = aMastar.j, Xnew = X.predict.Y, sel=sel[j])
    YaMastar.samp.mat[j,] <- as.vector(YaMastar.j)
    
    end.time.temp <- proc.time() 
    if(j%%50==0) print(paste("iter", j, "time: ", round((end.time.temp - start.time)["elapsed"]/60,2),"min")) 
  }
  toreturn <- apply(YaMastar.samp.mat,1,mean)
  
  return(toreturn)
}



##########################################################
####           Estimate NDE/NIE for BKMR              ####
##########################################################


mediation.bkmr <- function(a, astar, m.quant=c(0.25,0.5,0.75), fit.m, fit.y, fit.y.TE, X.predict.M, X.predict.Y, sel, seed, K){
  
  toreturn <- list()
  
  TE <- TE.bkmr(a=a, astar=astar, fit.y.TE=fit.y.TE, X.predict=X.predict.Y, sel=sel, seed=(seed+100))
  
  Ya     <- TE$Ya.samp
  Yastar <- TE$Yastar.samp
  
  
  YaMastar <- YaMastar.SamplePred(a=a, astar=astar, fit.m=fit.m, fit.y=fit.y,
                                        X.predict.M=X.predict.M, X.predict.Y=X.predict.Y, sel=sel, seed=seed, K=K)
  
  NDE <- YaMastar - Yastar
  NIE <- Ya - YaMastar
    

  toreturn$TE.samp <- TE$TE.samp
  toreturn$NDE.samp <- NDE
  toreturn$NIE.samp <- NIE
  
  toreturn$est <- matrix(NA, nrow=3, ncol=4, dimnames=list(c("TE","NDE","NIE"), c("mean","median","lower","upper")))
  toreturn$est[c("TE","NDE","NIE"),] <- rbind(postresults(TE$TE.samp) [c("mean","median","lower","upper")],
                                              postresults(NDE)[c("mean","median","lower","upper")],
                                              postresults(NIE)[c("mean","median","lower","upper")])
  
  return(toreturn)
}

