#################################################################
###   Bayesian kernel machine causal mediation analysis code  ###
###   developed by: Katrina Devick                            ###
###                                                           ###
###   To implement methods presented in arXiv: 1811.10453     ###       
##    last updated - 3/5/2019                                 ###
#################################################################


## load required library
library(bkmr)

## load the source file
## it needs to be in your working directory for it to load as is
source("source_BKMR_CMA.R")



##########################################################
###                 fit BKMR models                    ###
##########################################################


## let Z an n-by-L matrix containing the exposure mixture, y a vector of outcome data,  
## m a vector of mediator data, and Zm = cbind(Z,m) 

## *** NOTE: the last column of the Zm matrix MUST be your mediator in order for the functions to work properly! ***


## example from paper
# m <- dat.1$birthlength_z
# y <- dat.1$ccs_z
# Z <- cbind(dat.1$as_ln_z,dat.1$mn_ln_z,dat.1$pb_ln_z)
# X <- cbind(dat.1$momage_z,dat.1$age_z,dat.1$agesq,dat.1$smokenv,dat.1$sex,dat.1$momedud,dat.1$momIQ_z,dat.1$homescore_z,dat.1$proteintert)
# Zm <- cbind(Z,m)
# 
# colnames(Z) <- c("As","Mn","Pb")
# colnames(Zm) <- c("As","Mn","Pb","BL")



## fit BKMR models for the outcome, TE (outcome without the mediator), and mediator
## can use any features of BKMR (e.g. varsel=TRUE)
set.seed(1)
fit.y <- kmbayes(y=y, Z=Zm, X=X, iter=100000, verbose=TRUE, varsel=FALSE) 
save(fit.y,file="bkmr_y.RData")

set.seed(2)
fit.y.TE <- kmbayes(y=y, Z=Z, X=X, iter=100000, verbose=TRUE, varsel=FALSE) 
save(fit.y.TE,file="bkmr_y_TE.RData")

set.seed(3)
fit.m <- kmbayes(y=m, Z=Z, X=X, iter=100000, verbose=TRUE, varsel=FALSE) 
save(fit.m,file="bkmr_m.RData")


##### load models 
load("bkmr_y.RData")
load("bkmr_y_TE.RData")
load("bkmr_m.RData")




##################################################
### values at which to predict counterfactuals ###
##################################################


## mean level of confounders
X.predict <- matrix(colMeans(X),nrow=1)


## the change in exposure for which you want to estimate the mediation effects

## for the purpose of example, we will consider a change in all exposures from 
## their 25th to 75th percentiles. However, this contrast can be anything. 
astar <- apply(Z, 2, quantile, probs=0.25)
a     <- apply(Z, 2, quantile, probs=0.75)

## the index of the MCMC iterations to be used for inference 
sel <- seq(25001,100000,by=75)



##########################################################
####             Estimate TE for BKMR                 ####
##########################################################

## estimate the TE for a change in the exposures from astar to a
TE <- TE.bkmr(a=a, astar=astar, fit.y.TE=fit.y.TE, X.predict=X.predict, sel=sel, seed=122)

## look at the posterior mean, median, and 95% CI for TE
TE$est



##########################################################
####             Estimate CDE for BKMR                ####
##########################################################

## estimate the CDE for a change in the exposures from astar to a,
## fixing the mediator at its 25th, 50th, and 75th percentile
CDE <- CDE.bkmr(a=a, astar=astar, m.quant=c(0.25,0.5,0.75), fit.y=fit.y, sel=sel, seed=777)

## look at the posterior mean, median, and 95% CI for the CDEs 
CDE$est



##########################################################
####             Estimate NDE/NIE for BKMR            ####
##########################################################

## estimate the TE, NDE and NIE for a change in the exposures from astar to a

## *** NOTE: if the same confounders are used in both the mediation and outcome models, 
## X.predict.M and X.predict.Y are the same

## *** this step takes a while to run and will take longer for more complex bkmr fits, longer sel vectors and larger K
mediationeffects <- mediation.bkmr(a=a, astar=astar, fit.m=fit.m, fit.y=fit.y, fit.y.TE=fit.y.TE,
                                      X.predict.M=X.predict, X.predict.Y=X.predict, sel=sel, seed=22, K=1000)
## save this object
save(mediationeffects, file="mediationeffects.RData")

## look at the posterior mean, median, and 95% CI for the TE, NDE, and NIE
mediationeffects$est

