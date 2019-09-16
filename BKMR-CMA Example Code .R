################################################
###   Bayesian kernel machine regression--   ###
###   causal mediation analysis code         ###
###   developed by: Katrina Devick           ###
###                                          ###
###   To implement methods presented in      ###
###   arXiv: 1811.10453                      ###       
###                                          ###
###   last updated - 16 September 2019       ###
################################################


## load required library
library(bkmr)

## load the source file
## it needs to be in your working directory for it to load as is
source("source_BKMR_CMA.R")




##########################################################
###                 fit BKMR models                    ###
##########################################################


## let A be a n-by-L matrix containing an exposure mixture comprised of L elements,
## E.M and E.Y be effect modifiers of exposure-mediator and exposure-outcome relationship respectively,
## y, a vector of outcome data, and  
## m, a vector of mediator data

## Now, let Z.M be the exposures and effect modifers E.M and 
## let Z.Y be the exposures and effect modifers E.Y
## Z.M <- cbind(A,E.M) 
## Z.Y <- cbind(A,E.Y) 

## create one more matrix containing the exposures, effect modifier Z.Y and mediator, 
## precisely in that order (!!! This is very important !!!)
## Zm.Y <- cbind(Z.Y,m)
## *** NOTE: the last column of the Zm.Y matrix MUST be your mediator in order for the functions to work properly! ***


############################
##   example from paper   ##
############################

A <- cbind(dat$as_ln_z,dat$mn_ln_z,dat$pb_ln_z)
m <- dat$birthlength_z
y <- dat$ccs_z

## additional covariates assumed to have a linear effect with the mediator and outcome
X <- cbind(dat$momage_z,dat$smokenv,dat$sex,dat$momedud,dat$momIQ_z,dat$homescore_z,dat$proteintert)

E.M <- NULL
E.Y <- dat$age_z

Z.M <- cbind(A,E.M) 
Z.Y <- cbind(A,E.Y) 
Zm.Y <- cbind(Z.Y,m)

colnames(Z.M)  <- colnames(A)  <- c("As","Mn","Pb")
colnames(Z.Y)  <- c("As","Mn","Pb","age")
colnames(Zm.Y) <- c("As","Mn","Pb","age","BL")



## fit BKMR models for the outcome, TE (outcome without the mediator), and mediator
## can use any features of BKMR (e.g. varsel=TRUE)
set.seed(1)
fit.y <- kmbayes(y=y, Z=Zm.Y, X=X, iter=20000, verbose=TRUE, varsel=FALSE) 
save(fit.y,file="bkmr_y.RData")

set.seed(2)
fit.y.TE <- kmbayes(y=y, Z=Z.Y, X=X, iter=20000, verbose=TRUE, varsel=FALSE) 
save(fit.y.TE,file="bkmr_y_TE.RData")

set.seed(3)
fit.m <- kmbayes(y=m, Z=Z.M, X=X, iter=20000, verbose=TRUE, varsel=FALSE) 
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

## We will consider a change in all exposures from their 25th to 75th percentiles 
## fixing age at testing to its 10th and 90th percentiles 
## However, this contrast can be anything. 

## if modifiers are considered, you should fix the levels of the modifiers 

astar.age10 <- c(apply(A, 2, quantile, probs=0.25), quantile(E.Y, probs=0.1))
astar.age90 <- c(apply(A, 2, quantile, probs=0.25), quantile(E.Y, probs=0.9))

a.age10 <- c(apply(A, 2, quantile, probs=0.75), quantile(E.Y, probs=0.1))
a.age90 <- c(apply(A, 2, quantile, probs=0.75), quantile(E.Y, probs=0.9))


## the index of the MCMC iterations to be used for inference 
sel<-seq(5001,20000,by=15)



##########################################################
####             Estimate TE for BKMR                 ####
##########################################################

## estimate the TE for a change in the exposures from astar to a
## fixing age at testing to its 10th percentile 
TE.age10 <- TE.bkmr(a=a.age10, astar=astar.age10, fit.y.TE=fit.y.TE, X.predict=X.predict, sel=sel, seed=122)

## look at the posterior mean, median, and 95% CI for TE
TE.age10$est



## repeat, now fixing age at testing to its 90th percentile 
TE.age90 <- TE.bkmr(a=a.age90, astar=astar.age90, fit.y.TE=fit.y.TE, X.predict=X.predict, sel=sel, seed=122)

## look at the posterior mean, median, and 95% CI for TE
TE.age90$est




##########################################################
####             Estimate CDE for BKMR                ####
##########################################################

## estimate the CDE for a change in the exposures from astar to a,
## fixing the mediator at its 10th, 50th, and 75th percentile and
## age at testing at its 10th percentile 
CDE.age10 <- CDE.bkmr(a=a.age10, astar=astar.age10, m.quant=c(0.1,0.5,0.75), fit.y=fit.y, sel=sel, seed=777)

## look at the posterior mean, median, and 95% CI for the CDEs 
CDE.age10$est



##########################################################
####             Estimate NDE/NIE for BKMR            ####
##########################################################

## estimate the TE, NDE and NIE for a change in the exposures from astar to a
## fixing age at testing to its 90th percentile

## *** NOTE: if the same confounders are used in both the mediation and outcome models, 
## X.predict.M and X.predict.Y are the same

## *** this step takes a while to run and will take longer for more complex bkmr fits, longer sel vectors and larger K
mediationeffects.age90 <- mediation.bkmr(a=a.age90, astar=astar.age90, fit.m=fit.m, fit.y=fit.y, fit.y.TE=fit.y.TE,
                                      X.predict.M=X.predict, X.predict.Y=X.predict, sel=sel, seed=22, K=1000)
## save this object
save(mediationeffects.age90, file="mediationeffects_age90.RData")

## look at the posterior mean, median, and 95% CI for the TE, NDE, and NIE
mediationeffects.age90$est

