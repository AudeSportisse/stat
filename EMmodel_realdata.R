source("crossval.R")

library(missForest)
library(randomForest)
library(softImpute)

######
# Name: EMMNAR_Multivariate_realdata
# Date: 24/11/2019
# Description: 
# The function output is a list containing, in a list, the imputed datasets using our method and the imputation by the mean and in a vector, the names of the methods. 
# Arguments: 
#XNA: the incomplete matrix
#a, b: the logistic regression parameters.
#noise: sigma^2, the added noise to the parameter matrix.
#Ns: number of Monte Carlo simulations in the EM algorithm.
#nbcol: number of missing variables, the missing variables are then the first nbcol ones. 
#scale: boolean indicating if the variables are scaled or not for the EM MNAR method. 
#nbit: maximal number of iteration to be performed.
#####

EMMNAR_Multivariate_realdata <- function(XNA,a,b,noise,Ns,nbcol,scale,nbit){
  
  XNA <- as.matrix(XNA)
  p<-ncol(XNA)
  n<-nrow(XNA)
  colnamesXNA <- colnames(XNA)
  indexNA <- which(is.na(XNA))
  
  
  ###############
  ####### Missing-data pattern: 
  ###############
  
  M = 1-is.na(XNA)
  
  
  ###############
  ####### Concatenating the data matrix and the mask:
  ###############
  
  
  Y = cbind.data.frame(XNA,M[,1:nbcol])
  
  
  ###############
  ####### Optimization: Mean imputation and estimation
  ###############
  
  X.mean <- as.matrix(ImputeMean0(XNA))
  colnames(X.mean) <- colnamesXNA
  
  ThetaNew=X.mean
  ThetabisNew=ThetaNew
  aNew=a
  bNew=b
  diff=100
  ccompt<-0
  if(scale==TRUE){
    mean.thet <- apply(ThetaNew, 2, mean) 
    ThetaNew <- t(t(ThetaNew) - mean.thet)
    XNA <- t(t(XNA) - mean.thet)
    et.thet <- apply(ThetaNew, 2, sd)
    ThetaNew <- t(t(ThetaNew)/et.thet)
    XNA <- t(t(XNA)/et.thet)
  }else{
    mean.thet <- NULL
    et.thet <- NULL
  }
  while(diff>10^-5 & ccompt<nbit ){
    ParamNew <- IterEM_realdata(XNA,ThetaNew,aNew,bNew,M,Ns,noise,algo="soft",lam="Pred",nbcol=nbcol,scale=scale,mean.thet,et.thet)
    et.thet = ParamNew$et.fin
    mean.thet = ParamNew$mean.fin
    diff=ParamNew$diff
    ThetaNew=ParamNew$ThetaNew
    XNA_remp=ParamNew$XNA_remp
    aNew=ParamNew$a_initNew
    bNew=ParamNew$b_initNew
    XNA=ParamNew$XNA
    ccompt=ccompt+1
  }
  if (scale==TRUE){ 
    ThetaNew <- t(t(ThetaNew) * et.thet)
    ThetaNew <- t(t(ThetaNew) + mean.thet)
    XNA <- t(t(XNA) * et.thet)
    XNA <- t(t(XNA) + mean.thet)
  }
  print(paste("ccompt", ccompt))
  
  
  ThetaNew_imp=XNA
  ThetaNew_imp[indexNA]=ThetaNew[indexNA]
  
  colnames(ThetaNew) <- colnamesXNA
  colnames(ThetaNew_imp) <- colnamesXNA
  
  imputed=list(X.mean,ThetaNew_imp)
  names=c("Imputemean","modelsoftPred")
  results.list = list(data_imp=imputed,names=names)
  
  return(results.list)
}


#Imputation by the mean
ImputeMean0 <- function(tab){
  m <- apply(tab, 2, mean, na.rm = TRUE)
  tab <- sapply(1:ncol(tab), function(x) ifelse(is.na(tab[,x]), m[x], tab[,x]))
  tab <- as.data.frame(tab)
  return(tab)
}


######
# Name: IterEM_realdata
# Date: 24/11/2019
# Description: It computes one iteration of the EM algorithm for the univariate and multivariate cases. In particular, the mechanism parameters are the same ones for all the variables and the missing variables are the first nbcol variables.
# This function returns: the difference between the last matrix iteration and the new one, the new matrix, the new mechanisms parameters and if the GLM algorithm has converged.
# Arguments: 
#XNA: data matrix X containing missing values.  
#Thet: old matrix iteration.
#a, b: old mechanism parameters. 
#M: missing-data matrix.
#Ns: number of Monte Carlo simulations.
#noise: sigma^2.
#algo: "soft" or "FISTA" depending on the algorithm used for M-step. 
#lam: "Tot" or "Pred" depending on how the lambda is chosen. 
#nbcol: number of misisng variables.
#####


IterEM_realdata <- function(XNA,
                            Thet,
                            a,
                            b,
                            M,
                            Ns,
                            noise,
                            algo,
                            lam,
                            nbcol,
                            scale=TRUE, 
                            mean.init, 
                            et.init) {
  a_initOld = a
  b_initOld = b
  
  ThetaOld = Thet
  if (scale==TRUE){
    ThetaOld <- t(t(ThetaOld) * et.init)
    ThetaOld <- t(t(ThetaOld) + mean.init)
    XNA <- t(t(XNA) * et.init)
    XNA <- t(t(XNA) + mean.init)
    mean.fin <- apply(ThetaOld, 2, mean)
    ThetaOld <- t(t(ThetaOld) - mean.fin)
    XNA <- t(t(XNA) - mean.fin)
    et.fin <- apply(ThetaOld, 2, sd)
    ThetaOld <- t(t(ThetaOld)/et.fin)
    XNA <- t(t(XNA)/et.fin)
  }else{
    et.fin <- rep(1,ncol(XNA))
    mean.fin <- 1
  }
  
  SIR <- function(Theta, XNA, a, b, sigma, i, jind) {
    Nn = Ns * 1000
    sim = rnorm(Nn, mean = Theta[i, jind], sd = sigma)
    
    RapportRej <- function(x,a,b,i,jind) {
      g <- function(x) {
        res = ((1 / (1 + exp(-a * (
          x - b
        )))) ^ (1 - M[i, jind])) * ((1 - (1 / (
          1 + exp(-a * (x - b))
        ))) ^ M[i, jind])
        return(res)
      }
      
      num = g(x)
      return(num)
    }
    poids = RapportRej(sim,a,b,i,jind)
    res <- sample(sim, Ns, prob = poids, replace = TRUE)
    return(res)
  }
  
  sigma = sqrt(noise)
  XNA_remp = XNA
  data_remp <- array(data=NA,dim=c(Ns,nrow(XNA),nbcol))
  for (i in 1:nrow(XNA)) {
    for (j in 1:ncol(XNA)) {
      if (is.na(XNA[i, j])) {
        ech = SIR(ThetaOld, XNA, a_initOld[j], b_initOld[j], sigma, i, j)
        XNA_remp[i, j] = mean(ech)
        data_remp[,i,j] = ech
      }
    }
  }
  
  gridlambda1 = cv_soft(XNA_remp,len=100)
  gridParam = gridlambda1
  if (algo == "soft") {
    fit1 = softImpute(
      as.matrix(XNA_remp),
      rank = min(dim(XNA_remp)) - 1,
      lambda = gridlambda1,
      maxit = 10000,
      type = "svd"
    )
    if (fit1$d[1] == 0) {
      ThetaNew = as.matrix(ImputeMean(XNA_remp))
    } else if (length(fit1$d) == 1) {
      ThetaNew = (fit1$u * fit1$d) %*% t(fit1$v)
    } else{
      ThetaNew = (fit1$u %*% diag(fit1$d)) %*% t(fit1$v)
    }
  } else{
    ThetaNew = FISTA0(XNA_remp, gridParam, noise)
  }
  
  
  a_initNew <- numeric(length(a))
  b_initNew <- numeric(length(b))
  for (j in 1:nbcol){
    M_concat = rep(as.vector(M[,j]), Ns )
    X_concat = matrix(0, nrow = Ns * nrow(XNA) , ncol = 1)
    for (i in 1:nrow(XNA)) {
      for (k in (1:Ns)) {
        if (is.na(XNA[i,j]) == TRUE) {
          X_concat[nrow(XNA) * (k - 1) + i, 1] = data_remp[k,i,j]
        } else{
          X_concat[nrow(XNA) * (k - 1) + i, 1] = XNA[i, j]
        }
      }
    }
    
    GLM =  glm(as.factor(M_concat) ~ as.vector(X_concat), family = binomial(logit))
    a_initNew[j] = -GLM$coefficients[2]
    b_initNew[j] = -GLM$coefficients[1] / GLM$coefficients[2]
  }
  
  diff = norm(ThetaNew - as.matrix(ThetaOld), type = "F") / (norm(as.matrix(ThetaOld), type =
                                                                    "F") + 10 ^ -3)
  
  return(
    list(
      diff = diff,
      ThetaNew = ThetaNew,
      XNA_remp = XNA_remp,
      a_initNew = a_initNew,
      b_initNew = b_initNew,
      mean.fin = mean.fin,
      et.fin = et.fin,
      XNA = XNA
    )
  )
  
}



######
# Name: Initialize_theta
# Date: 27/12/2018
# Description: Initialization of the EM algorithm: it performs the SVD of the matrix X and keeps r dimensions. This function returns a matrix.
# Arguments: 
#X: matrix.
#r: the rank of the output matrix.
#####


Initialize_theta <- function(X, r) {
  res = svd(X)
  result = (res$u[, 1:r] * res$d[r]) %*% (t(res$v[, 1:r]))
  return(result)
}

###
# Functions for scale
###

moy.p <- function(V, poids) {
  res <- sum(V * poids, na.rm = TRUE)/sum(poids[!is.na(V)])
}
ec <- function(V, poids) {
  res <- sqrt(sum(V^2 * poids, na.rm = TRUE)/sum(poids[!is.na(V)]))
}