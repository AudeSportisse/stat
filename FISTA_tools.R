library(corpcor)
source("General_tools.R",local=TRUE)


prox_l1 <- function(x,alg,tau,ncp){ #computes the projection on the set M.*X = X_na
  if (alg=="hard"){res=(abs(x)>=sort(abs(x),decreasing=TRUE)[ncp])*x}else{res=max(0,1-tau/max(1e-15,abs(x)))*x}
  return(res)
}


soft_thresh_svd <- function(X,alg,gamma,ncp){ #performs soft thresholding in the value decomposition
  res=fast.svd(X)
  U=res$u
  V=res$v
  S=diag(res$d)
  if(alg=="hard"){
    s=prox_l1(diag(S),alg,gamma,ncp) 
  }else{
    s=sapply(diag(S),prox_l1,alg,gamma,ncp) 
  }
  if(alg=="hard"){
    if(ncp==0){
      X_res=as.matrix(ImputeMean(X))
    }else{
      S[1:length(s),1:length(s)] = diag(s)
      X_res=U%*%S%*%t(V) 
    }
  }else{
    S[1:length(s),1:length(s)] = diag(s)
    X_res=U%*%S%*%t(V) 
  }
  return(X_res)
}


######
# Name: FISTA0
# Date: 27/12/2018
# Description: FISTA algorithm without missing values.
# Arguments: 
  #X: matrix.
  #param: the regularization parameter of the optimization problem. 
  #noise: sigma^2. 
#####


FISTA0 <- function(X,param,noise){
  X=as.matrix(X) 
  niter=200
  lambdaNew=0.1
  xNew=yNew=matrix(0,ncol=ncol(X),nrow=nrow(X))
  diff=1
  while (diff>10^(-6)){
    xOld=xNew
    yOld=yNew
    lambdaOld=lambdaNew
    xNew=soft_thresh_svd(as.matrix(yOld-(yOld-X)),alg="soft",param,ncp=NULL) #soft-thresholding by default. 
    lambdaNew=(1+sqrt(1+4*lambdaOld**2))/2
    yNew=xNew+((lambdaOld-1)/lambdaNew)*(xNew-xOld)
    diff=norm(xNew-xOld,type="F")
  }
  return(xNew)
}


######
# Name: FISTANA
# Date: 27/12/2018
# Description: FISTA algorithm with missing values.
# Arguments: 
  #X: matrix.
  #M: the missing-data matrix.
  #param: the regularization parameter of the optimization problem. 
  #noise: sigma^2. 
  #alg: "soft" for soft-thresholding, "hard" for hard-thresholding, "thresh" for soft-thresholding and for each iteration of the FISTA algorithm, the missing value is imputed by the maximum between the threshold (assumed to be known) and the imputed value. 
#####


FISTANA <- function(X,M,param,noise,alg,ncp=NULL){ 
  X=as.matrix(X) 
  missing=which(is.na(X))
  X=ImputeMean0(X)
  lambdaNew=0.1
  xNew=yNew=matrix(0,ncol=ncol(X),nrow=nrow(X))
  diff=1
  while(diff>10^(-6)){
    xOld=xNew
    yOld=yNew
    lambdaOld=lambdaNew
    xNew=soft_thresh_svd(as.matrix(yOld-(M*(yOld-X))),alg,param,ncp)
    if(alg=="thresh"){
      thresh=max(xNew[missing]) 
      for( ind in missing){
        xNew[ind]=max(thresh,xNew[ind])
      }
    }
    lambdaNew=(1+sqrt(1+4*lambdaOld**2))/2
    yNew=xNew+((lambdaOld-1)/lambdaNew)*(xNew-xOld)
    diff=norm(xNew-xOld,type="F")
  }
  
  return(xNew)
}


######
# Name: FISTA
# Date: 27/12/2018
# Description: FISTA algorithm by modelling the mechanism with an intuitive idea if the missing values are introduced in the parameter matrix. 
# Arguments: 
  #X: matrix.
  #M: the missing-data matrix.
  #param: the regularization parameter of the optimization problem. 
  #a,b: mechanism parameters
  #noise: sigma^2. 
######

FISTA <- function(X,M,param,a,b,noise){ 
  X=as.matrix(X)
  X=ImputeMean0(X) 
  lambdaNew=0.1
  xNew=yNew=matrix(0,ncol=ncol(X),nrow=nrow(X))
  diff=1
  while(diff>10^(-6)){
    xOld=xNew
    yOld=yNew
    lambdaOld=lambdaNew
    grad=matrix(0,nrow=nrow(X),ncol=ncol(X))
    for(i in 1:nrow(X)){
      for(j in 1:ncol(X)){
        grad[i,j]=(a*exp(a*yOld[i,j])/(exp(a*b)+exp(a*yOld[i,j])))
      }
    }
    xNew=soft_thresh_svd(as.matrix(yOld-(1/((1/noise)+a^2/4))*((M*(1/noise)*(yOld-X)+M*grad))),param)
    lambdaNew=(1+sqrt(1+4*lambdaOld**2))/2
    yNew=xNew+((lambdaOld-1)/lambdaNew)*(xNew-xOld)
    diff=norm(xNew-xOld,type="F")
  }
  return(xNew)
}

