source("EM_tools.R",local=TRUE)
source("General_tools.R",local=TRUE)
source("PCA_tools.R",local=TRUE)
source("FISTA_tools.R",local=TRUE)
source("CodeMimi.R",local=TRUE)

library(missForest)
library(randomForest)
library(doParallel)
library(doSNOW)
cl <- makeCluster(20, outfile="") #here, change the number of clusters
registerDoSNOW(cl)

######
# Name: ComparMNAR_Multivariate
# Date: 27/12/2018
# Description: For the multivariate case, this function allow to compare different algorithms and methods to impute and estimate matrices which contain MNAR missing values.
# The function output is a list containing, for each simulation, the mean squared errors (the prediction error and the total error) for the different algorithms and methods.
# Arguments: 
  #Xtrue: the parameter matrix.
  #a, b: the logistic regression parameters.
  #r: rank of the parameter matrix.
  #noise: sigma^2, the added noise to the parameter matrix.
  #Ns: number of Monte Carlo simulations in the EM algorithm.
  #nbcol: number of missing variables, the missing variables are then the first nbcol ones. 
  #nbsim: number of simulations.
#####

ComparMNAR_Multivariate <- function(Xtrue,a,b,r,noise,Ns,nbcol,nbsim){
  
  p<-ncol(Xtrue)
  n<-nrow(Xtrue)
  
  results.list = foreach (ksim = 1:nbsim, .combine = "rbind")  %dopar% {
    source("EM_tools.R",local=TRUE)
    source("General_tools.R",local=TRUE)
    source("PCA_tools.R",local=TRUE)
    source("FISTA_tools.R",local=TRUE)
    source("CodeMimi.R",local=TRUE)
    
    print(paste("it globale",ksim))
    
    conv=FALSE
    cconv=0
    while(!conv & cconv<2){
      
      X=Xtrue+matrix(data=rnorm(n*p,0,sqrt(noise)),ncol=p)
      
      #MNAR with logistic regression
      select_prob <- function(x,a,b){ #probability of selecting coordinate Xij
        res=1/(1+exp(-a*(x-b)))
        return(res)
      }
      
      XNA=X
      for (i in 1:n){
        vec<-c()
        for(j in 1:nbcol){
          vec<-c(vec,select_prob(X[i,j],a,b))
        }
        
        proba<-combn(c(vec,1-vec),nbcol,prod)
        combi<-combn(c((1:nbcol),((1:nbcol)+nbcol)),nbcol)
        probanew<-c()
        combinew<-matrix(nrow=nbcol,ncol=1)
        for (k in 1:dim(combi)[2]){
          j=1
          vecbool=FALSE
          while(vecbool==FALSE & j<=nbcol){
            vecbool=is.element(j,combi[,k]) & is.element(j+nbcol,combi[,k])
            j=j+1
          }
          if(vecbool==FALSE){
            combinew=cbind(combinew,combi[,k])
            probanew=c(probanew,proba[k])
          }
        }
        ind<-sample(2:(2^nbcol+1),1,prob=probanew)
        XNA[i,combinew[,ind][which(combinew[,ind]<=nbcol)]]=NA
      }
      print(paste("fin simu", ksim))
      
      
      
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
      
      X.mean <- as.matrix(ImputeMean(XNA))
      
      
      ###############
      ####### EM with modell
      ###############
      
      #Fonction 
      ThetaNew=Initialize_theta(ImputeMean0(XNA),r)
      ThetabisNew=ThetaNew
      aNew=a-1
      bNew=b-1
      abisNew=a-1
      bbisNew=b-1
      diff=100
      ccompt<-0
      d1<-Sys.time()
      while(diff>10^-2 & ccompt<25){
        ParamNew <- IterEM(Xtrue,X,XNA,ThetaNew,aNew,bNew,M,Ns,noise,algo="soft",lam="Pred",nbcol=10)
        diff=ParamNew$diff
        print(diff)
        ThetaNew=ParamNew$ThetaNew
        aNew=ParamNew$a_initNew
        bNew=ParamNew$b_initNew
        conv1=ParamNew$conv
        ccompt=ccompt+1
      }
      d2<-Sys.time()
      print(paste("ccompt", ccompt, ksim))
      
      diff=100
      ccompt2<-0
      while(diff>10^-2 & ccompt2<25){
        ParamNew <- IterEM(Xtrue,X,XNA,ThetabisNew,abisNew,bbisNew,M,Ns,noise,algo="FISTA",lam="Pred",nbcol=10)
        diff=ParamNew$diff
        ThetabisNew=ParamNew$ThetaNew
        abisNew=ParamNew$a_initNew
        bbisNew=ParamNew$b_initNew
        conv2=ParamNew$conv
        ccompt2=ccompt2+1
      }
      print(paste("ccompt2", ccompt2, ksim))
      
      ThetaNewTot=Initialize_theta(ImputeMean0(XNA),r)
      ThetabisNewTot=ThetaNewTot
      aNew=a-1
      bNew=b-1
      abisNew=a-1
      bbisNew=b-1
      
      
      diff=100
      ccompt3=0
      while(diff>10^-2 & ccompt3<25){
        ParamNew <- IterEM(Xtrue,X,XNA,ThetaNewTot,aNew,bNew,M,Ns,noise,algo="soft",lam="Tot",nbcol=10)
        diff=ParamNew$diff
        print(diff)
        print( MSE(ThetaNewTot,Xtrue))
        ThetaNewTot=ParamNew$ThetaNew
        aNew=ParamNew$a_initNew
        bNew=ParamNew$b_initNew
        conv5=ParamNew$conv
        ccompt3=ccompt3+1
      }
      print(paste("ccompt3", ccompt3, ksim))
      
      
      diff=100
      ccompt4=0
      while(diff>10^-2 & ccompt4<25){
        ParamNew <- IterEM(Xtrue,X,XNA,ThetabisNewTot,abisNew,bbisNew,M,Ns,noise,algo="FISTA",lam="Tot",nbcol=10)
        diff=ParamNew$diff
        ThetabisNewTot=ParamNew$ThetaNew
        abisNew=ParamNew$a_initNew
        bbisNew=ParamNew$b_initNew
        conv6=ParamNew$conv
        ccompt4=ccompt4+1
      }
      print(paste("ccompt4", ccompt4, ksim))
      
      paste(c(conv1,conv2,conv5,conv6))
      if(conv1 & conv2  & conv5 & conv6 ){conv=TRUE}else(print("Pas de convergence"))
      cconv=cconv+1
    }
    
    
    
    ###############
    ####### Optimization: softImpute (SVD & softhresholding)
    ###############
    
    print(paste("softImpute",ksim))
    
    ####Without the mask
    RES=NULL
    RES2=NULL
    gridlambda1=seq(0, lambda0(X)*0.8, length = 300)
    for (i in 1:length(gridlambda1)){
      fit1=softImpute(as.matrix(XNA),rank=min(dim(as.matrix(XNA)))-1,lambda=gridlambda1[i],maxit = 10000,type="svd")
      if (fit1$d[1]==0){
        X1.soft=as.matrix(ImputeMean(XNA))
      }else if(length(fit1$d)==1){
        X1.soft= (fit1$u*fit1$d)%*%t(fit1$v)
      }else{
        X1.soft= (fit1$u%*%diag(fit1$d))%*%t(fit1$v)
      }
      RES2[i]=MSE(X1.soft*(1-M),X*(1-M))
      RES[i]=MSE(X1.soft,Xtrue)
    }
    fit1=softImpute(as.matrix(XNA),rank=min(dim(as.matrix(XNA)))-1,lambda=gridlambda1[which.min(RES)],maxit = 10000,type="svd")
    if (fit1$d[1]==0){
      X1.soft=as.matrix(ImputeMean(XNA))
    }else if(length(fit1$d)==1){
      X1.soft= (fit1$u*fit1$d)%*%t(fit1$v)
    }else{
      X1.soft= (fit1$u%*%diag(fit1$d))%*%t(fit1$v)
    }
    fit2=softImpute(as.matrix(XNA),rank=min(dim(as.matrix(XNA)))-1,lambda=gridlambda1[which.min(RES2)],maxit = 10000,type="svd")
    if (fit2$d[1]==0){
      X1.soft.bis=as.matrix(ImputeMean(XNA))
    }else if(length(fit2$d)==1){
      X1.soft.bis= (fit2$u*fit2$d)%*%t(fit2$v)
    }else{
      X1.soft.bis= (fit2$u%*%diag(fit2$d))%*%t(fit2$v)
    }
    
    ####With the mask
    RES=NULL
    RES2=NULL
    gridlambda1=seq(0, lambda0(X)*0.8, length = 300)
    for (i in 1:length(gridlambda1)){
      fit1=softImpute(as.matrix(Y),rank=min(dim(as.matrix(Y)))-1,lambda=gridlambda1[i],maxit = 10000,type="svd")
      if (fit1$d[1]==0){
        X1.soft.2=as.matrix(ImputeMean(Y))
      }else if(length(fit1$d)==1){ 
        X1.soft.2= (fit1$u*fit1$d)%*%t(fit1$v)
      }else{
        X1.soft.2= (fit1$u%*%diag(fit1$d))%*%t(fit1$v)
      }
      RES2[i]=MSE(X1.soft.2[,1:p]*(1-M),X*(1-M))
      RES[i]=MSE(X1.soft.2[,1:p],Xtrue)
    }
    fit1=softImpute(as.matrix(Y),rank=min(dim(as.matrix(Y)))-1,lambda=gridlambda1[which.min(RES)],maxit = 10000,type="svd")
    if (fit1$d[1]==0){
      X1.soft.2=as.matrix(ImputeMean(Y))
    }else if(length(fit1$d)==1){
      X1.soft.2= (fit1$u*fit1$d)%*%t(fit1$v)
    }else{
      X1.soft.2= (fit1$u%*%diag(fit1$d))%*%t(fit1$v)
    }
    fit2=softImpute(as.matrix(Y),rank=min(dim(as.matrix(Y)))-1,lambda=gridlambda1[which.min(RES2)],maxit = 10000,type="svd")
    if (fit2$d[1]==0){
      X1.soft.2.bis=as.matrix(ImputeMean(Y))
    }else if(length(fit2$d)==1){
      X1.soft.2.bis= (fit2$u*fit2$d)%*%t(fit2$v)
    }else{
      X1.soft.2.bis= (fit2$u%*%diag(fit2$d))%*%t(fit2$v)
    }
    
    ###############
    ####### Optimization: imputePCA (regularized iterative PCA)
    ###############
    
    print(paste("imputePCA",ksim))
    
    ####Without mask
    list=c(0,1,2,3,5,6,7,8,9,10)
    RES=NULL
    RES2=NULL
    for (i in 1:length(list)){
      X.pca.true=imputePCA(XNA,ncp=list[i],maxiter=10000)$fittedX
      RES2[i]=MSE(X.pca.true*(1-M),X*(1-M))
      RES[i]=MSE(X.pca.true,Xtrue)
    }
    X.pca.true=imputePCA(XNA,ncp=list[which.min(RES)],maxiter=10000)$fittedX
    X.pca.true.bis=imputePCA(XNA,ncp=list[which.min(RES2)],maxiter=10000)$fittedX
    
    ####With mask
    list=c(0,1,2,3,5,6,7,8,9,10,11,12)
    RES=NULL
    RES2=NULL
    for (i in 1:length(list)){
      X.pca.true.2=imputePCA(Y,ncp=list[i],maxiter=10000)$fittedX
      RES2[i]=MSE(X.pca.true.2[,1:p]*(1-M),X*(1-M))
      RES[i]=MSE(X.pca.true.2[,1:p],Xtrue)
    }
    X.pca.true.2=imputePCA(Y,ncp=list[which.min(RES)],maxiter=10000)$fittedX
    X.pca.true.2.bis=imputePCA(Y,ncp=list[which.min(RES2)],maxiter=10000)$fittedX
    
    ###############
    ####### Mimi Algorithm
    ###############
    
    print(paste("mimi",ksim))

    var.type = c(rep("gaussian", p),rep("binary", nbcol))
    test <- mimi(Y, model = "low-rank", var.type = var.type, lambda1 = 0.05, trace.it = F, max.rank = min(dim(X))-1)
    RES=NULL
    RES2=NULL
    for (i in 1:length(test$list.theta)){
      X.mimi <- test$list.theta[[i]]
      RES2[[i]]=MSE(X.mimi[,1:p]*(1-M),X*(1-M))
      RES[[i]]=MSE(X.mimi[,1:p],Xtrue)
    }
    X.mimi <- test$list.theta[[which.min(RES)]]
    X.mimi.bis <- test$list.theta[[which.min(RES2)]]
    
    ###############
    ####### FISTA Algorithm
    ###############
    
    print(paste("FISTA",ksim))
    
    RES=NULL
    RES2=NULL
    gridParam=seq(0,lambda0(X)*1.1, length = 100) 
    for(i in 1:length(gridParam)){
      X.FISTA2<- FISTANA(as.matrix(XNA),as.matrix(M),gridParam[i],noise,alg="other")
      RES2[i]=MSE(X.FISTA2[,1:p]*(1-M),X*(1-M))
      RES[i]=MSE(X.FISTA2[,1:p],Xtrue)
    }
    X.FISTA2<- FISTANA(XNA,M,gridParam[which.min(RES)],noise,alg="other")
    X.FISTA2.bis<- FISTANA(XNA,M,gridParam[which.min(RES2)],noise,alg="other")
    
    RES=NULL
    RES2=NULL
    gridParam=seq(0, lambda0(X)*1.1, length = 100) 
    for(i in 1:length(gridParam)){
      X.FISTA3<- FISTANA(Y,cbind.data.frame(M,M[,1:nbcol]),gridParam[i],noise,alg="other")
      RES2[i]=MSE(X.FISTA3[,1:p]*(1-M),X*(1-M))
      RES[i]=MSE(X.FISTA3[,1:p],Xtrue)
    }
    X.FISTA3<- FISTANA(Y,cbind.data.frame(M,M[,1:nbcol]),gridParam[which.min(RES)],noise,alg="other")
    X.FISTA3.bis<- FISTANA(Y,cbind.data.frame(M,M[,1:nbcol]),gridParam[which.min(RES2)],noise,alg="other")
    
    ###############
    ####### Random Forest
    ###############
    
    print(paste("rf",ksim))
    
    X.rf <- missForest(XNA,ntree=100)$ximp
    
    Yrf=Y
    Yrf[Yrf==0]=-10^6
    Yrf[Yrf==1]=10^6
    X.rf.2 <- missForest(Yrf,ntree=100)$ximp
    
    
    ###############
    ####### MSE
    ###############
    
    msePred=sapply(list(X.mean[,1:p]*(1-M),ThetaNew[,1:p]*(1-M),ThetabisNew[,1:p]*(1-M),ThetaNewTot[,1:p]*(1-M),ThetabisNewTot[,1:p]*(1-M),X1.soft[,1:p]*(1-M),X1.soft.bis[,1:p]*(1-M),X1.soft.2[,1:p]*(1-M),X1.soft.2.bis[,1:p]*(1-M),X.pca.true[,1:p]*(1-M),X.pca.true.bis[,1:p]*(1-M),X.pca.true.2[,1:p]*(1-M),X.pca.true.2.bis[,1:p]*(1-M),X.mimi[,1:p]*(1-M),X.mimi.bis[,1:p]*(1-M),X.FISTA2[,1:p]*(1-M),X.FISTA2.bis[,1:p]*(1-M),X.FISTA3[,1:p]*(1-M),X.FISTA3.bis[,1:p]*(1-M),X.rf[,1:p]*(1-M),X.rf.2[,1:p]*(1-M)),MSE,X2=X*(1-M))
    mseTrue=sapply(list(X.mean[,1:p],ThetaNew[,1:p],ThetabisNew[,1:p],ThetaNewTot[,1:p],ThetabisNewTot[,1:p],X1.soft[,1:p],X1.soft.bis[,1:p],X1.soft.2[,1:p],X1.soft.2.bis[,1:p],X.pca.true[,1:p],X.pca.true.bis[,1:p],X.pca.true.2[,1:p],X.pca.true.2.bis[,1:p],X.mimi[,1:p],X.mimi.bis[,1:p],X.FISTA2[,1:p],X.FISTA2.bis[,1:p],X.FISTA3[,1:p],X.FISTA3.bis[,1:p],X.rf[,1:p],X.rf.2[,1:p]),MSE,X2=Xtrue)
    
    names=c("Imputemean","modelsoftPred","modelFISTAPred","modelsoftTot","modelFISTATot","softTot","softPred","softmaskTot","softmaskPred","PCATot","PCAPred","PCAmaskTot","PCAmaskPred","mimiTot","mimiPred","FISTA","FISTAPred","FISTAmask","FISTAmaskPred","randomforest","randomforestmask")
    cbind(msePred,mseTrue,names)
  }
  
  
  return(results.list)
}
