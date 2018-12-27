source("EM_tools.R",local=TRUE)
source("General_tools.R",local=TRUE)
source("PCA_tools.R",local=TRUE)
source("FISTA_tools.R",local=TRUE)
source("CodeMimi.R",local=TRUE)

library(randomForest)
library(missForest)
library(doParallel)
library(doSNOW)
cl <- makeCluster(20, outfile="")
registerDoSNOW(cl)

ComparMNAR_Univariate <- function(Xtrue,a,b,r,bruit,Ns,modmecha,mecha,nbsim){
  
  nbcol=1
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
    while(!conv){
      
      set.seed(ksim)
      X=Xtrue+matrix(data=rnorm(n*p,0,sqrt(bruit)),ncol=p)
      
      #Logit or probit distribution
      select_prob <- function(x,modmecha){ #probability of selecting coordinate Xij
        if(modmecha=="logit"){
        res=1/(1+exp(-a*(x-b)))}else{ res=pnorm(x)}
        return(res)
      }
      
      #MNAR or MAR mechanism
      if(mecha=="MNAR"){
      prob <- sapply(X[,1],select_prob,modmecha)}else{prob <- sapply(X[,2],select_prob,modmecha)}
      compt=0
      missing=c()
      for (i in 1:n){
        u<-runif(1)
        compt=compt+(prob[i]>u)
        if(prob[i]>u){missing=c(missing,i)}
      }
      print(compt)
      XNA=X
      XNA[missing,1]=NA
      
      
      ###############
      ####### Missing-data pattern
      ###############
      
      M = 1-is.na(XNA)
      
      ###############
      ####### Concatenating the data matrix and the mask
      ###############
      
      Y = cbind.data.frame(XNA,M[,1])
      
      
      ###############
      ####### Optimization: Mean imputation and estimation
      ###############
      
      X.mean <- as.matrix(ImputeMean(XNA))
      
      
      ###############
      ####### EM with modell
      ###############
      
      
      #Fonction 
      ThetaNew=Initialize_theta(ImputeMean0(XNA),r)#on initialise Theta avec une matrice en rang inférieur
      ThetabisNew=ThetaNew
      aNew=a-1
      bNew=b-1
      abisNew=a-1
      bbisNew=b-1
      
      diff=100
      ccompt<-0
      while(ccompt<20){
        ParamNew <- IterEM(Xtrue,X,XNA,ThetaNew,aNew,bNew,M,Ns,bruit,algo="soft",lam="Pred",nbcol=1)
        diff=ParamNew$diff
        print(diff)
        ThetaNew=ParamNew$ThetaNew
        aNew=ParamNew$a_initNew
        bNew=ParamNew$b_initNew
        conv1=ParamNew$conv
        ccompt=ccompt+1
      }
      
      
      diff=100
      ccompt2<-0
      while(ccompt2 < 20){
        ParamNew <- IterEM(Xtrue,X,XNA,ThetabisNew,abisNew,bbisNew,M,Ns,bruit,algo="FISTA",lam="Pred",nbcol=1)
        diff=ParamNew$diff
        ThetabisNew=ParamNew$ThetaNew
        abisNew=ParamNew$a_initNew
        bbisNew=ParamNew$b_initNew
        conv2=ParamNew$conv
        ccompt2=ccompt2+1
      }
      
      Theta=list()
      Theta[[1]]=ThetaNew
      Thetabis=list()
      Thetabis[[1]]=ThetabisNew
      aListNew=list()
      bListNew=list()
      aListNew[[1]]=aNew
      bListNew[[1]]=bNew
      abisListNew=list()
      bbisListNew=list()
      abisListNew[[1]]=abisNew
      bbisListNew[[1]]=bbisNew
      Tt=10
      conv3=c()
      conv4=c()
      pb <- txtProgressBar(min = 0, max = 50, title = 'Blabrf', label = 'fgrtgr')
      for (t in 1:Tt){
        ParamNew <- IterEM(Xtrue,X,XNA,Theta[[t]],aListNew[[t]],bListNew[[t]],M,Ns,bruit,algo="soft",lam="Pred",nbcol=1)
        Theta[[t+1]]=ParamNew$ThetaNew
        aListNew[[t+1]]=ParamNew$a_initNew
        bListNew[[t+1]]=ParamNew$b_initNew
        conv3=c(conv3,ParamNew$conv)
        
        ParamNew <- IterEM(Xtrue,X,XNA,Thetabis[[t]],abisListNew[[t]],bbisListNew[[t]],M,Ns,bruit,algo="FISTA",lam="Pred",nbcol=1)
        Thetabis[[t+1]]=ParamNew$ThetaNew
        abisListNew[[t+1]]=ParamNew$a_initNew
        bbisListNew[[t+1]]=ParamNew$b_initNew
        conv4=c(conv4,ParamNew$conv)
        
        setTxtProgressBar(pb = pb, value = t)
      }
      if(sum(conv3)==0){conv3=FALSE}else{conv3=TRUE}
      if(sum(conv4)==0){conv4=FALSE}else{conv4=TRUE}
      
      
      ThetaNew=Theta[[Tt]]
      ThetabisNew=Thetabis[[Tt]]
      
      ThetaNewTot=Initialize_theta(ImputeMean0(XNA),r)
      ThetabisNewTot=ThetaNewTot
      aNew=a-1
      bNew=b-1
      abisNew=a-1
      bbisNew=b-1
      
      diff=100
      ccompt3=0
      while(ccompt3<20){
        ParamNew <- IterEM(Xtrue,X,XNA,ThetaNewTot,aNew,bNew,M,Ns,bruit,algo="soft",lam="Tot",nbcol=1)
        diff=ParamNew$diff
        ThetaNewTot=ParamNew$ThetaNew
        aNew=ParamNew$a_initNew
        bNew=ParamNew$b_initNew
        conv5=ParamNew$conv
        ccompt3=ccompt3+1
      }
      
      
      diff=100
      ccompt4=0
      while(ccompt4<20){
        ParamNew <- IterEM(Xtrue,X,XNA,ThetabisNewTot,abisNew,bbisNew,M,Ns,bruit,algo="FISTA",lam="Tot",nbcol=1)
        diff=ParamNew$diff
        ThetabisNewTot=ParamNew$ThetaNew
        abisNew=ParamNew$a_initNew
        bbisNew=ParamNew$b_initNew
        conv6=ParamNew$conv
        ccompt4=ccompt4+1
      }
      
      
      Theta=list()
      Theta[[1]]=ThetaNewTot
      Thetabis=list()
      Thetabis[[1]]=ThetabisNewTot
      aListNew=list()
      bListNew=list()
      aListNew[[1]]=aNew
      bListNew[[1]]=bNew
      abisListNew=list()
      bbisListNew=list()
      abisListNew[[1]]=abisNew
      bbisListNew[[1]]=bbisNew
      Tt=10
      conv7=c()
      conv8=c()
      pb <- txtProgressBar(min = 0, max = 50, title = 'Blabrf', label = 'fgrtgr')
      for (t in 1:Tt){
        ParamNew <- IterEM(Xtrue,X,XNA,Theta[[t]],aListNew[[t]],bListNew[[t]],M,Ns,bruit,algo="soft",lam="Tot",nbcol=1)
        Theta[[t+1]]=ParamNew$ThetaNew
        aListNew[[t+1]]=ParamNew$a_initNew
        bListNew[[t+1]]=ParamNew$b_initNew
        conv7=c(conv3,ParamNew$conv)
        
        ParamNew <- IterEM(Xtrue,X,XNA,Thetabis[[t]],abisListNew[[t]],bbisListNew[[t]],M,Ns,bruit,algo="FISTA",lam="Tot",nbcol=1)
        Thetabis[[t+1]]=ParamNew$ThetaNew
        abisListNew[[t+1]]=ParamNew$a_initNew
        bbisListNew[[t+1]]=ParamNew$b_initNew
        conv8=c(conv4,ParamNew$conv)
        
        setTxtProgressBar(pb = pb, value = t)
      }
      if(sum(conv7)==0){conv7=FALSE}else{conv7=TRUE}
      if(sum(conv8)==0){conv8=FALSE}else{conv8=TRUE}
      
      ThetaNewTot=Theta[[Tt]]
      ThetabisNewTot=Thetabis[[Tt]]
      
      paste(c(conv1,conv2,conv3,conv4,conv5,conv6,conv7,conv8))
      if(conv1 & conv2 & conv3 & conv4 & conv5 & conv6 & conv7 & conv8){conv=TRUE}else(print("Pas de convergence"))
      
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
    list=c(0,1,2,3)
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
    list=c(0,1,2,3)
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
    ####### Algorithm Genevieve
    ###############
    
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
    
    ####Without mask
    RES=NULL
    RES2=NULL
    gridParam=seq(0,lambda0(X)*1.1, length = 100) 
    for(i in 1:length(gridParam)){
      X.FISTA2<- FISTANA(as.matrix(XNA),as.matrix(M),gridParam[i],bruit,alg="other")
      RES2[i]=MSE(X.FISTA2[,1:p]*(1-M),X*(1-M))
      RES[i]=MSE(X.FISTA2[,1:p],Xtrue)
    }
    X.FISTA2<- FISTANA(XNA,M,gridParam[which.min(RES)],bruit,alg="other")
    X.FISTA2.bis<- FISTANA(XNA,M,gridParam[which.min(RES2)],bruit,alg="other")
    
    ####With mask
    RES=NULL
    RES2=NULL
    gridParam=seq(0, lambda0(X)*1.1, length = 100) 
    for(i in 1:length(gridParam)){
      X.FISTA3<- FISTANA(Y,cbind.data.frame(M,M[,1:nbcol]),gridParam[i],bruit,alg="other")
      RES2[i]=MSE(X.FISTA3[,1:p]*(1-M),X*(1-M))
      RES[i]=MSE(X.FISTA3[,1:p],Xtrue)
    }
    X.FISTA3<- FISTANA(Y,cbind.data.frame(M,M[,1:nbcol]),gridParam[which.min(RES)],bruit,alg="other")
    X.FISTA3.bis<- FISTANA(Y,cbind.data.frame(M,M[,1:nbcol]),gridParam[which.min(RES2)],bruit,alg="other")
    
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
    ####### Calcul of MSE
    ###############
    
    #Premier mécanisme
    
    msePred=sapply(list(X.mean[,1:p]*(1-M),ThetaNew[,1:p]*(1-M),ThetabisNew[,1:p]*(1-M),ThetaNewTot[,1:p]*(1-M),ThetabisNewTot[,1:p]*(1-M),X1.soft[,1:p]*(1-M),X1.soft.bis[,1:p]*(1-M),X1.soft.2[,1:p]*(1-M),X1.soft.2.bis[,1:p]*(1-M),X.pca.true[,1:p]*(1-M),X.pca.true.bis[,1:p]*(1-M),X.pca.true.2[,1:p]*(1-M),X.pca.true.2.bis[,1:p]*(1-M),X.mimi[,1:p]*(1-M),X.mimi.bis[,1:p]*(1-M),X.FISTA2[,1:p]*(1-M),X.FISTA2.bis[,1:p]*(1-M),X.FISTA3[,1:p]*(1-M),X.FISTA3.bis[,1:p]*(1-M),X.rf[,1:p]*(1-M),X.rf.2[,1:p]*(1-M)),MSE,X2=X*(1-M))
    mseTrue=sapply(list(X.mean[,1:p],ThetaNew[,1:p],ThetabisNew[,1:p],ThetaNewTot[,1:p],ThetabisNewTot[,1:p],X1.soft[,1:p],X1.soft.bis[,1:p],X1.soft.2[,1:p],X1.soft.2.bis[,1:p],X.pca.true[,1:p],X.pca.true.bis[,1:p],X.pca.true.2[,1:p],X.pca.true.2.bis[,1:p],X.mimi[,1:p],X.mimi.bis[,1:p],X.FISTA2[,1:p],X.FISTA2.bis[,1:p],X.FISTA3[,1:p],X.FISTA3.bis[,1:p],X.rf[,1:p],X.rf.2[,1:p]),MSE,X2=Xtrue)
    
    names=c("Imputemean","modelsoftPred","modelFISTAPred","modelsoftTot","modelFISTATot","softTot","softPred","softmaskTot","softmaskPred","PCATot","PCAPred","PCAmaskTot","PCAmaskPred","mimiTot","mimiPred","FISTA","FISTAPred","FISTAmask","FISTAmaskPred","randomforest","randomforestmask")
    cbind(msePred,mseTrue,names)
  }
  
  return(results.list)
}







