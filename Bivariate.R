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


ComparMNAR_Bivariate <- function(Xtrue,a,b,a2,b2,r,bruit,Ns,colbis,m1,m2,nbsim){
  
  nbcol<-2
  p<-ncol(Xtrue)
  n<-nrow(Xtrue)
  
  results.list = foreach (ksim = 1:nbsim, .combine = "rbind")  %dopar% {
    source("EM_tools.R",local=TRUE)
    source("General_tools.R",local=TRUE)
    source("PCA_tools.R",local=TRUE)
    source("FISTA_tools.R",local=TRUE)
    source("CodeMimi.R",local=TRUE)
    
    print(paste("it globale",ksim))
    
    ###############
    ####### Introduction of missing values: Logistic Regression
    ###############
    
    conv=FALSE
    while(!conv){
      
      set.seed(ksim)
      X=Xtrue+matrix(data=rnorm(n*p,0,sqrt(bruit)),ncol=p)
      
      #MNAR with logistic regression
      select_prob <- function(x,a,b){ #probability of selecting coordinate Xij
        res=1/(1+exp(-a*(x-b)))
        return(res)
      }
      
      XNA=X
      for (i in 1:n){
        p1 <- select_prob(X[i,m1],a,b)*select_prob(X[i,m2],a2,b2)  
        p2 <- select_prob(X[i,m1],a,b)*(1-select_prob(X[i,m2],a2,b2))  
        p3 <- (1-select_prob(X[i,m1],a,b))*select_prob(X[i,m2],a2,b2)  
        p4 <- (1-select_prob(X[i,m1],a,b))*(1-select_prob(X[i,m2],a2,b2)) 
        ind<-sample(1:4,1,prob=c(p1,p2,p3,p4))
        if(ind==1){
          XNA[i,1]=NA
          XNA[i,colbis]=NA
        }else if(ind==2){
          XNA[i,1]=NA
        }else if(ind==3){
          XNA[i,colbis]=NA
        }
      }
      print(paste("fin simu",ksim))
      
      
      
      ###############
      ####### Missing-data pattern: MECHANISM 2
      ###############
      
      M = 1-is.na(XNA)
      
      ###############
      ####### Concatenating the data matrix and the mask: MECHANISM 2
      ###############
      
      Y = cbind.data.frame(XNA,M[,c(1,colbis)])
      
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
      aNew2=a2-1
      bNew2=b2-1
      abisNew2=a2-1
      bbisNew2=b2-1
      
      diff=100
      ccompt<-0
      while(diff>10^-2 & ccompt<50){
        ParamNew <- IterEM_Bivariate(Xtrue,X,XNA,ThetaNew,aNew,bNew,aNew2,bNew2,M,Ns,bruit,algo="soft",lam="Pred",colbis)
        diff=ParamNew$diff
        ThetaNew=ParamNew$ThetaNew
        aNew=ParamNew$a_initNew
        bNew=ParamNew$b_initNew
        aNew2=ParamNew$a_initNew2
        bNew2=ParamNew$b_initNew2
        conv1=ParamNew$conv
        ccompt=ccompt+1
      }
      print(paste("ccompt",ccompt,ksim))
      
      diff=100
      ccompt2<-0
      while(diff>10^-2 & ccompt2<50){
        ParamNew <- IterEM_Bivariate(Xtrue,X,XNA,ThetabisNew,abisNew,bbisNew,abisNew2,bbisNew2,M,Ns,bruit,algo="FISTA",lam="Pred",colbis)
        diff=ParamNew$diff
        ThetabisNew=ParamNew$ThetaNew
        abisNew=ParamNew$a_initNew
        bbisNew=ParamNew$b_initNew
        abisNew2=ParamNew$a_initNew2
        bbisNew2=ParamNew$b_initNew2
        conv2=ParamNew$conv
        ccompt2=ccompt2+1
      }
      print(paste("ccompt2",ccompt2,ksim))
      
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
      aListNew2=list()
      bListNew2=list()
      aListNew2[[1]]=aNew2
      bListNew2[[1]]=bNew2
      abisListNew2=list()
      bbisListNew2=list()
      abisListNew2[[1]]=abisNew2
      bbisListNew2[[1]]=bbisNew2
      Tt=10
      conv3=c()
      conv4=c()
      for (t in 1:Tt){
        ParamNew <- IterEM_Bivariate(Xtrue,X,XNA,Theta[[t]],aListNew[[t]],bListNew[[t]],aListNew2[[t]],bListNew2[[t]],M,Ns,bruit,algo="soft",lam="Pred",colbis)
        Theta[[t+1]]=ParamNew$ThetaNew
        aListNew[[t+1]]=ParamNew$a_initNew
        bListNew[[t+1]]=ParamNew$b_initNew
        aListNew2[[t+1]]=ParamNew$a_initNew2
        bListNew2[[t+1]]=ParamNew$b_initNew2
        conv3=c(conv3,ParamNew$conv)
        
        ParamNew <- IterEM_Bivariate(Xtrue,X,XNA,Thetabis[[t]],abisListNew[[t]],bbisListNew[[t]],abisListNew2[[t]],bbisListNew[[t]],M,Ns,bruit,algo="FISTA",lam="Pred",colbis)
        Thetabis[[t+1]]=ParamNew$ThetaNew
        abisListNew[[t+1]]=ParamNew$a_initNew
        bbisListNew[[t+1]]=ParamNew$b_initNew
        abisListNew2[[t+1]]=ParamNew$a_initNew2
        bbisListNew2[[t+1]]=ParamNew$b_initNew2
        conv4=c(conv4,ParamNew$conv)
    
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
      aNew2=a2-1
      bNew2=b2-1
      abisNew2=a2-1
      bbisNew2=b2-1
      
      diff=100
      ccompt3=0
      while(diff>10^-2 & ccompt3<50){
        ParamNew <- IterEM_Bivariate(Xtrue,X,XNA,ThetaNewTot,aNew,bNew,aNew2,bNew2,M,Ns,bruit,algo="soft",lam="Tot",colbis)
        diff=ParamNew$diff
        ThetaNewTot=ParamNew$ThetaNew
        aNew=ParamNew$a_initNew
        bNew=ParamNew$b_initNew
        aNew2=ParamNew$a_initNew2
        bNew2=ParamNew$b_initNew2
        conv5=ParamNew$conv
        ccompt3=ccompt3+1
      }
      print(paste("ccompt3",ccompt3,ksim))
      
      diff=100
      ccompt4=0
      while(diff>10^-2 & ccompt4<50){
        ParamNew <- IterEM_Bivariate(Xtrue,X,XNA,ThetabisNewTot,abisNew,bbisNew,abisNew2,bbisNew2,M,Ns,bruit,algo="FISTA",lam="Tot",colbis)
        diff=ParamNew$diff
        ThetabisNewTot=ParamNew$ThetaNew
        abisNew=ParamNew$a_initNew
        bbisNew=ParamNew$b_initNew
        abisNew2=ParamNew$a_initNew2
        bbisNew2=ParamNew$b_initNew2
        conv6=ParamNew$conv
        ccompt4=ccompt4+1
      }
      print(paste("ccompt4",ccompt4,ksim))
      
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
      aListNew2=list()
      bListNew2=list()
      aListNew2[[1]]=aNew2
      bListNew2[[1]]=bNew2
      abisListNew2=list()
      bbisListNew2=list()
      abisListNew2[[1]]=abisNew2
      bbisListNew2[[1]]=bbisNew2
      Tt=10
      conv7=c()
      conv8=c()
      for (t in 1:Tt){
        ParamNew <- IterEM_Bivariate(Xtrue,X,XNA,Theta[[t]],aListNew[[t]],bListNew[[t]],aListNew2[[t]],bListNew2[[t]],M,Ns,bruit,algo="soft",lam="Tot",colbis)
        Theta[[t+1]]=ParamNew$ThetaNew
        aListNew[[t+1]]=ParamNew$a_initNew
        bListNew[[t+1]]=ParamNew$b_initNew
        aListNew2[[t+1]]=ParamNew$a_initNew2
        bListNew2[[t+1]]=ParamNew$b_initNew2
        conv7=c(conv3,ParamNew$conv)
        
        ParamNew <- IterEM_Bivariate(Xtrue,X,XNA,Thetabis[[t]],abisListNew[[t]],bbisListNew[[t]],abisListNew2[[t]],bbisListNew2[[t]],M,Ns,bruit,algo="FISTA",lam="Tot",colbis)
        Thetabis[[t+1]]=ParamNew$ThetaNew
        abisListNew[[t+1]]=ParamNew$a_initNew
        bbisListNew[[t+1]]=ParamNew$b_initNew
        abisListNew2[[t+1]]=ParamNew$a_initNew2
        bbisListNew2[[t+1]]=ParamNew$b_initNew2
        conv8=c(conv4,ParamNew$conv)
        
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
    list=c(0,1,2,3,4,5,6,7,8,9,10)
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
    list=c(0,1,2,3,4,5,6,7,8,9,10)
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
    ####### Mimi algorithm 
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
    ####### MSE
    ###############
    
    msePred=sapply(list(X.mean[,1:p]*(1-M),ThetaNew[,1:p]*(1-M),ThetabisNew[,1:p]*(1-M),ThetaNewTot[,1:p]*(1-M),ThetabisNewTot[,1:p]*(1-M),X1.soft[,1:p]*(1-M),X1.soft.bis[,1:p]*(1-M),X1.soft.2[,1:p]*(1-M),X1.soft.2.bis[,1:p]*(1-M),X.pca.true[,1:p]*(1-M),X.pca.true.bis[,1:p]*(1-M),X.pca.true.2[,1:p]*(1-M),X.pca.true.2.bis[,1:p]*(1-M),X.mimi[,1:p]*(1-M),X.mimi.bis[,1:p]*(1-M),X.FISTA2[,1:p]*(1-M),X.FISTA2.bis[,1:p]*(1-M),X.FISTA3[,1:p]*(1-M),X.FISTA3.bis[,1:p]*(1-M),X.rf[,1:p]*(1-M),X.rf.2[,1:p]*(1-M)),MSE,X2=X*(1-M))
    mseTrue=sapply(list(X.mean[,1:p],ThetaNew[,1:p],ThetabisNew[,1:p],ThetaNewTot[,1:p],ThetabisNewTot[,1:p],X1.soft[,1:p],X1.soft.bis[,1:p],X1.soft.2[,1:p],X1.soft.2.bis[,1:p],X.pca.true[,1:p],X.pca.true.bis[,1:p],X.pca.true.2[,1:p],X.pca.true.2.bis[,1:p],X.mimi[,1:p],X.mimi.bis[,1:p],X.FISTA2[,1:p],X.FISTA2.bis[,1:p],X.FISTA3[,1:p],X.FISTA3.bis[,1:p],X.rf[,1:p],X.rf.2[,1:p]),MSE,X2=Xtrue)
    
    names=c("Imputemean","modelsoftPred","modelFISTAPred","modelsoftTot","modelFISTATot","softTot","softPred","softmaskTot","softmaskPred","PCATot","PCAPred","PCAmaskTot","PCAmaskPred","mimiTot","mimiPred","FISTA","FISTAPred","FISTAmask","FISTAmaskPred","randomforest","randomforestmask")
    cbind(msePred,mseTrue,names)
}
  
  return(results.list)
}




