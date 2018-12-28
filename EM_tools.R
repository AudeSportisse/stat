library(softImpute)
library(randomForest)
library(missForest)
source("FISTA_tools.R", local = TRUE)
source("General_tools.R", local = TRUE)


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


######
# Name: IterEM
# Date: 27/12/2018
# Description: It computes one iteration of the EM algorithm for the univariate and multivariate cases. In particular, the mechanism parameters are the same ones for all the variables and the missing variables are the first nbcol variables.
# This function returns: the difference between the last matrix iteration and the new one, the new matrix, the new mechanisms parameters and if the GLM algorithm has converged.
# Arguments: 
  #Xtrue: parameter matrix. 
  #X: data matrix.
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


IterEM <- function(Xtrue,
                   X,
                   XNA,
                   Thet,
                   a,
                   b,
                   M,
                   Ns,
                   noise,
                   algo,
                   lam,
                   nbcol) {
  a_initOld = a
  b_initOld = b
  ThetaOld = Thet
  
  SIR <- function(Theta, XNA, a, b, sigma, i, jind) {
    Nn = Ns * 1000
    sim = rnorm(Nn, mean = Theta[i, jind], sd = sigma)
    
    RapportRej <- function(x) {
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
    poids = RapportRej(sim)
    res <- sample(sim, Ns, prob = poids, replace = TRUE)
    return(res)
  }
  
  sigma = sqrt(noise)
  XNA_remp = XNA
  data_remp <- list()
  nbNA = 0
  for (i in 1:nrow(XNA)) {
    nbNAcol = 0
    for (j in 1:ncol(XNA)) {
      if (is.na(XNA[i, j])) {
        ech = SIR(ThetaOld, XNA, a_initOld, b_initOld, sigma, i, j)
        XNA_remp[i, j] = mean(ech)
        if (nbNAcol == 0) {
          data_remp[[i]] = ech
        } else{
          data_remp[[i]] = append(data_remp[[i]], ech)
        }
        nbNAcol = nbNAcol + 1
      }
    }
    nbNA = nbNA + nbNAcol
  }
  
  if (algo == "soft") {
    RES = NULL
    gridlambda1 = seq(0, lambda0(X) * 1.1, length = 300)
    for (i in 1:length(gridlambda1)) {
      fit1 = softImpute(
        as.matrix(XNA_remp),
        rank = min(dim(XNA_remp)) - 1,
        lambda = gridlambda1[i],
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
      if (lam == "Tot") {
        RES[i] = MSE(ThetaNew, Xtrue)
      } else{
        RES[i] = MSE(ThetaNew * (1 - M), X * (1 - M))
      }
    }
    fit1 = softImpute(
      as.matrix(XNA_remp),
      rank = min(dim(XNA_remp)) - 1,
      lambda = gridlambda1[which.min(RES)],
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
    RES = NULL
    gridParam = seq(0, lambda0(X) * 1.1, length = 100)
    for (i in 1:length(gridParam)) {
      X.FISTA2 = FISTA0(XNA_remp, gridParam[i], noise)
      if (lam == "Tot") {
        RES[i] = MSE(X.FISTA2, Xtrue)
      } else{
        RES[i] = MSE(X.FISTA2 * (1 - M), X * (1 - M))
      }
    }
    ThetaNew = FISTA0(XNA_remp, gridParam[which.min(RES)], noise)
  }
  
  M_concat = rep(as.vector(M[,1]), Ns )
  X_concat = matrix(0, nrow = Ns * n , ncol = 1)
  for (i in 1:nrow(XNA)) {
    for (k in (1:Ns)) {
      if (is.na(XNA[i, 1]) == TRUE) {
        X_concat[n * (k - 1) + i, 1] = data_remp[[i]][k]
      } else{
        X_concat[n * (k - 1) + i, 1] = XNA[i, 1]
      }
    }
  }
  
  GLM = glm(M_concat ~ X_concat, family = binomial)
  a_initNew = -GLM$coeff[2]
  b_initNew = -GLM$coeff[1] / GLM$coeff[2]
  conv = GLM$converged
  
  diff = norm(ThetaNew - ThetaOld, type = "F") / (norm(ThetaOld, type =
                                                         "F") + 10 ^ -3)
  
  return(
    list(
      diff = diff,
      ThetaNew = ThetaNew,
      a_initNew = a_initNew,
      b_initNew = b_initNew,
      conv = conv
    )
  )
  
}


######
# Name: IterEM_Bivariate
# Date: 27/12/2018
# Description: It computes one iteration of the EM algorithm for the bivariate case. In particular, the missing variables are the first one and the colbis one.
# This function returns: the difference between the last matrix iteration and the new one, the new matrix, the new mechanisms parameters and if the GLM algorithm has converged.
# Arguments: 
  #Xtrue: parameter matrix. 
  #X: data matrix.
  #XNA: data matrix X containing missing values.  
  #Thet: old matrix iteration.
  #a, b, a2, b2: old mechanism parameters. 
  #M: missing-data matrix.
  #Ns: number of Monte Carlo simulations.
  #noise: sigma^2.
  #algo: "soft" or "FISTA" depending on the algorithm used for M-step. 
  #lam: "Tot" or "Pred" depending on how the lambda is chosen. 
  #colbis: the second missing variable.
#####


IterEM_Bivariate <-
  function(Xtrue,
           X,
           XNA,
           Thet,
           a,
           b,
           a2,
           b2,
           M,
           Ns,
           noise,
           algo,
           lam,
           colbis) {
    
    a_initOld = a
    b_initOld = b
    a_initOld2 = a2
    b_initOld2 = b2
    ThetaOld = Thet
    
    SIR = function(Theta, XNA, a, b, sigma, i, jind) {
      Nn = Ns * 1000
      sim = rnorm(Nn, mean = Theta[i, jind], sd = sigma)
      
      RapportRej <- function(x) {
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
      poids = RapportRej(sim)
      res = sample(sim, Ns, prob = poids, replace = TRUE)
      return(res)
    }
    
    sigma = sqrt(noise)
    XNA_remp = XNA
    data_remp = list()
    data_remp2 = list()
    nbNA = 0
    for (i in 1:nrow(XNA)) {
      if (is.na(XNA[i, 1]) & is.na(XNA[i, colbis])) {
        ech = SIR(ThetaOld, XNA, a_initOld, b_initOld, sigma, i, 1)
        ech2 =
          SIR(ThetaOld, XNA, a_initOld2, b_initOld2, sigma, i, colbis)
        XNA_remp[i, 1] = mean(ech)
        XNA_remp[i, colbis] = mean(ech2)
      } else if (is.na(XNA[i, 1])) {
        ech = SIR(ThetaOld, XNA, a_initOld, b_initOld, sigma, i, 1)
        ech2 = c()
        XNA_remp[i, 1] = mean(ech)
      } else if (is.na(XNA[i, colbis])) {
        ech = c()
        ech2 =
          SIR(ThetaOld, XNA, a_initOld2, b_initOld2, sigma, i, colbis)
        XNA_remp[i, colbis] = mean(ech2)
      } else{
        ech = c()
        ech2 = c()
      }
      data_remp[[i]] = ech
      data_remp2[[i]] = ech2
    }
    
    ###Minimisation pour Theta avec softImpute
    if (algo == "soft") {
      RES = NULL
      gridlambda1 = seq(0, lambda0(X) * 1.1, length = 300)
      for (i in 1:length(gridlambda1)) {
        fit1 = softImpute(
          as.matrix(XNA_remp),
          rank = min(dim(XNA_remp)) - 1,
          lambda = gridlambda1[i],
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
        if (lam == "Tot") {
          RES[i] = MSE(ThetaNew, Xtrue)
        } else{
          RES[i] = MSE(ThetaNew * (1 - M), X * (1 - M))
        }
      }
      fit1 = softImpute(
        as.matrix(XNA_remp),
        rank = min(dim(XNA_remp)) - 1,
        lambda = gridlambda1[which.min(RES)],
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
      RES = NULL
      gridParam = seq(0, lambda0(X) * 1.1, length = 100)
      for (i in 1:length(gridParam)) {
        X.FISTA2 = FISTA0(XNA_remp, gridParam[i], noise)
        if (lam == "Tot") {
          RES[i] = MSE(X.FISTA2, Xtrue)
        } else{
          RES[i] = MSE(X.FISTA2 * (1 - M), X * (1 - M))
        }
      }
      ThetaNew = FISTA0(XNA_remp, gridParam[which.min(RES)], noise)
    }
    
    
    M_concat = rep(M[, 1], Ns)
    X_concat = matrix(0, nrow = Ns * n, ncol = 1)
    for (k in 1:Ns) {
      for (i in 1:nrow(XNA)) {
        if (is.na(XNA[i, 1]) == TRUE) {
          X_concat[n * (k - 1) + i, 1] = data_remp[[i]][k]
        } else{
          X_concat[n * (k - 1) + i, 1] = XNA[i, 1]
        }
      }
    }
    
    GLM = glm(M_concat ~ X_concat, family = binomial)
    a_initNew = -GLM$coeff[2]
    b_initNew = -GLM$coeff[1] / GLM$coeff[2]
    conv = GLM$converged
    
    M_concat2 = rep(M[, colbis], Ns)
    X_concat2 = matrix(0, nrow = n * Ns, ncol = 1)
    for (k in 1:Ns) {
      for (i in 1:nrow(XNA)) {
        if (is.na(XNA[i, colbis]) == TRUE) {
          X_concat2[n * (k - 1) + i, 1] = data_remp2[[i]][k]
        } else{
          X_concat2[n * (k - 1) + i, 1] = XNA[i, colbis]
        }
      }
    }
    
    GLM = glm(M_concat2 ~ X_concat2, family = binomial)
    a_initNew2 = -GLM$coeff[2]
    b_initNew2 = -GLM$coeff[1] / GLM$coeff[2]
    convbis = GLM$converged
    
    if (convbis == FALSE) {
      conv = convbis
    }
    
    diff = norm(ThetaNew - ThetaOld, type = "F") / (norm(ThetaOld, type =
                                                           "F") + 10 ^ -3)

    return(
      list(
        diff = diff,
        ThetaNew = ThetaNew,
        a_initNew = a_initNew,
        b_initNew = b_initNew,
        a_initNew2 = a_initNew2,
        b_initNew2 = b_initNew2,
        conv = conv
      )
    )
    
  }
