mimi.lr<-function (y, var.type = c("gaussian", "binary", "categorical", 
                          "poisson"), lambda1, maxit = 100, theta0 = NULL, thresh = 1e-05, 
          trace.it = F, lambda1.max = NULL, length = 40, upper = 12, 
          lower = -12, offset = F, scale = F, max.rank = 20, wmax = NULL) 
{
  yy <- y
  nlevel = rep(1, ncol(y))
  if (sum(var.type == "binary") > 0) {
    for (j in 1:sum(var.type == "binary")) {
      y <- data.frame(y)
      y[, which(var.type == "binary")[j]] <- as.factor(y[, 
                                                         which(var.type == "binary")[j]])
      y[, which(var.type == "binary")[j]] <- sapply(as.character(y[, 
                                                                   which(var.type == "binary")[j]]), function(t) if (t %in% 
                                                                                                                     c(levels(y[, which(var.type == "binary")[j]])[1])) 
                                                                     1
                                                    else if (is.na(t)) 
                                                      NA
                                                    else 0)
    }
  }
  if (sum(var.type == "categorical") > 0) {
    y <- data.frame(y)
    for (j in 1:sum(var.type == "categorical")) {
      y[, which(var.type == "categorical")[j]] <- as.factor(y[, 
                                                              which(var.type == "categorical")[j]])
    }
    tab <- tab.disjonctif.prop(y[, var.type == "categorical"])
    tab[((tab > 0) * (tab < 1) == 1)] <- NA
    nlevel[var.type == "categorical"] <- sapply(which(var.type == 
                                                        "categorical"), function(t) nlevels(y[, t]))
    colnames(tab) <- paste(rep(colnames(y[, var.type == "categorical"]), 
                               nlevel[var.type == "categorical"]), ".", c(sapply(nlevel[var.type == 
                                                                                          "categorical"], function(k) 1:k)), sep = "")
    y2 <- rep(0, n)
    x2 <- rep(0, ncol(x))
    vt2 <- NULL
    count <- 1
    for (j in 1:ncol(y)) {
      if (!(var.type[j] == "categorical")) {
        y2 <- cbind(y2, y[j])
        vt2 <- c(vt2, var.type[j])
      }
      else {
        y2 <- cbind(y2, tab[, count:(count + nlevels(y[, 
                                                       j]) - 1)])
        vt2 <- c(vt2, rep("categorical", nlevels(y[, 
                                                   j])))
      }
    }
    y <- y2[, 2:ncol(y2)]
  }
  else vt2 <- var.type
  y <- as.matrix(y)
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  max.rank <- min(max.rank, min(n, p) - 1)
  y <- as.matrix(y)
  y <- matrix(as.numeric(y), nrow = n)
  if (is.null(wmax)) 
    wmax = 2 * max(y, na.rm = T)
  omega <- !is.na(y)
  if (is.null(lambda1.max)) 
    lambda1.max <- 1000 * lambda1
  lambda1.grid.log <- seq(log(lambda1.max), log(lambda1), length.out = length)
  if (is.null(theta0)) 
    theta0 <- matrix(rep(0, n * p), nrow = n)
  theta <- theta0
  iter <- 1
  list.theta <- list()
  list.theta[[1]] <- theta0
  list.res <- list()
  list.res[[1]] <- list(y = y, theta = theta0, objective = 0, 
                        y.imputed = apply(y, c(1, 2), function(xx) {
                          if (is.na(xx)) 0 else xx
                        }))
  for (i in 2:length) {
    res <- rwls_l1_nuc.lr(y, var.type = var.type, nlevel = nlevel, 
                          lambda1 = exp(lambda1.grid.log[i]), maxit = maxit, 
                          upper = upper, lower = lower, theta0 = list.theta[[iter]], 
                          thresh = thresh, trace.it = trace.it, offset = offset, 
                          max.rank = max.rank, vt2 = vt2, scale = scale, wmax = wmax)
    iter <- iter + 1
    list.res[[iter]] <- res
    list.theta[[iter]] <- res$theta
  }
  return(list(lambda1.grid = lambda1.grid.log, list.res = list.res, 
              list.theta = list.theta))
}



rwls_l1_nuc.lr<-function (y, var.type, lambda1, nlevel = NULL, maxit = 100, upper = NULL, 
                          lower = NULL, theta0 = NULL, thresh = 1e-05, trace.it = F, 
                          offset = F, scale = F, max.rank = 20, vt2 = NULL, wmax = NULL) 
{
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  y <- as.matrix(y)
  y <- matrix(as.numeric(y), nrow = n)
  omega <- !is.na(y)
  if (is.null(upper)) 
    upper <- 2 * max(abs(y), na.rm = T)
  if (is.null(lower)) 
    lower <- -2 * max(abs(y), na.rm = T)
  if (is.null(theta0)) 
    theta0 <- matrix(rep(0, n * p), nrow = n)
  theta <- theta0
  param <- theta
  objective <- NULL
  y0 <- y
  error <- 1
  iter <- 0
  while ((error > thresh) && (iter < maxit)) {
    iter <- iter + 1
    theta.tmp <- theta
    param.tmp <- param
    yv <- quad_approx(y0, param, var.type, nlevel, wmax)
    ytilde <- yv$ytilde
    vtilde2 <- yv$vtilde2
    if (scale == T) {
      scaling <- apply(y, 2, sd, na.rm = T)
      ytilde <- sweep(ytilde, 2, scaling, "/")
    }
    lambda1w <- lambda1/max(vtilde2)
    vtilde2 <- vtilde2/max(vtilde2)
    svd_theta <- wlra(x = ytilde, w = vtilde2, lambda = lambda1w, 
                      x0 = NULL, thresh = 0.1 * thresh, rank.max = max.rank, 
                      maxit = 1000)
    u <- svd_theta$u
    d <- svd_theta$d
    v <- svd_theta$v
    if (is.null(dim(u))) {
      theta <- d * u %*% t(v)
    }
    else {
      theta <- u %*% diag(d) %*% t(v)
    }
    if (sum(var.type == "categorical")) 
      theta[, vt2 == "categorical"] <- sweep(theta[, vt2 == 
                                                     "categorical"], 1, rowMeans(theta[, vt2 == "categorical"]))
    if (scale) {
      theta <- sweep(theta, 2, scaling, "*")
    }
    res <- bls.lr(y0, theta, theta.tmp, b = 0.5, lambda1, 
                  vt2, thresh)
    theta <- res$theta
    param <- theta
    objective <- c(objective, min(.Machine$double.xmax, res$objective))
    if (iter == 1) {
      error <- 1
    }
    else {
      if (objective[iter] == .Machine$double.xmax) {
        error <- 1
      }
      else error <- abs(objective[iter] - objective[iter - 
                                                      1])/abs(objective[iter])
    }
    if (all(var.type == "gaussian")) 
      error <- 0
    if (trace.it) {
      print(paste("iter ", iter, ": error ", error, " - objective: ", 
                  objective[iter]))
    }
  }
  if (error < thresh) 
    cvg = T
  else cvg = F
  y <- theta
  y[, vt2 == "poisson"] <- exp(y[, vt2 == "poisson"])
  y[, vt2 == "binary"] <- exp(y[, vt2 == "binary"])/(1 + exp(y[, 
                                                               vt2 == "binary"]))
  if (sum(var.type == "categorical") > 0) {
    for (j in 1:sum(var.type == "categorical")) {
      count <- sum(nlevel[1:(which(var.type == "categorical")[j] - 
                               1)])
      y[, (count + 1):(count + nlevel[which(vt2 == "categorical")[j]])] <- t(sapply(1:n, 
                                                                                    function(i) exp(y[i, (count + 1):(count + nlevel[which(vt2 == 
                                                                                                                                             "categorical")[j]])])/sum(exp(y[i, (count + 
                                                                                                                                                                                   1):(count + nlevel[which(vt2 == "categorical")[j]])]))))
    }
  }
  y.imputed <- y0
  y.imputed[is.na(y.imputed)] <- y[is.na(y.imputed)]
  return(list(y = y0, theta = theta, objective = objective, 
              y.imputed = y.imputed))
}


wlra<-function (x, w = NULL, lambda = 0, x0 = NULL, thresh = 1e-05, 
          maxit = 1000, rank.max = NULL) 
{
  d <- dim(x)
  n <- d[1]
  p <- d[2]
  if (is.null(w)) 
    w = matrix(rep(1, n * p), nrow = n)
  if (is.null(rank.max)) 
    rank.max <- min(n, p) - 1
  else rank.max <- min(rank.max, min(n, p) - 1)
  if (is.null(x0)) 
    x0 <- matrix(rep(0, n * p), nrow = n)
  xnas <- is.na(x)
  omega <- 1 * (!xnas)
  nz = n * p - sum(xnas)
  xfill <- x
  xfill[xnas] <- 0
  xfill <- omega * w * xfill + (1 - omega * w) * x0
  xhat <- x0
  x0 <- xfill
  iter <- 0
  error <- 100
  svd.xfill = svd(xfill)
  while ((error > thresh) && (iter < maxit)) {
    iter <- iter + 1
    svd.old = svd(xhat)
    xhat.old <- xhat
    d = svd.xfill$d
    d = pmax(d - lambda, 0)
    J <- min(rank.max, length(d))
    xhat <- svd.xfill$u[, seq(J)] %*% (d[seq(J)] * t(svd.xfill$v[, 
                                                                 seq(J)]))
    xfill <- omega * w * x0 + (1 - omega * w) * xhat
    svd.xfill = svd(xfill)
    denom <- sum(w * (x - xhat.old)^2, na.rm = T) + lambda * 
      sum(svd.old$d[seq(J)])
    error <- abs(sum(w * (x - xhat)^2, na.rm = T) + lambda * 
                   sum(d[seq(J)]) - denom)/denom
  }
  if (error < thresh) 
    cvg = T
  else cvg = F
  d <- svd.xfill$d[seq(J)]
  d = pmax(svd.xfill$d[seq(J)] - lambda, 0)
  J = min(sum(d > 0) + 1, J)
  svd.xfill = list(u = svd.xfill$u[, seq(J)], d = d[seq(J)], 
                   v = svd.xfill$v[, seq(J)])
  if (iter == maxit) 
    warning(paste("Convergence not achieved by", maxit, "iterations"))
  return(list(d = svd.xfill$d, u = svd.xfill$u, v = svd.xfill$v, 
              cvg = cvg, iter = iter))
}


mimi<-function (y, var.type = c("gaussian", "binary", "categorical", 
                          "poisson"), lambda1, maxit = 100, theta0 = NULL, thresh = 1e-05, 
          trace.it = F, lambda1.max = NULL, length = 40, upper = 12, 
          lower = -12, offset = F, scale = F, max.rank = 20, wmax = NULL) 
{
  yy <- y
  nlevel = rep(1, ncol(y))
  if (sum(var.type == "binary") > 0) {
    for (j in 1:sum(var.type == "binary")) {
      y <- data.frame(y)
      y[, which(var.type == "binary")[j]] <- as.factor(y[, 
                                                         which(var.type == "binary")[j]])
      y[, which(var.type == "binary")[j]] <- sapply(as.character(y[, 
                                                                   which(var.type == "binary")[j]]), function(t) if (t %in% 
                                                                                                                     c(levels(y[, which(var.type == "binary")[j]])[1])) 
                                                                     1
                                                    else if (is.na(t)) 
                                                      NA
                                                    else 0)
    }
  }
  if (sum(var.type == "categorical") > 0) {
    y <- data.frame(y)
    for (j in 1:sum(var.type == "categorical")) {
      y[, which(var.type == "categorical")[j]] <- as.factor(y[, 
                                                              which(var.type == "categorical")[j]])
    }
    tab <- tab.disjonctif.prop(y[, var.type == "categorical"])
    tab[((tab > 0) * (tab < 1) == 1)] <- NA
    nlevel[var.type == "categorical"] <- sapply(which(var.type == 
                                                        "categorical"), function(t) nlevels(y[, t]))
    colnames(tab) <- paste(rep(colnames(y[, var.type == "categorical"]), 
                               nlevel[var.type == "categorical"]), ".", c(sapply(nlevel[var.type == 
                                                                                          "categorical"], function(k) 1:k)), sep = "")
    y2 <- rep(0, n)
    x2 <- rep(0, ncol(x))
    vt2 <- NULL
    count <- 1
    for (j in 1:ncol(y)) {
      if (!(var.type[j] == "categorical")) {
        y2 <- cbind(y2, y[j])
        vt2 <- c(vt2, var.type[j])
      }
      else {
        y2 <- cbind(y2, tab[, count:(count + nlevels(y[, 
                                                       j]) - 1)])
        vt2 <- c(vt2, rep("categorical", nlevels(y[, 
                                                   j])))
      }
    }
    y <- y2[, 2:ncol(y2)]
  }
  else vt2 <- var.type
  y <- as.matrix(y)
  d <- dim(y)
  n <- d[1]
  p <- d[2]
  max.rank <- min(max.rank, min(n, p) - 1)
  y <- as.matrix(y)
  y <- matrix(as.numeric(y), nrow = n)
  if (is.null(wmax)) 
    wmax = 2 * max(y, na.rm = T)
  omega <- !is.na(y)
  if (is.null(lambda1.max)) 
    lambda1.max <- 1000 * lambda1
  lambda1.grid.log <- seq(log(lambda1.max), log(lambda1), length.out = length)
  if (is.null(theta0)) 
    theta0 <- matrix(rep(0, n * p), nrow = n)
  theta <- theta0
  iter <- 1
  list.theta <- list()
  list.theta[[1]] <- theta0
  list.res <- list()
  list.res[[1]] <- list(y = y, theta = theta0, objective = 0, 
                        y.imputed = apply(y, c(1, 2), function(xx) {
                          if (is.na(xx)) 0 else xx
                        }))
  for (i in 2:length) {
    res <- rwls_l1_nuc.lr(y, var.type = var.type, nlevel = nlevel, 
                          lambda1 = exp(lambda1.grid.log[i]), maxit = maxit, 
                          upper = upper, lower = lower, theta0 = list.theta[[iter]], 
                          thresh = thresh, trace.it = trace.it, offset = offset, 
                          max.rank = max.rank, vt2 = vt2, scale = scale, wmax = wmax)
    iter <- iter + 1
    list.res[[iter]] <- res
    list.theta[[iter]] <- res$theta
  }
  return(list(lambda1.grid = lambda1.grid.log, list.res = list.res, 
              list.theta = list.theta))
}

mimi<-function (y, model = c("groups", "covariates", "low-rank"), x = NULL, 
          groups = NULL, var.type = c("gaussian", "binary", "categorical", 
                                      "poisson"), lambda1, lambda2, maxit = 100, mu0 = NULL, 
          alpha0 = NULL, theta0 = NULL, thresh = 1e-06, trace.it = F, 
          lambda1.max = NULL, lambda2.max = NULL, length = 20, upper = 12, 
          lower = -12, offset = T, scale = T, max.rank = 20, wmax = NULL) 
{
  if (model == "groups") {
    return(mimi.multi(y = y, groups = groups, var.type = var.type, 
                      lambda1 = lambda1, lambda2 = lambda2, maxit = maxit, 
                      mu0 = mu0, alpha0 = alpha0, theta0 = theta0, thresh = thresh, 
                      trace.it = trace.it, lambda1.max = lambda1.max, lambda2.max = lambda2.max, 
                      length = length, upper = upper, lower = lower, offset = offset, 
                      scale = scale, max.rank = max.rank, wmax = wmax))
  }
  else if (model == "covariates") {
    return(mimi.cov(y = y, x = x, var.type = var.type, lambda1 = lambda1, 
                    lambda2 = lambda2, maxit = maxit, mu0 = mu0, alpha0 = alpha0, 
                    theta0 = theta0, thresh = thresh, trace.it = trace.it, 
                    lambda1.max = lambda1.max, lambda2.max = lambda2.max, 
                    length = length, upper = upper, lower = lower, offset = offset, 
                    scale = scale, max.rank = max.rank, wmax = wmax))
  }
  else {
    return(mimi.lr(y = y, var.type = var.type, lambda1 = lambda1, 
                   maxit = maxit, theta0 = theta0, thresh = thresh, 
                   trace.it = trace.it, lambda1.max = lambda1.max, length = length, 
                   upper = upper, lower = lower, offset = offset, scale = scale, 
                   max.rank = max.rank, wmax = wmax))
  }
}

quad_approx<-function (y, param, var.type, nlevel, wmax) 
{
  n <- nrow(y)
  p <- length(var.type)
  yt <- rep(0, n)
  vt <- rep(0, n)
  count <- 1
  for (j in 1:p) {
    w <- lapply(1:n, function(i) wght(y = y[i, count:(count + 
                                                        nlevel[j] - 1)], param = param[i, count:(count + 
                                                                                                   nlevel[j] - 1)], var.type = var.type[j], nlevel = nlevel[j], 
                                      wmax))
    count <- count + nlevel[j]
    if (nlevel[j] == 1) {
      ytilde <- sapply(1:n, function(i) w[[i]]$ytilde)
      vtilde2 <- sapply(1:n, function(i) w[[i]]$vtilde2)
    }
    else {
      ytilde <- t(sapply(1:n, function(i) w[[i]]$ytilde))
      vtilde2 <- t(sapply(1:n, function(i) w[[i]]$vtilde2))
    }
    yt <- cbind(yt, as.matrix(ytilde))
    vt <- cbind(vt, as.matrix(vtilde2))
  }
  return(list(ytilde = yt[, 2:ncol(yt)], vtilde2 = vt[, 2:ncol(vt)]))
}


wght <- function (y, param, var.type, nlevel = NULL, wmax) 
{
  if (var.type == "gaussian") {
    ytilde <- y
    vtilde2 <- 1
  }
  else if (var.type == "binary") {
    vtilde2 <- 1/4
    ytilde <- y/vtilde2 - (exp(param)/(1 + exp(param)))/vtilde2 + 
      param
  }
  else if (var.type == "poisson") {
    vtilde2 <- exp(param)
    ytilde <- (y - exp(param))/vtilde2 + param
  }
  else if (var.type == "categorical") {
    vtilde2 <- rep(0.25, nlevel)
    ytilde <- y/vtilde2 - (exp(param)/sum(exp(param)))/vtilde2 + 
      param
  }
  else {
    print(var.type)
    stop("Incorrect type of variable. Should be 'gaussian', 'binary', 'poisson' or 'categorical'.")
  }
  return(list(ytilde = ytilde, vtilde2 = vtilde2))
}


bls.lr<-function (y0, theta, theta.tmp, b = 0.5, lambda1, var.type, thresh = 1e-05) 
{
  d <- dim(y0)
  n <- d[1]
  p <- d[2]
  omega <- !is.na(y0)
  param <- theta
  param.tmp <- theta.tmp
  gaus.tmp <- (1/2) * sum((y0[, var.type == "gaussian"] - param.tmp[, 
                                                                    var.type == "gaussian"])^2, na.rm = T)
  pois.tmp <- sum(-(y0[, var.type == "poisson"] * param.tmp[, 
                                                            var.type == "poisson"]) + exp(param.tmp[, var.type == 
                                                                                                      "poisson"]), na.rm = T)
  binom.tmp <- sum(-(y0[, var.type == "binary"] * param.tmp[, 
                                                            var.type == "binary"]) + log(1 + exp(param.tmp[, var.type == 
                                                                                                             "binary"])), na.rm = T)
  truc <- rep(0, n)
  if (sum(var.type == "categorical") > 0) {
    for (j in 1:sum(var.type == "categorical")) {
      tt <- rowSums(exp(param[, which(var.type == "categorical")[j]:(which(var.type == 
                                                                             "categorical")[j] + nlevel[which(var.type == 
                                                                                                                "categorical")[j]] - 1)]))
      truc <- cbind(truc, matrix(rep(tt, nlevel[which(var.type == 
                                                        "categorical")[j]]), nrow = n))
    }
    truc <- truc[, 2:ncol(truc)]
    cat.tmp <- sum(-(y0[, vt2 == "categorical"] * (param.tmp[, 
                                                             vt2 == "categorical"]) - log(truc)), na.rm = T)
  }
  else cat.tmp <- 0
  d.tmp <- svd(theta.tmp)$d
  gaus <- (1/2) * sum((y0[, var.type == "gaussian"] - param[, 
                                                            var.type == "gaussian"])^2, na.rm = T)
  pois <- sum(-(y0[, var.type == "poisson"] * param[, var.type == 
                                                      "poisson"]) + exp(param[, var.type == "poisson"]), na.rm = T)
  binom <- sum(-(y0[, var.type == "binary"] * param[, var.type == 
                                                      "binary"]) + log(1 + exp(param[, var.type == "binary"])), 
               na.rm = T)
  truc <- rep(0, n)
  if (sum(var.type == "categorical") > 0) {
    for (j in 1:sum(var.type == "categorical")) {
      tt <- rowSums(exp(param[, which(var.type == "categorical")[j]:(which(var.type == 
                                                                             "categorical")[j] + nlevel[which(var.type == 
                                                                                                                "categorical")[j]] - 1)]))
      truc <- cbind(truc, matrix(rep(tt, nlevel[which(var.type == 
                                                        "categorical")[j]]), nrow = n))
    }
    truc <- truc[, 2:ncol(truc)]
    cat <- sum(-(y0[, vt2 == "categorical"] * (param[, vt2 == 
                                                       "categorical"]) - log(truc)), na.rm = T)
  }
  else cat <- 0
  d <- svd(theta)$d
  t <- 1
  theta2 <- (1 - t) * theta.tmp + t * theta
  param2 <- (1 - t) * param.tmp + t * param
  diff <- pois.tmp - pois + gaus.tmp - gaus + binom.tmp - binom + 
    lambda1 * (sum(d.tmp) - sum(d))
  number <- gaus.tmp + pois.tmp + binom.tmp + lambda1 * sum(d.tmp)
  while (diff < -abs(number) * thresh/2) {
    t <- b * t
    theta2 <- (1 - t) * theta.tmp + t * theta
    param2 <- (1 - t) * param.tmp + t * param
    gaus <- (1/2) * sum((y0[, var.type == "gaussian"] - param2[, 
                                                               var.type == "gaussian"])^2, na.rm = T)
    pois <- sum(-(y0[, var.type == "poisson"] * param2[, 
                                                       var.type == "poisson"]) + exp(param2[, var.type == 
                                                                                              "poisson"]), na.rm = T)
    binom <- sum(-(y0[, var.type == "binary"] * param2[, 
                                                       var.type == "binary"]) + log(1 + exp(param2[, var.type == 
                                                                                                     "binary"])), na.rm = T)
    truc <- rep(0, n)
    if (sum(var.type == "categorical") > 0) {
      for (j in 1:sum(var.type == "categorical")) {
        tt <- rowSums(exp(param[, which(var.type == "categorical")[j]:(which(var.type == 
                                                                               "categorical")[j] + nlevel[which(var.type == 
                                                                                                                  "categorical")[j]] - 1)]))
        truc <- cbind(truc, matrix(rep(tt, nlevel[which(var.type == 
                                                          "categorical")[j]]), nrow = n))
      }
      truc <- truc[, 2:ncol(truc)]
      cat <- sum(-(y0[, vt2 == "categorical"] * (param[, 
                                                       vt2 == "categorical"]) - log(truc)), na.rm = T)
    }
    else cat <- 0
    d <- svd(theta2)$d
    diff <- pois.tmp - pois + gaus.tmp - gaus + binom.tmp - 
      binom + lambda1 * (sum(d.tmp) - sum(d))
  }
  obj <- pois + gaus + binom + lambda1 * d
  return(list(theta = theta2, objective = obj, t = t))
}