function (model.m, model.y, sims = 1000, boot = FALSE, boot.ci.type = "perc", 
          treat = "treat.name", mediator = "med.name", covariates = NULL, 
          outcome = NULL, control = NULL, conf.level = 0.95, control.value = 0, 
          treat.value = 1, long = TRUE, dropobs = FALSE, robustSE = FALSE, 
          cluster = NULL, group.out = NULL, ...) 
{
  cl <- match.call()
  if (match("INT", names(cl), 0L)) {
    warning("'INT' is deprecated - existence of interaction terms is now automatically detected from model formulas")
  }
  if (robustSE && boot) {
    warning("'robustSE' is ignored for nonparametric bootstrap")
  }
  if (!is.null(cluster) && boot) {
    warning("'cluster' is ignored for nonparametric bootstrap")
  }
  if (robustSE & !is.null(cluster)) {
    stop("choose either `robustSE' or `cluster' option, not both")
  }
  if (boot.ci.type != "bca" & boot.ci.type != "perc") {
    stop("choose either `bca' or `perc' for boot.ci.type")
  }

  isGam.y <- inherits(model.y, "gam")
  isGam.m <- inherits(model.m, "gam")
  isGlm.y <- inherits(model.y, "glm")
  isGlm.m <- inherits(model.m, "glm")
  isLm.y <- inherits(model.y, "lm")
  isLm.m <- inherits(model.m, "lm")
  isVglm.y <- inherits(model.y, "vglm")
  isRq.y <- inherits(model.y, "rq")
  isRq.m <- inherits(model.m, "rq")
  isOrdered.y <- inherits(model.y, "polr")
  isOrdered.m <- inherits(model.m, "polr")
  isSurvreg.y <- inherits(model.y, "survreg")
  isSurvreg.m <- inherits(model.m, "survreg")
  isMer.y <- inherits(model.y, "merMod")
  isMer.m <- inherits(model.m, "merMod")


 
  m.data <- model.frame(model.m)
  y.data <- model.frame(model.y)
  if (!is.null(cluster)) {
    row.names(m.data) <- 1:nrow(m.data)
    row.names(y.data) <- 1:nrow(y.data)
    if (!is.null(model.m$weights)) {
      m.weights <- as.data.frame(model.m$weights)
      m.name <- as.character(model.m$call$weights)
      names(m.weights) <- m.name
      m.data <- cbind(m.data, m.weights)
    }
    if (!is.null(model.y$weights)) {
      y.weights <- as.data.frame(model.y$weights)
      y.name <- as.character(model.y$call$weights)
      names(y.weights) <- y.name
      y.data <- cbind(y.data, y.weights)
    }
  }


 
  
    group.m <- NULL
    group.y <- NULL
    group.out <- NULL
  
    group.id.m <- NULL
  
    group.id.y <- NULL
    
    group.id <- NULL
    group.name <- NULL
    
  n.m <- nrow(m.data)
  n.y <- nrow(y.data)
 
  m <- length(sort(unique(model.frame(model.m)[, 1])))
  
  weights.m <- model.weights(m.data)
  weights.y <- model.weights(y.data)


  weights.m <- rep(1, nrow(m.data))

  weights.y <- rep(1, nrow(y.data))
  
  weights <- weights.m
 
  if (is.character(m.data[, treat])) {
    m.data[, treat] <- factor(m.data[, treat])
  }
  
  if (is.character(y.data[, treat])) {
    y.data[, treat] <- factor(y.data[, treat])
  }
  
  if (is.character(y.data[, mediator])) {
    y.data[, mediator] <- factor(y.data[, mediator])
  }
  
  isFactorT.m <- is.factor(m.data[, treat])
  isFactorT.y <- is.factor(y.data[, treat])
  
  
  isFactorT <- isFactorT.y

  cat.0 <- control.value
  cat.1 <- treat.value
  
  isFactorM <- is.factor(y.data[, mediator])
  
  
  indexmax <- function(x) {
    order(x)[length(x)]
  }
  
  getvcov <- function(dat, fm, cluster) {
    cluster <- factor(cluster)
    M <- nlevels(cluster)
    N <- sum(!is.na(cluster))
    K <- fm$rank
    dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
    uj <- apply(estfun(fm), 2, function(x) tapply(x, cluster, 
                                                  sum))
    dfc * sandwich(fm, meat. = crossprod(uj)/N)
  }
  
  if (!isOrdered.y) {
    if (!boot) {
      
        MModel.coef <- coef(model.m)
        scalesim.m <- FALSE
        
        MModel.var.cov <- vcov(model.m)

        YModel.coef <- coef(model.y)
        scalesim.y <- FALSE
     
        
        
        YModel.var.cov <- vcov(model.y)
        
      
      se.ranef.new <- function(object) {
        se.bygroup <- lme4::ranef(object, condVar = TRUE)
        n.groupings <- length(se.bygroup)
        for (m in 1:n.groupings) {
          vars.m <- attr(se.bygroup[[m]], "postVar")
          K <- dim(vars.m)[1]
          J <- dim(vars.m)[3]
          names.full <- dimnames(se.bygroup[[m]])
          se.bygroup[[m]] <- array(NA, c(J, K))
          for (j in 1:J) {
            se.bygroup[[m]][j, ] <- sqrt(diag(as.matrix(vars.m[, 
                                                               , j])))
          }
          dimnames(se.bygroup[[m]]) <- list(names.full[[1]], 
                                            names.full[[2]])
        }
        return(se.bygroup)
      }
      
      
      MModel <- rmvnorm(sims, mean = MModel.coef, 
                          sigma = MModel.var.cov)
 
      YModel <- rmvnorm(sims, mean = YModel.coef, 
                          sigma = YModel.var.cov)
      
     
     
      pred.data.t <- pred.data.c <- m.data
      
      pred.data.t[, treat] <- factor(cat.1, levels = t.levels)
      pred.data.c[, treat] <- factor(cat.0, levels = t.levels)
      
     
      mmat.t <- model.matrix(terms(model.m), data = pred.data.t)
      mmat.c <- model.matrix(terms(model.m), data = pred.data.c)
    
    
      
      sigma <- summary(model.m)$sigma
      error <- rnorm(sims * n, mean = 0, sd = sigma)
      muM1 <- tcrossprod(MModel, mmat.t)
      muM0 <- tcrossprod(MModel, mmat.c)
      PredictM1 <- muM1 + matrix(error, nrow = sims)
      PredictM0 <- muM0 + matrix(error, nrow = sims)
      rm(error)
     
     
      rm(mmat.t, mmat.c)
    
      effects.tmp <- array(NA, dim = c(n, sims, 4))
     
      for (e in 1:4) {
        tt <- switch(e, c(1, 1, 1, 0), c(0, 0, 1, 0), 
                     c(1, 0, 1, 1), c(1, 0, 0, 0))
        Pr1 <- matrix(, nrow = n, ncol = sims)
        Pr0 <- matrix(, nrow = n, ncol = sims)
        for (j in 1:sims) {
          pred.data.t <- pred.data.c <- y.data
        
          cat.t <- ifelse(tt[1], cat.1, cat.0)
          cat.c <- ifelse(tt[2], cat.1, cat.0)
          cat.t.ctrl <- ifelse(tt[1], cat.0, cat.1)
          cat.c.ctrl <- ifelse(tt[2], cat.0, cat.1)

            pred.data.t[, treat] <- factor(cat.t, levels = t.levels)
            pred.data.c[, treat] <- factor(cat.c, levels = t.levels)

            red.data.t[, control] <- factor(cat.t.ctrl, levels = t.levels)
            pred.data.c[, control] <- factor(cat.c.ctrl, levels = t.levels)
            
          
         
          PredictMt <- PredictM1[j, ] * tt[3] + PredictM0[j, ] * (1 - tt[3])
          PredictMc <- PredictM1[j, ] * tt[4] + PredictM0[j, ] * (1 - tt[4])
        
            pred.data.t[, mediator] <- PredictMt
            pred.data.c[, mediator] <- PredictMc
          
          ymat.t <- model.matrix(terms(model.y), data = pred.data.t)
          ymat.c <- model.matrix(terms(model.y), data = pred.data.c)
          
          
         
          
            Pr1[, j] <- t(as.matrix(YModel[j, ])) %*% t(ymat.t)
            Pr0[, j] <- t(as.matrix(YModel[j, ])) %*% t(ymat.c)
          
          rm(ymat.t, ymat.c, pred.data.t, pred.data.c)
        }

          Pr1 <- apply(Pr1, 2, model.y$family$linkinv)
          Pr0 <- apply(Pr0, 2, model.y$family$linkinv)
        
        
       
        effects.tmp[, , e] <- Pr1 - Pr0
        rm(Pr1, Pr0)
      }

      rm(PredictM1, PredictM0, YModel, MModel)
     
      et1 <- effects.tmp[, , 1]
      et2 <- effects.tmp[, , 2]
      et3 <- effects.tmp[, , 3]
      et4 <- effects.tmp[, , 4]
      
      delta.1 <- t(as.matrix(apply(et1, 2, weighted.mean, w = weights)))
      delta.0 <- t(as.matrix(apply(et2, 2, weighted.mean, w = weights)))
      zeta.1 <- t(as.matrix(apply(et3, 2, weighted.mean,  w = weights)))
      zeta.0 <- t(as.matrix(apply(et4, 2, weighted.mean,  w = weights)))
      rm(effects.tmp)
      
      tau <- (zeta.1 + delta.0 + zeta.0 + delta.1)/2
      nu.0 <- delta.0/tau
      nu.1 <- delta.1/tau
      delta.avg <- (delta.1 + delta.0)/2
      zeta.avg <- (zeta.1 + zeta.0)/2
      nu.avg <- (nu.1 + nu.0)/2
      d0 <- mean(delta.0)
      d1 <- mean(delta.1)
      z1 <- mean(zeta.1)
      z0 <- mean(zeta.0)
      tau.coef <- mean(tau)
      n0 <- median(nu.0)
      n1 <- median(nu.1)
      d.avg <- (d0 + d1)/2
      z.avg <- (z0 + z1)/2
      n.avg <- (n0 + n1)/2
     
    }
    else {
     
      Call.M <- getCall(model.m)
      Call.Y <- getCall(model.y)
      delta.1 <- matrix(NA, sims, 1)
      delta.0 <- matrix(NA, sims, 1)
      zeta.1 <- matrix(NA, sims, 1)
      zeta.0 <- matrix(NA, sims, 1)
      tau <- matrix(NA, sims, 1)
      for (b in 1:(sims + 1)) {
        index <- sample(1:n, n, replace = TRUE)
        if (b == sims + 1) {
          index <- 1:n
        }
       
        
          Call.M$data <- m.data[index, ]
        
       
        
          Call.Y$data <- y.data[index, ]
        Call.M$weights <- m.data[index, "(weights)"]
        Call.Y$weights <- y.data[index, "(weights)"]
        
        new.fit.M <- eval.parent(Call.M)
        new.fit.Y <- eval.parent(Call.Y)
        pred.data.t <- pred.data.c <- m.data
        if (isFactorT) {
          pred.data.t[, treat] <- factor(cat.1, levels = t.levels)
          pred.data.c[, treat] <- factor(cat.0, levels = t.levels)
        }
        else {
          pred.data.t[, treat] <- cat.1
          pred.data.c[, treat] <- cat.0
        }
        if (!is.null(covariates)) {
          for (p in 1:length(covariates)) {
            vl <- names(covariates[p])
            if (is.factor(pred.data.t[, vl])) {
              pred.data.t[, vl] <- pred.data.c[, vl] <- factor(covariates[[p]], 
                                                               levels = levels(m.data[, vl]))
            }
            else {
              pred.data.t[, vl] <- pred.data.c[, vl] <- covariates[[p]]
            }
          }
        }
       
        
        if (isLm.m) {
          sigma <- summary(new.fit.M)$sigma
          error <- rnorm(n, mean = 0, sd = sigma)
          PredictM1 <- predict(new.fit.M, type = "response", 
                               newdata = pred.data.t) + error
          PredictM0 <- predict(new.fit.M, type = "response", 
                               newdata = pred.data.c) + error
          rm(error)
        }
        else {
          stop("mediator model is not yet implemented")
        }
        effects.tmp <- matrix(NA, nrow = n, ncol = 4)
        for (e in 1:4) {
          tt <- switch(e, c(1, 1, 1, 0), c(0, 0, 1, 
                                           0), c(1, 0, 1, 1), c(1, 0, 0, 0))
          pred.data.t <- pred.data.c <- y.data
          if (!is.null(covariates)) {
            for (p in 1:length(covariates)) {
              vl <- names(covariates[p])
              if (is.factor(pred.data.t[, vl])) {
                pred.data.t[, vl] <- pred.data.c[, vl] <- factor(covariates[[p]], 
                                                                 levels = levels(y.data[, vl]))
              }
              else {
                pred.data.t[, vl] <- pred.data.c[, vl] <- covariates[[p]]
              }
            }
          }
          cat.t <- ifelse(tt[1], cat.1, cat.0)
          cat.c <- ifelse(tt[2], cat.1, cat.0)
          cat.t.ctrl <- ifelse(tt[1], cat.0, cat.1)
          cat.c.ctrl <- ifelse(tt[2], cat.0, cat.1)
          if (isFactorT) {
            pred.data.t[, treat] <- factor(cat.t, levels = t.levels)
            pred.data.c[, treat] <- factor(cat.c, levels = t.levels)
            if (!is.null(control)) {
              pred.data.t[, control] <- factor(cat.t.ctrl, 
                                               levels = t.levels)
              pred.data.c[, control] <- factor(cat.c.ctrl, 
                                               levels = t.levels)
            }
          }
          else {
            pred.data.t[, treat] <- cat.t
            pred.data.c[, treat] <- cat.c
            if (!is.null(control)) {
              pred.data.t[, control] <- cat.t.ctrl
              pred.data.c[, control] <- cat.c.ctrl
            }
          }
          PredictM1.tmp <- PredictM1
          PredictM0.tmp <- PredictM0
          PredictMt <- PredictM1 * tt[3] + PredictM0 * 
            (1 - tt[3])
          PredictMc <- PredictM1 * tt[4] + PredictM0 * 
            (1 - tt[4])
          if (isFactorM) {
            pred.data.t[, mediator] <- factor(PredictMt, 
                                              levels = 1:m, labels = m.levels)
            pred.data.c[, mediator] <- factor(PredictMc, 
                                              levels = 1:m, labels = m.levels)
          }
          else {
            pred.data.t[, mediator] <- PredictMt
            pred.data.c[, mediator] <- PredictMc
          }
         
          
            pr.1 <- predict(new.fit.Y, type = "response", 
                            newdata = pred.data.t)
            pr.0 <- predict(new.fit.Y, type = "response", 
                            newdata = pred.data.c)
          
          pr.mat <- as.matrix(cbind(pr.1, pr.0))
          effects.tmp[, e] <- pr.mat[, 1] - pr.mat[, 
                                                   2]
          rm(pred.data.t, pred.data.c, pr.1, pr.0, pr.mat)
        }
        if (b == sims + 1) {
          d1 <- weighted.mean(effects.tmp[, 1], weights)
          d0 <- weighted.mean(effects.tmp[, 2], weights)
          z1 <- weighted.mean(effects.tmp[, 3], weights)
          z0 <- weighted.mean(effects.tmp[, 4], weights)
        }
        else {
          delta.1[b] <- weighted.mean(effects.tmp[, 
                                                  1], weights)
          delta.0[b] <- weighted.mean(effects.tmp[, 
                                                  2], weights)
          zeta.1[b] <- weighted.mean(effects.tmp[, 3], 
                                     weights)
          zeta.0[b] <- weighted.mean(effects.tmp[, 4], 
                                     weights)
        }
      }
      tau.coef <- (d1 + d0 + z1 + z0)/2
      n0 <- d0/tau.coef
      n1 <- d1/tau.coef
      d.avg <- (d1 + d0)/2
      z.avg <- (z1 + z0)/2
      n.avg <- (n0 + n1)/2
      tau <- (delta.1 + delta.0 + zeta.1 + zeta.0)/2
      nu.0 <- delta.0/tau
      nu.1 <- delta.1/tau
      delta.avg <- (delta.0 + delta.1)/2
      zeta.avg <- (zeta.0 + zeta.1)/2
      nu.avg <- (nu.0 + nu.1)/2
    }
    low <- (1 - conf.level)/2
    high <- 1 - low
    if (boot & boot.ci.type == "bca") {
      BC.CI <- function(theta) {
        z.inv <- length(theta[theta < mean(theta)])/sims
        z <- qnorm(z.inv)
        U <- (sims - 1) * (mean(theta) - theta)
        top <- sum(U^3)
        under <- 6 * (sum(U^2))^{
          3/2
        }
        a <- top/under
        lower.inv <- pnorm(z + (z + qnorm(low))/(1 - 
                                                   a * (z + qnorm(low))))
        lower2 <- lower <- quantile(theta, lower.inv)
        upper.inv <- pnorm(z + (z + qnorm(high))/(1 - 
                                                    a * (z + qnorm(high))))
        upper2 <- upper <- quantile(theta, upper.inv)
        return(c(lower, upper))
      }
      d0.ci <- BC.CI(delta.0)
      d1.ci <- BC.CI(delta.1)
      tau.ci <- BC.CI(tau)
      z1.ci <- BC.CI(zeta.1)
      z0.ci <- BC.CI(zeta.0)
      n0.ci <- BC.CI(nu.0)
      n1.ci <- BC.CI(nu.1)
      d.avg.ci <- BC.CI(delta.avg)
      z.avg.ci <- BC.CI(zeta.avg)
      n.avg.ci <- BC.CI(nu.avg)
    }
    else {
      d0.ci <- quantile(delta.0, c(low, high), na.rm = TRUE)
      d1.ci <- quantile(delta.1, c(low, high), na.rm = TRUE)
      tau.ci <- quantile(tau, c(low, high), na.rm = TRUE)
      z1.ci <- quantile(zeta.1, c(low, high), na.rm = TRUE)
      z0.ci <- quantile(zeta.0, c(low, high), na.rm = TRUE)
      n0.ci <- quantile(nu.0, c(low, high), na.rm = TRUE)
      n1.ci <- quantile(nu.1, c(low, high), na.rm = TRUE)
      d.avg.ci <- quantile(delta.avg, c(low, high), na.rm = TRUE)
      z.avg.ci <- quantile(zeta.avg, c(low, high), na.rm = TRUE)
      n.avg.ci <- quantile(nu.avg, c(low, high), na.rm = TRUE)
    }
    d0.p <- pval(delta.0, d0)
    d1.p <- pval(delta.1, d1)
    d.avg.p <- pval(delta.avg, d.avg)
    z0.p <- pval(zeta.0, z0)
    z1.p <- pval(zeta.1, z1)
    z.avg.p <- pval(zeta.avg, z.avg)
    n0.p <- pval(nu.0, n0)
    n1.p <- pval(nu.1, n1)
    n.avg.p <- pval(nu.avg, n.avg)
    tau.p <- pval(tau, tau.coef)
 
    INT <- paste(treat, mediator, sep = ":") %in% attr(terms(model.y), 
                                                       "term.labels") | paste(mediator, treat, sep = ":") %in% 
      attr(terms(model.y), "term.labels")
   
    if (long && !isMer.y && !isMer.m) {
      out <- list(d0 = d0, d1 = d1, d0.ci = d0.ci, d1.ci = d1.ci, 
                  d0.p = d0.p, d1.p = d1.p, d0.sims = delta.0, 
                  d1.sims = delta.1, z0 = z0, z1 = z1, z0.ci = z0.ci, 
                  z1.ci = z1.ci, z0.p = z0.p, z1.p = z1.p, z0.sims = zeta.0, 
                  z1.sims = zeta.1, n0 = n0, n1 = n1, n0.ci = n0.ci, 
                  n1.ci = n1.ci, n0.p = n0.p, n1.p = n1.p, n0.sims = nu.0, 
                  n1.sims = nu.1, tau.coef = tau.coef, tau.ci = tau.ci, 
                  tau.p = tau.p, tau.sims = tau, d.avg = d.avg, 
                  d.avg.p = d.avg.p, d.avg.ci = d.avg.ci, d.avg.sims = delta.avg, 
                  z.avg = z.avg, z.avg.p = z.avg.p, z.avg.ci = z.avg.ci, 
                  z.avg.sims = zeta.avg, n.avg = n.avg, n.avg.p = n.avg.p, 
                  n.avg.ci = n.avg.ci, n.avg.sims = nu.avg, boot = boot, 
                  boot.ci.type = boot.ci.type, treat = treat, 
                  mediator = mediator, covariates = covariates, 
                  INT = INT, conf.level = conf.level, model.y = model.y, 
                  model.m = model.m, control.value = control.value, 
                  treat.value = treat.value, nobs = n, sims = sims, 
                  call = cl, robustSE = robustSE, cluster = cluster)
      class(out) <- "mediate"
    }
    if (!long && !isMer.y && !isMer.m) {
      out <- list(d0 = d0, d1 = d1, d0.ci = d0.ci, d1.ci = d1.ci, 
                  d0.p = d0.p, d1.p = d1.p, z0 = z0, z1 = z1, 
                  z0.ci = z0.ci, z1.ci = z1.ci, z0.p = z0.p, z1.p = z1.p, 
                  n0 = n0, n1 = n1, n0.ci = n0.ci, n1.ci = n1.ci, 
                  n0.p = n0.p, n1.p = n1.p, tau.coef = tau.coef, 
                  tau.ci = tau.ci, tau.p = tau.p, d.avg = d.avg, 
                  d.avg.p = d.avg.p, d.avg.ci = d.avg.ci, z.avg = z.avg, 
                  z.avg.p = z.avg.p, z.avg.ci = z.avg.ci, n.avg = n.avg, 
                  n.avg.p = n.avg.p, n.avg.ci = n.avg.ci, boot = boot, 
                  boot.ci.type = boot.ci.type, treat = treat, 
                  mediator = mediator, covariates = covariates, 
                  INT = INT, conf.level = conf.level, model.y = model.y, 
                  model.m = model.m, control.value = control.value, 
                  treat.value = treat.value, nobs = n, sims = sims, 
                  call = cl, robustSE = robustSE, cluster = cluster)
      class(out) <- "mediate"
    }
    
    
    out
  }
  else {
    if (boot != TRUE) {
      warning("ordered outcome model can only be used with nonparametric bootstrap - option forced")
      boot <- TRUE
    }
    
    n.ycat <- length(unique(model.response(y.data)))
    delta.1 <- matrix(NA, sims, n.ycat)
    delta.0 <- matrix(NA, sims, n.ycat)
    zeta.1 <- matrix(NA, sims, n.ycat)
    zeta.0 <- matrix(NA, sims, n.ycat)
    tau <- matrix(NA, sims, n.ycat)
    for (b in 1:(sims + 1)) {
      index <- sample(1:n, n, replace = TRUE)
      if (b == sims + 1) {
        index <- 1:n
      }
      Call.M <- model.m$call
      Call.Y <- model.y$call
  
        Call.M$data <- m.data[index, ]
  
      Call.Y$data <- y.data[index, ]
      Call.M$weights <- m.data[index, "(weights)"]
      Call.Y$weights <- y.data[index, "(weights)"]
      new.fit.M <- eval.parent(Call.M)
      new.fit.Y <- eval.parent(Call.Y)
      
      pred.data.t <- pred.data.c <- m.data
      if (isFactorT) {
        pred.data.t[, treat] <- factor(cat.1, levels = t.levels)
        pred.data.c[, treat] <- factor(cat.0, levels = t.levels)
      }
      else {
        pred.data.t[, treat] <- cat.1
        pred.data.c[, treat] <- cat.0
      }
      if (!is.null(covariates)) {
        for (p in 1:length(covariates)) {
          vl <- names(covariates[p])
          if (is.factor(pred.data.t[, vl])) {
            pred.data.t[, vl] <- pred.data.c[, vl] <- factor(covariates[[p]], 
                                                             levels = levels(m.data[, vl]))
          }
          else {
            pred.data.t[, vl] <- pred.data.c[, vl] <- covariates[[p]]
          }
        }
      }
    
      if (isLm.m) {
        sigma <- summary(new.fit.M)$sigma
        error <- rnorm(n, mean = 0, sd = sigma)
        PredictM1 <- predict(new.fit.M, type = "response", 
                             newdata = pred.data.t) + error
        PredictM0 <- predict(new.fit.M, type = "response", 
                             newdata = pred.data.c) + error
        rm(error)
      }
      else {
        stop("mediator model is not yet implemented")
      }
      effects.tmp <- array(NA, dim = c(n, n.ycat, 4))
      for (e in 1:4) {
        tt <- switch(e, c(1, 1, 1, 0), c(0, 0, 1, 0), 
                     c(1, 0, 1, 1), c(1, 0, 0, 0))
        pred.data.t <- pred.data.c <- y.data
        if (!is.null(covariates)) {
          for (p in 1:length(covariates)) {
            vl <- names(covariates[p])
            if (is.factor(pred.data.t[, vl])) {
              pred.data.t[, vl] <- pred.data.c[, vl] <- factor(covariates[[p]], 
                                                               levels = levels(y.data[, vl]))
            }
            else {
              pred.data.t[, vl] <- pred.data.c[, vl] <- covariates[[p]]
            }
          }
        }
        cat.t <- ifelse(tt[1], cat.1, cat.0)
        cat.c <- ifelse(tt[2], cat.1, cat.0)
        cat.t.ctrl <- ifelse(tt[1], cat.0, cat.1)
        cat.c.ctrl <- ifelse(tt[2], cat.0, cat.1)
        if (isFactorT) {
          pred.data.t[, treat] <- factor(cat.t, levels = t.levels)
          pred.data.c[, treat] <- factor(cat.c, levels = t.levels)
          if (!is.null(control)) {
            pred.data.t[, control] <- factor(cat.t.ctrl, 
                                             levels = t.levels)
            pred.data.c[, control] <- factor(cat.c.ctrl, 
                                             levels = t.levels)
          }
        }
        else {
          pred.data.t[, treat] <- cat.t
          pred.data.c[, treat] <- cat.c
          if (!is.null(control)) {
            pred.data.t[, control] <- cat.t.ctrl
            pred.data.c[, control] <- cat.c.ctrl
          }
        }
        PredictM1.tmp <- PredictM1
        PredictM0.tmp <- PredictM0
        PredictMt <- PredictM1 * tt[3] + PredictM0 * 
          (1 - tt[3])
        PredictMc <- PredictM1 * tt[4] + PredictM0 * 
          (1 - tt[4])
        if (isFactorM) {
          pred.data.t[, mediator] <- factor(PredictMt, 
                                            levels = 1:m, labels = m.levels)
          pred.data.c[, mediator] <- factor(PredictMc, 
                                            levels = 1:m, labels = m.levels)
        }
        else {
          pred.data.t[, mediator] <- PredictMt
          pred.data.c[, mediator] <- PredictMc
        }
        probs_p1 <- predict(new.fit.Y, newdata = pred.data.t, 
                            type = "probs")
        probs_p0 <- predict(new.fit.Y, newdata = pred.data.c, 
                            type = "probs")
        effects.tmp[, , e] <- probs_p1 - probs_p0
        rm(pred.data.t, pred.data.c, probs_p1, probs_p0)
      }
      if (b == sims + 1) {
        d1 <- apply(effects.tmp[, , 1], 2, weighted.mean, 
                    w = weights)
        d0 <- apply(effects.tmp[, , 2], 2, weighted.mean, 
                    w = weights)
        z1 <- apply(effects.tmp[, , 3], 2, weighted.mean, 
                    w = weights)
        z0 <- apply(effects.tmp[, , 4], 2, weighted.mean, 
                    w = weights)
      }
      else {
        delta.1[b, ] <- apply(effects.tmp[, , 1], 2, 
                              weighted.mean, w = weights)
        delta.0[b, ] <- apply(effects.tmp[, , 2], 2, 
                              weighted.mean, w = weights)
        zeta.1[b, ] <- apply(effects.tmp[, , 3], 2, 
                             weighted.mean, w = weights)
        zeta.0[b, ] <- apply(effects.tmp[, , 4], 2, 
                             weighted.mean, w = weights)
      }
    }
    tau.coef <- (d1 + d0 + z1 + z0)/2
    tau <- (zeta.1 + zeta.0 + delta.0 + delta.1)/2
    low <- (1 - conf.level)/2
    high <- 1 - low
    if (boot.ci.type == "bca") {
      BC.CI <- function(theta) {
        z.inv <- length(theta[theta < mean(theta)])/sims
        z <- qnorm(z.inv)
        U <- (sims - 1) * (mean(theta) - theta)
        top <- sum(U^3)
        under <- 6 * (sum(U^2))^{
          3/2
        }
        a <- top/under
        lower.inv <- pnorm(z + (z + qnorm(low))/(1 - 
                                                   a * (z + qnorm(low))))
        lower2 <- lower <- quantile(theta, lower.inv)
        upper.inv <- pnorm(z + (z + qnorm(high))/(1 - 
                                                    a * (z + qnorm(high))))
        upper2 <- upper <- quantile(theta, upper.inv)
        return(c(lower, upper))
      }
      d0.ci <- BC.CI(delta.0)
      d1.ci <- BC.CI(delta.1)
      tau.ci <- BC.CI(tau)
      z1.ci <- BC.CI(zeta.1)
      z0.ci <- BC.CI(zeta.0)
    }
    else {
      CI <- function(theta) {
        return(quantile(theta, c(low, high), na.rm = TRUE))
      }
      d0.ci <- apply(delta.0, 2, CI)
      d1.ci <- apply(delta.1, 2, CI)
      tau.ci <- apply(tau, 2, CI)
      z1.ci <- apply(zeta.1, 2, CI)
      z0.ci <- apply(zeta.0, 2, CI)
    }
    d0.p <- d1.p <- z0.p <- z1.p <- tau.p <- rep(NA, n.ycat)
    for (i in 1:n.ycat) {
      d0.p[i] <- pval(delta.0[, i], d0[i])
      d1.p[i] <- pval(delta.1[, i], d1[i])
      z0.p[i] <- pval(zeta.0[, i], z0[i])
      z1.p[i] <- pval(zeta.1[, i], z1[i])
      tau.p[i] <- pval(tau[, i], tau.coef[i])
    }
    INT <- paste(treat, mediator, sep = ":") %in% attr(model.y$terms, 
                                                       "term.labels") | paste(mediator, treat, sep = ":") %in% 
      attr(model.y$terms, "term.labels")
    if (long) {
      out <- list(d0 = d0, d1 = d1, d0.ci = d0.ci, d1.ci = d1.ci, 
                  d0.p = d0.p, d1.p = d1.p, d0.sims = delta.0, 
                  d1.sims = delta.1, tau.coef = tau.coef, tau.ci = tau.ci, 
                  tau.p = tau.p, z0 = z0, z1 = z1, z0.ci = z0.ci, 
                  z1.ci = z1.ci, z0.p = z0.p, z1.p = z1.p, z1.sims = zeta.1, 
                  z0.sims = zeta.0, tau.sims = tau, boot = boot, 
                  boot.ci.type = boot.ci.type, treat = treat, 
                  mediator = mediator, covariates = covariates, 
                  INT = INT, conf.level = conf.level, model.y = model.y, 
                  model.m = model.m, control.value = control.value, 
                  treat.value = treat.value, nobs = n, sims = sims, 
                  call = cl, robustSE = robustSE, cluster = cluster)
    }
    else {
      out <- list(d0 = d0, d1 = d1, d0.ci = d0.ci, d1.ci = d1.ci, 
                  d0.p = d0.p, d1.p = d1.p, tau.coef = tau.coef, 
                  tau.ci = tau.ci, tau.p = tau.p, z0 = z0, z1 = z1, 
                  z0.ci = z0.ci, z1.ci = z1.ci, z0.p = z0.p, z1.p = z1.p, 
                  boot = boot, boot.ci.type = boot.ci.type, treat = treat, 
                  mediator = mediator, covariates = covariates, 
                  INT = INT, conf.level = conf.level, model.y = model.y, 
                  model.m = model.m, control.value = control.value, 
                  treat.value = treat.value, nobs = n, sims = sims, 
                  call = cl, robustSE = robustSE, cluster = cluster)
    }
    class(out) <- "mediate.order"
    out
  }
}
