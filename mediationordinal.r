data <- d7
data$mRS <- ordered(as.factor(data$mRS)) %>% forcats::fct_rev()


stat.fcn <- function(data, indices) {
  covs <- ' + age + NIHSS + vol0'
  med <- 'PC1'
  d.temp <- data[indices,] # allows boot to select sample
  mod.c.star <- polr(as.formula(paste0('mRS ~ ', med, ' + treatment', covs)), data = d.temp, Hess=TRUE,  na.action = na.omit, method = 'logistic')
  mod.c <- polr(as.formula(paste0('mRS ~ treatment', covs)), data = d.temp, Hess=TRUE,  na.action = na.omit, method = 'logistic')
  
  OR.c <- summary(mod.c)$coefficients['treatmentrtPA','Value'] %>% exp()
  OR.c.star <- summary(mod.c.star)$coefficients['treatmentrtPA','Value'] %>% exp()
  prop.med <- log(OR.c/OR.c.star)/log(OR.c)

  return(prop.med)
}

res <- boot::boot(data = data, statistic = stat.fcn, R = 1e4)
plot(res)
res2 <- res
res2$t <- lapply(res2$t, function(x)(max(x,-100))) %>% unlist() %>% as.matrix()
boot::boot.ci(res, conf = .95)
res$t[res$t<0] %>% sort()
