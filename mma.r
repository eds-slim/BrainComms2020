d.mma.pre <- ddd %>% dplyr::select(tail(names(.),43), treatment, mRS, NIHSS, age) %$% .[complete.cases(.),]


require(mma)

x <- d.mma.pre %>% dplyr::select(head(names(.),43), NIHSS, age) 
x <- x[,c(names(out$vars),'NIHSS','age')] %>% as.data.frame()
y <- as.numeric(d.mma.pre$mRS <= 2)
pred <- d.mma.pre$treatment

dataorg <- data.org(x = x
         , y = y
         , pred = pred
         , mediator = names(out$vars)
         , contmed = names(out$vars)
         , jointm = list(n=1, j1=names(out$vars))
         , cova = d.mma.pre %>% dplyr::select(NIHSS)
         , alpha = 1, alpha2 = 1)

med(dataorg)

res <- boot.med(dataorg, n2 = 50)


res$estimation
res$model
res$bootsresults$te %>% .[abs(.)<1e3] %>% mean()
res$bootsresults$te %>% .[abs(.)<1e3]%>% sd()


res$bootsresults$de %>% .[abs(.)<1e10] %>% mean()
res$bootsresults$de %>% .[abs(.)<1e10]%>% sd()

bind_cols(te=res$bootsresults$te,de=res$bootsresults$de) %>% 
  mutate(excess = te-de) %>% 
  summarise(mean(excess>0))

summary(res, RE=F, alpha = 0.05)
