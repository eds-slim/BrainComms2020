nlme.exp<-function(d, response, preds = NULL){
  
  cols <- preds %>% unlist() %>% as.character()
  cols <- cols[cols != 1] %>% unique()
  d <- d %>% dplyr::select(c('numID','tp',response,cols)) %>% droplevels() %>% distinct()

  deriv_mod <- deriv(~a+del*(1-exp(-b*tp)), c('tp','a','b','del'),function(tp,a,b,del){})
  
  x <- d %>% filter(tp==0) %>% pull(response) %>% mean(na.rm=TRUE)
  y <- d %>% filter(tp==1) %>% pull(response) %>% mean(na.rm=TRUE)
  z <- d %>% filter(tp==12) %>% pull(response) %>% mean(na.rm=TRUE)
  params.a <- x
  params.del <- z - x
  params.b <- log(1-(y-x)/(z-x))
  
    
  tryCatch({
  m<-nls(formula = as.formula(paste0(response,'~a+del*(1-exp(-b*tp))'))
         #, iter = c(5,5,5) 
         #, start_lower  = list(a=params.a/2, b=params.b/2, del=params.del/2)
         #, start_upper  = list(a=params.a*2, b=params.b*2, del=params.del*2)
         , start = list(a=params.a, b=params.b, del=params.del)
         , data = d
         , control = nls.control(minFactor = 1e-5, maxiter = 1e3, tol = 1e-1)
         #, control = nls.lm.control(maxiter = 1e3)
         , na.action = na.exclude
  )
  
  coefs<-coefficients(m)
  }, finally = {coefs=c(a=1,b=1,del=1)}, error = function(e)e
)
  
  if(is.null(preds)){
    fixed <- list(a+b+del~1)
    start <- coefs
  } else {
    fixed <- foreach (p = preds, n = names(preds)) %do% {as.formula(paste0(n,'~',p))} %>% unlist()
    start <- foreach (p = preds, n = names(preds)) %do% {
        if(p == 1){
          c(coefs[n])
        }
        else if(class(d[[p]]) == 'factor'){
          c(coefs[n], rep(0,nlevels(d[[p]])-1))
        } else if(class(d[[p]]) == 'numeric'){
          c(coefs[n], 0)
        }
      } %>% unlist() %>% as.numeric()
  }
  nlme(model=as.formula(paste0(response,'~deriv_mod(tp,a,b,del)'))
       , fixed=fixed
       , random=a~1|numID
       , start=start
       , data=d
       #, correlation = corCAR1(form = ~tp|numID/relpos)
       , na.action = na.exclude
       , control = nlmeControl(minScale = 1e-10
                               , maxIter=1e3
                               , msMaxIter = 1e2
                               , pnlsMaxIter = 1e2
                               , opt = 'nlminb'
                               , eval.max = 1e2
                               , pnlsTol = 5e-1
                               , msVerbose = F
                               , tolerance = 1e-2
       )
       , method = "ML"
       , verbose = F
  )
}


