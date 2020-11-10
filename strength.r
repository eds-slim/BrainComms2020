
## lesion volume
stats.LGP.vol <- d %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='NIHSS') %>% 
  group_by(node,LGP.name) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme)(LGP.value~deriv_mod(tp,a,b,del)
                                    , fixed=list(a+b~1,del~logvol*relpos)
                                    , random=a~1|numID
                                    , data=.
                                    , start=c(1,1,1,1,1,1)
                                    , na.action = na.omit
  )$result
  )
  ) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))


stats.LGP.vol %>% 
  unnest(tidy,.drop = T) %>%
  filter(node %in% nodes.strength) %>% 
  select(-c(effect,group)) %T>%
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4)) %>% 
  write.table(file = './../Word/Floats/Tables/10-LGP_logvol.tab', sep = ",", quote = FALSE, row.names = F)



stats.LGP.vol.relpos <- d %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='NIHSS') %>% 
  group_by(node,relpos, LGP.name) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme.exp)(., response = 'LGP.value', preds=list(a='1', b='1',del='logvol'))$result)) %>% 
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))


stats.LGP.vol.relpos %>% 
  arrange(LGP.name,relpos) %>% 
  unnest(tidy,.drop = T) %>% 
  filter(node %in% nodes.strict) %>% 
  select(-c(effect,group,df,statistic)) %T>%
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4)) %>% 
  write.table(file = './../Word/Floats/Tables/11-LGP_logvol_by_relpos.tab', sep = ",", quote = FALSE, row.names = F)  


stats.LGP.vol.relpos <- d.LGP %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='NIHSS') %>% 
  group_by(node,relpos, LGP.name) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme.exp)(., response = 'LGP.value', preds=list(a='1', b='1',del='recovery'))$result)) %>% 
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))


# clinical

## mass univariate
stats.LGP.clinical <- d.LGP %>% 
  merge(data.clinical_long) %>% 
  group_by(node,LGP.name,clin.meas.name) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme)(LGP.value~deriv_mod(tp,a,b,del)
                                    , fixed=list(a+b~1,del~clin.meas.value*relpos)
                                    , random=a~1|numID
                                    , data=.
                                    , start=c(1,1,1,1,1,1)
                                    , na.action = na.omit
  )$result
  )
  ) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))

stats.LGP.clinical %>% 
  unnest(tidy,.drop = T) %>%
  select(-c(effect,group)) %T>%
  print(n=Inf) %>% 
  mutate_if(is.numeric, ~signif(.,digits = 4)) %>% 
  write.table(file = './../Word/Floats/Tables/13-LGP_clinical.tab', sep = ",", quote = FALSE, row.names = F)  

stats.LGP.clinical.relpos <- d.LGP %>% 
  merge(data.clinical_long) %>% 
  filter(relpos=='ipsi') %>% 
  group_by(node,relpos, LGP.name, clin.meas.name) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme.exp)(., response = 'LGP.value', preds=list(a='1', b='1',del='clin.meas.value'))$result)) %>% 
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))

stats.LGP.clinical.relpos %>% 
  unnest(tidy,.drop = T) %>%
  select(-c(effect,group)) %>% 
  filter(term=='del.clin.meas.value' & p.value<0.05/3/3) %>% 
  arrange(LGP.name) %>% 
  print(n=Inf)

# mixed effects of node
deriv_mod <- deriv(~a+del*(1-exp(-b*tp)), c('tp','a','b','del'),function(tp,a,b,del){})
test<-d.LGP
test$Dummy <- factor(1)
test$numID<-as.factor(test$numID)
test<-subset(test,is.finite(LGP.value))
stats.LGP.clinical.mixed<-test %>% as_tibble() %>%
  merge(data.clinical_long) %>% 
  filter(relpos=='ipsi') %>%
  group_by(clin.meas.name, LGP.name) %>%
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme)(LGP.value~deriv_mod(tp,a,b,del)
                              , fixed=list(a~1,b~1,del~clin.meas.value)
                              , random = list(Dummy=pdIdent(a ~ node)
                                              ,numID=pdIdent(a ~ 1))
                              #, random = c~1|numID
                              #, random = c~1|node/numID
                              , data=subset(., node %in% (nodes.liberal %>% filter(LGP.name==unique(.$LGP.name)) %>% pull('node') %>% as.character()))
                              , start=c(1,1,1,1)
                              #, correlation = corCAR1(form = ~tp|Dummy/numID/relpos/node)
                              , na.action = na.omit
                              , control = nlmeControl(tolerance = 1e-2
                                                      ,minScale = 1e-7
                                                      , maxIter = 100
                                                      , pnlsTol = 1e-4
                                                      ,opt = 'nlminb'),
                              method = 'ML')$result
  )
  ) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))


  
  
stats.LGP.clinical.mixed %>% 
  unnest(tidy, .drop = T) %>% 
  select(-c(effect,group,df)) %>% 
  filter(term=='del.clin.meas.value' & p.value<=0.05)
    
  

m %>% summary()



## plot nonlinear fit ###


d.plt<-d.LGP %>% 
  as_tibble() %>%
  merge(data.clinical_long) %>% 
  filter(relpos=='ipsi' & clin.meas.name=='NIHSS' & LGP.name == 'efficiency') %>% 
  as_tibble() 


rhs_fcn<-function(data,coefs){
  d<-c(NA,coef(m)[['a']]) %>% enframe(name=NULL,value='a') %>% 
    cbind(unique(d.plt$numID) %>% enframe(name=NULL,value='numID')) %>% 
    as_tibble() %>% 
    join(data) %>% 
    as_tibble()%>% 
    mutate(rhs=coefs[['a']]+(coefs[['del.(Intercept)']]+coefs[['del.clin.meas.value']]*clin.meas.value)*(1-exp(-coefs[['b']]*tp)))
  return(d)
}

plt_fcn<-function(data,node){
  my.rmc<-rmcorr(participant = numID, measure2 = LGP.value, measure1 = rhs, dataset = data)
  if(my.rmc$p>0.05) return(nullGrob())
  data %>% 
    ggplot(aes(x=LGP.value,y=rhs, color=as.factor(tp)))+
    geom_point()+
    geom_smooth(aes(group=numID), method = 'lm', se = F) +
    geom_smooth(aes(group=1), method = 'lm', se = T, color='green')+
    ggtitle(paste0(node,'\np.value in rmcorr = ', sprintf('%f',my.rmc$p), '\n r = ' ,sprintf('%f',my.rmc$r))) +
    guides(color=F)+
    theme(title = element_text(size = 6))
}

d.plt %>% 
  group_by(node) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme.exp)(., response = 'LGP.value', preds=list(a='1', b='1',del='clin.meas.value'))$result)) %>% 
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe())) %>% 
  mutate(coefs=map(mdl,fixef)) %>% 
  mutate(data=pmap(list(data,coefs),~rhs_fcn(.x,.y))) %>% 
  mutate(plt=pmap(list(data,node),~plt_fcn(.x,.y))) %$% 
  do.call("grid.arrange", list(grobs = .$plt, nrow = 9))

m<-nlme(model=LGP.value~deriv_mod(tp,a,b,del)
                                    , fixed=list(a~1,b~1,del~clin.meas.value)
                                    , random = list(#Dummy=pdIdent(a ~ node)
                                      numID=pdIdent(a ~ 1))
                                    #, random = c~1|numID
                                    #, random = c~1|node/numID
                                    , data=subset(d.plt, node %in% c('thalamus'))
                                    , start=c(1,1,1,1)
                                    #, correlation = corCAR1(form = ~tp|Dummy/numID/relpos/node)
                                    , na.action = na.omit
                                    , control = nlmeControl(tolerance = 1e-2
                                                            ,minScale = 1e-7
                                                            , maxIter = 100
                                                            , pnlsTol = 1e-4
                                                            ,opt = 'nlminb'),
                                    method = 'ML')
m %>% tidy()
coefs<-fixef(m)
dd.plt<-c(NA,coef(m)[['a']]) %>% enframe(name=NULL,value='a') %>% 
  cbind(unique(d.plt$numID) %>% enframe(name=NULL,value='numID')) %>% 
  as_tibble() %>% 
  join(d.plt) %>% 
  as_tibble()%>% 
  mutate(rhs=coefs[['a']]+(coefs[['del.(Intercept)']]+coefs[['del.clin.meas.value']]*clin.meas.value)*(1-exp(-coefs[['b']]*tp)))
  
my.rmc<-dd.plt %>% 
  rmcorr(participant = numID, measure2 = LGP.value, measure1 = rhs, dataset = .)
plot(my.rmc)

dd.plt %>% 
  ggplot(aes(x=LGP.value,y=rhs, color=as.factor(tp)))+
  geom_point()+
  geom_smooth(aes(group=numID), method = 'lm', se = F) +
  geom_smooth(aes(group=1), method = 'lm', se = T, color='green')

dd.plt %>% 
  lme(rhs~LGP.value, random = ~1|numID, data=., na.action = na.omit) %>% 
  tidy()

dd.plt %>% 
  lm(rhs~LGP.value+numID, data=., na.action = na.omit) %>% 
  drop1(., ~., test = "F")


## simplify relationship clinical / degeneration

str(d)
dd<-d.LGP %>%
  merge(data.clinical_long) %>% 
  filter(tp %in% c(0,1,3,12) & relpos == 'ipsi') %>%
  filter(LGP.name == 'efficiency' & clin.meas.name == 'rGS' & node %in% nodes.strict$node) %>%
  group_by(numID,node) %>%
  arrange(tp, .by_group = TRUE) %>%
  mutate(diff.LGP.value = LGP.value - lag(LGP.value, n=3, default = NA))%>%
  mutate(diff.clin.meas.value = clin.meas.value - lag(clin.meas.value, default = first(clin.meas.value))) %>% 
  filter(tp==12  & diff.LGP.value>-5)

dd %>% 
  ggplot(aes(y=clin.meas.value,x=diff.LGP.value, color=node))+
  geom_point()+
  geom_smooth(aes(group=node), method = 'lm', se = F)

dd %>% 
  group_by(node) %>% 
  nest() %>% 
  mutate(mdl=map(data,~lm(clin.meas.value~diff.LGP.value, data=., na.action = na.omit))) %>%
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term=='diff.LGP.value' & p.value<.05)

##




################# experimental ########
plotheatmap<-function(d,cmn,nd,LGPn)
{
  dd<-d %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name==cmn
         & relpos=='ipsi' 
         & node == nd 
         & LGP.name==LGPn
         & tp>=0) %>% as_tibble()

d3.lag1<-dd %>% 
  group_by(numID,LGP.name,clin.meas.name,node,relpos) %>% 
  arrange(tp) %>% 
  dplyr::mutate(LGP.value = LGP.value - dplyr::lag(LGP.value, n = 1, default = NA)
                ,clin.meas.value = clin.meas.value - dplyr::lag(clin.meas.value, n = 1, default = NA)
                ,tp = (tp + dplyr::lag(tp, n = 1, default = NA))/2)
d3.lag2<-dd %>% 
  group_by(numID,LGP.name,clin.meas.name,node,relpos) %>% 
  arrange(tp) %>% 
  dplyr::mutate(LGP.value = LGP.value - dplyr::lag(LGP.value, n = 2, default = NA)
                ,clin.meas.value = clin.meas.value - dplyr::lag(clin.meas.value, n = 2, default = NA)
                ,tp = (tp + dplyr::lag(tp, n = 2, default = NA))/2)
d3.lag3<-dd %>% 
  group_by(numID,LGP.name,clin.meas.name,node,relpos) %>% 
  arrange(tp) %>% 
  dplyr::mutate(LGP.value = LGP.value - dplyr::lag(LGP.value, n = 3, default = NA)
                ,clin.meas.value = clin.meas.value - dplyr::lag(clin.meas.value, n = 3, default = NA)
                ,tp = (tp + dplyr::lag(tp, n = 3, default = NA))/2)

d3<-bind_rows(d3.lag1,d3.lag2, d3.lag3) %>% 
  filter(!is.na(tp)) %>% 
  bind_rows(dd) %>% 
  arrange(numID,LGP.name,tp)


p.list=c()
for( t1 in unique(d3$tp)){ #LGP
  for( t2 in unique(d3$tp)){ # clinical
    if(t1>t2){p.list<-c(p.list,NA); next()}
    x<-d3 %>% filter(tp==t1) %>% select(-c(tp,clin.meas.value))
    y<-d3 %>% filter(tp==t2) %>% select(-c(tp,LGP.value))
    df<-merge(x,y) %>% as_tibble()
    m<-df %>% 
      group_by(clin.meas.name,LGP.name,relpos,node) %>% 
      nest() %>% 
      mutate(mdl=map(data,~lm(clin.meas.value~LGP.value,data=., na.action = na.omit))) %>% 
      mutate(tidy=map(mdl,tidy))
    p<-m %>% 
      unnest(tidy) %>% 
      filter(term=='LGP.value') %>% 
      pull('p.value')
    p.list<-c(p.list,p)
    
    next()
    m<-lm(y~x,df,na.action = na.omit) %>% tidy
    p<-m%>% filter(term=='x') %>% pull('p.value')
    beta<-m%>% filter(term=='x') %>% pull('estimate')
    if(p<0.05){sprintf('t1=%f, t2=%f; p=%f, beta=%f',t1,t2,p,beta) %>% print()}
  } 
}

plt<-matrix(p.list, ncol = length(unique(d3$tp))) %>% 
  apply(2, function(y) as.numeric(y > 0.01)) %>% 
  heatmap.2(Rowv = NA, Colv = NA, dendrogram = 'none', trace='none', key=F
            , labRow = unique(d3$tp), labCol = unique(d3$tp)
            , symm=F,symkey=F,symbreaks=T, scale="none")
return(plt)
}

plt<-plotheatmap(d,'NIHSS','thalamus','efficiency')
