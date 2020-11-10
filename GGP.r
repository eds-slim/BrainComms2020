## relation clinical

stats.GGP.clinical<-d.GGP %>% 
  merge(data.clinical_long) %>% 
  group_by(clin.meas.name, lab, threshold) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme)(GGP~deriv_mod(tp,a,b,del)
                                    , fixed=list(a+b~1,del~clin.meas.value*relpos)
                                    , random=a~1|numID
                                    , data=.
                                    , start=c(1,1,1,1,1,1)
                                    , na.action = na.omit
  )$result
  )
  ) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))




tbl.GGP.clinical<-stats.GGP.clinical %>%
  filter(threshold==1) %>% 
  unnest(tidy) %>% 
  select(-c(effect,group))  %T>% 
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4))
tbl.GGP.clinical %>% 
  write.table(file = './../Word/Floats/Tables/13-GGP_clinical.tab', sep = ",", quote = FALSE, row.names = F)


stats.GGP.clinical.relpos<-d.GGP %>% 
  merge(data.clinical_long) %>% 
  group_by(clin.meas.name,relpos, lab, threshold) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme)(GGP~deriv_mod(tp,a,b,del)
                                    , fixed=list(a~1,b~1,del~clin.meas.value)
                                    , random=a~1|numID
                                    , data=.
                                    , start=c(1,1,1,1)
                                    , na.action = na.omit
  )$result
  )
  ) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))

tbl.GGP.relpos <- stats.GGP.clinical.relpos %>% 
  filter(threshold==1 & relpos == 'ipsi') %>% 
  unnest(tidy) %>% 
  dplyr::select(-c(effect,group)) %T>% 
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4))

tbl.GGP.relpos %>% 
  write.table(file = './../Word/Floats/Tables/14-GGP_clinical_by_relpos.tab', sep = ",", quote = FALSE, row.names = F)
