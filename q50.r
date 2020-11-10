## relation clinical

stats.GC.clinical<-d.conn %>% 
  merge(data.clinical_long) %>% 
  group_by(clin.meas.name) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme)(conn.meas.value~deriv_mod(tp,a,b,del)
                                    , fixed=list(a+b~1,del~clin.meas.value*relpos)
                                    , random=a~1|numID
                                    , data=.
                                    , start=c(1,1,1,1,1,1)
                                    , na.action = na.omit
  )$result
  )
  ) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))

stats.GC.clinical %>% 
  unnest(tidy) %>% 
  dplyr::select(-c(effect,group))  %T>% 
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4)) %>% 
  write.table(file = './../Word/Floats/Tables/11-GC_clinical.tab', sep = ",", quote = FALSE, row.names = F)


stats.GC.clinical.relpos<-d.conn %>% 
  merge(data.clinical_long) %>% 
  group_by(clin.meas.name,relpos) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme)(conn.meas.value~deriv_mod(tp,a,b,del)
                                    , fixed=list(a+b~1,del~clin.meas.value)
                                    , random=a~1|numID
                                    , data=.
                                    , start=c(1,1,1,1)
                                    , na.action = na.omit
  )$result
  )
  ) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe()))

tbl.GC.clinical.relpos<-stats.GC.clinical.relpos %>% 
  unnest(tidy) %>% 
  dplyr::select(-c(effect,group)) %T>% 
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4))
tbl.GC.clinical.relpos %>% 
  write.table(file = './../Word/Floats/Tables/12-GC_clinical_by_relpos.tab', sep = ",", quote = FALSE, row.names = F)



