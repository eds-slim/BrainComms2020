d<-data.clinical_long %>% 
  mutate(clin.meas.value=if_else(clin.meas.name=='FM',66-clin.meas.value,clin.meas.value)) %>% 
  merge(d.LGP) %>% 
  filter(clin.meas.name %in% c('NIHSS','FM','rGS') & relpos=='ipsi' & tp>=1)

d.LGP.clinical.joint <- d %>% 
  filter(node %in% unique(nodes$node) & is.finite(LGP.value.diff)) %>% 
  group_by(node, LGP.name,clin.meas.name) %>% 
  nest() %>% 
  mutate(plt= pmap(list(data,node,clin.meas.name)
                   ,~joint_fit_CV('LGP.value.diff',..1,..2,..3
                                  , fam = if_else(..3 %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
                   )
  )
  ) %>% 
  mutate(mdl= pmap(list(data,node,clin.meas.name)
                   ,~joint_fit_CV_mdl('LGP.value.diff',..1,..2,..3
                                      , fam = if_else(..3 %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
                                      , covs = ''
                   )
  )
  )

d.LGP.clinical.vol.joint <- d %>% 
  filter(node %in% unique(nodes$node) & is.finite(LGP.value.diff)) %>% 
  group_by(node, LGP.name,clin.meas.name) %>% 
  nest() %>% 
  mutate(plt= pmap(list(data,node,clin.meas.name)
                   ,~joint_fit_CV('LGP.value.diff',..1,..2,..3
                                  , fam = if_else(..3 %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
                   )
  )
  ) %>% 
  mutate(mdl= pmap(list(data,node,clin.meas.name)
                   ,~joint_fit_CV_mdl('LGP.value.diff',..1,..2,..3
                                      , fam = if_else(..3 %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
                                      , covs = ' + logvol'
                   )
  )
  )

tbl.LGP.joint.preformat <- d.LGP.clinical.joint %>% 
  #bind_rows(d.LGP.clinical.vol.joint %>% mutate(clin.meas.name = paste0('vol_', clin.meas.name))) %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term == 'fitted(fm)' & LGP.name %in% c('strength','efficiency','clustering') & p.value < 0.05) %>% 
  dplyr::select(-c(term, statistic)) %>% 
  gather(statsname,value,estimate,std.error,p.value) %>% 
  mutate(value = as.numeric(value)) %>% 
  unite(temp,statsname,clin.meas.name) %>% 
  spread(temp,value,fill = NA) %>% 
  dplyr::select(node,LGP.name, estimate_NIHSS, std.error_NIHSS, p.value_NIHSS
                             , estimate_FM, std.error_FM, p.value_FM
                            , estimate_rGS, std.error_rGS, p.value_rGS) %>% 
  arrange(node)

tbl.LGP.vol.joint.preformat <- d.LGP.clinical.vol.joint %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term == 'fitted(fm)' & LGP.name %in% c('strength','efficiency','clustering') & p.value < 0.05) %>% 
  dplyr::select(-c(term, statistic)) %>% 
  gather(statsname,value,estimate,std.error,p.value) %>% 
  mutate(value = as.numeric(value)) %>% 
  unite(temp,statsname,clin.meas.name) %>% 
  spread(temp,value,fill = NA) %>% 
  dplyr::select(node,LGP.name, estimate_NIHSS, std.error_NIHSS, p.value_NIHSS
                , estimate_FM, std.error_FM, p.value_FM
                , estimate_rGS, std.error_rGS, p.value_rGS) %>% 
  arrange(node)

tbl.LGP.vol.joint.preformat
