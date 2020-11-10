stats.LGP.vol <- d.LGP %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='NIHSS'& LGP.name %in% c('strength','efficiency','clustering')) %>% 
  group_by(node, LGP.name) %>% 
  arrange(node) %>% 
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
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=F)$result))

tbl.LGP.vol.preformat <- stats.LGP.vol %>% select(tidy) %>%  unnest(tidy) %>% 
  filter(term=='del.logvol:relposcontra' & p.value<0.05 & node %in% unique(nodes$node)) %>% 
  dplyr::select(-c(effect,group)) %>% 
  group_by(node) %>%
  nest() %>% 
  mutate(p.min = map_dbl(data, ~ min(.$p.value))) %>% 
  unnest() %>% 
  dplyr::select(-c(statistic,df)) %>% 
  arrange(p.min) %>% dplyr::select(-p.min,-term) %>% arrange(node,LGP.name)
tbl.LGP.vol.preformat %>% kable(format = 'html') %>% kable_styling() %>% collapse_rows()


stats.LGP.vol.relpos <- d.LGP %>% 
  merge(data.clinical_long) %>% 
  filter(clin.meas.name=='NIHSS' & relpos == 'ipsi') %>% 
  group_by(node, LGP.name, relpos) %>% 
  arrange(node) %>% 
  nest() %>% 
  mutate(mdl=map(data,~safely(nlme)(LGP.value~deriv_mod(tp,a,b,del)
                                    , fixed=list(a+b~1,del~logvol)
                                    , random=a~1|numID
                                    , data=.
                                    , start=c(1,1,1,1)
                                    , na.action = na.omit
  )$result
  )
  ) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=F)$result))

tbl.LGP.vol.ipsi.preformat <- stats.LGP.vol.relpos %>% select(tidy) %>%  unnest(tidy) %>% 
  filter(term=='del.logvol' & p.value<0.05 & node %in% unique(nodes$node) & LGP.name %in% c('strength','efficiency','clustering')) %>% 
  dplyr::select(-c(effect,group)) %>% 
  group_by(node) %>%
  nest() %>% 
  mutate(p.min = map_dbl(data, ~ min(.$p.value))) %>% 
  unnest() %>% 
  dplyr::select(-c(relpos,term,statistic,df)) %>% 
  arrange(p.min) %>% dplyr::select(-p.min)%>% arrange(node,LGP.name)

tbl.LGP.vol.ipsi <- tbl.LGP.vol.ipsi.preformat %>%   
  kable(format = 'html', escape = F) %>% 
  kable_styling(bootstrap_options = c('striped','bordered','condensed')) %>% collapse_rows()
tbl.LGP.vol.ipsi
