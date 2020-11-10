
stats.LGP<-d.LGP %>%
  group_by(LGP.name, node) %>%
  nest() %>%
  mutate(mdl=map(data,~safely(nlme.exp)(., response = 'LGP.value', preds=list(a='relpos', b='1',del='relpos'))$result)) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result), AIC = map(mdl,~safely(AIC)(.)$result %>% enframe())) %>% 
  select(tidy) %>% unnest(tidy,.drop = T)

stats.LGP.ipsi<-d.LGP %>%
  filter(relpos=='ipsi') %>% 
  group_by(LGP.name, node, relpos) %>%
  nest() %>%
  mutate(mdl=map(data,~safely(nlme.exp)(., response = 'LGP.value', preds=list(a='1', b='1',del='1'))$result)) %>%
  mutate(tidy = map(mdl, ~safely(tidy)(., conf.int=T)$result)) %>% 
  select(tidy) %>% unnest(tidy,.drop = T)

nodes <-stats.LGP.ipsi %>%  filter(term=='del' & p.value<0.01 & estimate<0) %>% ungroup() %>% dplyr::select(c(LGP.name,node))


tbl.LGP.preformat<-stats.LGP %>% 
  dplyr::select(-c(group, effect, df, statistic, conf.low, conf.high)) %>% 
  filter(LGP.name %in% c('strength','efficiency','clustering') & term %in% c('a.relposcontra', 'del.relposcontra') & node != 'lentiform') %>% 
  group_by(node,LGP.name) %>% 
  nest() %>% 
  mutate(p.min = map_dbl(data,~min(.$p.value))) %>% 
  filter(p.min<0.05) %>% dplyr::select(-p.min) %>% 
  unnest(data) %>%
  group_by(node) %>% 
  nest() %>% 
  mutate(p.min = map_dbl(data,~min(.$p.value))) %>% 
  unnest(data) %>% 
  arrange(p.min) %>% dplyr::select(-p.min) %>% dplyr::select(node,LGP.name,everything()) %>% arrange(node,LGP.name)

#tbl.LGP.preformat %>% kable(format = 'html') %>% kable_styling() %>% collapse_rows()

tbl.LGP.ipsi.preformat <- stats.LGP.ipsi %>% 
  dplyr::select(-c(group, effect, df, statistic, conf.low, conf.high)) %>% 
  filter(LGP.name %in% c('strength','efficiency','clustering') & term %in% c('del') & p.value<0.05 & node != 'lentiform') %>% 
  group_by(node) %>% 
  nest() %>% 
  mutate(p.min = map_dbl(data,~min(.$p.value))) %>% 
  unnest(data) %>% 
  arrange(p.min) %>% dplyr::select(-p.min,-relpos,-term)

tbl.LGP.ipsi.preformat %>% 
  kable(format='html', escape=F, caption = 'Stats for modelling of time course of LGPs in ipsilesional hemispheres-') %>% 
  kable_styling(bootstrap_options = c('striped','bordered')) %>% 
  collapse_rows(columns = c(1))


