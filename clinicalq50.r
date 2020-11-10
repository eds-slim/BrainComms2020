d<-data.clinical_long %>% 
  mutate(clin.meas.value=if_else(clin.meas.name=='FM',66-clin.meas.value,clin.meas.value)) %>% 
  merge(d.conn) %>% 
  filter(clin.meas.name %in% c('NIHSS','FM','rGS') & relpos=='ipsi' & tp>=1 & numID != 40)
  
d.GC.clinical<-d %>% 
  group_by(conn.meas.name,tp, clin.meas.name) %>% 
  arrange(clin.meas.name,tp) %>% 
  nest() %>% 
  mutate(mdl=pmap(list(data,clin.meas.name),~glm(I(clin.meas.value)~1+conn.meas.value.diff
                           , data=.x
                           , family = if_else(.y %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
                           , na.action = na.exclude)))

tbl.GC.clinical <- d.GC.clinical %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term == 'conn.meas.value.diff') %>% 
  mutate(stats=cell_spec(sprintf('%1.2f +/- %1.2f, p=%1.4f', estimate, std.error, p.value), format = 'html', bold = p.value<0.05)) %>% 
  dplyr::select(-c(term, estimate,std.error,statistic,p.value)) %>%
  unite(temp,conn.meas.name,tp) %>% 
  spread(temp,stats)


p.GC.clinical.list <- d.GC.clinical %>% 
  mutate(data=pmap(list(data,mdl),~add_preds(.x,.y))) %>% 
  ungroup() %>% 
  group_by(conn.meas.name,clin.meas.name) %>% 
  nest() %>% 
  mutate(mdl=map(data,~.$mdl[[1]])) %>% 
  mutate(data=map(data,function(d){d %>% dplyr::select(-mdl) %>%  unnest(data)})) %>% 
  mutate(plt=pmap(list(data,conn.meas.name,clin.meas.name, mdl),~plot_preds('conn.meas.value.diff',..1,..2,..3,..4)))


p.GC.clinical.time<-do.call(arrangeGrob,c(p.GC.clinical.list$plt, ncol=3, as.table=T)) 
p.GC.clinical.time %>% 
  grid::grid.draw()


d.GC.clinical.joint <- d %>% 
  group_by(conn.meas.name,clin.meas.name) %>% 
  nest() %>% 
  mutate(plt= pmap(list(data,conn.meas.name,clin.meas.name)
                   ,~joint_fit_CV('conn.meas.value.diff',..1,..2,..3
                                  , fam = if_else(..3 %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
                                  )
                   )
         ) %>% 
  mutate(mdl= pmap(list(data,conn.meas.name,clin.meas.name)
                   ,~joint_fit_CV_mdl('conn.meas.value.diff',..1,..2,..3
                                      , fam = if_else(..3 %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
                                      )
                   )
         )
  


p.GC.clinical.joint <- do.call(arrangeGrob,c(d.GC.clinical.joint$plt
                                             , ncol=length(unique(d.GC.clinical.joint$clin.meas.name))
                                             , as.table=F)
                               )

p.GC.clinical.joint %>% 
  grid::grid.draw()


tbl.GC.clinical <- d.GC.clinical.joint %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term == 'fitted(fm)') %>% 
  mutate(stats=cell_spec(sprintf('%1.2f +/- %1.2f, p=%1.4f', estimate, std.error, p.value), format = 'html', bold = p.value<0.05)) %>% 
  dplyr::select(-c(term, estimate,std.error,statistic,p.value)) %>% 
  spread(conn.meas.name,stats) %>% 
  merge(tbl.GC.clinical) %>% 
  kable(format = 'html', escape = F) %>% 
  kable_styling(bootstrap_options = c('striped','bordered'))
tbl.GC.clinical
