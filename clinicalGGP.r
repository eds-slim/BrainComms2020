d<-data.clinical_long %>% 
  mutate(clin.meas.value=if_else(clin.meas.name=='FM',66-clin.meas.value,clin.meas.value)) %>% 
  merge(d.GGP) %>% 
  filter(clin.meas.name %in% c('NIHSS','FM','rGS') & relpos=='ipsi' & tp>=1  & threshold == 1)


d.GGP.clinical<-d %>% 
  group_by(lab,tp, clin.meas.name) %>% 
  arrange(clin.meas.name,tp) %>% 
  nest() %>% 
  mutate(mdl=pmap(list(data,clin.meas.name),~glm(I(clin.meas.value)~1+GGP
                                                 , data=.x
                                                 , family = if_else(.y %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
                                                 , na.action = na.exclude)))

tbl.GGP.clinical.time <- d.GGP.clinical %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  select(tidy) %>% unnest(tidy) %>% 
  filter(term == 'GGP') %>% 
  ungroup() %>% 
  mutate(lab=as_labeller(GGP_names)(lab) %>% as.character()) %>% 
  dplyr::select(-c(term,statistic)) %>% 
  gather(statsname,value,estimate,std.error,p.value) %>% 
  unite(temp,statsname,tp) %>% 
  spread(temp,value) %>% 
  dplyr::select(c(clin.meas.name,lab,estimate_1,std.error_1,p.value_1,estimate_3,std.error_3,p.value_3,estimate_12,std.error_12,p.value_12)) %>% 
  arrange(clin.meas.name,lab)


tbl.GGP.clinical.time %>% kable(format = 'html') %>% kable_styling(bootstrap_options = 'bordered') %>% collapse_rows()

p.GGP.clinical.list <- d.GGP.clinical %>% 
  mutate(data=pmap(list(data,mdl),~add_preds(.x,.y))) %>% 
  ungroup() %>% 
  mutate(clin.meas.name = forcats::fct_relevel(as.factor(clin.meas.name), 'NIHSS', 'rGS')) %>% 
  group_by(lab,clin.meas.name) %>% 
  arrange(lab, clin.meas.name) %>% 
  nest() %>% 
  mutate(mdl=map(data,~.$mdl[[1]])) %>% 
  mutate(data=map(data,function(d){d %>% dplyr::select(-mdl) %>%  unnest(data)})) %>% 
  mutate(plt=pmap(list(data,lab,clin.meas.name,mdl),~plot_preds('GGP',..1,as_labeller(GGP_names)(..2),as.character(..3),..4)))


p.GGP.clinical.time<-do.call(arrangeGrob,c(p.GGP.clinical.list$plt, ncol=3, as.table=F)) 
p.GGP.clinical.time %>% 
  grid::grid.draw()
ggsave('./../Word/Floats/Figures/S7-GGPclin-time-cs.tiff', plot = p.GGP.clinical.time %>% as_ggplot()
       , device = 'tiff', width = 18.5, height = 18.5, units = 'cm', dpi = 1200)

d.GGP.clinical.joint <- d %>% 
  group_by(lab,clin.meas.name) %>%
  nest() %>% 
  mutate(plt= pmap(list(data,lab,clin.meas.name)
                   ,~joint_fit_CV('GGP.diff',..1,as_labeller(GGP_names)(..2),..3
                                  , fam = if_else(..3 %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
                   )
  )
  ) %>% 
  mutate(mdl= pmap(list(data,lab,clin.meas.name)
                   ,~joint_fit_CV_mdl('GGP.diff',..1,as_labeller(GGP_names)(..2),..3
                                      , fam = if_else(..3 %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
                   )
  )
  ) %>% 
  ungroup() %>% 
  mutate(clin.meas.name = forcats::fct_relevel(as.factor(clin.meas.name), 'NIHSS','rGS','FM')) %>% 
  arrange(clin.meas.name,lab)
  


p.GGP.clinical.joint <- do.call(arrangeGrob,c(d.GGP.clinical.joint$plt
                                             , ncol=length(unique(d.GGP.clinical.joint$clin.meas.name))
                                             , as.table=T)
)

p.GGP.clinical.joint %>% 
  grid::grid.draw()
ggsave('./../Word/Floats/Figures/S8-GGPclin.tiff', plot = p.GGP.clinical.joint %>% as_ggplot()
       , device = 'tiff', width = 18.5, height = 18.5, units = 'cm', dpi = 1200)




tbl.GGP.clinical.preformat <- d.GGP.clinical.joint %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  select(clin.meas.name,lab,tidy) %>% unnest(tidy) %>% 
  filter(term == 'fitted(fm)') %>% 
  dplyr::select(-c(statistic,term)) %>% 
  ungroup() %>% 
  mutate(lab=as_labeller(GGP_names)(lab) %>% as.character()) %>% 
  dplyr::select(clin.meas.name,lab,everything()) %>% 
  arrange(clin.meas.name,lab) %>% 
  bind_cols(tbl.GGP.clinical.time %>% dplyr::select(estimate_1:p.value_12))
#  merge(tbl.GGP.clinical.time)
tbl.GGP.clinical.preformat %>% kable(format = 'html') %>% kable_styling(bootstrap_options = 'striped') %>% collapse_rows()



######################################
## volume corrected tables ###########
######################################

d.GGP.clinical.vol<-d %>% 
  group_by(lab,tp, clin.meas.name) %>% 
  arrange(clin.meas.name,tp) %>% 
  nest() %>% 
  mutate(mdl=pmap(list(data,clin.meas.name),~glm(
                                                #as.formula(paste0('I(clin.meas.value)~ GGP.diff + initial', clin.meas.name))    
                                                 I(clin.meas.value)~ GGP.diff + logvol
                                                 , data=.x
                                                 , family = if_else(.y %in% c('2NIHSS','2FM'), 'quasipoisson','gaussian')
                                                 , na.action = na.exclude)))

tbl.GGP.clinical.vol.time <- d.GGP.clinical.vol %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  select(tidy) %>% unnest(tidy) %>% 
  filter(term == 'GGP.diff') %>% 
  ungroup() %>% 
  mutate(lab=as_labeller(GGP_names)(lab) %>% as.character()) %>% 
  dplyr::select(-c(term,statistic)) %>% 
  gather(statsname,value,estimate,std.error,p.value) %>% 
  unite(temp,statsname,tp) %>% 
  spread(temp,value) %>% 
  dplyr::select(c(clin.meas.name,lab,estimate_1,std.error_1,p.value_1,estimate_3,std.error_3,p.value_3,estimate_12,std.error_12,p.value_12)) %>% 
  arrange(clin.meas.name,lab)


tbl.GGP.clinical.vol.time %>% kable(format = 'html') %>% kable_styling(bootstrap_options = 'bordered') %>% collapse_rows()


d.GGP.clinical.vol.joint <- d %>% 
  group_by(lab,clin.meas.name) %>% 
  nest() %>% 
  mutate(mdl= pmap(list(data,lab,clin.meas.name)
                   ,~joint_fit_CV_mdl('GGP.diff',..1,as_labeller(GGP_names)(..2),..3
                                      , fam = if_else(..3 %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
                                      , covs = ' + tp + logvol'
                                      #, covs = paste0(' + tp + initial', clin.meas.name)
                                      
                                      )
                   )
  )



tbl.GGP.clinical.vol.preformat <- d.GGP.clinical.vol.joint %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  select(tidy) %>% unnest(tidy) %>% 
  filter(term == 'fitted(fm)') %>% 
  dplyr::select(-c(statistic,term)) %>% 
  ungroup() %>% 
  mutate(lab=as_labeller(GGP_names)(lab) %>% as.character()) %>% 
  dplyr::select(clin.meas.name,lab,everything()) %>% 
  arrange(clin.meas.name,lab) %>% 
  bind_cols(tbl.GGP.clinical.vol.time %>% dplyr::select(estimate_1:p.value_12))
#  merge(tbl.GGP.clinical.time)
tbl.GGP.clinical.vol.preformat %>% kable(format = 'html') %>% kable_styling(bootstrap_options = 'striped') %>% collapse_rows()


