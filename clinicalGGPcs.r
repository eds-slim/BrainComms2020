d.temp<-data.clinical_long %>% 
  mutate(clin.meas.value=if_else(clin.meas.name=='FM',66-clin.meas.value,clin.meas.value)) %>% 
  merge(d.GGP) %>% 
  filter(clin.meas.name %in% c('NIHSS','FM','rGS') & relpos=='ipsi' & tp>=1  & threshold == 1)

d.GGP.clinical.temp<-d.temp %>% 
  group_by(lab,tp, clin.meas.name) %>% 
  arrange(clin.meas.name,tp) %>% 
  nest() %>% 
  mutate(mdl=pmap(list(data,clin.meas.name),~glm(I(clin.meas.value)~1+GGP
                                                 , data=.x
                                                 , family = if_else(.y %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
                                                 , na.action = na.exclude)))

p.GGP.clinical.list.temp <- d.GGP.clinical.temp %>% 
  mutate(data=pmap(list(data,mdl),~add_preds(.x,.y))) %>% 
  ungroup() %>% 
  mutate(clin.meas.name = forcats::fct_relevel(as.factor(clin.meas.name), 'NIHSS', 'rGS')) %>% 
  group_by(lab,clin.meas.name) %>% 
  arrange(lab, clin.meas.name) %>% 
  nest() %>% 
  mutate(mdl=map(data,~.$mdl[[1]])) %>% 
  mutate(data=map(data,function(d){d %>% dplyr::select(-mdl) %>%  unnest(data)})) %>% 
  mutate(plt=pmap(list(data,lab,clin.meas.name,mdl),~plot_preds('GGP',..1,as_labeller(GGP_names)(..2),as.character(..3),..4)))


p.GGP.clinical.time.temp<-do.call(arrangeGrob,c(p.GGP.clinical.list.temp$plt, ncol=3, as.table=F)) 
p.GGP.clinical.time.temp %>% 
  grid::grid.draw()

ggsave('./../Word/Floats/Figures/S7-GGPclin-time-cs.tiff', plot = p.GGP.clinical.time.temp %>% as_ggplot()
       , device = 'tiff', width = 9, height = 9, units = 'cm', dpi = 1200)
