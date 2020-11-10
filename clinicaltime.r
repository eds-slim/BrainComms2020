pd <- position_dodge(0.4)
col.idx<-swatch()[3]
p.clinical.timecourse<-data.clinical_long %>%
  ungroup() %>% 
  filter(clin.meas.name %in% c('NIHSS','FM','rGS')) %>% 
  mutate(clin.meas.name = forcats::fct_relevel(clin.meas.name,'NIHSS','rGS')) %>% 
  ggplot(aes(x=tp, y=clin.meas.value, group=numID, shape=Stroke.side)) +
  geom_line(size=.1,position=pd, color=col.idx, alpha=0.32)+
  geom_jitter(size=.1, position=pd, color=col.idx, alpha=0.32)+
  #stat_summary(aes(group = recovery, color=recovery), geom = "point", fun.y = mean, shape = 20, size = 3,)+
  #stat_summary(aes(group = recovery, color=recovery, width=.5), size=1, fun.data = mean_se, geom = "errorbar")+
  geom_smooth(aes(group=1)
              , formula = 'y~a+d*(1-exp(-b*x))'
              , method.args = list(start=c(a=10,b=.5, d=5)
                                   , control=nls.control(minFactor=1e-10, maxiter = 100))
              , method = 'nls'
              , se = F
              , color=col.idx
              , size = .5
  )+
  stat_summary(aes(group = 1, width=.5), size=.25, fun.data = mean_se, geom = "errorbar", color=col.idx)+
  stat_summary(aes(group = 1), geom = "point", fun.y = mean, shape = 21, color='black', size = .5, fill=col.idx, stroke=.5)+
  facet_wrap(.~clin.meas.name
             , scales = 'free', ncol = 3# strip.position = "left"
             , labeller = as_labeller(c('NIHSS'='NIHSS score','rGS'='Relative grip strength','FM'='Fugl-Meyer score')))  +
  ylab(NULL) +
  sc_x_tp+
  scale_shape_manual(values = c(17,18))+
  guides(shape=F,color=F, linetype=F)+
  #thm
  theme(strip.background = element_blank(), strip.placement = "inside"
        , text = element_text(size = 8)
        , axis.text = element_text(size = 4)
        , axis.text.x = element_text(hjust = c(.75,.25,.5,.5))
        , axis.line = element_line(size = .25, colour = "black")
        , axis.ticks = element_line(size = .25, colour = 'black'))
p.clinical.timecourse
ggsave('./../Word/Floats/Figures/1-clinicaltimecourse.tiff'
       , plot = p.clinical.timecourse
       , device = 'tiff'
       , width = 9, height = 4.5, units = 'cm', dpi = 1200)


stats.clinicaltime.exp<-data.clinical_long %>%
  group_by(clin.meas.name) %>%
  nest() %>%
  mutate(mdl=map(data,~safely(nlme.exp)(.,'clin.meas.value', preds=list(a=1,b=1,del='1'))$result)) %>% 
  mutate(tidy = map(mdl, tidy), AIC = map(mdl,AIC))

AIC.exp<-stats.clinicaltime.exp %>% 
  unnest(AIC,.drop = T)

tbl.clinicaltime.exp <- stats.clinicaltime.exp %>% 
  select(tidy) %>% unnest(tidy,.drop = T) %>% 
  dplyr::select(-c(effect,group)) %T>%
  print(n=Inf) %>%
  mutate_if(is.numeric, ~signif(.,digits = 4))
tbl.clinicaltime.exp %>% 
  write.table(file = './../Word/Floats/Tables/01-clinicaltimecourse.tab', sep = ",", quote = FALSE, row.names = F)

  
stats.clinicaltime.lin<-data.clinical_long %>%
    group_by(clin.meas.name) %>%
    nest() %>%
    mutate(mdl=map(data,~lme(clin.meas.value~1+tp, random=~1|numID, data=., na.action = na.omit, method = 'ML'))) %>% 
    mutate(tidy = map(mdl, tidy), AIC = map(mdl,AIC))
  
AIC.lin<-stats.clinicaltime.lin %>% 
    unnest(AIC,.drop = T)
  
stats.clinicaltime.exp %>% 
    select(tidy) %>% unnest(tidy,.drop = T) %>% 
    dplyr::select(-c(effect,group)) %T>%
    print(n=Inf) %>%
    mutate_if(is.numeric, ~signif(.,digits = 4)) %>% 
    write.table(file = './../Word/Floats/Tables/01-clinicaltimecourse.tab', sep = ",", quote = FALSE, row.names = F)  
  
  