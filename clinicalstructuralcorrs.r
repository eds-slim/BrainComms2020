## q50
d<-data.clinical_long %>% 
  merge(d.conn) %>%
  filter(relpos=='ipsi' & numID !=40 & tp>=3) %>% ## conn.meas.value.diff at t=12 way too small (outlier)
  group_by(clin.meas.name,conn.meas.name,relpos) %>% 
  nest()
tbl.GC.clinical.chronic.gaussian<-d %>% 
  filter(clin.meas.name=='rGS') %>% 
  mutate(mdl=map(data,~glm(clin.meas.value~conn.meas.value.diff+tp, data=., family = 'gaussian'))) %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term=='conn.meas.value.diff') %>% arrange(clin.meas.name)
tbl.GC.clinical.chronic.quasipoisson<-d %>% 
  filter(clin.meas.name %in% c('NIHSS','FM')) %>% 
  mutate(mdl=map(data,~glm(clin.meas.value~conn.meas.value.diff+tp, data=., family = 'quasipoisson'))) %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term=='conn.meas.value.diff') %>% arrange(clin.meas.name)


tbl.GC.clinical.vol.chronic.gaussian<-d %>% 
  filter(clin.meas.name=='rGS') %>% 
  mutate(mdl=map(data,~glm(clin.meas.value~conn.meas.value.diff+tp+logvol, data=., family = 'gaussian'))) %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term=='conn.meas.value.diff') %>% arrange(clin.meas.name)
tbl.GC.clinical.vol.chronic.quasipoisson<-d %>% 
  filter(clin.meas.name %in% c('NIHSS','FM')) %>% 
  mutate(mdl=map(data,~glm(clin.meas.value~conn.meas.value.diff+tp+logvol, data=., family = 'quasipoisson'))) %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term=='conn.meas.value.diff') %>% arrange(clin.meas.name)

## GGP

d<-data.clinical_long %>% 
  merge(d.GGP) %>%
  filter(relpos=='ipsi' & numID !=40 & tp==3) %>% ## conn.meas.value.diff at t=12 way too small (outlier)
  group_by(clin.meas.name,lab, threshold,relpos) %>% 
  nest()
tbl.GGP.clinical.chronic.gaussian<-d %>% 
  filter(clin.meas.name=='rGS') %>% 
  mutate(mdl=map(data,~glm(clin.meas.value~GGP.diff+tp, data=., family = 'gaussian'))) %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term=='GGP.diff') %>% arrange(clin.meas.name)
tbl.GGP.clinical.chronic.quasipoisson<-d %>% 
  filter(clin.meas.name %in% c('NIHSS','FM')) %>% 
  mutate(mdl=map(data,~glm(clin.meas.value~GGP.diff+tp, data=., family = 'quasipoisson'))) %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term=='GGP.diff') %>% arrange(clin.meas.name)

tbl.GGP.clinical.vol.chronic.gaussian<-d %>% 
  filter(clin.meas.name=='rGS') %>% 
  mutate(mdl=map(data,~glm(clin.meas.value~GGP.diff+tp+logvol, data=., family = 'gaussian'))) %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term=='GGP.diff') %>% arrange(clin.meas.name)
tbl.GGP.clinical.vol.chronic.quasipoisson<-d %>% 
  filter(clin.meas.name %in% c('NIHSS','FM')) %>% 
  mutate(mdl=map(data,~glm(clin.meas.value~GGP.diff+tp+logvol, data=., family = 'quasipoisson'))) %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term=='GGP.diff') %>% arrange(clin.meas.name)

#### LGP ####
data.clinical_long %>% 
  merge(d.LGP) %>%
  filter(relpos == 'ipsi' & numID != 40 & is.finite(LGP.value.diff) & is.finite(clin.meas.value) & LGP.name %in% c('strength','efficiency','clustering')) %>% 
  filter(node %in% (nodes.strict$node %>% as.character())) %>% 
  group_by(clin.meas.name,LGP.name, node, relpos, tp) %>% 
  nest() %>% 
  mutate(mdl=map(data,~glm(clin.meas.value~1+LGP.value.diff+logvol,data=., family = "quasipoisson"))) %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(term=='LGP.value.diff' & p.value<0.05  & estimate < 1000) %>% 
  arrange(clin.meas.name,LGP.name,tp,node) %>% 
  print(n=Inf)


d<-data.clinical_long %>% 
  merge(d.LGP) %>%
  filter(relpos == 'ipsi' & numID != 40 & is.finite(LGP.value.diff) & is.finite(clin.meas.value) & LGP.name %in% c('strength','efficiency','clustering')) %>% 
  filter(clin.meas.name=='NIHSS' & LGP.name=='strength' & node=='thalamus' & relpos=='ipsi' & tp==12)
m1<-d %>% 
  glm(clin.meas.value~1+LGP.value.diff,data=., family = "poisson")
m2<-d %>% 
  glm(clin.meas.value~1+LGP.value.diff,data=., family = "quasipoisson")

AIC(m1)
AIC(m2)
with(m1, cbind(res.deviance = deviance, df = df.residual,
               p = pchisq(deviance, df.residual, lower.tail=FALSE)))
AER::dispersiontest(m1)


d<-data.clinical_long %>% 
  merge(d.GGP) %>% 
  filter(clin.meas.name=='NIHSS' & lab==4 & threshold==1 & relpos=='ipsi' & tp==12 & numID != 40)


add_preds<-function(d,mdl){
  d$pred<-predict(mdl, type = 'response', se.fit = T)$fit
  d$pred.se<-predict(mdl, type = 'response', se.fit = T)$se.fit
  return(d)
}
plot_preds<-function(d,lab,cmn){
  if(lab>=1){
    s<-scale_y_continuous(cmn)
    xlab<-paste0('Change in ',GGP_names[lab])
  }else
  {
    #s<-scale_y_continuous('',breaks = NULL)
    xlab<-paste0('... ',GGP_names[lab])
  }
  d %>% 
    ggplot(aes(x=GGP.diff, y=clin.meas.value, color=as.factor(tp)))+
    geom_point()+
    geom_line(aes(y=pred))+
    geom_line(aes(y=pred+pred.se), linetype=2, alpha=.32)+
    geom_line(aes(y=pred-pred.se), linetype=2, alpha=.32)+
    s+
    scale_x_continuous(xlab)+
    scale_alpha_manual(values=c(.5,1))+
    #ggtitle(lab)+
    guides(alpha=F, color=F)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())%>% 
    return()
}

d<-data.clinical_long %>% 
  merge(d.conn) %>% 
  filter(clin.meas.name=='NIHSS' & relpos=='ipsi' & tp>=1 & numID != 40) %>% 
  group_by(tp,clin.meas.name) %>% 
  arrange(tp) %>% 
  nest() %>% 
  mutate(mdl=map(data,~glm(I(66-clin.meas.value)~1+GGP.diff,data=., family = "quasipoisson", na.action = na.exclude))) %>% 
  mutate(data=pmap(list(data,mdl),~add_preds(.x,.y))) %>% 
  mutate(plt=pmap(list(data,tp,clin.meas.name),~plot_preds(..1,..2,..3)))

d %>% mutate(tidy=map(mdl,tidy)) %>% unnest(tidy)
d$plt

do.call(arrangeGrob,c(d$plt, ncol=3, as.table=F)) %>% 
  grid::grid.draw()


d<-data.clinical_long %>% 
  merge(d.GGP) %>% 
  mutate(clin.meas.value=if_else(clin.meas.name=='FM',66-clin.meas.value,clin.meas.value)) %>% 
  filter(clin.meas.name %in% c('NIHSS','FM'), threshold==1 & relpos=='ipsi' & lab %in% c(1,2,3) & tp>=1 & numID != 40) %>% 
  group_by(lab,tp, clin.meas.name) %>% 
  arrange(clin.meas.name,tp, lab) %>% 
  nest() %>% 
  mutate(mdl=map(data,~glm(I(clin.meas.value)~1+GGP.diff,data=., family = "quasipoisson", na.action = na.exclude))) %>% 
  mutate(data=pmap(list(data,mdl),~add_preds(.x,.y))) %>% 
  ungroup() %>% 
  group_by(lab,clin.meas.name) %>% 
  nest() %>% 
  mutate(data=map(data,function(d){d %>% select(-mdl) %>%  unnest(data)})) %>% 
  mutate(plt=pmap(list(data, lab, clin.meas.name),~plot_preds(..1,..2,..3)))


d %>% mutate(tidy=map(mdl,tidy)) %>% unnest(tidy)
d$plt

do.call(arrangeGrob,c(d$plt, ncol=3, as.table=T)) %>% 
  grid::grid.draw()



d<-data.clinical_long %>% 
  merge(d.LGP) %>% 
  filter(clin.meas.name=='rGS'& relpos=='ipsi' & tp>=1 & numID != 40 & is.finite(LGP.value) & node %in% nodes.liberal$node) %>% 
  group_by(node,LGP.name,tp) %>% 
  arrange(node,tp) %>% 
  nest() %>% 
  mutate(mdl=map(data,~glm(I(clin.meas.value)~1+LGP.value.diff,data=., family = "gaussian", na.action = na.exclude))) %>% 
  mutate(data=pmap(list(data,mdl),~add_preds(.x,.y))) %>% 
  mutate(plt=pmap(list(data,node),~plot_preds(.x,.y)))

d %>% mutate(tidy=map(mdl,tidy)) %>% unnest(tidy) %>% 
  filter(p.value<0.05 & term=='LGP.value.diff')
d$plt

do.call(arrangeGrob,c(d$plt, ncol=3, as.table=F)) %>% 
  grid::grid.draw()



data.clinical_long %>% 
  filter(tp==12) %>%
  select(-tp) %>% 
  merge(d.GGP %>% filter(tp==3) %>% select(-tp)) %>% as_tibble() %>% 
  filter(clin.meas.name=='NIHSS'& relpos=='ipsi' & numID != 40 & lab == 3 & threshold == 1) %>% 
  glm(clin.meas.value~GGP.diff, data=., family = 'quasipoisson') %>% 
  tidy()



d<-data.clinical_long %>% 
  merge(d.GGP) %>% 
  filter(clin.meas.name=='FM'& relpos=='ipsi' & tp>=1 & !(numID == 40 & tp%in%c(1,3,12))  & lab==1 & threshold==1)
  mutate(clin.meas.value = clin.meas.value/67)

d %>% 
  ggplot(aes(x=tp,y=GGP.diff))+
  geom_point()+
  geom_smooth(aes(group=numID), method = 'lm', se = F)

fm <- lme(GGP.diff ~ tp, data = d, random = ~ 1 | numID, na.action = na.exclude)
mdl<-glm(clin.meas.value ~ fitted(fm)+logvol, data = d, family = 'quasipoisson', na.action = na.exclude)
mdl<-glm(clin.meas.value>5 ~ fitted(fm)+tp, family='binomial', data = d, na.action = na.omit)

summary(mdl)
d$pred<-fitted(fm)

preds<-predict(mdl, type = 'response', se.fit = T)
d$.fitted<-preds
d$.fitted<-preds$fit
d$.se<-preds$se.fit

d %>% 
  ggplot(aes(x=pred,y=clin.meas.value, color=as.factor(tp)))+
  geom_point()+
  geom_line(aes(y=.fitted), color='black')
  #geom_line(aes(y=.fitted+.se), color='black', linetype=2)+
  #geom_line(aes(y=.fitted-.se), color='black', linetype=2)
