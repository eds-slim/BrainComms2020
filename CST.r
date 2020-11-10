CST <- read.csv('./../SFB_structural_connectomes/CST_FA_55_w0001_onlya.csv')

d.GGP.fa <- CST %>% 
  as_tibble() %>% 
  mutate(ID = substr(ID,2,3) %>% as.numeric()) %>% 
  mutate(status = recode(status, 'healthy'='contra','stroke'='ipsi')) %>% 
  mutate(tp = recode(tp, '1'='0','2'='1','3'='3','4'='12') %>% as.double()) %>% 
  dplyr::rename(relpos = status) %>% 
  merge(data.patients) %>% 
  mutate(threshold=1)

d.GGP.fa <-  d.GGP.fa %>%   
  complete(numID,tp,relpos) %>% 
  group_by(numID,relpos) %>% 
  arrange(tp) %>% 
  mutate(fa.diff=fa-fa[tp==0]) %>% 
  merge(d.GGP)


data.clinical_long <- data.clinical_long %>%   
  group_by(numID, clin.meas.name) %>% 
  arrange(tp) %>% 
  mutate(clin.meas.value.diff=clin.meas.value.diff-clin.meas.value.diff[tp==0])


d<-data.clinical_long %>% 
  mutate(clin.meas.value=if_else(clin.meas.name=='FM',66-clin.meas.value,clin.meas.value)) %>% 
  merge(d.GGP.fa) %>% 
  filter(clin.meas.name %in% c('NIHSS','FM','rGS') & relpos=='ipsi' & tp>=0  & threshold == 1)


d.GGP.clinical<-d %>% 
  filter(clin.meas.name=='rGS') %>% 
  group_by(lab,tp, clin.meas.name) %>% 
  arrange(clin.meas.name,tp) %>% 
  nest() %>% 
  mutate(mdl=pmap(list(data,clin.meas.name),~glm(I(clin.meas.value.diff)~fa.diff
                                                 , data=.x
                                                 , family = if_else(.y %in% c('NIHSS','FM'), 'quasipoisson','gaussian')
                                                 , na.action = na.exclude)))

d.GGP.clinical %>% 
  mutate(tidy=map(mdl,tidy)) %>%
  unnest(tidy) %>% 
  filter(term!='(Intercept)') %>% 
  arrange(lab) %>% 
  print(n=Inf)
  

d %>% 
  filter(relpos=='ipsi' & tp>=1) %>% 
ggplot(aes(x=fa.diff, y=clin.meas.value.diff, color=as.factor(tp)))+
  geom_point() +
  geom_smooth(method = 'lm')+
  facet_wrap(lab~clin.meas.name, scales = 'free')
