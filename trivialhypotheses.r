## clinical lag correlations
data.clinical_long %>% 
  group_by(numID,clin.meas.name) %>% 
  arrange(tp) %>% 
  dplyr::mutate(clin.meas.value.lag = dplyr::lag(clin.meas.value, n = 1, default = NA)) %>% 
  filter(tp>0) %>% 
  ggplot(aes(x=clin.meas.value, y=clin.meas.value.lag, color=as.factor(tp)))+
  geom_point()+
  geom_smooth(aes(group=tp), method = 'lm', se = F)+
  geom_smooth(aes(group=numID), method = 'lm', se = F, color='grey', size=.5)+
  facet_wrap(.~clin.meas.name, scales = 'free')


## q50 - clinical time course 1-1
data.clinical_long %>% 
  merge(d.conn) %>% 
  filter(relpos=='ipsi') %>% 
  ggplot(aes(x=conn.meas.value, y=clin.meas.value, color=as.factor(tp)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F)+
  geom_smooth(aes(group=numID), method = 'lm', se = F, color='grey', size=.5)+
  facet_wrap(conn.meas.name~clin.meas.name, scales = 'free')


## q50 - clinical time course linear delta
d<-data.clinical_long %>% 
  merge(d.conn) %>% 
  group_by(numID,clin.meas.name,conn.meas.name,relpos) %>% 
  arrange(tp) %>% 
  dplyr::mutate(conn.meas.value = conn.meas.value - dplyr::lag(conn.meas.value, n = 1, default = NA)
                ,clin.meas.value = clin.meas.value - dplyr::lag(clin.meas.value, n = 1, default = NA)
                ,tp = (tp + dplyr::lag(tp, n = 1, default = NA))/2) %>% 
  ungroup(numID) %>% 
  filter(relpos=='ipsi') %>% 
  filter(tp %in% c(0.5,2,7.5))

d %>% 
  group_by(clin.meas.name,conn.meas.name,relpos) %>% 
  nest() %>% 
  mutate(mdl=map(data,~lme(clin.meas.value~conn.meas.value, random=~1|numID, data=., na.action = na.omit))) %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(effect=='fixed')
d %>% 
  ggplot(aes(x=conn.meas.value, y=clin.meas.value, color=as.factor(tp)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F)+
  geom_smooth(aes(group=numID), method = 'lm', se = F, color='grey', size=.5)+
  facet_wrap(conn.meas.name~clin.meas.name, scales = 'free')
  

## q50 - clinical time course linear delta delay adjusted
d<-data.clinical_long %>% 
  merge(d.conn) %>% 
  group_by(numID,clin.meas.name,conn.meas.name,relpos) %>% 
  arrange(tp) %>% 
  dplyr::mutate(tp.lag = dplyr::lag(tp, n = 1, default = NA)
                ,tp.delta = tp-tp.lag
                ,tp.mean = (tp+tp.lag)/2
                ,conn.meas.value.lag = (conn.meas.value - dplyr::lag(conn.meas.value, n = 1, default = NA))/tp.delta
                ,clin.meas.value.lag = (clin.meas.value - dplyr::lag(clin.meas.value, n = 1, default = NA))/tp.delta
                ) %>% 
  ungroup(numID) %>% 
  filter(relpos=='ipsi') %>% 
  select(-c(ID,Stroke.side,age,gender,dominant_affected,vol,logvol))

d %>% 
  ggplot(aes(x=conn.meas.value.lag, y=clin.meas.value.lag, color=as.factor(tp.mean)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F)+
  geom_smooth(aes(group=numID), method = 'lm', se = F, color='grey', size=.5)+
  facet_wrap(conn.meas.name~clin.meas.name, scales = 'free')


################################################
######## same for GGPs st threshold 1 ##########
################################################



## GGP - clinical time course 1-1
d.GGP %>% 
  merge(data.clinical_long) %>% 
  merge(d.conn) %>% 
  filter(relpos=='ipsi' & threshold == 1) %>% 
  ggplot(aes(x=GGP, y=clin.meas.value, color=as.factor(tp)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F)+
  geom_smooth(aes(group=numID), method = 'lm', se = F, color='grey', size=.5)+
  facet_wrap(lab~clin.meas.name, scales = 'free')


## GGP - clinical time course linear delta
d<-data.clinical_long %>% 
  merge(d.GGP) %>% 
  group_by(numID,clin.meas.name,lab,relpos) %>% 
  arrange(tp) %>% 
  dplyr::mutate(clin.meas.value = clin.meas.value - dplyr::lag(clin.meas.value, n = 1, default = NA)
                ,GGP = GGP - dplyr::lag(GGP, n = 1, default = NA)
                ,tp = (tp + dplyr::lag(tp, n = 1, default = NA))/2) %>% 
  ungroup(numID) %>% 
  filter(relpos=='ipsi') %>% 
  filter(tp %in% c(0.5,2,7.5))

d %>% 
  group_by(clin.meas.name,lab,relpos) %>% 
  nest() %>% 
  mutate(mdl=map(data,~lme(clin.meas.value~GGP, random=~1|numID, data=., na.action = na.omit))) %>% 
  mutate(tidy=map(mdl,tidy)) %>% 
  unnest(tidy) %>% 
  filter(effect=='fixed')
d %>% 
  ggplot(aes(x=GGP, y=clin.meas.value, color=as.factor(tp)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F)+
  geom_smooth(aes(group=numID), method = 'lm', se = F, color='grey', size=.5)+
  facet_wrap(lab~clin.meas.name, scales = 'free')


## GGP - clinical time course linear delta delay adjusted
d<-data.clinical_long %>% 
  merge(d.GGP) %>% 
  group_by(numID,clin.meas.name,lab,relpos) %>% 
  arrange(tp) %>% 
  dplyr::mutate(tp.lag = dplyr::lag(tp, n = 1, default = NA)
                ,tp.delta = tp-tp.lag
                ,tp.mean = (tp+tp.lag)/2
                ,clin.meas.value.lag = (clin.meas.value - dplyr::lag(clin.meas.value, n = 1, default = NA))/tp.delta
                ,GGP.delta = (GGP - dplyr::lag(GGP, n = 1, default = NA))/tp.delta
  ) %>% 
  ungroup(numID) %>% 
  filter(relpos=='ipsi') %>% 
  select(-c(ID,Stroke.side,age,gender,dominant_affected,vol,logvol))

d %>% 
  ggplot(aes(x=GGP.delta, y=clin.meas.value.lag, color=as.factor(tp.mean)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F)+
  geom_smooth(aes(group=numID), method = 'lm', se = F, color='grey', size=.5)+
  facet_wrap(lab~clin.meas.name, scales = 'free')
        

##################
### LGP ##########
##################
## GGP - clinical time course 1-1
d<-d.LGP %>% 
  merge(data.clinical_long) %>% 
  merge(d.conn) %>% 
  filter(relpos=='ipsi') 

d %>% 
  filter(node %in% (nodes.liberal$node %>% as.character())) %>% 
  filter(LGP.name=='efficiency' & clin.meas.name == 'NIHSS') %>% 
  ggplot(aes(y=LGP.value, x=clin.meas.value, color=as.factor(tp)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F)+
  geom_smooth(aes(group=numID), method = 'lm', se = F, color='grey', size=.5)+
  geom_smooth(aes(group=1), method='lm', se = F, color='green')+
  facet_wrap(.~node, scales = 'free', ncol = 4)
  
d %>% ggplot(aes(y=LGP.value, x=clin.meas.value, color=as.factor(tp)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F)+
  geom_smooth(aes(group=numID), method = 'lm', se = F, color='grey', size=.5)+
  geom_smooth(aes(group=1), method='lm', se = F, color='green')+
  facet_wrap(LGP.name~clin.meas.name, scales = 'free', ncol = 4)

d %>% 
  filter(clin.meas.name=='NIHSS' & LGP.name == 'efficiency' & tp %in% c(12)) %>%
  filter(node %in% (nodes.strict$node %>% as.character())) %>% 
  lmer(LGP.value~clin.meas.value+node + (1|numID), data=., na.action=na.omit) %>% 
  tidy()

d %>% 
  filter(clin.meas.name=='rGS' & LGP.name == 'efficiency' & tp %in% c(3,12)) %>% 
  filter(node=='thalamus') %>% 
  filter(node %in% (nodes.liberal$node %>% as.character())) %>% 
  rmcorr(participant = tp, measure1 = LGP.value, measure2 = clin.meas.value) %>% 
  plot()


#### linear delta

d<-data.clinical_long %>% 
  merge(d.LGP) %>% 
  group_by(numID,clin.meas.name,LGP.name,node,relpos) %>% 
  arrange(tp) %>% 
  dplyr::mutate(clin.meas.value = clin.meas.value - dplyr::lag(clin.meas.value, n = 1, default = NA)
                ,LGP.value = LGP.value - dplyr::lag(LGP.value, n = 1, default = NA)
                ,tp = (tp + dplyr::lag(tp, n = 1, default = NA))/2) %>% 
  ungroup() %>% 
  filter(relpos=='ipsi') %>% 
  filter(tp %in% c(0.5,2,7.5))
